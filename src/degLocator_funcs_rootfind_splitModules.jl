using NelderMeadTest

struct DegSimplices
	params::DegParams;
	lnSimp::Int64;
	lnSimpExtra::Int64;
	lnSimpAll::Int64;
	
	shLstSimp::Array{Int64,2};
	shLstSelf::Vector{Int64};
	
	eigMatsThr::ThrStruct{EigTmpMats};
	
	thresBndRatio::Float64;
	hashIds::Array{Vector{Vector{Float64}}};
	
	# ElstThr = ThrArray{Float64,1};
	# vLstThr = ThrArray{ComplexF64,2};
	# HlstThr = ThrArray{ComplexF64,2};
end

function DegSimplices( params::DegParams; HlstThr=nothing, thresBndRatio = 10 )
	lnSimp = 2^(params.nDim-1);
	lnSimpExtra = (params.nDim>=3) ? 1 : 0;
	lnSimpAll = lnSimp + lnSimpExtra;
	 
	shLstSimp = ones(Int64, params.nDim, lnSimpAll);
	id = 1;
	for iDim1 = 1:params.nDim, iDim2 = iDim1+1:params.nDim
		id+=1;
		shLstSimp[iDim1,id] = -1;
		shLstSimp[iDim2,id] = -1;
	end
	
	shLstSelf = zeros(Int64, lnSimpAll);
	shLstSelf[end] = 1;
	
	eigMatsThr = thrStructCopy( eigTmpMatsInit( params.N ) );
	
	hashIds = [
		Vector{Vector{Float64}}(undef,0)
		for pos in params.posLst];
	# ElstThr = threaded_zeros(params.N);
	# vLstThr = threaded_zeros(ComplexF64,params.N,params.N);
	# if isnothing(HlstThr)
		# HlstThr = threaded_zeros(ComplexF64,params.N,params.N);
	# end
	
	return DegSimplices( params, lnSimp, lnSimpExtra, lnSimpAll, shLstSimp, shLstSelf, eigMatsThr, thresBndRatio, hashIds );
end

function clearHashIds( degSmplx::DegSimplices )
	for pos in degSmplx.params.posLst
		empty!( degSmplx.hashIds[pos] );
	end
end

function tmpArrsDegRootFind( params::DegParams )
	HmatThr = threaded_zeros( ComplexF64, params.N, params.N );
	degMats = matsGridHThreaded( params, HmatThr );
	nmArrsThr = thrStructCopy( nmArrsConstruct( degMats.params.nDim ) );
	degSmplx = DegSimplices( params );
	
	return [degMats, degSmplx, nmArrsThr];
end

function locateRootFind( degMats::DegMatsOnGrid, degSmplx::DegSimplices, nmArrsThr::ThrStruct{NelderMeadArrs}; HmatFun, thresVal = 1e-9, thresSz = 1e-9, thresRelaxRatio = 30 )
	@info("prelim Eigen:")
	startNextEigen( degMats );
	Utils.@timeInfo eigenAll( degMats; HmatFun = HmatFun );
	
	locCount = zeros(Int64, degMats.params.N);
	locLst = Vector{Array{Float64,2}}(undef,degMats.params.N);
	
	locLstRaw = zeros(degMats.params.nDim,degSmplx.lnSimpAll,degMats.params.divLst...,degMats.params.N-1);
	gapLstRaw = zeros(degSmplx.lnSimpAll,degMats.params.divLst...,degMats.params.N-1);
	
	@info("rootFind: ")
	Utils.@timeInfo begin
	for n = 1 : degMats.params.N - 1
		gapFunN = (xLst) -> gapFunThr( xLst, HmatFun, degSmplx.eigMatsThr, n );
		Threads.@threads for pos in degMats.params.posLst
			for iSimp = 1 : degSmplx.lnSimpAll
				degSimpPtsInit!( 
				getThrInst(nmArrsThr), 
					degSmplx, pos, iSimp);
				nmOpt!( getThrInst(nmArrsThr), gapFunN; thresVal = thresVal, thresSz = thresSz );
				@view(locLstRaw[:,iSimp,pos,n]) .= 
					@view( getThrInst(nmArrsThr).ptLst[:,1] );
				gapLstRaw[iSimp,pos,n] = 
					getThrInst(nmArrsThr).valLst[1];
			end
		end
	end
	end
	
	@info( "pick gap 0 points" )
	Utils.@timeInfo locLst0 = pick0gapLocs( locLstRaw, gapLstRaw, degSmplx, thresRelaxRatio * thresVal );
	
	@info( "remove duplicates" )
	Utils.@timeInfo locLst0Distilled = locLstCollisionRemove( locLst0, degSmplx, thresRelaxRatio * thresSz );
	Nlst = [ size(locLst0Distilled[n],1) for n = 1:degSmplx.params.N-1 ];
	@infiltrate
	
	return Nlst, Nlst, locLst0Distilled, locLst0Distilled;
end

function gapFunThr( xLst, HmatFun, eigMatsThr::ThrStruct{EigTmpMats}, n )
	return gapFun( xLst, HmatFun, getThrInst(eigMatsThr), n );
end

function gapFun( xLst, HmatFun, eigMats::EigTmpMats, n::Int64 )
	HmatFun( eigMats.Hmat, xLst );
	eigOnTmpMats!( eigMats );
	
	return abs( eigMats.Elst[n+1] - eigMats.Elst[n] );
end

function degSimpPtsInit!( nmArrs::NelderMeadArrs, degSmplx::DegSimplices, pos, iSimp )
	for iDim = 1 : degSmplx.params.nDim
		@view(nmArrs.ptLstOrg[iDim,:]) .= pos[iDim];
		nmArrs.ptLstOrg[iDim,end] += degSmplx.shLstSelf[iSimp];
	end
	for iDim2 = 1 : degSmplx.params.nDim
		nmArrs.ptLstOrg[iDim2,iDim2] +=  degSmplx.shLstSimp[iDim2,iSimp];
	end
	nmArrs.ptLst .= nmArrs.ptLstOrg;
end

function pick0gapLocsNonStruct( locLstRaw, gapLstRaw, mSz, nDim, thresVal )
	locLst = Vector{Array{Float64,2}}(undef,mSz-1);
	for n = 1:mSz-1
		zerosIdLst = findall( x->x<thresVal, selectdim(gapLstRaw,1+nDim+1,n) );
		locLst[n] = zeros(nDim, length(zerosIdLst));
		for ii = 1 : length(zerosIdLst), iDim = 1 : nDim
			locLst[n][iDim,ii] = 
				locLstRaw[iDim, zerosIdLst[ii], n];
		end
	end
	return locLst;
end

# function pick0gapLocs( locLstRaw, gapLstRaw, degSmplx::DegSimplices, thresVal )
	# locLst = Vector{Array{Float64,2}}(undef,degSmplx.params.N-1);
	# for n = 1:degSmplx.params.N-1
		# zerosIdLst = findall( x->x<thresVal, selectdim(gapLstRaw,1+degSmplx.params.nDim+1,n) );
		# locLst[n] = zeros(degSmplx.params.nDim, length(zerosIdLst));
		# for ii = 1 : length(zerosIdLst), iDim = 1 : degSmplx.params.nDim
			# locLst[n][iDim,ii] = 
				# locLstRaw[iDim, zerosIdLst[ii], n];
		# end
	# end
	# return locLst;
# end

function pick0gapLocs( locLstRaw, gapLstRaw, degSmplx::DegSimplices, thresVal )
	return pick0gapLocsNonStruct( locLstRaw, gapLstRaw, degSmplx.params.N, degSmplx.params.nDim, thresVal );
end

function pick0gapLocsFromFile( mSz, divLst, itNum, seed; thresVal, thresSz, fMod = "", dim = 3, fExt = jld2Type, thresValRatio = 30 )
	attrMoreLst = ["thresVal", "thresSz"];
	valMoreLst = [thresVal, thresSz];
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; attrMoreLst = attrMoreLst, valMoreLst = valMoreLst, dim = dim );
	fMain = "locsRaw";
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	locLstRawLst = load( fName, "locLstRawLst" );
	gapLstRawLst = load( fName, "gapLstRawLst" );
	locLstNdim = ndims(locLstRawLst);
	gapLstNdim = ndims(gapLstRawLst);
	
	locLst0Lst = Vector{Vector{Array{Float64,2}}}(undef,itNum);
	
	Threads.@threads for it = 1 : itNum
		locLstRaw = selectdim( locLstRawLst, locLstNdim, it );
		gapLstRaw = selectdim( gapLstRawLst, gapLstNdim, it );
		locLst0Lst[it] = pick0gapLocsNonStruct(locLstRaw, gapLstRaw, mSz, dim, thresValRatio * thresVal );
	end
	
	# locLst0 = pick0gapLocsNonStruct(locLstRaw, gapLstRaw, mSz, dim, thresValRatio * thresVal );
	
	push!( attrLst, "thresValRatio" );
	push!( valLst, thresValRatio );
	# @infiltrate
	fOutMain = "locLst0Lst";
	fOutName = fNameFunc( fOutMain, attrLst, valLst, fExt; fMod = fMod );
	save( fOutName, "locLst0Lst", locLst0Lst );
end

function locLstCollisionRemove( locLst, degSmplx, thresSz )
	idInt = zeros(Int64, degSmplx.params.nDim);
	idRemain = zeros(degSmplx.params.nDim);
	idDiff = zeros(degSmplx.params.nDim);
	shBndLst = zeros(Int64, degSmplx.params.nDim, 2^degSmplx.params.nDim);
	
	locLstOfLstDistilled = [
		Vector{Vector{Float64}}(undef,0)
		for n = 1 : degSmplx.params.N-1];
	
	for n = 1 : degSmplx.params.N-1
		clearHashIds( degSmplx );
		for ii = 1 : size(locLst[n],2)
			shBndLst .= 0;
			loc = @view( locLst[n][:,ii] );
			getThrInst( degSmplx.params.locItThr ) .= 
			# idInt .= 
				floor.( loc ./ degSmplx.params.stepLst );
			idRemain .= loc .- getThrInst( degSmplx.params.locItThr ) .* degSmplx.params.stepLst;
			getThrInst( degSmplx.params.locItThr ) .= mod.( getThrInst( degSmplx.params.locItThr ), degSmplx.params.divLst ) .+ 1;
			
			isCollided = false;
			idLst = degSmplx.hashIds[ linItCurrent( degSmplx.params ) ];
			for iCol = 1 : length(idLst)
				idDiff .= mod.(loc .- idLst[iCol] .+ degSmplx.params.maxLst ./ 2, degSmplx.params.maxLst) .- degSmplx.params.maxLst ./ 2;
				idDiff 
				if norm(idDiff) < thresSz
					isCollided = true;
					break;
				end
			end
			if !isCollided
				push!( locLstOfLstDistilled[n], loc );
				nBnd = 1;
				shBnd = 0;
				for iDim = 1 : degSmplx.params.nDim
					if idRemain[iDim] < degSmplx.thresBndRatio * thresSz
						shBnd = -1;
					elseif idRemain[iDim] > degSmplx.params.stepLst[iDim] - degSmplx.thresBndRatio * thresSz
						shBnd = 1;
					end
					if shBnd != 0
						for iBnd = 1 : nBnd
							@view(shBndLst[:,iBnd+nBnd]) .= @view(shBndLst[:,iBnd]);
							shBndLst[iDim,iBnd+nBnd] = shBnd;
						end
						nBnd *= 2;
					end
				end
				for iBnd = 1 : nBnd
					locBnd = @view(shBndLst[:,iBnd]);
					locBnd .+= getThrInst( degSmplx.params.locItThr );
					wrapIdVec!(locBnd, degSmplx.params);
					idLinBnd = linIdFromIdVec( locBnd, degSmplx.params );
					push!( degSmplx.hashIds[idLinBnd], copy(loc) );
				end
			end
		end
	end
	
	locLstDistilled = [
		[ locLstOfLstDistilled[n][iId][iDim] for iId = 1 : length(locLstOfLstDistilled[n]), iDim = 1 : degSmplx.params.nDim ]
		for n = 1 : degSmplx.params.N - 1]; 
	return locLstDistilled;
end

function locLstCollisionRemoveFromFile( mSz, divLst, itNum, seed; fMod = "", fExt = jld2Type, thresVal, thresSz, thresValRatio, thresSzRatio, dim = 3 )
	attrMoreLst = ["thresVal", "thresSz", "thresValRatio"];
	valMoreLst = Any[thresVal, thresSz, thresValRatio];
	
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim = dim, attrMoreLst = attrMoreLst, valMoreLst = valMoreLst );
	fMain = "locLst0Lst";
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	locLst0Lst = load(fName, "locLst0Lst");
	
	minNum = 0;
	maxNum = 2*pi;
	paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, dim );
	degSmplx = DegSimplices( paramsFull );
	degSmplxThr = thrStructCopy( degSmplx );
	
	locLstDistilledLst = Vector{Vector{Array{Float64,2}}}(undef,itNum);
	NLstDistilledLst = zeros( Int64, itNum );
	
	
	Threads.@threads for it = 1 : itNum
		locLstDistilledLst[it] = locLstCollisionRemove( locLst0Lst[it], getThrInst( degSmplxThr ), thresSzRatio * thresSz );
	end
	
	Threads.@threads for it = 1 : itNum
		NLstDistilledLst[it] = sum( (x->size(x,1)).( locLstDistilledLst[it] ) );
	end
	
	# locLstDistilled = locLstCollisionRemoveFromFile( locLst0, degSmplx, thresSzRatio * thresSz );
	
	fOutMain = "locLstRootDistilled";
	push!( attrLst, "thresSzRatio" );
	push!( valLst, thresSzRatio );
	
	fName = fNameFunc( fOutMain, attrLst, valLst, fExt; fMod = fMod );
	save( fName, "locLstDistilledLst", locLstDistilledLst );
	
	fOutMain = "NLstRootDistilled";
	fName = fNameFunc( fOutMain, attrLst, valLst, fExt; fMod = fMod );
	save( fName, "NLstDistilledLst", NLstDistilledLst );
end

function locateRootFindRaw( degMats::DegMatsOnGrid, degSmplx::DegSimplices, nmArrsThr::ThrStruct{NelderMeadArrs}; HmatFun, thresVal = 1e-9, thresSz = 1e-9 )
	@info("prelim Eigen:")
	startNextEigen( degMats );
	Utils.@timeInfo eigenAll( degMats; HmatFun = HmatFun );
	
	locCount = zeros(Int64, degMats.params.N);
	locLst = Vector{Array{Float64,2}}(undef,degMats.params.N);
	
	locLstRaw = zeros(degMats.params.nDim,degSmplx.lnSimpAll,degMats.params.divLst...,degMats.params.N-1);
	gapLstRaw = zeros(degSmplx.lnSimpAll,degMats.params.divLst...,degMats.params.N-1);
	
	@info("rootFind: ")
	Utils.@timeInfo begin
	for n = 1 : degMats.params.N - 1
		gapFunN = (xLst) -> gapFunThr( xLst, HmatFun, degSmplx.eigMatsThr, n );
		Threads.@threads for pos in degMats.params.posLst
			for iSimp = 1 : degSmplx.lnSimpAll
				# print("\r",pos,", iSimp = $iSimp")
				degSimpPtsInit!( 
				getThrInst(nmArrsThr), 
					degSmplx, pos, iSimp);
				nmOpt!( getThrInst(nmArrsThr), gapFunN; thresVal = thresVal, thresSz = thresSz );
				@view(locLstRaw[:,iSimp,pos,n]) .= 
					@view( getThrInst(nmArrsThr).ptLst[:,1] );
				gapLstRaw[iSimp,pos,n] = 
					getThrInst(nmArrsThr).valLst[1];
			end
		end
	end
	end
	
	return locLstRaw, gapLstRaw;
end
