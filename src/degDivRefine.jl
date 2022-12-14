
function degDivRefineFromFile( mSz, divLst, itNum, seed; fMod = "", fOutMod = "", dim = 3, fExt = jld2Type, numRes = 3, thres = 0.1 )
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim=dim );
	fMain = fDeg;
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	posLocLst = load( fName, varNamePosLoc );
	negLocLst = load( fName, varNameNegLoc );
	HLstLst = load( fName, varNameHlst );
	HLstArr = Array{Matrix{ComplexF64},3}(undef,2,dim,itNum);
	HLstArr = [
		HLstLst[it][iCos,iDim]
		for iCos = 1:2, iDim = 1 : dim, it = 1 : itNum];
	
	locLstPol = [posLocLst,negLocLst];
	
	param_min_num = 0;
	param_max_num = 2*pi;
	params = degParamsPeriodic( mSz, divLst, param_min_num, param_max_num, dim );
	
	HmatThr = threaded_zeros( ComplexF64, params.N, params.N );
	degMatsGrid = matsGridHThreaded( params, HmatThr );
	
	divBLst = [ [ 
		[ zeros( size(locLstPol[iPol][it][lev],1), numRes ) for lev = 1:mSz ] 
			for it = 1:itNum ] 
				for iPol = 1:2 ];
	BfieldLst = [
		[ [ zeros( 2, 3, size( locLstPol[iPol][it][n], 1 ), numRes )
			for n = 1 : mSz ]
				for it = 1 : itNum ] 
					for iPol = 1 : 2];
	divBErrCntLst = zeros( Int64, mSz, itNum, 2 );
	
	divLstCube = fill(2,params.nDim);
	
	matsGridCubeLst = Vector{DegMatsOnGrid}(undef, numRes);
	degBerrysLst = Vector{DegBerrys}(undef,numRes);
	
	divLstCube = ones(Int64,params.nDim);
	for iRes = 1 : numRes
		divLstCube .= 2^(iRes-1);
		matsGridCubeLst[iRes] = 
			matsGridHThreaded( 
				degParamsNonInit( params.N, divLstCube, params.nDim; isNonPeriodic = true ), 
				HmatThr
			);
		degBerrysLst[iRes] = 
			degBerrysInit( 
				matsGridCubeLst[iRes].params, 
				matsGridCubeLst[iRes] );
	end
	
	numSh = 2^params.nDim;
	minLstTmp = zeros(params.nDim);
	maxLstTmp = zeros(params.nDim);
	iShTmp = zeros(Int64,params.nDim);
	for it = 1 : itNum
		# println( "it = ", it, "/", itNum );
		HLst = HLstLst[it];
		# HLst = selectdim( HLstArr, 3, it );
		HmatFun = ( (H,xLst) -> Hmat_3comb!(H,xLst,HLst) );
		# HmatGridFun reevaluate
		HmatGridFun = (xLst) -> HmatFun( getHLoc( iLst, degMatsGrid ), 
			degMatsGrid.mesh[
				linIdFromIdVec( iLst, params )] 
		);
		startNextEigen( degMatsGrid );
		for n = 1 : mSz, iPol = 1 : 2
			print( "it = ", it, " / ", itNum, ", ", "n = ", n, " / ", mSz, "              \r" );
			# println( "n = ", n, );
			locLst = locLstPol[iPol][it][n];
			# Threads.@threads 
			@info( "re eigen corners" )
			Utils.@timeInfo for iLoc = 1 : size(locLst,1)
				loc = @view(locLst[iLoc,:]);
				
				for iSh = 1 : numSh
					setCurrentLoc( degMatsGrid.params, loc );
					shLocIt!( degMatsGrid.params, iSh );
					
					eigenAtCurrentLoc( degMatsGrid; HmatFun = HmatFun );
				end
			end
			@info( "finer eigen each cell" )
			# Threads.@threads 
			Utils.@timeInfo for iLoc = 1 : size(locLst, 1)
				loc = @view(locLst[iLoc,:]);
				linId = linIdFromIdVec( loc, params );
				maxLstTmp .= params.mesh[linId] .+ params.stepLst;
				for iRes = 1 : numRes
					matsCube = matsGridCubeLst[iRes];
					updateParamsRange( params.mesh[linId], maxLstTmp, matsCube.params );
					
					startNextEigen( matsCube );
					if iRes == 1
						matsGridTransfer!( matsCube, degMatsGrid; locS = loc );
					else
						matsGridTransferSurfaceDouble!(matsCube, matsGridCubeLst[iRes-1]);
					end
					
					divBSurfaceOutput( degBerrysLst[iRes], HmatFun; transferredResults = true );
					divBLst[iPol][it][n][iLoc,iRes] = 
					real(degBerrysLst[iRes].divBSurface[n]);
					for iDim = 1 : params.nDim, iBnd = 1 : 2
						BfieldLst[iPol][it][n][iBnd,iDim,iLoc,iRes] = real( degBerrysLst[iRes].BfieldLstSurface[iBnd,iDim,n] );
					end
				end
				if abs(divBLst[iPol][it][n][iLoc,numRes]) < 1-thres
					divBErrCntLst[n,it,iPol] += 1;
				end
			end
		end
	end
	 
	attrLstRes = deepcopy( attrLst );
	valLstRes = deepcopy( valLst );
	push!( attrLstRes, "nRes" );
	push!( valLstRes, numRes );
	fResMain = "divBrefined";
	fOutModFull = fMod;
	if fOutMod != ""
		fOutModFull = fMod * fOutMod;
	end
	fResName = fNameFunc( fResMain, attrLstRes, valLstRes, fExt; fMod = fOutModFull );
	save( fResName, "divBLst", divBLst, "BfieldLst", BfieldLst ); 
	# , "divBErrNeighborLst", divBErrNeighborLst );
end

function divBErrNeighbor( mSz, divLst, itNum, seed; fModLocs = "", fMod = "", dim = 3, numRes = 3, fModOut = "", fExt = jld2Type, thres = 0.1 )
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim=dim, attrMoreLst = ["nRes"], valMoreLst = [numRes] );
	fMain = fDeg;
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fModLocs );
	
	posLocLst = load( fName, varNamePosLoc );
	negLocLst = load( fName, varNameNegLoc );
	HLstLst = load( fName, varNameHlst );
	
	locLstPol = [posLocLst,negLocLst];
	
	param_min_num = 0;
	param_max_num = 2*pi;
	params = degParamsPeriodic( mSz, divLst, param_min_num, param_max_num, dim );
	
	HmatThr = threaded_zeros( ComplexF64, params.N, params.N );
	degMatsGrid = matsGridHThreaded( params, HmatThr );
	
	divLstCube = fill(2,params.nDim);
	
	matsGridCubeLst = Vector{DegMatsOnGrid}(undef, numRes);
	degBerrysLst = Vector{DegBerrys}(undef,numRes);
	
	divLstCube = ones(Int64,params.nDim);
	for iRes = 1 : numRes
		divLstCube .= 2^(iRes-1);
		matsGridCubeLst[iRes] = 
			matsGridHThreaded( 
				degParamsNonInit( params.N, divLstCube, params.nDim; isNonPeriodic = true ), 
				HmatThr
			);
		degBerrysLst[iRes] = 
			degBerrysInit( 
				matsGridCubeLst[iRes].params, 
				matsGridCubeLst[iRes] );
	end
	
	numSh = 2^params.nDim;
	minLstTmp = zeros(params.nDim);
	maxLstTmp = zeros(params.nDim);

	# attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim = dim, attrMoreLst = ["nRes"], valMoreLst = [numRes] );
	fMain = "divBrefined";
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	divBLst = load( fName, "divBLst" );
	
	divBErrCntLst = zeros( Int64, mSz, itNum, 2 );
	
	for iPol = 1:2, it = 1:itNum, n = 1:mSz
		for iLoc = 1 : size(locLstPol[iPol][it][n],1)
			if abs(divBLst[iPol][it][n][iLoc,numRes]) < 1-thres
				divBErrCntLst[n,it,iPol] += 1;
			end			
		end
	end
	
	locLstErrPol = [[[
		ones(Int64, dim, divBErrCntLst[n,it,iPol])
		for n = 1:mSz]
		for it = 1:itNum]
		for iPol = 1:2];
	locIdLstErrPol = [[[
		ones(Int64, divBErrCntLst[n,it,iPol])
		for n = 1 : mSz]
		for it = 1:itNum]
		for iPol = 1:2];
	divBErrNeighborLst = [[[
		zeros( 2, params.nDim, divBErrCntLst[n,it,iPol], 2 )
		for n = 1:mSz]
		for it = 1:itNum]
		for iPol = 1:2];
	for iPol = 1:2, it = 1:itNum, n = 1:mSz
		locLst = locLstPol[iPol][it][n];
		HLst = HLstLst[it];
		HmatFun = ( (H,xLst) -> Hmat_3comb!(H,xLst,HLst) );
		idErr = 0;
		print( "it = ", it, " / ", itNum, ", ", "n = ", n, " / ", mSz, "              \r" );
		for iLoc = 1 : size( locLst, 1 )
			if abs(divBLst[iPol][it][n][iLoc,numRes]) < 1-thres
				idErr += 1;
				loc = @view( locLst[iLoc,:] );
				@views locLstErrPol[iPol][it][n][:,idErr] .= loc;
				locIdLstErrPol[iPol][it][n][idErr] = iLoc;
				for iBnd = 1:2, iDim = 1:params.nDim, iRes = 1 : 2
					if iRes == 1
						nRes = 1;
					else
						nRes = numRes;
					end
					setCurrentLoc( params, loc );
					minLstTmp .= getCurrentArrSh( params.mesh, params; dimSh = iDim, iSh = (-1)^iBnd );
					maxLstTmp .= minLstTmp .+ params.stepLst;
					matsCube = matsGridCubeLst[nRes];
					updateParamsRange( minLstTmp, maxLstTmp, matsCube.params );
					divBSurfaceOutput( degBerrysLst[nRes], HmatFun );
					divBErrNeighborLst[iPol][it][n][iBnd,iDim,idErr,iRes] = real( degBerrysLst[nRes].divBSurface[n] );
				end
			end
		end
	end
	
	# for iPol = 1:2, it = 1:itNum, n = 1:mSz
		# for iErrLoc = 1 : divBErrCntLst[n,it,iPol]
			
		# end
	# end
	
	divBErrNeighborSumLst = [[[
		sumDropDims( divBErrNeighborLst[iPol][it][n][:,:,:,2]-divBErrNeighborLst[iPol][it][n][:,:,:,1]; dims = [1,2] )
		for n = 1 : mSz]
		for it = 1 : itNum]
		for iPol = 1:2];
	
	fOutMain = "divBErrNeighbors";
	fOutName = fNameFunc( fOutMain, attrLst, valLst, fExt; fMod = [fMod, fModOut] );
	save( fOutName, "divBErrNeighborLst", divBErrNeighborLst, "locLstErrPol", locLstErrPol,  "locIdLstErrPol", locIdLstErrPol, "divBErrCntLst", divBErrCntLst, "divBErrNeighborSumLst", divBErrNeighborSumLst );
end

function divBErrClassify( mSz, divLst, itNum, seed; fMod = "", numRes = 3, dim = 3, fExt = jld2Type, fModOut = "" )
	fMain = "divBErrNeighbors";
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim = dim, attrMoreLst = ["nRes"], valMoreLst = [numRes] );
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	divBErrCntLst = load( fName, "divBErrCntLst" );
	locIdLstErrPol = load( fName, "locIdLstErrPol" );
	divBErrNeighborLst = load( fName, "divBErrNeighborLst" );
	
	locIdLstJustMoved = [[[
		zeros(Int64,0)
		for n = 1:mSz] for it = 1:itNum]
		for iPol = 1:2];
	locIdErrLstJustMoved = [[[
		zeros(Int64,0)
		for n = 1:mSz] for it = 1:itNum]
		for iPol = 1:2];
	locIdLstNotJustMoved = deepcopy( locIdLstJustMoved );
	locIdErrLstNotJustMoved = deepcopy(locIdErrLstJustMoved);
	locIdLstElse = deepcopy(locIdLstJustMoved);
	locIdErrLstElse = deepcopy(locIdErrLstJustMoved);
	
	fluxOrgNeighborSum = 0;
	fluxDiffNeighborSum = 0;
	
	for iPol = 1:2, it = 1:itNum, n = 1:mSz
		for iLoc = 1 : divBErrCntLst[iPol,it,n]
			fluxOrgNeighborSum = sum( abs.(divBNeighborLst[iPol][it][n][:,:,iLoc,1]) );
			fluxDiffNeighborSum = sum( abs.(divBNeighborLst[iPol][it][n][:,:,iLoc,end]) ) - 
			fluxOrgNeighborSum;
			if fluxOrgNeighborSum < thres
				if fluxDiffNeighborSum > 1-thres && fluxDiffNeighborSum < 1+thres
					push!( locIdErrLstJustMoved[iPol][it][n], iLoc );
					push!( locIdLstJustMoved[iPol][it][n], locIdLstErrPol[iPol][it][n][iLoc] );
				else
					push!( locIdErrLstNotJustMoved[iPol][it][n], iLoc );
					push!( locIdLstNotJustMoved[iPol][it][n], locIdLstErrPol[iPol][it][n][iLoc] );
				end
			else
				push!( locIdLstElse[iPol][it][n], iLoc );
				push!( locIdErrLstElse[iPol][it][n], locIdLstErrPol[iPol][it][n][iLoc] );
			end
		end
	end
	
	fOutMain = "degDivRefineClassify"
	fOutName = fNameFunc( fOutMain, attrLst, valLst, fExt; fMod = [fMod , fModOut]);
	
	save( fOutName, "locIdLstJustMoved", locIdLstJustMoved, "locIdErrLstJustMoved", locIdErrLstJustMoved, "locIdLstNotJustMoved", locIdLstNotJustMoved, "locIdErrLstNotJustMoved", locIdErrLstNotJustMoved );
end

function degDivRefineClassifyBfields( mSz, divLst, itNum, seed; fMod = "", dim = 3, fExt = jld2 )
	attrMoreLst = ["nRes"]; 
	valMoreLst = [numRes];
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim = dim, attrMoreLst = attrMoreLst, valMoreLst = valMoreLst );
	fDBfieldDivBMain = "divBRefined";
	fName = fNameFunc( fDBfieldDivBMain, attrLst, valLst, fExt; fMod = fMod );
	
	BfieldLst = load( fName, "BfieldLst" );
	
	fNeighborMain = "divBErrNeighborLst";
	fNeighborName = fNameFunc( fNeighborMain, attrLst, valLst, fExt; fMod = fMod );
	divBErrNeighborLst = load( fNeighborName, "divBErrNeighborLst" );
	
	fClassifyMain = "degDivRefineClassify";
	fClassifyName = fNameFunc( fClassifyMain, attrLst, valLst, fExt; fMod = fMod, attrMoreLst = attrMoreLst, valMoreLst = valMoreLst );
	
	locIdLstJustMore = load( fClassifyName, "locIdLstJustMore" );
	locIdErrLstJustMoved = load( fClassifyName, "locIdErrLstJustMore" );
	locIdLstNotJustMore = load( fClassifyName, "locIdLstNotJustMore" );
	locIdLstNotJustMore = load( fClassifyName, "locIdLstNotJustMore" );
	locIdErrLstNotJustMore = load( fClassifyName, "locIdErrLstNotJustMore" );
	
	BDiffMaxErrLstJustMore = [[[
		zeros( length( locIdLstJustMore[iPol][it][n] ) )
		for n = 1:mSz] for it = 1:itNum]
		for iPol = 1:2];
	BMaxErrLstJustMore = [[[
		zeros( length( locIdLstJustMore[iPol][it][n] ) )
		for n = 1:mSz] for it = 1:itNum]
		for iPol = 1:2];
	BDiffMaxErrLstNotJustMore = [[[
		zeros( length( locIdLstJustMore[iPol][it][n] ) )
		for n = 1:mSz] for it = 1:itNum]
		for iPol = 1:2];
	BMaxErrLstNotJustMore = [[[
		zeros( length( locIdLstJustMore[iPol][it][n] ) )
		for n = 1:mSz] for it = 1:itNum]
		for iPol = 1:2];
	for iPol = 1:2, it = 1:itNum, n = 1:mSz
		for iLoc = 1:length( locIdLstJustMore[iPol][it][n] )
			iErr = locIdErrLstJustMore[iPol][it][n][iLoc];
			iOrg = locIdLstJustMore[iPol][it][n][iLoc];
			idMax = argmax( abs.( @view( divBErrNeighborLst[iPol][it][n][:,:,iErr,end] ) ) );
			BDiffMaxErrLstJustMore[iPol][it][n][iLoc] = BfieldLst[iPol][it][n][idMax,iOrg,end] - BfieldLst[iPol][it][n][idMax,iOrg,1];
			BMaxErrLstJustMore[iPol][it][n][iLoc] = BfieldLst[iPol][it][n][idMax,iOrg,1];
		end
	end
	for iPol = 1:2, it = 1:itNum, n = 1:mSz
		for iLoc = 1:length( locIdLstNotJustMore[iPol][it][n] )
			iErr = locIdErrLstNotJustMore[iPol][it][n][iLoc];
			iOrg = locIdLstNotJustMore[iPol][it][n][iLoc];
			idMax = argmax( abs.( @view( divBErrNeighborLst[iPol][it][n][:,:,iErr,end] ) ) );
			BDiffMaxErrLstJustMore[iPol][it][n][iLoc] = BfieldLst[iPol][it][n][idMax,iOrg,end] - BfieldLst[iPol][it][n][idMax,iOrg,1];
			BMaxErrLstJustMore[iPol][it][n][iLoc] = BfieldLst[iPol][it][n][idMax,iOrg,1];
		end
	end
end

function divBErrNeighborSumCalc( mSz, divLst, itNum, seed; fMod = "", numRes = 3, dim = 3, fExt = jld2Type, fModOut = "" )
	fMain = "divBErrNeighbors";
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim = dim, attrMoreLst = ["nRes"], valMoreLst = [numRes] );
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	divBErrNeighborLst = load( fName, "divBErrNeighborLst" );
	divBErrNeighborAbsSumLst = [[[
		sumDropDims( abs.(divBErrNeighborLst[iPol][it][n][:,:,:,2]-divBErrNeighborLst[iPol][it][n][:,:,:,1]); dims = [1,2] )
		for n = 1 : mSz]
		for it = 1 : itNum]
		for iPol = 1:2];
	divBErrNeighborAbsSumVcat = [
		vcat( (x->vcat(x...)).(divBErrNeighborAbsSumLst[iPol]) ... )
		for iPol = 1 :2];
	
	fOutMain = "divBErrNeighborsVcat";
	fOutName = fNameFunc( fOutMain, attrLst, valLst, fExt; fMod = [fMod, fModOut] );
	save( fOutName, "divBErrNeighborAbsSumVcat", divBErrNeighborAbsSumVcat );
end

function divBRefinedStats( mSz, divLst, itNum, seed; fMod = "", dim = 3, numRes = 3, fModOut = "", fExt = jld2Type, thres=0.1 )
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim = dim, attrMoreLst = ["nRes"], valMoreLst = [numRes] );
	fMain = "divBrefined";
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	divBLst = load( fName, "divBLst" );
	BfieldLst = load( fName, "BfieldLst" );
	
	divBLstVcat = [ 
		vcat( ((x->vcat(x...)).( divBLst[iPol]) )... )
		for iPol = 1 : 2 ];
	
	
	BMaxLst = [ [ [ 
		maxDropDims( BfieldLst[iPol][it][n]; dims = [1,2] )
		for n = 1 : mSz ]
		for it = 1 : itNum ]
		for iPol = 1 : 2 ];
	BRemainLst = [ [ [ 
		zeros( 2*dim-1, size( divBLst[iPol][it][n], 1 ), numRes )
		for n = 1 : mSz ]
		for it = 1 : itNum ]
		for iPol = 1 : 2 ];
	BDiffLst = [[[
		BfieldLst[iPol][it][n][:,:,:,end] - BfieldLst[iPol][it][n][:,:,:,1]
		for n = 1 : mSz]
		for it = 1 : itNum]
		for iPol = 1 : 2];
	BDiffLstCat = [ 
		cat( (x->cat(x...; dims=3)).(BDiffLst[iPol]) ...; dims = 3 )
		for iPol = 1 : 2 ];
	BDiffLstVcat = [ reshape(BDiffLstCat[iPol],length(BDiffLstCat[iPol])) for iPol = 1 : 2 ];
	
	for iPol = 1:2, it = 1:itNum, n = 1:mSz
		for iRes = 1 : numRes, iLoc = 1 : size( divBLst[iPol][it][n], 1 )
			isFnd = false;
			valMax = BMaxLst;
			id = 1;
			for iBnd = 1:2, iDim = 1 : dim
				if isFnd || BfieldLst[iPol][it][n][iBnd,iDim,iLoc,iRes] != BMaxLst[iPol][it][n][iLoc,iRes]
					BRemainLst[iPol][it][n][id,iLoc,iRes] = BfieldLst[iPol][it][n][iBnd,iDim,iLoc,iRes];
					id += 1;
				else
					isFnd = true;
				end
			end
		end
	end
	
	BMaxLstVcat = [ 
		vcat( (x->vcat(x...)).( BMaxLst[iPol] )... )
		for iPol = 1:2 ];
	BRemainLstCat2 = [
		cat( (x->cat(x...; dims=2)).( BRemainLst[iPol] ) ...; dims=2 )
		for iPol = 1:2 ];
	BRemainLstVcat = [
		reshape( 
		cat( (x->cat(x...; dims=2)).( BRemainLst[iPol] ) ...; dims=2 ), :, numRes )
		for iPol = 1:2 ];
	BfieldLstVcat = [
		cat( (x->cat(x...; dims=3)).( BfieldLst[iPol] ) ...; dims=3 )
		for iPol = 1:2 ];
	BfieldLstErrVcat = [
		BfieldLstVcat[iPol][:,:,abs.(@view(divBLstVcat[iPol][:,3])).<(1-thres),:]
		for iPol = 1:2];
		
	BMaxLstErrVcat = [ 
		BMaxLstVcat[iPol][
			abs.(@view(divBLstVcat[iPol][:,3])).<(1-thres)
			,:] 
		for iPol = 1:2 ];
	# @infiltrate
	BMaxLstNonErrVcat = [ 
		BMaxLstVcat[iPol][
			abs.(@view(divBLstVcat[iPol][:,3])).>(1-thres)
			,:] 
		for iPol = 1:2 ];
	BRemainLstErrCat2 = [ 
		BRemainLstCat2[iPol][:,
			abs.(@view(divBLstVcat[iPol][:,3])).<(1-thres)
			,:] 
		for iPol = 1:2 ];
	BRemainLstNonErrCat2 = [ 
		BRemainLstCat2[iPol][:,
			abs.(@view(divBLstVcat[iPol][:,3])).>(1-thres)
			,:] 
		for iPol = 1:2 ];
	BRemainLstErrVcat = [
		reshape( BRemainLstErrCat2[iPol], :, 3 )
		for iPol = 1:2];
	BRemainLstNonErrVcat = [
		reshape( BRemainLstNonErrCat2[iPol], :, 3 )
		for iPol = 1:2];
	BDiffLstErrCat = [ 
		BDiffLstCat[iPol][:,:,abs.(@view(divBLstVcat[iPol][:,3])).<(1-thres)]
		for iPol = 1 : 2 ];
	BDiffLstNonErrCat = [ 
		BDiffLstCat[iPol][:,:,abs.(@view(divBLstVcat[iPol][:,3])).>(1-thres)]
		for iPol = 1 : 2 ];
	BMaxDiffLstErrCat = [ 
		(-1)^(iPol-1) *minDropDims( (-1)^(iPol-1) * BDiffLstErrCat[iPol]; dims = [1,2] )
		for iPol = 1 : 2 ];
	
	BRemainDiffLstErrCat = [
		zeros(2*dim-1,length(BMaxDiffLstErrCat[iPol]))
		for iPol = 1 : 2];
	
	# @infiltrate
	for iPol = 1 : 2, ii = 1 : length(BMaxDiffLstErrCat[iPol])
		isFnd = false;
		id = 1;
		for iBnd = 1:2, iDim = 1 : dim
			if isFnd || BDiffLstErrCat[iPol][iBnd,iDim,ii] != BMaxDiffLstErrCat[iPol][ii]
				BRemainDiffLstErrCat[iPol][id,ii] = BDiffLstErrCat[iPol][iBnd,iDim,ii];
				id += 1;
			else
				isFnd = true;
			end
		end
	end
	# @infiltrate
	
	BRemainDiffLstErrVcat = [
		reshape( BRemainDiffLstErrCat[iPol], length(BRemainDiffLstErrCat[iPol]) )
		for iPol = 1 : 2];
	
	fOutMain = "divBRefinedVcat";
	
	fOutName = fNameFunc( fOutMain, attrLst, valLst, fExt; fMod = [fMod, fModOut] );
	
	save( fOutName, "divBLstVcat", divBLstVcat , "BMaxLst", BMaxLst, "BRemainLst", BRemainLst, "BMaxLstVcat", BMaxLstVcat, "BRemainLstVcat", BRemainLstVcat, "BMaxLstErrVcat", BMaxLstErrVcat, "BMaxLstNonErrVcat", BMaxLstNonErrVcat, "BRemainLstErrVcat", BRemainLstErrVcat, "BRemainLstNonErrVcat", BRemainLstNonErrVcat, "BDiffLstVcat", BDiffLstVcat, "BMaxDiffLstErrCat", BMaxDiffLstErrCat, "BRemainDiffLstErrCat", BRemainDiffLstErrCat, "BRemainDiffLstErrVcat", BRemainDiffLstErrVcat);
end
