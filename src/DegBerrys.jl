struct DegBerrys
	params::DegParams;
	degMats::DegMatsOnGrid;
	BfieldLn::Int64;
	dimLstRev::Vector{Int64};
	
	enumSaveMem::EnumSaveMem;
	
	linkLst::Vector{Array{ Vector{ComplexF64} }};
	BfieldLst::Vector{Array{ Vector{ComplexF64} }};
	divBLst::Array{ Vector{Complex{Float64}} };
	
	divBSurface::Vector{ComplexF64};
	BfieldLstSurface::Array{ComplexF64};
	
	linkRatioThr::Vector{ThrArray{ComplexF64,1}};
end

function calcBfieldDims( nDim::Int64 )
	BfieldLn = Int64( nDim * ( nDim -1 ) / 2 );
	dimLstRev = [nDim:-1:1;];
	return BfieldLn, dimLstRev;
end

function degBerrysInit( params::DegParams, degMats::DegMatsOnGrid; isFullInit = false, enumSaveMem = memNone )
	BfieldLn, dimLstRev, linkLst, BfieldLst,  linkRatioThr = initLinkBfield( params );
	
	divBLst = Array{Vector{Float64},params.nDim}(undef,params.divLst...);
	divBSurface = zeros(ComplexF64, params.N);
	BfieldLstSurface = zeros(ComplexF64, 2, params.nDim, params.N);
	
	degBerrys = DegBerrys( params, degMats, BfieldLn, dimLstRev, enumSaveMem, linkLst, BfieldLst, divBLst, divBSurface, BfieldLstSurface, linkRatioThr );
	
	if isFullInit
		degBerrysArrsInitFull( degBerrys );
	end
	
	return degBerrys;
end

function initLinkBfield( params::DegParams )
	BfieldLn, dimLstRev = calcBfieldDims( params.nDim );
	linkLst = Vector{Array{ Vector{ComplexF64} , params.nDim}}(undef,params.nDim);
	BfieldLst = Vector{Array{ Vector{ComplexF64}, params.nDim}}(undef, BfieldLn);
	
	szLst = ones(Int64,params.nDim);
	for iDim = 1 : params.nDim
		szLst .= params.divLst .+ params.nonPeriodicLst;
		szLst[iDim] = params.divLst[iDim];
		linkLst[iDim] = 
				Array{ Vector{ComplexF64}, params.nDim}(undef,szLst...);
	end
	iB = 1;
	for iDim1 = 1 : params.nDim
		for iDim2 = iDim1+1 : params.nDim
			szLst .= params.divLst .+ params.nonPeriodicLst;
			szLst[iDim1] = params.divLst[iDim1];
			szLst[iDim2] = params.divLst[iDim2];
			BfieldLst[iB] = Array{Vector{ComplexF64}, params.nDim }(undef,szLst...);
			iB += 1;
		end
	end
	
	linkRatioThr = [ threaded_zeros( ComplexF64, params.N ) for iDim = 1 : 2 ];
	
	return BfieldLn, dimLstRev, linkLst, BfieldLst, linkRatioThr;
end

function degBerrysEigLayered( params::DegParams )
	enumSaveMem = memEig;
	
	divLstSlice = deepcopy( params.divLst );
	divLstSlice[end] = 3;
	
	# maxLstSlice = deepcopy( params.minLst );
	# maxLstSlice[end] = params.minLst[end] + 2*params.stepLst[end];
	
	paramsSlice = degParamsNonInit( params.N, divLstSlice, params.nDim );
	
	HmatLst = threaded_zeros( ComplexF64, params.N, params.N );
	
	matsSlice = matsGridHThreaded( paramsSlice, HmatLst );
	matsGridInitAll( matsSlice );
	
	degBerrysInit( params, matsSlice; isFullInit = true, enumSaveMem = enumSaveMem );
end

function degBerrysArrsInitFull( degBerrys::DegBerrys )
	Threads.@threads for iLin = 1 : length(degBerrys.params.posLst)
		if !isassigned( degBerrys.divBLst, iLin )
			degBerrys.divBLst[iLin] = zeros(ComplexF64, degBerrys.params.N);
		end	
		iB = 1;
		for iDim = 1 : degBerrys.params.nDim
			if !isassigned( degBerrys.linkLst[iDim], iLin )
				degBerrys.linkLst[iDim][iLin] = zeros(ComplexF64, degBerrys.params.N);
			end
			for iDim2 = iDim + 1 : degBerrys.params.nDim
				if !isassigned( degBerrys.BfieldLst[iB], iLin )
					degBerrys.BfieldLst[iB][iLin] = zeros(ComplexF64, degBerrys.params.N);
				end
				iB += 1;
			end
		end
	end
end

function linksCalcSurface( degBerrys::DegBerrys )
	for iDim = 1 : degBerrys.params.nDim
		for pos in degBerrys.params.posLst
			isSurface = false;
			if pos[iDim] == size( degBerrys.params.posLst, iDim )
				continue;
			else
				for iSf = 1 : degBerrys.params.nDim
					if iSf != iDim
						isSurface = (isSurface || (pos[iSf] == 1 ) ||
						(pos[iSf] == degBerrys.params.divLst[iSf]+1) );
						if isSurface
							break;
						end
					end
				end
			end
			if !isSurface
				continue;
			end
			setCurrentLoc( degBerrys.params, pos );
			initLinkLst( degBerrys, iDim, pos );
			# @infiltrate
			dotEachCol!( 
				degBerrys.linkLst[iDim][pos], 
				getNextVLst( degBerrys.degMats, iDim ),
				getCurrentVLst(degBerrys.degMats)
				);
			degBerrys.linkLst[iDim][pos] .= 
				degBerrys.linkLst[iDim][pos ] ./ abs.(degBerrys.linkLst[iDim][pos] );
		end
	end
end

function linksCalcAll( degBerrys::DegBerrys )
	for iDim = 1 : degBerrys.params.nDim
		vLstSh = ShiftedArrays.circshift( 
			degBerrys.degMats.vLst, 
			ntuple( i -> i== iDim ? -1 : 0, 
				degBerrys.params.nDim));
		Threads.@threads for pos in degBerrys.params.posLst
			setCurrentLoc( degBerrys.params, pos );
			# initLinkLst( degBerrys, iDim, pos );
			dotEachCol!( 
				degBerrys.linkLst[iDim][pos], 
				# getNextVLst( degBerrys.degMats, iDim ),
				vLstSh[pos],
				# getCurrentVLst(degBerrys.degMats)
				degBerrys.degMats.vLst[pos]
				);
			degBerrys.linkLst[iDim][pos] .= 
				degBerrys.linkLst[iDim][pos] ./ abs.(degBerrys.linkLst[iDim][pos] );
		end
	end
end

function linksCalcAllLayered( degBerrys::DegBerrys, HmatFun )
	i3Slc = 1;
	# minLstSlc = deepcopy( degBerrys.params.minLst );
	# stepLstSlc = deepcopy( degBerrys.params.stepLst );
	# stepLstSlc[end] = 0;
	degBerrys.degMats.params.minLst .= degBerrys.params.minLst;
	degBerrys.degMats.params.stepLst .=degBerrys.params.stepLst;
	degBerrys.degMats.params.stepLst[end] = 0;
	i3Slc = 1;
	for i3 = 1 : degBerrys.params.divLst[end]
		# minLstSlc[end] = 
		degBerrys.degMats.params.minLst[end] = 
		degBerrys.params.minLst[end] + degBerrys.params.stepLst[end] * (i3-1);
		i3Slc = mod(i3,2)+1;
		i3Prev = mod(i3-1,2)+1;
		updateParamsRangeSteps( degBerrys.degMats.params.minLst, degBerrys.degMats.params.stepLst, degBerrys.degMats.params );
		# updateParamsRangeSteps( minLstSlc, stepLstSlc, degBerrys.degMats.params );
		# @infiltrate
		startNextEigen( degBerrys.degMats );
		eigenLayer( degBerrys.degMats, i3Slc; HmatFun = HmatFun );
		
		vLstSlc = selectdim( degBerrys.degMats.vLst, degBerrys.params.nDim, i3Slc );
		if i3 == 1 
			vLstEnd = selectdim( degBerrys.degMats.vLst, degBerrys.params.nDim, degBerrys.degMats.params.divLst[end] );
			Threads.@threads for linId = 1 : length(vLstSlc)
				vLstEnd[linId] .= vLstSlc[linId];
			end
		end
		for iDim = 1 : degBerrys.params.nDim-1
			linkLstSlc = selectdim( degBerrys.linkLst[iDim], degBerrys.params.nDim, i3 );
			vLstSh = ShiftedArrays.circshift( 
			vLstSlc, 
			ntuple( (i -> (i == iDim ? -1 : 0)), degBerrys.params.nDim ) );
			Threads.@threads for linId = 1 : length(linkLstSlc)
				dotEachCol!( 
				linkLstSlc[linId],
				vLstSh[linId],
				vLstSlc[linId]
				);
			end
		end
		# if i3Link != degBerrys.params.divLst[end]
			
		# else
			# vLstSh3 = selectdim( degBerrys.degMats.vLst, degBerrys.nDim, degBerrys.degMats.params.divLst[end] );
		# end
		if i3 != 1
			vLstSh3 = selectdim( degBerrys.degMats.vLst, degBerrys.params.nDim, i3Slc );
			vLstSlc = selectdim( degBerrys.degMats.vLst, degBerrys.params.nDim, i3Prev );
			linkLstSlc = selectdim( degBerrys.linkLst[degBerrys.params.nDim], degBerrys.params.nDim, i3-1 );
			Threads.@threads for linId = 1 : length(vLstSlc)
			dotEachCol!( 
				linkLstSlc[linId],
				vLstSh3[linId],
				vLstSlc[linId]
				);
			end
		end
		if i3 == degBerrys.params.divLst[end]
			vLstSh = selectdim( degBerrys.degMats.vLst, degBerrys.params.nDim, degBerrys.degMats.params.divLst[end] );
			vLstSlc = selectdim( degBerrys.degMats.vLst, degBerrys.params.nDim, i3Slc );
			linkLstSlc = selectdim( degBerrys.linkLst[degBerrys.params.nDim], degBerrys.params.nDim, i3 );
			Threads.@threads for linId = 1 : length(vLstSlc)
			dotEachCol!( 
				linkLstSlc[linId],
				vLstSh[linId],
				vLstSlc[linId]
				);
			end
		end
	end
	
end

function BfieldCalcSurface( degBerrys::DegBerrys )
	iB = 1;
	for iDim1 = 1 : degBerrys.params.nDim
		for iDim2 = iDim1+1 : degBerrys.params.nDim
			iDim3 = degBerrys.dimLstRev[iB];
			for i3 = 1 : size( degBerrys.BfieldLst[iB], iDim3 )-1 : size( degBerrys.params.posLst, iDim3 )
				degBerrys.params.locItThr[iDim3] = i3;
				for i1 = 1 : degBerrys.params.divLst[iDim1], i2 = 1 : degBerrys.params.divLst[iDim2]
					degBerrys.params.locItThr[iDim1] = i1;
					degBerrys.params.locItThr[iDim2] = i2;
					getThrInst( degBerrys.linkRatioThr[1] ).= 
					getLinkLst( degBerrys, iDim2; dimSh = iDim1, iSh = 1 ) ./ 
					getLinkLst( degBerrys, iDim2 );
					getThrInst( degBerrys.linkRatioThr[2] ).= getLinkLst( degBerrys, iDim1; dimSh = iDim2, iSh = 1 ) ./ 
					getLinkLst( degBerrys, iDim1 );
					
					initBfieldLst( degBerrys, iB, getThrInst( degBerrys.params.locItThr ) );
					
					parity = (-1)^(iDim1+iDim2-1);
					getBfieldLst( degBerrys, iB ) .= parity ./ 1im .* 
					log.( getThrInst( degBerrys.linkRatioThr[1] )./ getThrInst( degBerrys.linkRatioThr[2] ) );
				end
			end
			iB += 1;
		end
	end
end

function BfieldCalcAll( degBerrys::DegBerrys )
	iB = 1;
	for iDim1 = 1 : degBerrys.params.nDim
		for iDim2 = iDim1+1 : degBerrys.params.nDim
			parity = (-1)^(iDim1+iDim2-1);
			# @time begin
			linkLstShTmp1 = ShiftedArrays.circshift( degBerrys.linkLst[iDim1], 
			ntuple( i -> i==iDim2 ? -1 : 0, degBerrys.params.nDim ) );
			linkLstShTmp2 = ShiftedArrays.circshift( degBerrys.linkLst[iDim2], 
			ntuple( i -> i==iDim1 ? -1 : 0, degBerrys.params.nDim ) );
			Threads.@threads for pos in degBerrys.params.posLst
				# @time 
				setCurrentLoc( degBerrys.params, pos );
				# @time 
				# getThrInst( degBerrys.linkRatioThr[1] ).= 
				# getLinkLst( degBerrys, iDim2; dimSh = iDim1, iSh = 1 ) ./ 
				# getLinkLst( degBerrys, iDim2 );
				getThrInst( degBerrys.linkRatioThr[1] ).= 
				linkLstShTmp2[pos] ./ 
				degBerrys.linkLst[iDim2][pos];
				# @time 
				# getThrInst( degBerrys.linkRatioThr[2] ).= getLinkLst( degBerrys, iDim1; dimSh = iDim2, iSh = 1 ) ./ 
				# getLinkLst( degBerrys, iDim1 );
				getThrInst( degBerrys.linkRatioThr[2] ).= linkLstShTmp1[pos] ./ 
				degBerrys.linkLst[iDim1][pos];
				
				# @time 
				# initBfieldLst( degBerrys, iB, getThrInst( degBerrys.params.locItThr ) );
				
				# @time 
				# getBfieldLst( degBerrys, iB ) .= parity ./ 1im .* 
				# log.( getThrInst( degBerrys.linkRatioThr[1] )./ getThrInst( degBerrys.linkRatioThr[2] ) );
				degBerrys.BfieldLst[iB][pos] .= parity ./ 1im .* 
				log.( getThrInst( degBerrys.linkRatioThr[1] )./ getThrInst( degBerrys.linkRatioThr[2] ) );
				# @infiltrate 
			end
			# end
			# @infiltrate
			iB += 1;
		end
	end
end

function divBCalcAll( degBerrys::DegBerrys )
	Threads.@threads for pos in degBerrys.params.posLst
		# initDivBLst( degBerrys, pos );
		degBerrys.divBLst[pos] .= 0;
	end
	for iB = 1 : degBerrys.params.nDim
		iDimB = degBerrys.dimLstRev[iB];
		BfieldSh = ShiftedArrays.circshift( 
			degBerrys.BfieldLst[iB], ntuple( i -> i==iDimB ? -1 : 0, degBerrys.params.nDim ) );
		Threads.@threads for pos in degBerrys.params.posLst
			setCurrentLoc( degBerrys.params, pos );
			# degBerrys.divBLst[pos] .+= getBfieldLst( degBerrys, iB; dimSh = iDimB, iSh = 1 );
			degBerrys.divBLst[pos] .+= BfieldSh[pos];
			# degBerrys.divBLst[pos] .-= getBfieldLst( degBerrys, iB );
			degBerrys.divBLst[pos] .-= degBerrys.BfieldLst[iB][pos];
		end
	end
end

function divBCalcSurface( degBerrys::DegBerrys )
	degBerrys.divBSurface .= 0;
	iB = 1;
	for iB = 1 : degBerrys.params.nDim		
		iDimB = degBerrys.dimLstRev[iB];
		@view(degBerrys.BfieldLstSurface[1,iDimB,:]) .= 
			-sum( selectdim( degBerrys.BfieldLst[iB], iDimB, 1 ) );
		@view(degBerrys.BfieldLstSurface[2,iDimB,:]) .= 
			sum( selectdim( degBerrys.BfieldLst[iB], iDimB, size( degBerrys.BfieldLst[iB], iDimB ) ) );
		degBerrys.divBSurface .+= @view(degBerrys.BfieldLstSurface[2,iDimB,:]);
		degBerrys.divBSurface .+= @view( degBerrys.BfieldLstSurface[1,iDimB,:] );
		# degBerrys.divBSurface .+= sum( selectdim( degBerrys.BfieldLst[iB], iDimB, size( degBerrys.BfieldLst[iB], iDimB ) ) );
		# degBerrys.divBSurface .-= sum( selectdim( degBerrys.BfieldLst[iB], iDimB, 1 ) );
	end
	# degBerrys.divBSurface .= sum( degBerrys.BfieldLstSurface );
	degBerrys.divBSurface ./= 2*pi;
end

function divBSurfaceOutput( degBerrys::DegBerrys, HmatFun; transferredResults = false )
	if !transferredResults
		startNextEigen( degBerrys.degMats );
	end
	eigenOnSurface( degBerrys.degMats; HmatFun = HmatFun );
	linksCalcSurface( degBerrys );
	BfieldCalcSurface( degBerrys );
	divBCalcSurface( degBerrys );
end

function initLinkLst( degBerrys::DegBerrys, iDim, pos )
	linId = linIdFromIdVecArr( pos, degBerrys.linkLst[iDim] );
	if !isassigned( degBerrys.linkLst[iDim], linId )
		degBerrys.linkLst[iDim][linId] = zeros( ComplexF64, degBerrys.params.N );
	end
end

function initBfieldLst( degBerrys::DegBerrys, iB, pos )
	isInitted, linId = needInitArr( degBerrys.BfieldLst[iB], pos );
	if !isInitted
		degBerrys.BfieldLst[iB][linId] = zeros(ComplexF64, degBerrys.params.N);
	end
end

function initDivBLst( degBerrys::DegBerrys, pos )
	isInitted, linId = needInitArr( degBerrys.divBLst, pos );
	if !isInitted
		degBerrys.divBLst[linId] = zeros(ComplexF64, degBerrys.params.N);
	end
end

function getLinkLst( degBerrys::DegBerrys, iDim; dimSh = 0, iSh = 0 )
	getCurrentArrSh( degBerrys.linkLst[iDim], degBerrys.params; dimSh = dimSh, iSh = iSh );
end

function getBfieldLst( degBerrys::DegBerrys, iB; dimSh = 0, iSh = 0 )
	getCurrentArrSh( degBerrys.BfieldLst[iB], degBerrys.params; dimSh = dimSh, iSh = iSh );
end
