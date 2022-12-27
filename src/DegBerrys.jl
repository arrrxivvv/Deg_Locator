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
	BfieldLstSurface; # ::Array{ComplexF64};
	
	linkRatioThr::Vector{ThrArray{ComplexF64,1}};
end

function calcBfieldDims( nDim::Int64 )
	BfieldLn = Int64( nDim * ( nDim -1 ) / 2 );
	dimLstRev = [nDim:-1:1;];
	return BfieldLn, dimLstRev;
end

function degBerrysGen( params::DegParams; isFullInit = false, enumSaveMem = memNone, typeElm = ComplexF64 )
	if enumSaveMem == memNone
		paramsEig = params;
	else
		divLstSlice = deepcopy( params.divLst );
		divLstSlice[end] = 3;
		
		paramsEig = degParamsNonInit( params.N, divLstSlice, params.nDim );
	end
	paramsLst = Vector{DegParams}(undef,params.nDim);
	
	HmatLst = threaded_zeros( typeElm, params.N, params.N );
	degMats = matsGridHThreaded( paramsEig, HmatLst; typeElm = typeElm );
	matsGridInitAll( degMats );
	
	iFullStart = 1;
	if enumSaveMem == memEigLink
		iFullStart = 2;
	elseif enumSaveMem == memEigBfield
		iFullStart = 3;
	end
	
	for iP = 1 : iFullStart-1
		paramsLst[iP] = paramsEig;
	end
	for iP = iFullStart : params.nDim
		paramsLst[iP] = params;
	end
	# paramsLst[1:iFullStart-1] .= degMats.params;
	# paramsLst[iFullStart:end] .= paramsEig;
	
	BfieldLn, dimLstRev, linkLst, BfieldLst,  linkRatioThr = initLinkBfield( paramsLst[2] );
	
	divBLst = Array{Vector{Float64},params.nDim}(undef,params.divLst...);
	divBSurface = zeros(ComplexF64, params.N);
	BfieldLstSurface = zeros(ComplexF64, 2, params.nDim, params.N);
	
	degBerrys = DegBerrys( params, degMats, BfieldLn, dimLstRev, enumSaveMem, linkLst, BfieldLst, divBLst, divBSurface, BfieldLstSurface, linkRatioThr );
	
	degBerrysArrsInit( degBerrys );
	
	return degBerrys;
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
	linkLst = Vector{Array{ Vector{ComplexF64} , params.nDim}}(undef,params.nDim);
	
	szLst = ones(Int64,params.nDim);
	for iDim = 1 : params.nDim
		szLst .= params.divLst .+ params.nonPeriodicLst;
		szLst[iDim] = params.divLst[iDim];
		linkLst[iDim] = 
			Array{ Vector{ComplexF64}, params.nDim}(undef,szLst...);
	end
	BfieldLst, BfieldLn, dimLstRev = initBfield( params );
	
	linkRatioThr = [ threaded_zeros( ComplexF64, params.N ) for iDim = 1 : 2 ];
	
	return BfieldLn, dimLstRev, linkLst, BfieldLst, linkRatioThr;
end

function initBfield( params::DegParams )
	BfieldLn, dimLstRev = calcBfieldDims( params.nDim );
	BfieldLst = Vector{Array{ Vector{ComplexF64}, params.nDim}}(undef, BfieldLn);
	
	szLst = ones(Int64,params.nDim);
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
	
	return BfieldLst, BfieldLn, dimLstRev;
end

function degBerrysEigLayered( params::DegParams; typeElm = ComplexF64 )
	enumSaveMem = memEig;
	
	divLstSlice = deepcopy( params.divLst );
	divLstSlice[end] = 3;
	
	paramsSlice = degParamsNonInit( params.N, divLstSlice, params.nDim );
	
	HmatLst = threaded_zeros( typeElm, params.N, params.N );
	
	matsSlice = matsGridHThreaded( paramsSlice, HmatLst; typeElm = typeElm );
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

function degBerrysArrsInit( degBerrys::DegBerrys )
	Threads.@threads for pos in degBerrys.params.posLst
		degBerrys.divBLst[pos] = zeros(ComplexF64, degBerrys.params.N);
	end
	iB = 1;
	for iDim = 1 : degBerrys.params.nDim
		for pos in eachindex(  degBerrys.linkLst[iDim])
			degBerrys.linkLst[iDim][pos] = zeros(ComplexF64, degBerrys.params.N);
		end
		for iDim2 = iDim+1 : degBerrys.params.nDim
			for pos in eachindex( degBerrys.BfieldLst[iB] )
				degBerrys.BfieldLst[iB][pos] = zeros(ComplexF64, degBerrys.params.N);
			end
			iB += 1;
		end
	end
	
	GC.gc();
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
	degBerrys.degMats.params.minLst .= degBerrys.params.minLst;
	degBerrys.degMats.params.stepLst .=degBerrys.params.stepLst;
	degBerrys.degMats.params.stepLst[end] = 0;
	i3Slc = 1;
	for i3 = 1 : degBerrys.params.divLst[end]
		degBerrys.degMats.params.minLst[end] = 
		degBerrys.params.minLst[end] + degBerrys.params.stepLst[end] * (i3-1);
		i3Slc = mod(i3,2)+1;
		i3Prev = mod(i3-1,2)+1;
		updateParamsRangeSteps( degBerrys.degMats.params.minLst, degBerrys.degMats.params.stepLst, degBerrys.degMats.params );
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
	for iDim = 1 : degBerrys.params.nDim
		Threads.@threads for pos in degBerrys.params.posLst
			degBerrys.linkLst[iDim][pos] .= 
				degBerrys.linkLst[iDim][pos] ./ abs.(degBerrys.linkLst[iDim][pos] );
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
				getThrInst( degBerrys.linkRatioThr[1] ).= 
				linkLstShTmp2[pos] ./ 
				degBerrys.linkLst[iDim2][pos];
				getThrInst( degBerrys.linkRatioThr[2] ).= linkLstShTmp1[pos] ./ 
				degBerrys.linkLst[iDim1][pos];
				
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

function divBOutput( degBerrys::DegBerrys, HmatFun )
	if degBerrys.params.nDim == 3 && degBerrys.enumSaveMem == memEigBfield
		@info("Eigen, link, Bfield and divB layered:")
		Utils.@timeInfo divBCalcAllLayered( degBerrys, HmatFun );
	else 
		if degBerrys.enumSaveMem >= memEig 
			@info("Eigen and Link layered:")
			Utils.@timeInfo linksCalcAllLayered( degBerrys, HmatFun );
		else
			@info("Eigen:")
			startNextEigen( degBerrys.degMats );
			Utils.@timeInfo eigenAll( degBerrys.degMats; HmatFun = HmatFun );

			@info("Link:")
			Utils.@timeInfo linksCalcAll( degBerrys );
		end
		@info("Bfield:")
		Utils.@timeInfo BfieldCalcAll( degBerrys );
		if degBerrys.params.nDim >= 3
			@info("DivB:")
			Utils.@timeInfo divBCalcAll( degBerrys );
		end
	end
end

function divBCalcAll( degBerrys::DegBerrys )
	divBCalcAll( degBerrys, degBerrys.BfieldLst );
end

function divBCalcAll( degBerrys::DegBerrys, BfieldLstIn )
	Threads.@threads for pos in degBerrys.params.posLst
		# initDivBLst( degBerrys, pos );
		degBerrys.divBLst[pos] .= 0;
	end
	for iB = 1 : degBerrys.params.nDim
		iDimB = degBerrys.dimLstRev[iB];
		BfieldSh = ShiftedArrays.circshift( 
			BfieldLstIn[iB], ntuple( i -> i==iDimB ? -1 : 0, degBerrys.params.nDim ) );
		Threads.@threads for pos in degBerrys.params.posLst
			setCurrentLoc( degBerrys.params, pos );
			# degBerrys.divBLst[pos] .+= getBfieldLst( degBerrys, iB; dimSh = iDimB, iSh = 1 );
			degBerrys.divBLst[pos] .+= BfieldSh[pos];
			# degBerrys.divBLst[pos] .-= getBfieldLst( degBerrys, iB );
			degBerrys.divBLst[pos] .-= BfieldLstIn[iB][pos];
		end
	end
end

function divBCalcAllLayered( degBerrys::DegBerrys, HmatFun )
	BfieldLnLayer = Int64( (degBerrys.params.nDim-2) * (degBerrys.params.nDim-1) / 2 );
	posLstSlc = selectdim( degBerrys.params.posLst, degBerrys.params.nDim, 1 );
	
	degBerrys.degMats.params.minLst .= degBerrys.params.minLst;
	degBerrys.degMats.params.stepLst .=degBerrys.params.stepLst;
	degBerrys.degMats.params.stepLst[end] = 0;
	
	arrLayers = ( x -> arrSlcLst( x, 
		degBerrys.degMats.params.nDim, degBerrys.degMats.params.divLst[end] ) );
	arrLayerShs = ( x -> 
		[arrShAllLst( x[i3], degBerrys.degMats.params.nDim-1 ) 
			for i3 = 1 : degBerrys.degMats.params.divLst[end] ] );
	arrLayerBackup! = ( (arr, i3) -> assignArrOfArrs!( arr[end],  arr[i3]) );
	vLstSlcs = arrLayers( degBerrys.degMats.vLst );
	vLstShs = arrLayerShs( vLstSlcs );
	linkLstSlcs = [
		arrLayers( degBerrys.linkLst[iDim] )
		for iDim = 1 : degBerrys.params.nDim];
	BfieldLstSlcs= [
		arrLayers( degBerrys.BfieldLst[iB] )
		for iB = 1 : degBerrys.BfieldLn];
	linkLstShs = [ arrLayerShs( linkLstSlcs[iDim] ) 
		for iDim = 1 : degBerrys.params.nDim ];
	BfieldLstShs = [ arrLayerShs( BfieldLstSlcs[iB] )
		for iB = 1 : degBerrys.BfieldLn];
	BLstForDiv = [[ BfieldLstSlcs[iB][i3] 
		for iB = 1 : degBerrys.BfieldLn ]
		for i3 = 1 : 2 ]; 
	BShLstForDiv = [ similar( BLstForDiv[i3], supertype(eltype(BLstForDiv[i3])) )
		for i3 = 1 : 2 ];
	# for i3 = 1 : 2
		# i3Next = mod(i3,2) + 1;
		# for iB = 1 : BfieldLnLayer
			# BShLstForDiv[i3][iB] = BLstForDiv[i3Next][iB];
		# end
	# end
	for i3 = 1 : 2
		for iB = BfieldLnLayer + 1 : degBerrys.params.nDim
			BShLstForDiv[i3][iB] = BfieldLstShs[iB][i3][degBerrys.dimLstRev[iB]];
		end
	end
	
	i3Slc = 1;
	for i3 = 1:degBerrys.params.divLst[end] + 1
		if i3 <= degBerrys.params.divLst[end]
			i3Slc = mod(i3-1,2)+1;
			i3Prev = mod(i3,2)+1;
		else
			i3Prev = i3Slc;
			i3Slc = degBerrys.degMats.params.divLst[end];
		end
		
		if i3 <= degBerrys.params.divLst[end]
			degBerrys.degMats.params.minLst[end] = degBerrys.params.gridLst[end][i3];
			updateParamsRangeSteps( degBerrys.degMats.params.minLst, degBerrys.degMats.params.stepLst, degBerrys.degMats.params );
			
			startNextEigen( degBerrys.degMats );
			eigenLayer( degBerrys.degMats, i3Slc; HmatFun = HmatFun );
			
			for iDim = 1 : degBerrys.params.nDim-1
				linkFromVLstSh!.( linkLstSlcs[iDim][i3Slc], vLstShs[i3Slc][iDim], vLstSlcs[i3Slc] );
			end
			iB = 1;
			for iDim1 = 1 : degBerrys.params.nDim-1
				for iDim2 = iDim1 + 1 : degBerrys.params.nDim-1
					parity = (-1)^(iDim1+iDim2-1);
					BfieldFromLinkSh!.( BfieldLstSlcs[iB][i3Slc], linkLstSlcs[iDim1][i3Slc], linkLstSlcs[iDim2][i3Slc], linkLstShs[iDim1][i3Slc][iDim2], linkLstShs[iDim2][i3Slc][iDim1], parity );
					iB += 1;
					# ( (b, lSh, l) -> ( b .= lSh ./ l ) ).( BfieldSlcs[iB][i3Slc], linkLstShs[iDim2][i3Slc][iDim1], linkLstSlc[iDim2][i3Slc] );
					# ( (b, lSh, l) -> ( b .*= l ./ lSh ) ).(BfieldSlcs[iB][i3Slc], linkLstSlc[iDim1][i3Slc], linkLstShs[iDim1][i3Slc][iDim2]);
					# iB += 1;
				end
			end
		end
		
		if i3 == 1 
			arrLayerBackup!( vLstSlcs, i3Slc );
			for iDim = 1 : degBerrys.params.nDim
				arrLayerBackup!( linkLstSlcs[iDim], i3Slc );
			end
			for iB = 1 : degBerrys.BfieldLn
				arrLayerBackup!( BfieldLstSlcs[iB], i3Slc );
			end
		end
		if i3 > 1
			iDim = degBerrys.params.nDim;
			linkFromVLstSh!.( linkLstSlcs[iDim][i3Prev], vLstSlcs[i3Slc], vLstSlcs[i3Prev] );
			iDim2 = degBerrys.params.nDim;
			iB = BfieldLnLayer + 1;
			for iDim1 = 1 : degBerrys.params.nDim-1
				parity = (-1)^(iDim1+iDim2-1);
				BfieldFromLinkSh!.( BfieldLstSlcs[iB][i3Prev], linkLstSlcs[iDim1][i3Prev], linkLstSlcs[iDim2][i3Prev], linkLstSlcs[iDim1][i3Slc], linkLstShs[iDim2][i3Prev][iDim1], parity );
				iB += 1;
			end
			divBSlc = selectdim( degBerrys.divBLst, degBerrys.params.nDim, i3-1 );
			
			for iB = 1 : BfieldLnLayer
				BShLstForDiv[i3Prev][iB] = BfieldLstSlcs[iB][i3Slc];
			end
			
			divBFromBSh!( divBSlc, BLstForDiv[i3Prev], BShLstForDiv[i3Prev], posLstSlc );
		end
		# if i3 == degBerrys.params.divLst[end]
			# i3Prev = i3Slc;
			# i3Slc = degBerrys.degMats.params.divLst[end];
			
			# iDim = degBerrys.params.nDim;
			# linkFromVLstSh!.( linkLstSlcs[iDim][i3Prev], vLstSlcs[i3Slc], vLstSlcs[i3Prev] );
			# iDim2 = degBerrys.params.nDim;
			# iB = BfieldLnLayer + 1;
			# for iDim1 = 1 : degBerrys.params.nDim-1
				# parity = (-1)^(iDim1+iDim2-1);
				# BfieldFromLinkSh!( BfieldSlcs[iB][i3Prev], linkLstSlcs[iDim1][i3Prev], linkLstSlcs[iDim2][i3Prev], linkLstSlcs[iDim1][i3Slc], linkLstShs[iDim2][i3Prev][iDim1], parity );
				# iB += 1;
			# end
			# divBSlc = selectdim( degBerrys.divBLst, degBerrys.params.nDim, i3 );
			
			# divBFromBSh!( divBSlc, BLstForDiv[i3Prev], BShLstForDiv[i3Prev] );
		# end
	end
end

function linkFromVLstSh!( link, vSh, v )
	dotEachCol!( link, vSh, v );
	link ./= abs.(link);
end

function BfieldFromLinkSh!( Bfield, link1, link2, link1Sh2, link2Sh1, parity )
	Bfield .= parity .* 1/1im .* log.( link2Sh1 ./ link2 ./ link1Sh2 .* link1 );
end

function divBFromBSh!( divB, BLst, BShLst, posLstSlc )
	for pos in posLstSlc
		divB[pos] .= 0;
		for iDim = 1 : 3
			divB[pos] .+= BShLst[iDim][pos];
			divB[pos] .-= BLst[iDim][pos];
		end
	end
	# divB .= 0;
	# for iDim = 1 : 3
		# divB .+= BShLst[iDim];
		# divB .-= BLst[iDim];
	# end
end

# function divBCalcAll( degBerrys::DegBerrys )
	# Threads.@threads for pos in degBerrys.params.posLst
		# # initDivBLst( degBerrys, pos );
		# degBerrys.divBLst[pos] .= 0;
	# end
	# for iB = 1 : degBerrys.params.nDim
		# iDimB = degBerrys.dimLstRev[iB];
		# BfieldSh = ShiftedArrays.circshift( 
			# degBerrys.BfieldLst[iB], ntuple( i -> i==iDimB ? -1 : 0, degBerrys.params.nDim ) );
		# Threads.@threads for pos in degBerrys.params.posLst
			# setCurrentLoc( degBerrys.params, pos );
			# # degBerrys.divBLst[pos] .+= getBfieldLst( degBerrys, iB; dimSh = iDimB, iSh = 1 );
			# degBerrys.divBLst[pos] .+= BfieldSh[pos];
			# # degBerrys.divBLst[pos] .-= getBfieldLst( degBerrys, iB );
			# degBerrys.divBLst[pos] .-= degBerrys.BfieldLst[iB][pos];
		# end
	# end
# end

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

function getBfieldLst( degBerrys::DegBerrys, iB, params::DegParams; dimSh = 0, iSh = 0 )
	getCurrentArrSh( degBerrys.BfieldLst[iB], params; dimSh = dimSh, iSh = iSh );
end
