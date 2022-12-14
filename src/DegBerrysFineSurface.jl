# struct DegBerrysFineSurface
	# params::DegParams;
	# paramsSurfSlice::DegParams;
	# degMats::DegMatsOnGrid;
	# BfieldLn::Int64;
	# dimLstRev::Vector{Int64};
	
	# enumSaveMem::EnumSaveMem;
	
	# linkLst::Vector{Array{ Vector{ComplexF64} }};
	# BfieldLst::Vector{Array{ Vector{ComplexF64} }};
	# divBLst::Array{ Vector{Complex{Float64}} };
	
	# linkRatioThr::Vector{ThrArray{ComplexF64,1}};
# end

function degBerrysFineSurfaceInit( params::DegParams, ratFine::Int64; enumSaveMem = memNone )
	divLstsurfSlice = params.divLst .* ratFine;
	divLstsurfSlice[end] = ratFine + 1;
	nonPeriodicLst = fill(false, params.nDim);
	nonPeriodicLst[end] = true;
	paramsSurfSlice = degParamsNonInit( params.N, divLstsurfSlice, params.nDim; isNonPeriodic = nonPeriodicLst );
	paramsSurfSlice.stepLst .= params.stepLst ./ ratFine;
	
	matsGrid = matsGridHThreaded( paramsSurfSlice, threaded_zeros( ComplexF64, params.N, params.N ) );
	
	BfieldLn, dimLstRev, linkLst, BfieldLst,  linkRatioThr = initLinkBfield( paramsSurfSlice );
	
	divBLst = Array{Vector{ComplexF64},params.nDim}(undef,params.divLst...);
	divBSurface = zeros(ComplexF64, params.N);
	BfieldLstSurface, _, _ = initBfield( params );
	# BfieldLstSurface = zeros(ComplexF64, 2, params.nDim, params.N);
	
	degBerrys = DegBerrys( params, matsGrid, BfieldLn, dimLstRev, enumSaveMem, linkLst, BfieldLst, divBLst, divBSurface, BfieldLstSurface, linkRatioThr );
	
	degBerrysInitCells( degBerrys, ratFine );
	
	return degBerrys;
end

function degBerrysInitCells( degBerrys::DegBerrys, gapLn::Int64 )
	matsInitCellsExtraLayer( degBerrys.degMats, gapLn );
	Threads.@threads for pos in degBerrys.degMats.params.posLst
		for iDim = 1:degBerrys.degMats.params.nDim
			if pos[degBerrys.degMats.params.nDim] > degBerrys.degMats.params.divLst[end] && iDim != degBerrys.degMats.params.nDim
				degBerrys.linkLst[iDim][pos] = zeros( eltype(eltype(degBerrys.degMats.vLst)), degBerrys.params.N );
			elseif iDim == degBerrys.degMats.params.nDim && pos[iDim] >= degBerrys.degMats.params.divLst[end]
				continue;
			end
			for iDimWall = 1:degBerrys.degMats.params.nDim
				if iDimWall == iDim
					continue;
				elseif mod( pos[iDimWall], gapLn ) == 1
					degBerrys.linkLst[iDim][pos] = zeros( eltype(eltype(degBerrys.degMats.vLst)), degBerrys.params.N );
				end
			end
		end
		iB = 0;
		for iDim = 1 : degBerrys.degMats.params.nDim
			for iDim2 = 1 + iDim : degBerrys.degMats.params.nDim
				iB += 1;
				if iDim2 == degBerrys.degMats.params.nDim && pos[iDim2] >= degBerrys.degMats.params.divLst[iDim2]
					continue;
				elseif pos[degBerrys.degMats.params.nDim] > degBerrys.degMats.params.divLst[end]
					degBerrys.BfieldLst[iB][pos] = zeros( eltype(eltype(degBerrys.degMats.vLst)), degBerrys.params.N );
				end
				for iDimWall = 1:degBerrys.degMats.params.nDim
					if iDimWall == iDim || iDimWall == iDim2
						continue;
					elseif mod( pos[iDimWall], gapLn ) == 1
						degBerrys.BfieldLst[iB][pos] = zeros( eltype(eltype(degBerrys.degMats.vLst)), degBerrys.params.N );
					end
				end
			end
		end
	end
	Threads.@threads for pos in degBerrys.params.posLst
		degBerrys.divBLst[pos] = zeros( eltype(eltype(degBerrys.degMats.vLst)), degBerrys.params.N );
		iB = 0;
		for iDim1 = 1 : degBerrys.params.nDim
			for iDim2 = iDim1+1 : degBerrys.params.nDim
				iB += 1;
				degBerrys.BfieldLstSurface[iB][pos] = zeros( eltype(eltype(degBerrys.degMats.vLst)), degBerrys.params.N );
			end
		end
	end
end

function divBCellLayered( degBerrys::DegBerrys, HmatFun, gapLn )
	degBerrys.degMats.params.minLst .= degBerrys.params.minLst;
	
	vLstShLst = [
		ShiftedArrays.circshift( degBerrys.degMats.vLst, ntuple(x-> x==iDim ? -1 : 0, degBerrys.degMats.params.nDim ) )
		for iDim = 1 : degBerrys.degMats.params.nDim];
	linkLstShLst = [
		ShiftedArrays.circshift( degBerrys.linkLst[iDim2], ntuple( x->x==iDim1 ? -1 : 0, degBerrys.degMats.params.nDim ) )
		for iDim1 = 1:degBerrys.degMats.params.nDim, iDim2 = 1:degBerrys.degMats.params.nDim];
	BfieldLstShLst = [
		ShiftedArrays.circshift( degBerrys.BfieldLst[iDim], ntuple( x->x==degBerrys.dimLstRev[iDim] ? -gapLn : 0, degBerrys.params.nDim ) )
		for iDim = 1 : degBerrys.params.nDim];
	arrCopyLst = [degBerrys.degMats.vLst, degBerrys.BfieldLst[1]];
	for iDim = 1 : degBerrys.degMats.params.nDim - 1
		push!( arrCopyLst, degBerrys.linkLst[iDim] );
	end

	for i3 = 0 : degBerrys.params.divLst[end]
		degBerrys.degMats.params.minLst[end] = (i3-1)*degBerrys.params.stepLst[end];
		refreshMeshGrid( degBerrys.degMats.params );
		
		if i3 == degBerrys.params.divLst[end]
			for arr in arrCopyLst # [degBerrys.degMats.vLst, degBerrys.BfieldLst[1]]
				( (x,y)->(x.=y) ).( selectdim( arr, degBerrys.degMats.params.nDim, degBerrys.degMats.params.divLst[end] ), selectdim( arr, degBerrys.degMats.params.nDim, degBerrys.degMats.params.divLst[end]+1 ) );
			end
		else
			startNextEigen( degBerrys.degMats );
			eigenLayer( degBerrys.degMats, degBerrys.degMats.params.divLst[end]; HmatFun );
			
			Threads.@threads for pos3 in selectdim( degBerrys.degMats.params.posLst, degBerrys.degMats.params.nDim, degBerrys.degMats.params.divLst[end] )
				for iDim = 1 : degBerrys.degMats.params.nDim-1
					dotEachCol!( degBerrys.linkLst[iDim][pos3], vLstShLst[iDim][pos3], degBerrys.degMats.vLst[pos3] );
					degBerrys.linkLst[iDim][pos3] ./= abs.(degBerrys.linkLst[iDim][pos3]);
				end
			end
			
			Threads.@threads for pos3 in selectdim( degBerrys.degMats.params.posLst, degBerrys.degMats.params.nDim, degBerrys.degMats.params.divLst[end] )
				iB = 0;
				for iDim1 = 1 : degBerrys.degMats.params.nDim-1
					for iDim2 = iDim1+1 : degBerrys.degMats.params.nDim-1
						iB += 1;
						parity = (-1)^(iDim1+iDim2-1);
						degBerrys.BfieldLst[iB][pos3] .= parity ./ 1im .* log.( linkLstShLst[iDim1,iDim2][pos3] ./ degBerrys.linkLst[iDim2][pos3] ./ linkLstShLst[iDim2,iDim1][pos3] .* degBerrys.linkLst[iDim1][pos3] );
					end
				end
			end

		end
		
		if i3 == 0
			for arr in arrCopyLst
				((x,y)->(x.=y)).( selectdim(arr, degBerrys.degMats.params.nDim, degBerrys.degMats.params.divLst[end]+1), selectdim(arr, degBerrys.degMats.params.nDim, degBerrys.degMats.params.divLst[end]) );
			end
		else
			Threads.@threads for pos in degBerrys.degMats.params.posLst
				if pos[degBerrys.degMats.params.nDim] == 1 || pos[degBerrys.degMats.params.nDim] >= degBerrys.degMats.params.divLst[end]
					continue;
				end
				for iDim = 1 : degBerrys.degMats.params.nDim-1
					if mod( pos[iDim], gapLn ) == 1
						eigenAtLocForced( pos, degBerrys.degMats, HmatFun );
					end
				end
			end
			Threads.@threads for pos in degBerrys.degMats.params.posLst
				if pos[degBerrys.degMats.params.nDim] >= degBerrys.degMats.params.divLst[end]
					continue;
				end
				for iDim = 1 : degBerrys.degMats.params.nDim-1
					if mod( pos[iDim], gapLn ) == 1
						dotEachCol!( degBerrys.linkLst[end][pos], vLstShLst[end][pos], degBerrys.degMats.vLst[pos] );
						degBerrys.linkLst[end][pos] ./= abs.(degBerrys.linkLst[end][pos]);
					end
				end
			end
			Threads.@threads for pos in degBerrys.degMats.params.posLst
				if pos[degBerrys.degMats.params.nDim] == 1 || pos[degBerrys.degMats.params.nDim] >= degBerrys.degMats.params.divLst[end]
					continue;
				end
				for iDim = 1 : degBerrys.degMats.params.nDim-1
					for iDimWall = 1 : degBerrys.degMats.params.nDim-1
						if iDimWall == iDim
							continue;
						end
						if mod(pos[iDimWall], gapLn) == 1
							dotEachCol!( degBerrys.linkLst[iDim][pos], vLstShLst[iDim][pos], degBerrys.degMats.vLst[pos] );
							degBerrys.linkLst[iDim][pos] ./= abs.(degBerrys.linkLst[iDim][pos]);
						end
					end
				end
			end
			Threads.@threads for pos in degBerrys.degMats.params.posLst
				if pos[degBerrys.degMats.params.nDim] >= degBerrys.degMats.params.divLst[end]
					continue;
				end
				iB = 0;
				for iDim1 = 1 : degBerrys.degMats.params.nDim
					for iDim2 = iDim1+1 : degBerrys.degMats.params.nDim
						iB += 1;
						if iB == 1
							continue;
						end
						for iDimWall = 1 : degBerrys.degMats.params.nDim
							if iDimWall == iDim1 || iDimWall == iDim2
								continue;
							elseif mod( pos[iDimWall], gapLn ) == 1
								parity = (-1)^(iDim1+iDim2-1);
								degBerrys.BfieldLst[iB][pos] .= parity ./ 1im .* log.( linkLstShLst[iDim1,iDim2][pos] ./ degBerrys.linkLst[iDim2][pos] ./ linkLstShLst[iDim2,iDim1][pos] .* degBerrys.linkLst[iDim1][pos] );
							end
						end
					end
				end
			end			
			Threads.@threads for posOrg in selectdim( degBerrys.params.posLst, degBerrys.params.nDim, i3 )
				degBerrys.divBLst[posOrg] .= 0;
				for iB = 1 : degBerrys.BfieldLn
					degBerrys.BfieldLstSurface[iB][posOrg] .= 0;
				end
				for iDim = 1 : degBerrys.params.nDim - 1
					degBerrys.degMats.params.locItThr[iDim] = (posOrg[iDim]-1)*gapLn + 1;
				end
				degBerrys.degMats.params.locItThr[end] = 1;
				wrapCurrentIdVec( degBerrys.degMats.params );
				iB = 0;
				for iDim1 = 1 : degBerrys.params.nDim
					for iDim2 = iDim1+1 : degBerrys.params.nDim
						iB += 1;
						# for iBnd = 1:2
							for i1 = 1:gapLn
								for i2 = 1:gapLn
									getCurrentLinId( degBerrys.degMats.params );
									linId = getThrInst( degBerrys.degMats.params.linIdThr );
									degBerrys.BfieldLstSurface[iB][posOrg] .+= degBerrys.BfieldLst[iB][linId];
									# degBerrys.divBLst[posOrg] .-= degBerrys.BfieldLst[iB][linId];
									# degBerrys.divBLst[posOrg] .+= BfieldLstShLst[iB][linId];
									# @infiltrate
									degBerrys.degMats.params.locItThr[iDim2] += 1;
									wrapCurrentIdVec( degBerrys.degMats.params );
								end
								degBerrys.degMats.params.locItThr[iDim2] -= gapLn;
								wrapCurrentIdVec( degBerrys.degMats.params );
								degBerrys.degMats.params.locItThr[iDim1] += 1;
								wrapCurrentIdVec( degBerrys.degMats.params );
							end
							degBerrys.degMats.params.locItThr[iDim1] -= gapLn;
						wrapCurrentIdVec( degBerrys.degMats.params );
					end
				end
			end
		end
		
		divBCalcAll( degBerrys, degBerrys.BfieldLstSurface );
		
		for arr in arrCopyLst
			( (x,y)->(x.=y) ).( selectdim(arr, degBerrys.degMats.params.nDim, 1), selectdim(arr, degBerrys.degMats.params.nDim, degBerrys.degMats.params.divLst[end]) );
		end
	end
end
