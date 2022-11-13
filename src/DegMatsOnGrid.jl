# using Infiltrator
@enum HgridStoreOpt HgridFull=1 Hlayered HthreadedOnly

struct DegMatsOnGrid
	params::DegParams
	Elst::Array{Vector{Float64}};
	vLst::Array{Matrix{ComplexF64}};
	
	HmatLst; # ::Array{Matrix{ComplexF64}};
	HgridOpt::HgridStoreOpt;
	eigWorkTh::ThrStruct{EigCustom.EigWork};
	
	eigenId::Base.RefValue{Int64};
	eigenIdGrid::Array{Int64};
	
	function DegMatsOnGrid( params, Elst, vLst, HmatLst, HgridOpt, eigWorkTh )
		eigenId = 0;
		eigenIdGrid = makeArrOverGrid( Int64, params );
		eigenIdGrid .= 0;
		new( params, Elst, vLst, HmatLst, HgridOpt, eigWorkTh, Ref(eigenId), eigenIdGrid );
	end
end

function matsGridBase( params::DegParams, HmatLst, HOpt::HgridStoreOpt )
	Elst = makeArrOverGrid( Vector{Float64}, params );
	vLst = makeArrOverGrid( Matrix{ComplexF64}, params );
	
	eigWorkTh = thrStructCopy( eigWorkStructFromNum!( params.N ) );
	
	return DegMatsOnGrid( params, Elst, vLst, HmatLst, HOpt, eigWorkTh );
end

function matsGridHinternal( params::DegParams, Hopt::HgridStoreOpt )
	if Hopt == HgridFull
		HmatLst = makeArrOverGrid( Matrix{ComplexF64}, params );
	elseif Hopt == Hlayered
		HmatLst = makeArrOverGrid( Matrix{	ComplexF64}, params );
	end
	
	return matsGridBase( params, HmatLst, Hopt );
end

function matsGridHFull( params::DegParams )
	matsGridHinternal( params, HgridFull );
end

function matsGridHExt( params::DegParams, HmatLst, Hopt::HgridStoreOpt )
	return matsGridBase( params, HmatLst, HthreadedOnly );
end

function matsGridHThreaded( params::DegParams, HmatLst )
	return matsGridHExt( params, HmatLst, HthreadedOnly );
end

function getHCurrent( matsGrid::DegMatsOnGrid )
	return getHLoc( getThrInst( matsGrid.params.locItThr ), matsGrid );
end

function getHLoc( loc, matsGrid::DegMatsOnGrid )
	if matsGrid.HgridOpt == HgridFull
		if isa(loc, Vector)
			idLin = linIdFromIdVec(loc);
		else
			idLin = loc;
		end
		if !isassigned( matsGrid.HmatLst, idLin )
			matsGrid.HmatLst[idLin] = zeros(ComplexF64, matsGrid.params.N, matsGrid.params.N);
		end
		return matsGrid.HmatLst[idLin];
	elseif matsGrid.HgridOpt == HthreadedOnly
		return getThrInst( matsGrid.HmatLst );
	end
end

function startNextEigen( matsGrid::DegMatsOnGrid )
	matsGrid.eigenId[] += 1;
end

function checkEigenDone( matsGrid::DegMatsOnGrid, loc )
	if isa(loc, Vector)
		idLin = linIdFromIdVec( loc, matsGrid.params );
	elseif isa(loc, CartesianIndex)
		for iDim = 1 : matsGrid.params.nDim
			matsGrid.params.locItThr[iDim] = loc[iDim];
		end
		idLin = linItCurrent( matsGrid.params );
	else
		idLin = loc;
	end
	if matsGrid.eigenIdGrid[idLin] ==  matsGrid.eigenId[]
		return true;
	else
		matsGrid.eigenIdGrid[idLin] = matsGrid.eigenId[];
		return false;
	end
end

function eigenAtCurrentLoc( matsGrid::DegMatsOnGrid; noReEigen = true, HmatFun = nothing )
	idLin = linItCurrent( matsGrid.params );
	if checkEigenDone( matsGrid, idLin )
		return;
	end
	Hmat = getHLoc( idLin, matsGrid );
	if !isnothing(HmatFun)
		HmatFun(Hmat, matsGrid.params.mesh[idLin]);
	end
	
	if !isassigned( matsGrid.Elst, idLin )
		matsGrid.Elst[idLin] = zeros(matsGrid.params.N);
	end
	if !isassigned( matsGrid.vLst, idLin )
		matsGrid.vLst[idLin] = zeros( ComplexF64, matsGrid.params.N, matsGrid.params.N );
	end
	
	eigenZheevrStruct!( Hmat, matsGrid.Elst[idLin], matsGrid.vLst[idLin], getThrInst(matsGrid.eigWorkTh) );
end

function eigenAtLoc( loc, matsGrid::DegMatsOnGrid; noReEigen = true )
	setCurrentLoc( matsGrid.params, loc );
	eigenAtCurrentLoc( matsGrid; noReEigen = noReEigen );
end

function setCurrentDone( matsGrid::DegMatsOnGrid )
	linId = linItCurrent( matsGrid.params );
	matsGrid.eigenIdGrid[linId] = matsGrid.eigenId[];
end

function getCurrentElst( matsGrid::DegMatsOnGrid )
	linId = linIdFromIdVec( 
		getThrInst( matsGrid.params.locItThr ), matsGrid.params );
	return matsGrid.Elst[linId];
end

function setCurrentElst( matsGrid::DegMatsOnGrid, mat )
	linId = linIdFromIdVec( 
		getThrInst( matsGrid.params.locItThr ), matsGrid.params );
	if !isassigned( matsGrid.Elst, linId )
		matsGrid.Elst[linId] = mat;
	else
		matsGrid.Elst[linId] .= mat;
	end
end

function getCurrentVLst( matsGrid::DegMatsOnGrid )
	linId = linIdFromIdVec( 
		getThrInst( matsGrid.params.locItThr ), matsGrid.params );
	return matsGrid.vLst[linId];
end

# function getNextVLst( matsGrid::DegMatsOnGrid, iDim::Int64 )
	# matsGrid.params.locItThr[iDim] += 1;
	# linId = linIdFromIdVec( 
		# getThrInst( matsGrid.params.locItThr ), matsGrid.params );
	# matsGrid.params.locItThr[iDim] -= 1;
	# return matsGrid.vLst[linId];
# end

function getNextVLst( matsGrid::DegMatsOnGrid, iDim::Int64 )
	return getCurrentArrSh( matsGrid.vLst, matsGrid.params; dimSh = iDim, iSh = 1 );
end

function setCurrentVLst( matsGrid::DegMatsOnGrid, mat )
	linId = linIdFromIdVec( 
		getThrInst( matsGrid.params.locItThr ), matsGrid.params );
	if !isassigned( matsGrid.vLst, linId )
		matsGrid.vLst[linId] = mat;
	else
		matsGrid.vLst[linId] .= mat;
	end
end

function matsGridTransfer!( matsGrT, matsGrS; locS = nothing )
	for iSh = 1 : 2^matsGrT.params.nDim
		if !isnothing(locS)
			shIdVec!( matsGrS.params, locS, iSh );
		else
			locItCorner!( matsGrS.params, iSh );
		end
		locItCorner!( matsGrT.params, iSh );
		# @infiltrate;
		setCurrentElst( matsGrT, 
			getCurrentElst( matsGrS ) );
		setCurrentVLst( matsGrT, 
			getCurrentVLst( matsGrS ) );
		setCurrentDone( matsGrT );
	end
end

function matsGridTransferSurfaceDouble!( matsGrT, matsGrS )
	for pos in matsGrS.params.posLst
		if isPosSurface( pos, matsGrS.params )
			setCurrentLoc( matsGrS.params, pos );
			setCurrentLoc( matsGrT.params, pos );
			setDoubleLoc( matsGrT.params );
			setCurrentElst( matsGrT, 
				getCurrentElst( matsGrS ) );
			setCurrentVLst( matsGrT, 
				getCurrentVLst( matsGrS ) );
			setCurrentDone( matsGrT );
		end
	end
end

function eigenOnSurface( matsGrid::DegMatsOnGrid; HmatFun = nothing )
	Threads.@threads for pos in matsGrid.params.posLst
		isSurface = false;
		for iDim = 1 : matsGrid.params.nDim
			isSurface = ( isSurface 
				|| ( pos[iDim] == 1 )
				|| ( pos[iDim] > matsGrid.params.divLst[iDim] ) );
			if isSurface
				break;
			end
		end
		if !isSurface
			continue;
		end
		if !isnothing(HmatFun)
			HmatFun( getHLoc( pos, matsGrid ), matsGrid.params.mesh[pos] );
		end
		eigenAtLoc( pos, matsGrid );
	end
end

function eigenAll( matsGrid::DegMatsOnGrid; HmatFun = nothing )
	Threads.@threads for pos in matsGrid.params.posLst
		# if checkEigenDone( matsGrid, pos );
			# continue;
		# end
		# setCurrentLoc( matsGrid.params, pos );
		if !isnothing(HmatFun)
			HmatFun( getHLoc( pos, matsGrid ), matsGrid.params.mesh[pos] );
		end
		eigenAtLoc( pos, matsGrid );
	end
end

# function eigenOnSurface( matsGrid::DegMatsOnGrid )
	# for iDim = 1 : matsGrid.params.nDim
		# broadcasrAssign!(matsGrid.params.stepsItThr, 1);
		# matsGrid.params.stepsItThr[iDim] = 
			# matsGrid.params.divLst[iDim];
		# for iDiv = 1 : matsGrid.params.nDim
			# matsGrid.params.divsItThr[iDiv] =  matsGrid.params.divLst[iDiv] + 1;
		# end
		# matsGrid.params.divItThr[iDim] = ;
		# for i1 = 1 : matsGrid.params.stepsItThr[1] : matsGrid.params.divLst[1]+1
			# for i2 = 1 : matsGrid.params.stepsItThr[2] : matsGrid.params.divLst[2] + 1
				# for i3 = 1 : matsGrid.params.stepsItThr[3] : matsGrid.params.divLst[3] + 1
					
				# end
			# end
		# end
	# end
	# for pos in matsGrid.params.posLst
		# isSurface = true;
		# for iDim = 1 : params.nDim
			# isSurface = ( isSurface 
				# && ( pos[iDim] > 1 )
				# && ( pos[iDim] <= matsGrid.params.divLst[iDim] ) );
			# if !isSurface
				# break;
			# end
		# end
		# if !isSurface
			# continue;
		# end
		# eigen
	# end
	# for iCoDim = nDim : -1 : 1
		# dimBndLst = ones(Int64, nDim - iCoDim);
		# for ii = 1 : nDim - iCoDim
			# dimBndLst[ii] = ii;
		# end
		# while
			# broadcastAssign!( matsGrid.params.stepsItThr, 1 );
			# for ii = 1 : nDim
				# matsGrid.params.divsItThr[ii] = matsGrid.params.divLst[ii]+1;
			# end
			# iDim = 1;
			# broadcastAssign!( matsGrid.params.locItThr, 2 );
			# for ii = 1 : nDim - iCoDim
				# matsGrid.params.divsItThr[dimBndLst[ii]] = 2;
				# matsGrid.params.stepsItThr[dimBndLst[ii]] = matsGrid.params.divLst[dimBndLst[ii]];
				# matsGrid.params.locItThr[dimBndLst[ii]] = 1;
			# end
			
			# iDim = 1;
			# broadcastAssign( matsGrid.params.locItThr, 2 );
			
			# while
				
			# end
			
			# iDim = 1;
			# dimBndLst[1] += 1;
			# while dimBndLst[iDim] > nDim - iDim +1 && iDim < nDim-iCoDim
				# iDim += 1;
				# dimBndLst[iDim] += 1;
			# end
			# if iDim == nDim - iCoDim
				# break;
			# else
				# while iDim < nDim - iCoDim
					# dimBndLst[iDim+1] = dimBndLst[iDim] + 1;
					# iDim += 1;
				# end
			# end
		# end
		# for 
			# ;
		# end
	# end
# end
