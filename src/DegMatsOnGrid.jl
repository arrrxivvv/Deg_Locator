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
		eigenIdGrid = zeros( Int64, params.divLst... );
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

function getHLoc( loc::Vector{Int64}, matsGrid::DegMatsOnGrid )
	if matsGrid.HgridOpt == HgridFull
		idLin = linIdFromIdVec(loc);
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

function checkEigenDone( matsGrid::DegMatsOnGrid, loc::Vector{Int64} )
	idLin = linIdFromIdVec( loc, matsGrid.params );
	if matsGrid.eigenIdGrid[idLin] ==  matsGrid.eigenId[]
		return true;
	else
		matsGrid.eigenIdGrid[idLin] = matsGrid.eigenId[];
		return false;
	end
end

function eigenAtLoc( loc::Vector{Int64}, matsGrid::DegMatsOnGrid; noReEigen = true )
	idLin = linIdFromIdVec( loc, matsGrid.params );
	if checkEigenDone( matsGrid, loc )
		return;
	end
	
	Hmat = getHLoc( loc, matsGrid );
	if !isassigned( matsGrid.Elst, idLin )
		matsGrid.Elst[idLin] = zeros(matsGrid.params.N);
	end
	if !isassigned( matsGrid.vLst, idLin )
		matsGrid.vLst[idLin] = zeros( ComplexF64, matsGrid.params.N, matsGrid.params.N );
	end
	
	eigenZheevrStruct!( Hmat, matsGrid.Elst[idLin], matsGrid.vLst[idLin], getThrInst(matsGrid.eigWorkTh) );
end

function setCurrentDone( matsGrid::DegMatsOnGrid )
	linId = linItCurrent( matsGrid.params );
	matsGrid.eigenIdGrid[linId] = matsGrid.eigenId;
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

function setCurrentVLst( matsGrid::DegMatsOnGrid, mat )
	linId = linIdFromIdVec( 
		getThrInst( matsGrid.params.locItThr ), matsGrid.params );
	if !isassigned( matsGrid.vLst, linId )
		matsGrid.vLst[linId] = mat;
	else
		matsGrid.vLst[linId] .= mat;
	end
end

function matsGridTransfer!( matsGsT, matsGrS, locS )
	for iSh = 1 : 2^matsGrT.params.nDim
		shIdVec!( matsGrS.params, locS, iSh );
		locItCorner!( matsGrT.params, iSh );
		setCurrentElst( matsGrT, 
			getCurrentElst( matsGrS ) );
		setCurrentVLst( matsGrT, 
			getCurrentVLst( matsGrS ) );
		setCurrentDone( matsGrT );
	end
end

function eigenOnSurface( matsGrid::DegMatsOnGrid )
	for iCoDim = nDim : -1 : 1
		dimBndLst = ones(Int64, nDim - iCoDim);
		for ii = 1 : nDim - iCoDim
			dimBndLst[ii] = ii;
		end
		while
			broadcastAssign!( matsGrid.params.stepsItThr, 1 );
			for ii = 1 : nDim
				matsGrid.params.divsItThr[ii] = matsGrid.params.divLst[ii]+1;
			end
			iDim = 1;
			broadcastAssign( matsGrid.params.locItThr, 2 );
			for ii = 1 : nDim - iCoDim
				matsGrid.params.divsItThr[dimBndLst[ii]] = 2;
				matsGrid.params.stepsItThr[dimBndLst[ii]] = matsGrid.params.divLst[dimBndLst[ii]];
				matsGrid.params.locItThr[dimBndLst[ii]] = 1;
			end
			
			iDim = 1;
			broadcastAssign( matsGrid.params.locItThr, 2 );
			
			while
				
			end
			
			iDim = 1;
			dimBndLst[1] += 1;
			while dimBndLst[iDim] > nDim - iDim +1 && iDim < nDim-iCoDim
				iDim += 1;
				dimBndLst[iDim] += 1;
			end
			if iDim == nDim - iCoDim
				break;
			else
				while iDim < nDim - iCoDim
					dimBndLst[iDim+1] = dimBndLst[iDim] + 1;
					iDim += 1;
				end
			end
		end
		for 
			;
		end
	end
end
