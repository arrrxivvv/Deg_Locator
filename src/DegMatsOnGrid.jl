@enum HgridStoreOpt HgridFull=1 Hlayered HthreadedOnly

struct DegMatsOnGrid
	params::DegParams
	Elst::Array{Vector{Float64}};
	vLst::Array{Matrix{ComplexF64}};
	
	HmatLst; # ::Array{Matrix{ComplexF64}};
	HgridOpt::HgridStoreOpt;
	eigWorkTh::ThrStruct{EigCustom.EigWork};
end

function matsGridBase( params::DegParams, HmatLst, HOpt::HgridStoreOpt )
	Elst = makeArrOverGrid( Vector{Float64}, params );
	vLst = makeArrOverGrid( Matrix{ComplexF64}, params );
	
	# if Hopt == HgridFull
		# HmatLst = makeArrOverGrid( Matrix{ComplexF64}, params );
	# elseif Hopt == HthreadedOnly
		# HmatLst = threaded_zeros( ComplexF64, params.N, params.N );
	# end
	
	eigWorkTh = thrStructCopy( eigWorkStructFromNum!( params.N ) );
	
	return DegMatsOnGrid( params, Elst, vLst, HmatLst, HOpt, eigWorkTh );
end

function matsGridHinternal( params::DegParams, Hopt::HgridStoreOpt )
	if Hopt == HgridFull
		HmatLst = makeArrOverGrid( Matrix{ComplexF64}, params );
	elseif Hopt == Hlayered
		HmatLst = makeArrOverGrid( Matrix{	ComplexF64}, params );
	# elseif Hopt == HthreadedOnly
		# HmatLst = threaded_zeros( ComplexF64, params.N, params.N );
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

function eigenAtLoc( loc::Vector{Int64}, matsGrid::DegMatsOnGrid )
	idLin = linIdFromIdVec( loc, matsGrid.params );
	Hmat = getHLoc( loc, matsGrid );
	if !isassigned( matsGrid.Elst, idLin )
		matsGrid.Elst[idLin] = zeros(matsGrid.params.N);
	end
	if !isassigned( matsGrid.vLst, idLin )
		matsGrid.vLst[idLin] = zeros( ComplexF64, matsGrid.params.N, matsGrid.params.N );
	end
	
	eigenZheevrStruct!( Hmat, matsGrid.Elst[idLin], matsGrid.vLst[idLin], getThrInst(matsGrid.eigWorkTh) );
end
