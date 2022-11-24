struct EigWork
	lwork::BlasInt;
	lrwork::BlasInt;
	liwork::BlasInt;
	
	workLst::Vector{ComplexF64};
	rworkLst::Vector{Float64};
	iworkLst::Vector{BlasInt};
end

function eigenPreworkStruct!( A::Array{ComplexF64,2}, w::Vector{Float64}, Z::Array{ComplexF64,2} )
	lwork, lrwork, liwork = eigenPrework!( A, w, Z );
	
	work = Vector{ComplexF64}(undef,Int64(lwork));
	rwork = Vector{Float64}(undef,lrwork);
	iwork = Vector{BlasInt}(undef,liwork);
	
	return EigWork( lwork, lrwork, liwork, work, rwork, iwork );
end

function eigWorkStructFromNum!( mSz::Int64 )
	A = zeros(ComplexF64, mSz, mSz);
	w = zeros(mSz);
	Z = zeros(ComplexF64, mSz, mSz);
	
	lwork, lrwork, liwork = eigenPrework!( A, w, Z );
	
	work = Vector{ComplexF64}(undef,Int64(lwork));
	rwork = Vector{Float64}(undef,lrwork);
	iwork = Vector{BlasInt}(undef,liwork);
	
	return EigWork( lwork, lrwork, liwork, work, rwork, iwork );
end
