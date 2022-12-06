struct EigWork{T<:Number}
	lwork::BlasInt;
	lrwork::BlasInt;
	liwork::BlasInt;
	
	workLst::Vector{T};
	rworkLst::Vector{Float64};
	iworkLst::Vector{BlasInt};
	isuppz::Vector{BlasInt};
	
	# m::Ref{BlasInt};
	# info::Ref{BlasInt};
end

function eigenPreworkStruct!( A::Array{T,2}, w::Vector{Float64}, Z::Array{T,2} ) where {T<:Number}
	lwork, lrwork, liwork = eigenPrework!( A, w, Z );
	
	work = Vector{T}(undef,Int64(lwork));
	rwork = Vector{Float64}(undef,lrwork);
	iwork = Vector{BlasInt}(undef,liwork);
	
	n = checksquare(A);
	isuppz = similar(A, BlasInt, 2*n);
	
	# m = Ref{BlasInt}(); 
	# info = Ref{BlasInt}();
	
	return EigWork( lwork, lrwork, liwork, work, rwork, iwork, isuppz );
	# , m, info );
end

function eigWorkStructFromNum!( mSz::Int64; typeElm = ComplexF64 )
	A = zeros(typeElm, mSz, mSz);
	w = zeros(mSz);
	Z = zeros(typeElm, mSz, mSz);
	
	return eigenPreworkStruct!( A, w, Z );
end
