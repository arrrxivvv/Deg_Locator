const liblapack = Base.liblapack_name

import ..LinearAlgebra.BLAS.@blasfunc

import ..LinearAlgebra: BlasFloat, BlasInt, LAPACKException,
    DimensionMismatch, SingularException, PosDefException, chkstride1, checksquare

using ..LinearAlgebra: triu, tril, dot

using Base: iszero, require_one_based_indexing

function callZheevr!( A::AbstractMatrix{ComplexF64}, w::AbstractMatrix{ComplexF64}, Z::AbstractMatrix{ComplexF64} )
	@time begin
	
	laRange = 'A';
	jobz = 'V';
	vl, vu, il, iu = 0, 0, 0, 0;
	abstol = -1;
	uplo = 'U';
	
	if laRange == 'I' && !(1 <= il <= iu <= n)
		throw(ArgumentError("illegal choice of eigenvalue indices (il = $il, iu=$iu), which must be between 1 and n = $n"));
	end
	if laRange == 'V' && vl >= vu
		throw(ArgumentError("lower boundary, $vl, must be less than upper boundary, $vu"));
	end
	
	n = checksquare(A);
	lda = max(1,stride(A,2));
	ldz = n;
	m = Ref{BlasInt}();
	# w = similar(A, Float64, n);
	# if jobz == 'N'
		# ldz = 1;
		# Z = similar(A, ComplexF64, ldz, 0);
	# elseif jobz == 'V'
		# ldz = n;
		# Z = similar(A, ComplexF64, ldz, n);
	# end
	isuppz = similar(A, BlasInt, 2*n);
	work   = Vector{ComplexF64}(undef, 1);
	lwork  = BlasInt(-1);
	rwork  = Vector{Float64}(undef, 1);
	lrwork = BlasInt(-1);
	iwork  = Vector{BlasInt}(undef, 1);
	liwork = BlasInt(-1);
	info   = Ref{BlasInt}();
	end
# end

	@debug("n: ($n)", "lda: ($lda)");

# function callZheevr()
	# @info(lwork);
	for i = 1:2  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
		@info("ccall: ")
		@time begin
		ccall((@blasfunc(zheevr_), liblapack), Cvoid,
			  (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
			   Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64}, Ref{ComplexF64},
			   Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{BlasInt},
			   Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
			   Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
			   Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
			   Clong, Clong, Clong),
			  jobz, laRange, uplo, n,
			  A, lda, vl, vu,
			  il, iu, abstol, m,
			  w, Z, ldz, isuppz,
			  work, lwork, rwork, lrwork,
			  iwork, liwork, info,
			  1, 1, 1)
		end
		# chklapackerror(info[])
		if i == 1
			# @info("resize")
			@info("lwork")
			@time begin
			lwork = BlasInt(real(work[1]));
			resize!(work, lwork);
			lrwork = BlasInt(rwork[1])
			resize!(rwork, lrwork);
			liwork = iwork[1];
			resize!(iwork, liwork);
			end
		end
	end
	@debug("m[]: ($m[])")
	# return w,Z;
end
