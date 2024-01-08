function syevr!(jobz::AbstractChar, laRange::AbstractChar, uplo::AbstractChar, A::AbstractMatrix{ComplexF64},
                        vl::AbstractFloat, vu::AbstractFloat, il::Integer, iu::Integer, abstol::AbstractFloat)
		chkstride1(A);
		n = checksquare(A);
		if laRange == 'I' && !(1 <= il <= iu <= n)
			throw(ArgumentError("illegal choice of eigenvalue indices (il = $il, iu=$iu), which must be between 1 and n = $n"))
		end
		if laRange == 'V' && vl >= vu
			throw(ArgumentError("lower boundary, $vl, must be less than upper boundary, $vu"))
		end
		lda = max(1,stride(A,2))
		m = Ref{BlasInt}()
		w = similar(A, Float64, n)
		if jobz == 'N'
			ldz = 1
			Z = similar(A, ComplexF64, ldz, 0)
		elseif jobz == 'V'
			ldz = n
			Z = similar(A, ComplexF64, ldz, n)
		end
		isuppz = similar(A, BlasInt, 2*n);
		work   = Vector{ComplexF64}(undef, 1);
		lwork  = BlasInt(-1);
		rwork  = Vector{Float64}(undef, 1);
		lrwork = BlasInt(-1);
		iwork  = Vector{BlasInt}(undef, 1);
		liwork = BlasInt(-1);
		info   = Ref{BlasInt}();
		for i = 1:2  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
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
			chklapackerror(info[])
			if i == 1
				lwork = BlasInt(real(work[1]))
				resize!(work, lwork)
				lrwork = BlasInt(rwork[1])
				resize!(rwork, lrwork)
				liwork = iwork[1]
				resize!(iwork, liwork)
			end
		end
		w[1:m[]], Z[:,1:(jobz == 'V' ? m[] : 0)]
	end
