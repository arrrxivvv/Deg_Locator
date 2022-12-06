module EigCustom
	# using Infiltrator
	
	const liblapack = Base.liblapack_name

	import LinearAlgebra.BLAS.@blasfunc

	import LinearAlgebra: BlasFloat, BlasInt, LAPACKException,
		DimensionMismatch, SingularException, PosDefException, chkstride1, checksquare

	using LinearAlgebra: triu, tril, dot
	
	using LinearAlgebra

	using Base: iszero, require_one_based_indexing
	
	export eigenZheevr!, eigenZheevrStruct!, eigenPrework!, eigenPreworkStruct!, eigWorkStructFromNum!, eigenWorkThrdInit!, testZheevTime, testZheevTime2, testZheevInside;
	
	include("EigCustom_preworks.jl")
	export EigWork, eigenPreworkStruct!, eigWorkStructFromNum!
	
	include("EigCustom_tmpMats.jl")
	export EigTmpMats, eigTmpMatsInit, eigOnTmpMats!
	
	function eigenZheevrStruct!( A::Array{ComplexF64,2}, w::Vector{Float64}, Z::Array{ComplexF64,2}, eigWork::EigWork; jobz = 'V' )
		eigenZheevr!( A, w, Z, eigWork.lwork,eigWork.lrwork, eigWork.liwork; work = eigWork.workLst, rwork = eigWork.rworkLst, iwork = eigWork.iworkLst, isuppz = eigWork.isuppz, jobz = jobz )
		# , m = eigWork.m, info = eigWork.info );
	end
	
	function eigenZheevrStruct!( A::Array{Float64,2}, w::Vector{Float64}, Z::Array{Float64,2}, eigWork::EigWork; jobz = 'V' )
		eigenDsyver!( A, w, Z, eigWork.lwork,eigWork.liwork; work = eigWork.workLst,iwork = eigWork.iworkLst, isuppz = eigWork.isuppz, jobz = jobz )
		# , m = eigWork.m, info = eigWork.info );
	end
	
	function eigenZheevr!( A::Array{ComplexF64,2}, w::Vector{Float64}, Z::Array{ComplexF64,2}, lwork = BlasInt(-1), lrwork = BlasInt(-1), liwork = BlasInt(-1); work= Vector{ComplexF64}(undef, 1), rwork = Vector{Float64}(undef, 1), iwork = Vector{BlasInt}(undef, 1), jobz = 'V', isuppz = zeros(BlasInt,0) )
	# , m = Ref{BlasInt}(), info = Ref{BlasInt}() )
		# @time begin
		laRange = 'A';
		vl, vu, il, iu = 0, 0, 0, 0;
		abstol = -1;
		uplo = 'U';
		
		n = checksquare(A);
		lda = max(1,stride(A,2));
		ldz = n;
		m = Ref{BlasInt}();
		if isempty(isuppz)
			isuppz = similar(A, BlasInt, 2*n);
		end
		info = Ref{BlasInt}();
	
		for i = 1:2  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
			if lwork != -1 && i == 1
				continue;
			end
			# @debug("ccall: ")
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
			if i == 1
				lwork = BlasInt(real(work[1]));
				resize!(work, lwork);
				lrwork = BlasInt(rwork[1])
				resize!(rwork, lrwork);
				liwork = iwork[1];
				resize!(iwork, liwork);
			end
		end
		# end
	end
	
	function eigenDsyver!( A::Array{Float64,2}, w::Vector{Float64}, Z::Array{Float64,2}, lwork = BlasInt(-1), liwork = BlasInt(-1); work= Vector{Float64}(undef, 1), iwork = Vector{BlasInt}(undef, 1), jobz = 'V', isuppz = zeros(BlasInt,0) )
		laRange = 'A';
		vl, vu, il, iu = 0, 0, 0, 0;
		abstol = -1;
		uplo = 'U';
		
		n = checksquare(A);
		lda = max(1,stride(A,2));
		ldz = n;
		m = Ref{BlasInt}();
		if isempty(isuppz)
			isuppz = similar(A, BlasInt, 2*n);
		end
		info = Ref{BlasInt}();
	
		for i = 1:2  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
			if lwork != -1 && i == 1
				continue;
			end
			ccall((@blasfunc(dsyevr_), liblapack), Cvoid,
				(Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
					Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
					Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{BlasInt},
					Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
					Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
					Ptr{BlasInt}, Clong, Clong, Clong),
				jobz, laRange, uplo, n,
				A, max(1,lda), vl, vu,
				il, iu, abstol, m,
				w, Z, max(1,ldz), isuppz,
				work, lwork, iwork, liwork,
				info, 1, 1, 1)
			if i == 1
				lwork = BlasInt(real(work[1]));
				resize!(work, lwork);
				liwork = iwork[1];
				resize!(iwork, liwork);
			end
		end
	end
	
	function eigenPrework!( A::Array{ComplexF64,2}, w::Vector{Float64}, Z::Array{ComplexF64,2} )
		laRange = 'A';
		jobz = 'V';
		vl, vu, il, iu = 0, 0, 0, 0;
		abstol = -1;
		uplo = 'U';
		
		n = checksquare(A);
		lda = max(1,stride(A,2));
		ldz = n;
		m = Ref{BlasInt}();
		# isuppz = similar(A, BlasInt, 2*n);
		isuppz = similar(A, BlasInt, 2*n);
		work   = Vector{ComplexF64}(undef, 1);
		lwork  = BlasInt(-1);
		rwork  = Vector{Float64}(undef, 1);
		lrwork = BlasInt(-1);
		iwork  = Vector{BlasInt}(undef, 1);
		liwork = BlasInt(-1);
		info   = Ref{BlasInt}();
		
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
		
		lwork = BlasInt(real(work[1]));
		lrwork = BlasInt(rwork[1]);
		liwork = iwork[1];
		
		return lwork, lrwork, liwork;
	end
	
	function eigenPrework!( A::Array{Float64,2}, w::Vector{Float64}, Z::Array{Float64,2} )
		laRange = 'A';
		jobz = 'V';
		vl, vu, il, iu = 0, 0, 0, 0;
		abstol = -1;
		uplo = 'U';
		
		n = checksquare(A);
		lda = max(1,stride(A,2));
		ldz = n;
		m = Ref{BlasInt}();
		isuppz = similar(A, BlasInt, 2*n);
		work   = Vector{Float64}(undef, 1);
		lwork  = BlasInt(-1);
		iwork  = Vector{BlasInt}(undef, 1);
		liwork = BlasInt(-1);
		info   = Ref{BlasInt}();
		
		ccall((@blasfunc(dsyevr_), liblapack), Cvoid,
			(Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
				Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
				Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{BlasInt},
				Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
				Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
				Ptr{BlasInt}, Clong, Clong, Clong),
			jobz, laRange, uplo, n,
			A, max(1,lda), vl, vu,
			il, iu, abstol, m,
			w, Z, max(1,ldz), isuppz,
			work, lwork, iwork, liwork,
			info, 1, 1, 1)
		
		lwork = BlasInt(real(work[1]));
		lrwork = 0;
		liwork = iwork[1];
		
		return lwork, lrwork, liwork;
	end
	
	
	function eigenWorkThrdInit!( HmatLst, Elst, vecLst )
		# @infiltrate
		nThreads = Threads.nthreads();
		lwork, lrwork, liwork = eigenPrework!( HmatLst, Elst, vecLst );
		workLst = Array{ComplexF64,2}(undef, Int64(lwork),nThreads);
		rworkLst = Array{Float64,2}(undef, lrwork,nThreads);
		iworkLst = Array{BlasInt,2}(undef, liwork,nThreads);
		
		return lwork, lrwork, liwork, workLst, rworkLst, iworkLst;
	end
	
	function testZheevTime( N, itLnNum; isPreWork = false )
		mat = rand(ComplexF64, N, N);
		mat = mat + mat';
		
		matTmpArr = [ similar(mat) for it = 1 : itLnNum ];
		ZLst = [ similar(mat) for it = 1 : itLnNum ];
		wLst = [ similar(mat,Float64,N) for it = 1 : itLnNum ];
		
		@time begin
		if isPreWork
			matTmpArr[1] .= mat;
			lwork, lrwork, liwork = eigenPrework!( matTmpArr[1], wLst[1], ZLst[1] );
			@debug("($lwork), ($lrwork), ($liwork)");
			workLst = [ Vector{ComplexF64}(undef, lwork) for ii = 1 : itLnNum ];
			rworkLst = [ Vector{Float64}(undef, lrwork) for ii = 1 : itLnNum ];
			iworkLst = [ Vector{BlasInt}(undef, liwork) for ii = 1 : itLnNum ];
		end
		end
		
		@time begin
			Threads.@threads for ii in 1:itLnNum
				matTmpArr[ii] .= mat;
				if isPreWork
					eigenZheevr!( matTmpArr[ii], wLst[ii], ZLst[ii], lwork, lrwork, liwork; work=workLst[ii], rwork=rworkLst[ii], iwork=iworkLst[ii] );
				else
					eigenZheevr!( matTmpArr[ii], wLst[ii], ZLst[ii] );
				end
			end
		end
		
		# @time eigenZheevr!( mat, wN, ZN );
	end
	
	function testZheevTime2( N, itLnNum; isPreWork = false, isEigen = false )
		mat = rand(ComplexF64, N, N);
		mat = mat + mat';
		
		matTmpArr = [ similar(mat) for it = 1 : itLnNum ];
		ZLst = [ similar(mat) for it = 1 : itLnNum ];
		wLst = [ similar(mat,Float64,N) for it = 1 : itLnNum ];
		
		@time begin
		if isPreWork
			matTmpArr[1] .= mat;
			lwork, lrwork, liwork = eigenPrework!( matTmpArr[1], wLst[1], ZLst[1] );
			@debug("($lwork), ($lrwork), ($liwork)");
			workLst = Array{ComplexF64,2}(undef, lwork,itLnNum);
			rworkLst = Array{Float64,2}(undef, lrwork,itLnNum);
			iworkLst = Array{BlasInt,2}(undef, liwork,itLnNum);
		end
		end
		
		@time begin
			Threads.@threads for ii in 1:itLnNum
				matTmpArr[ii] .= mat;
				if isPreWork
					eigenZheevr!( matTmpArr[ii], wLst[ii], ZLst[ii], lwork, lrwork, liwork; work=view(workLst,:,ii), rwork=view(rworkLst,:,ii), iwork=view(iworkLst,:,ii) );
				elseif isEigen
					wLst[ii], ZLst[ii] = eigen( matTmpArr[ii] );
				else
					eigenZheevr!( matTmpArr[ii], wLst[ii], ZLst[ii] );
				end
			end
		end
		
		# @time eigenZheevr!( mat, wN, ZN );
	end
	
	function testZheevInside( N, itLnNum; isPreWork = false, isEigen = false )
		mat = rand(ComplexF64, N, N);
		mat = mat + mat';
		
		matTmpArr = [ similar(mat) for it = 1 : itLnNum ];
		ZLst = [ similar(mat) for it = 1 : itLnNum ];
		wLst = [ similar(mat,Float64,N) for it = 1 : itLnNum ];
		
		matTmpArr[1] .= mat;
		A = matTmpArr[1];
		
		laRange = 'A';
		jobz = 'V';
		vl, vu, il, iu = 0, 0, 0, 0;
		abstol = -1;
		uplo = 'U';
		
		n = checksquare(A);
		lda = max(1,stride(A,2));
		ldz = n;
		m = Ref{BlasInt}();
		isuppz = similar(A, BlasInt, 2*n);
		info = Ref{BlasInt}();
		
		@time begin
		if isPreWork
			matTmpArr[1] .= mat;
			lwork, lrwork, liwork = eigenPrework!( matTmpArr[1], wLst[1], ZLst[1] );
			@debug("($lwork), ($lrwork), ($liwork)");
			workLst = Array{ComplexF64,2}(undef, lwork,itLnNum);
			rworkLst = Array{Float64,2}(undef, lrwork,itLnNum);
			iworkLst = Array{BlasInt,2}(undef, liwork,itLnNum);
		end
		end
		
		@info( "ccalls" );
		@time begin
		Threads.@threads for ii = 1:itLnNum  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
			@debug("ii: ($ii)");
			if lwork != -1 && ii == 1
				continue;
			end
			matTmpArr[ii] .= mat;
			if !isPreWork
				break
			end
			@debug("ccall: ")
			# @time begin
			ccall((@blasfunc(zheevr_), liblapack), Cvoid,
				  (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
				   Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64}, Ref{ComplexF64},
				   Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{BlasInt},
				   Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
				   Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
				   Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
				   Clong, Clong, Clong),
				  jobz, laRange, uplo, n,
				  matTmpArr[ii], lda, vl, vu,
				  il, iu, abstol, m,
				  wLst[ii], ZLst[ii], ldz, isuppz,
				  view(workLst,:,ii), lwork, view(rworkLst,:,ii), lrwork,
				  view(iworkLst,:,ii), liwork, info,
				  1, 1, 1)
			# end
		end
		end
	end
	
	function eigenZheevrRep!( itNum, A::Array{ComplexF64,2}, w::Vector{Float64}, Z::Array{ComplexF64,2}, lwork = BlasInt(-1), lrwork = BlasInt(-1), liwork = BlasInt(-1); work= Vector{ComplexF64}(undef, 1), rwork = Vector{Float64}(undef, 1), iwork = Vector{BlasInt}(undef, 1) )
		# @time begin
		
		Aorg = copy(A);
		
		laRange = 'A';
		jobz = 'V';
		vl, vu, il, iu = 0, 0, 0, 0;
		abstol = -1;
		uplo = 'U';
		
		n = checksquare(A);
		lda = max(1,stride(A,2));
		ldz = n;
		m = Ref{BlasInt}();
		isuppz = similar(A, BlasInt, 2*n);
		info = Ref{BlasInt}();
		
		# end
	# end

		@debug("n: ($n)", "lda: ($lda)");
		@debug( "lrwork: ($lrwork), rwork size: ($(size(rwork)))" );
	
		# @debug(lwork);
		for i = 1:itNum  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
			if lwork != -1 && i == 1
				continue;
			end
			A .= Aorg;
			@debug("ccall: ")
			# @time begin
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
			# end
			
			if i == 1
				@debug("resize")
				@debug("info: ($info)")
				# @time begin
				lwork = BlasInt(real(work[1]));
				@debug("lwork: ($lwork)")
				@debug("work size: ($(size(work)))")
				resize!(work, lwork);
				@debug("work size: ($(size(work)))")
				lrwork = BlasInt(rwork[1])
				resize!(rwork, lrwork);
				@debug("lrwork: ($lrwork)")
				liwork = iwork[1];
				@debug("liwork: ($liwork)")
				resize!(iwork, liwork);
				# end
			end
		end
	end
	
end
