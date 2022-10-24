using LinearAlgebra
using Infiltrator
using EllipsisNotation
using Random
using Logging
using LoggingExtras
using CUDA
using EigCustom
import LinearAlgebra: BlasFloat, BlasInt, LAPACKException,
	DimensionMismatch, SingularException, PosDefException, chkstride1, checksquare

function locToParam( loc, param_grids, param_dim )
	paramArr = zeros( param_dim );
	for iDim = 1 : param_dim
		paramArr[iDim] = param_grids[iDim][loc[iDim]]
	end
	
	return paramArr;
end

function locToParamArr( locArr, param_grids, param_dim )
	paramArr = zeros( size(locArr) );
	for it = 1 : size(locArr)[1]
		paramArr[it,:] = locToParam( locArr[it,:], param_grids, param_dim );
	end
	return paramArr;
end

function HlstFunc( Hfunc, param_dim, N, isOrtho = false )
	Hlst = fill( Array{Complex{Float64}}(undef,(N,N)), (2,param_dim) );
	# Hlst = [ zeros( Complex{Float64}, N,N ) for ii = 1:2, jj = 1:param_dim ];
	for it in eachindex(Hlst)
		Hlst[it] = Hfunc(N);
		if isOrtho
			for jj = 1 : it-1
				Hlst[it] -= Hlst[jj] * tr(Hlst[it] * Hlst[jj]) / tr(Hlst[jj] * Hlst[jj]);
			end
			Hlst[it] += Hlst[it]';
			Hlst[it] /= 2;
		end
	end
	return Hlst;
end

function paramToLoc( param, param_dim, param_min, param_step )
	loc = floor.( Int64, ( param .- param_min ) ./ param_step ) .+ 1;
	
	return loc;
end

function paramToLocArr( paramArr, param_dim, param_min, param_step )
	locArr = ones( Int64, size(paramArr) );
	for it = 1 : size( paramArr, 1 )
		locArr[it, :] = paramToLoc( paramArr[it,:], param_dim, param_min, param_step );
	end
	return locArr;
end

function locateNon0DivOld( divBlst, threshold, N, param_dim )
	divBlstReal = (p-> real.(p)).(divBlst);
	divBposInd = broadcast.( >, divBlstReal, threshold );
	divBnegInd = broadcast.( <, divBlstReal, -threshold );
	divBposArr = arrLstToArr( divBposInd );
	divBnegArr = arrLstToArr( divBnegInd );
	posLocLst = Vector{ Array{ Int64,2 } }(undef,N);
	negLocLst = Vector{ Array{ Int64,2 } }(undef,N);
	for n = 1:N
		posLocLst[n] = cartIndLstToArr( findall( divBposArr[..,n] ), param_dim );
		negLocLst[n] = cartIndLstToArr( findall( divBnegArr[..,n] ), param_dim );
	end	
	posCount = sum( divBposInd );
	negCount = sum( divBnegInd );
	
	return posCount, negCount, posLocLst, negLocLst;
end

function locateNon0Div( divBlst, threshold, N, param_dim; clipBnd = 0, isGOE3d=false, Bdim = 1 )
	divBlstReal = (p-> real.(p)).(divBlst);
	param_divide = size(divBlst, 1);
	if clipBnd > 0
		sliceLst = []
		for dim = 1 : param_dim
			ln = size( divBlst, dim );
			append!(sliceLst, [1:ln-clipBnd]);
		end
		divBlstReal = divBlstReal[sliceLst...];
	end
	divBposInd = broadcast.( >, divBlstReal, threshold );
	divBnegInd = broadcast.( <, divBlstReal, -threshold );
	divBposArr = arrLstToArr( divBposInd );
	divBnegArr = arrLstToArr( divBnegInd );
	if !isGOE3d 
		posCount = sum( divBposInd );
		negCount = sum( divBnegInd );
		posLocLst = Vector{ Array{ Int64,2 } }(undef,N);
		negLocLst = Vector{ Array{ Int64,2 } }(undef,N);
		for n = 1:N
			posLocLst[n] = cartIndLstToArr( findall( divBposArr[..,n] ), param_dim );
			negLocLst[n] = cartIndLstToArr( findall( divBnegArr[..,n] ), param_dim );
		end	
	else 
		dimLst = [1 : param_dim;];
		deleteat!(dimLst, Bdim);
		# @infiltrate
		posCount = dropdims( sum( divBposInd; dims=dimLst ); dims=Tuple(dimLst) );
		negCount = dropdims( sum( divBnegInd; dims=dimLst ); dims=Tuple(dimLst) );
		# @infiltrate
		posLocLst = Array{ Array{ Int64,2 }, 2 }(undef,param_divide,N);
		negLocLst = Array{ Array{ Int64,2 }, 2 }(undef,param_divide,N);
		locLstPol = [posLocLst, negLocLst];
		divBindPol = [divBposArr, divBnegArr];
		for n = 1:N, id = 1:param_divide, ip = 1 : 2
			# @infiltrate
			divLst = fill(param_divide, param_dim);
			divLst[Bdim] = 1;
			locLstPol[ip][id,n] = cartIndLstToArr( findall( reshape( selectdim( (@view divBindPol[ip][..,n]), Bdim, id ), Tuple(divLst) )), param_dim );
			locLstPol[ip][id,n][:,Bdim] .= id;
			# posLocLst[id,n] = cartIndLstToArr( findall( reshape( selectdim( (@view divBposArr[..,n]), Bdim, id ), Tuple(divLst) )), param_dim );
			# negLocLst[id,n] = cartIndLstToArr( findall( reshape( selectdim( (@view divBnegArr[..,n]), Bdim, id ), Tuple(divLst) ) ), param_dim );
		end	
	end
	
	return posCount, negCount, posLocLst, negLocLst;
end

function locLstPurify_old( posLocLst, negLocLst, posNlst )
	itNum = size(posLocLst, 1);
	N = size(posLocLst, 2);
	nLst = ones(Int64, itNum, N);
	pureNlst = posNlst[:,2:end] - posNlst[:,1:end-1];
	
	posLocLstPure = [ [ zeros( Int64, pureNlst[it,n], 3 ) for n = 1:N-1 ] for it = 1:itNum ];
	negLocLstPure = [ [ zeros( Int64, pureNlst[it,n], 3 ) for n = 1:N-1 ] for it = 1:itNum ];
	
	locLstPols = [posLocLst, negLocLst];
	locLstPolsPure = [posLocLstPure, negLocLstPure];
	polInvLst = [1,2];
	
	Threads.@threads for it = 1 : itNum
		for n = 2 : N
			for pol = 1:2
				if n == 1
					locLstPolsPure[pol][it][1,:] = locLstPols[pol][it][1,:];
					continue ;
				end
				polInv = polInvLst[pol];
				
				lnThis = posNlst[itNum,n];
				lnPrev = pureNlst[itNum,n-1];
				idPrev = 1;
				idPure = 1;
				for id = 1 : lnThis
					locTested = locLstPos[pol][it][n][id,:];
					if idPrev > lnPrev
						locLstPolsPure[pol][it][n][idPure,:] = locTested;
						idPure+=1;
					elseif locTested != locLstPolsPure[polInv][it][n][idPrev,:]
						locLstPolsPure[pol][it][n][idPure,:] = locTested;
						idPure+=1;
					else
						idPrev+=1;
					end
				end
			end
		end
	end
	
	return locLstPolsPure[1], locLstPolsPure[2], pureNlst;
end

function locLstDistill( N, param_dim, posLocLst, negLocLst )
	posLocLstPure = Vector{ Array{ Array{Int64} } }(undef,N);
	negLocLstPure = Vector{ Array{ Array{Int64} } }(undef,N);
	locLstPols = broadcast( (x->arrToArrLst.(x)), [posLocLst, negLocLst] );
	locLstPolsPure = [posLocLstPure, negLocLstPure];
	idPolsOppLst = [2,1];
	for ip = 1 : 2, n in [1,N]
		locLstPolsPure[ip][n] = locLstPols[ip][n];
	end
	for n = 2 : N-1
		for idPols = 1 : 2
			idPolsOpp = idPolsOppLst[idPols];
			locLstPolsPure[idPols][n] = filter( !(x->x in locLstPolsPure[idPolsOpp][n-1]), locLstPols[idPols][n] );
			if length( locLstPolsPure[idPols][n] ) != length( locLstPols[idPols][n] ) - length(locLstPolsPure[idPolsOpp][n-1])
				@warn( "imprecise exclusion: ", string( length( locLstPols[idPols][n] ) - length(locLstPolsPure[idPolsOpp][n-1]) - length( locLstPolsPure[idPols][n] ) ) );
			end
		end
	end
	locLstPolsPure = broadcast( (x->arrLstToArr.(x)), locLstPolsPure );
	
	return locLstPolsPure[1], locLstPolsPure[2];
end

function locator_Eigen_threaded!( degObj, HmatFun; isGC = false )
	isGC = degObj.isStrongGC;
	@info( "Create Hmats" )
	# @infiltrate
	tFull = @timed begin
		Threads.@threads for ind in degObj.posLst
			HmatFun( degObj.HmatLst[ind], degObj.param_mesh[ind] );
		end
	end
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	
	nThreads = Threads.nthreads();
	lwork, lrwork, liwork = eigenPrework!( degObj.HmatLst[1], degObj.Elst[1], degObj.vecLst[1] );
	workLst = Array{ComplexF64,2}(undef, Int64(lwork),nThreads);
	rworkLst = Array{Float64,2}(undef, lrwork,nThreads);
	iworkLst = Array{BlasInt,2}(undef, liwork,nThreads);
	
	@info( "Eigen: " )
	tFull = @timed begin
	Threads.@threads for ii in degObj.posLst
		eigenZheevr!( degObj.HmatLst[ii], degObj.Elst[ii], degObj.vecLst[ii], lwork, lrwork, liwork; work=view(workLst,:,Threads.threadid()), rwork=view(rworkLst,:,Threads.threadid()), iwork=view(iworkLst,:,Threads.threadid()) );
	end
	end
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	# @infiltrate
	if isGC
		@info( "GC: " )
		tFull = @timed GC.gc();
		# @infiltrate
		# asdlfja;sdkfj;adkja
		@info( timeMemStr( tFull.time, tFull.bytes ) )
	end
	
	return degObj.Elst, degObj.vecLst;
end

function locator_Eigen_threaded_comments!( degObj, HmatFun, )
	@info( "Create Hmats" )
	tFull = @timed begin
		# param_arr = zeros(degObj.param_dim);
		# for ind in degObj.posLst
			# # param_arr = [ degObj.param_grids[j][ind[j]] for j = 1 : degObj.param_dim ];
			# for j = 1 : degObj.param_dim
				# param_arr[j] = param_grids[j][ind[j]];
			# end
			# # HmatLst[ind] = zeros(Complex{Float64},N,N);
			# # HmatFun( HmatLst[ind], param_arr... );
			# HmatLst[ind] = HmatFun( param_arr... );
		# end
		# degObj.HmatLst .= HmatFun.( degObj.param_mesh );
		Threads.@threads for ind in degObj.posLst
			# degObj.HmatLst[ind] = HmatFun( degObj.param_mesh[ind] );
			HmatFun( degObj.HmatLst[ind], degObj.param_mesh[ind] );
		end
		@debug "HmatLst" degObj.HmatLst[1,2,1];
	end
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	
	nThreads = Threads.nthreads();
	lwork, lrwork, liwork = eigenPrework!( degObj.HmatLst[1], degObj.Elst[1], degObj.vecLst[1] );
	workLst = Array{ComplexF64,2}(undef, Int64(lwork),nThreads);
	rworkLst = Array{Float64,2}(undef, lrwork,nThreads);
	iworkLst = Array{BlasInt,2}(undef, liwork,nThreads);
	
	# lwork, lrwork, liwork = eigenPrework!( degObj.HmatLst[1], degObj.Elst[1], degObj.vecLst[1] );
	# workLst = Array{ComplexF64,4}(undef, Int64(lwork),degObj.param_divide...);
	# rworkLst = Array{Float64,4}(undef, lrwork,degObj.param_divide...);
	# iworkLst = Array{BlasInt,4}(undef, liwork,degObj.param_divide...);
	
	@info( "Eigen: " )
	@time begin
	Threads.@threads for ii in degObj.posLst
		# val, degObj.vecLst[ii] = eigen( degObj.HmatLst[ii] );
		# degObj.Elst[ii] = Real.(val);
		
		
		# eigenZheevr!( degObj.HmatLst[ii], degObj.Elst[ii], degObj.vecLst[ii], lwork, lrwork, liwork; work=view(workLst,:,ii), rwork=view(rworkLst,:,ii), iwork=view(iworkLst,:,ii) );
		
		
		eigenZheevr!( degObj.HmatLst[ii], degObj.Elst[ii], degObj.vecLst[ii], lwork, lrwork, liwork; work=view(workLst,:,Threads.threadid()), rwork=view(rworkLst,:,Threads.threadid()), iwork=view(iworkLst,:,Threads.threadid()) );
	end
	end
	@debug "Elst" degObj.Elst[1,2,1];
	@debug "vecLst" degObj.vecLst[1,2,1];
	
	return degObj.Elst, degObj.vecLst;
end

function locator_Eigen_CUDA( N, HmatFun, param_dim, param_divide, param_grids, dummyMat )
	param_grids_multi = Array{Array{Float64,param_dim},1}();
	# (undef,param_dim);
	@debug typeof(param_grids_multi);
	for iDim = 1 : param_dim
		shDim = ones(Int64,param_dim);
		shDim[iDim] = param_divide[iDim];
		push!(param_grids_multi, reshape( param_grids[iDim], shDim... ) );
	end
	
	HmatArr = zeros( Complex{Float32}, N, N, param_divide... );
	@info( "Create Hmats arrays" )
	@time begin
		@debug( size(HmatArr) )
		param_arr = zeros(param_dim);
		for ind in CartesianIndices(dummyMat)
			for j = 1 : param_dim
				param_arr[j] = param_grids[j][ind[j]];
			end
			HmatFun( view( HmatArr, :,:,ind ), param_arr... );
		end
	end
	
	@info( "cuArray Hmat" )
	@time HmatLstCu = cu( reshape( HmatArr, (N,N,:) ) );
	# paramMesh = zeros( param_divide..., param_dim );
	
	# @debug "Hmat" HmatLstCu[:,:,1]
	
	@info( "cuArray solve" )
	@time sols = CUDA.CUSOLVER.heevjBatched!('V','U',HmatLstCu);
	
	ElstCpu = Array( reshape( sols[1], N, param_divide... ) );
	# @debug "ElstCpu" ElstCpu[:,1:3,1,1]
	# readline();
	vecLstCpu = Array( reshape( sols[2], N, N, param_divide... ) );
	
	Elst = [ ElstCpu[:,ind] for ind in CartesianIndices(dummyMat) ];
	vecLst = [ vecLstCpu[:,:,ind] for ind in CartesianIndices(dummyMat) ];
	
	return Elst, vecLst;
end

function locator_div(degObj, HmatFun; isGOE3d = false)
	dimLst = [1:degObj.param_dim;];
	
	locator_Eigen_threaded!( degObj, HmatFun );
	
	@info( "link" );
	tFull = @timed linkLstFun!( degObj );
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	@info( "Bfield" )
	tFull = @timed BfieldFromLink!( degObj );
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	@debug "Bfield" degObj.BfieldLst[1][2,2,2];
	
	if isGOE3d
		non0Lst = ( q -> (p-> real.(p)).(q) ).(degObj.BfieldLst);
	elseif degObj.param_dim == 2
		non0Lst = (p-> real.(p)).(degObj.BfieldLst[1]);
	elseif degObj.param_dim >= 3 && !isGOE3d
		@info "divB"
		tFull = @timed divB!( degObj );
		@info( timeMemStr( tFull.time, tFull.bytes ) )
		non0Lst = degObj.divBLst;
	end
	
	threshold = 0.00001;
	divBlstReal = (p-> real.(p)).(degObj.divBLst);
	if !degObj.nonPeriodic
		clipBnd = 0;
	else
		clipBnd = 1;
	end
	@info "locateNon0"
	if !isGOE3d
		tFull = @timed posCount, negCount, posLocLst, negLocLst = locateNon0Div( non0Lst, threshold, degObj.N, degObj.param_dim, clipBnd = clipBnd );
		@info( timeMemStr( tFull.time, tFull.bytes ) )
		
		if degObj.param_dim == 2
			posLocLst = vcat.(posLocLst, negLocLst);
			posCount = posCount .+ negCount;
		end
	else
		posCount = []; negCount = [];
		posLocLst = []; negLocLst = [];
		BdimLst = [3,2,1];
		for iB = 1 : length(degObj.BfieldLst)
			tFull = @timed posCount1, negCount1, posLocLst1, negLocLst1 = locateNon0Div( non0Lst[iB], threshold, degObj.N, degObj.param_dim; clipBnd = clipBnd, isGOE3d = isGOE3d, Bdim = BdimLst[iB] );
			@info( timeMemStr( tFull.time, tFull.bytes ) )
			# @infiltrate
			posLocLst1 = vcat.(posLocLst1, negLocLst1);
			posCount1 = posCount1 .+ negCount1;
			append!( posCount, [posCount1] );
			append!( negCount, [negCount1] );
			append!( posLocLst, [posLocLst1] );
			append!( negLocLst, [negLocLst1] );
		end
	end
	
	return posCount, negCount, posLocLst, negLocLst, degObj.BfieldLst, divBlstReal;
end

function locator_div_fine_uneditted(N, HmatFun, param_dim, param_divide, param_max, param_min; maxDepth = 0, depth = 0, nSpecific = 0, dim_hold = 0)
	dimLst = collect(1:param_dim);
	param_divide = Int64.( param_divide );
	if !( param_divide isa Array{Int64} )
		param_divide_num = param_divide;
		param_divide = param_divide_num * ones(param_dim);
		param_divide = Int64.( param_divide );
	end
	param_step = ( param_max - param_min ) ./ param_divide;
	param_grids = [ collect( range( param_min[iDim], param_max[iDim] - param_step[iDim], length = param_divide[iDim] ) ) for iDim = 1:param_dim ];
	dummyMat = zeros( Tuple(param_divide) );
		
	@time HmatLst = [ zeros( Complex{Float64}, N, N ) for i1 in CartesianIndices(dummyMat) ];
	
	@time begin
		param_arr = zeros(param_dim);
		for ind in CartesianIndices(dummyMat)
			for j = 1 : param_dim
				param_arr[j] = param_grids[j][ind[j]];
			end
			HmatLst[ind] = HmatFun( param_arr... );
		end
	end
	
	@time eigenLst = eigen.( HmatLst );
	Elst = (p-> Real.( p.values )).( eigenLst );
	vecLst = (p->p.vectors).( eigenLst );
	if depth > 0
		Nlimited = 1;
		vecLst = (x->x[:,nSpecific]).(vecLst);
		Elst = (x->[x[nSpecific]]).(Elst);
		clipBnd = 1;
	else
		Nlimited = N;
		clipBnd = 0;
	end
	@time EigenSort!.( Elst, vecLst );
	
	HmatLst = nothing;
	GC.gc();

	@time linkLst = linkLstFun( param_dim, vecLst );
	@time Bfield_ln, BfieldLst = BfieldFromLink( param_dim, linkLst );
	if param_dim == 2
		non0Lst = BfieldLst[1];
		divBlst = [];
	elseif param_dim >= 3
		@time divBlst = divB( BfieldLst, Bfield_ln, N );
		non0Lst = divBlst;
	end
	
	threshold = 0.00001;
	divBlstReal = (p-> real.(p)).(divBlst);
	posCount, negCount, posLocLst, negLocLst = locateNon0Div( non0Lst, threshold, N, param_dim );
	
	if param_dim == 2
		posLocLst = vcat.(posLocLst, negLocLst);
		posCount = posCount .+ negCount;
	end
	
	return posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal;
end

function locator_div_fine_recursive(N, HmatFun, param_dim, param_divide, param_max, param_min; maxDepth = 0, depth = 0, nSpecific = 0, dim_hold = 0, isPeriodic = true, param_divide_next = nothing)
	# println(depth);
	if param_divide_next == nothing
		param_divide_next = 2 * ones(Int64, param_dim);
	end
	param_divide_fail_incre = ones(Int64, param_dim);
	for dim = 1 : param_dim
		if param_divide_next[dim] == 1
			param_divide_fail_incre[dim] = 0;
		end
	end
	# println("incre:", param_divide_fail_incre);
	dimLst = collect(1:param_dim);
	param_divide = Int64.( param_divide );
	if !( param_divide isa Array{Int64} )
		param_divide_num = param_divide;
		param_divide = param_divide_num * ones(Int64, param_dim);
		# param_divide = Int64.( param_divide );
	end
	param_step = ( param_max - param_min ) ./ param_divide;
	if isPeriodic
		param_divide_bnd = param_divide;
		grid_max = param_max .- param_step;
	else
		param_divide_bnd = param_divide .+ 1;
		grid_max = param_max;
	end
	param_grids = [ collect( range( param_min[iDim], grid_max[iDim], length = param_divide_bnd[iDim] ) ) for iDim = 1:param_dim ];
	dummyMat = zeros( Tuple(param_divide_bnd) );
		
	HmatLst = [ zeros( Complex{Float64}, N, N ) for i1 in CartesianIndices(dummyMat) ];
	
	
	param_arr = zeros(param_dim);
	for ind in CartesianIndices(dummyMat)
		for j = 1 : param_dim
			param_arr[j] = param_grids[j][ind[j]];
		end
		HmatLst[ind] = HmatFun( param_arr... );
	end
	
	
	eigenLst = eigen.( HmatLst );
	Elst = (p-> Real.( p.values )).( eigenLst );
	vecLst = (p->p.vectors).( eigenLst );
	if depth > 0
		Nlimited = 1;
		vecLst = (x->x[:,nSpecific]).(vecLst);
		Elst = (x->[x[nSpecific]]).(Elst);
		clipBnd = 1;
	else
		Nlimited = N;
		clipBnd = 0;
	end
	EigenSort!.( Elst, vecLst );
	
	HmatLst = nothing;
	if depth == 0
		GC.gc();
	end

	linkLst = linkLstFun( param_dim, vecLst );
	Bfield_ln, BfieldLst = BfieldFromLink( param_dim, linkLst );
	if param_dim == 2
		non0Lst = BfieldLst[1];
		divBlst = [];
	elseif param_dim >= 3
		divBlst = divB( BfieldLst, Bfield_ln, Nlimited );
		non0Lst = divBlst;
	end
	
	threshold = 0.00001;
	divBlstReal = (p-> real.(p)).(divBlst);
	posCount, negCount, posLocLst, negLocLst = locateNon0Div( non0Lst, threshold, Nlimited, param_dim; clipBnd = clipBnd );
	
	numPol = 2;
	if param_dim == 2
		posLocLst = vcat.(posLocLst, negLocLst);
		posCount = posCount .+ negCount;
		numPol = 1;
	end
	
	locToParamHere = x -> locToParamArr( x, param_grids, param_dim );
	
	posParamLst = locToParamHere.( posLocLst );
	negParamLst = locToParamHere.( negLocLst );
	
	paramLstPol = [posParamLst, negParamLst];
	locLstPol = [posLocLst, negLocLst];
	countPol = [posCount, negCount];
	
	# println( "recursive check: ", maxDepth > 0 && depth < maxDepth );
	# readline();
	
	# @time begin
	notFoundCnt = 0;
	finerFound = false;
	if maxDepth > 0 && depth < maxDepth
		# println("recursive start");
		for iPol = 1 : numPol
			paramLst = paramLstPol[iPol];
			locLst = locLstPol[iPol];
			degCount = countPol[iPol];
			for n = 1 : Nlimited
				itStop = degCount[n];
				if depth > 0
					itStop = min(1,itStop);
				end
				for it = 1 : itStop
					param = paramLst[n][it,:];
					param_min_next = param;
					param_max_next = param .+ param_step;
					# if dim_hold != 0
						# param_max_next[dim_hold] = param + param_step;
						# param_divide_next[dim_hold] = 1
					# end
					depthNext = depth + 1;
					if depth > 0
						nSpecificNext = nSpecific;
					else
						nSpecificNext = n;
					end
					trNum = 3;
					param_divide_next_more = param_divide_next;
					finerFound = false;
					for iTr = 1 : trNum
						# println("divqide_shallow", param_divide_next);
						# readline();
						posCountFine,negCountFine,posParamLstFine,negParamLstFine, BfieldLstFine, divBlstFine = locator_div_fine_recursive( N, HmatFun, param_dim, param_divide_next_more, param_max_next, param_min_next; maxDepth=maxDepth, depth=depthNext,nSpecific=nSpecificNext, isPeriodic = false, param_divide_next = param_divide_next );
						# readline();
						paramLstFinePol = [posParamLstFine, negParamLstFine];
						countFinePol = [posCountFine,negCountFine];
						countFine = countFinePol[iPol];
						paramLstFine = paramLstFinePol[iPol];
						nLimInd = 1;
						# println( "deeper_cnt:", countFine[nLimInd], " ", countFine[nLimInd] > 0 );
						if countFine[nLimInd] > 0
							paramLst[n][it,:] = paramLstFine[nLimInd][1,:];
							# @info("found");
							finerFound = true;
							break;
						else 
							# println("divide_deep_before", param_divide_next_more);
							param_divide_next_more = param_divide_next_more .+ param_divide_fail_incre;
							# @info( "not found" );
							# println( "original divB: ", non0Lst[posLocLst[n][it,:]...][n] );
							# println("divide_deep", param_divide_next_more);
							# readline();
						end
					end
				end
				if !finerFound
					notFoundCnt += 1;
				end
			end
		end
	end
	# end
	
	if depth == 0
		println("notFoundCnt: ", notFoundCnt);
	end
	
	return posCount, negCount, posParamLst, negParamLst, BfieldLst, divBlstReal, posLocLst, negLocLst;
end

function locator_div_fine(N, HmatFun, param_dim, param_divide_full, param_max, param_min; maxDepth = 0, param_divide_next = nothing)
	if param_divide_next == nothing
		param_divide_next = 2 * ones(Int64, param_dim);
	end
	
	param_divide = Int64.( param_divide_full ./ (param_divide_next.^maxDepth) );
	
	posCount, negCount, posParamLst, negParamLst, BfieldLst, divBlstReal, posLocLst, negLocLst = locator_div_fine_recursive( N, HmatFun, param_dim, param_divide, param_max, param_min; maxDepth = maxDepth, depth = 0, param_divide_next = param_divide_next ); 
	
	param_step = ( param_max .- param_min ) ./ param_divide_full;
	
	paramToLocArrHere = x -> paramToLocArr( x, param_dim, param_min, param_step );
	
	posLocLst = paramToLocArrHere.( posParamLst );
	negLocLst = paramToLocArrHere.( negParamLst );
	
	return posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal;
end

function locator_div_GUE( N, param_divide_num; Hmat_seed = 1234, polarity = +1 )
	param_dim = 3;
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	
	H_GUElst = HlstFunc( H_GUE, param_dim, N );
	HmatFun = (x,y,z)->Hmat_3comb([x,y,z],H_GUElst);
	
	posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min ); 
	return posCount, negCount, posLocLst, negLocLst, H_GUElst, BfieldLst, divBlstReal;
end

function locator_div_GUE_scale( N, param_divide_num, scale; Hmat_seed = 1234, polarity = +1 )
	param_dim = 3;
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	
	H_GUElst = HlstFunc( H_GUE, param_dim, N );
	HmatFun = (x,y,z)-> scale .* Hmat_3comb([x,y,z],H_GUElst);
	
	return locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min );
end

function locator_div_GUE_ratio( degObj; ratio=[1,1,1], Hmat_seed = 1234, polarity = +1, maxDepth = 0, param_divide_next = [1,2,2], isOrtho = false, alpha = 0 )
	if degObj.param_dim !=3
		@error "GUE should be dim 3"
	end
	
	H_GUElst = HlstFunc( H_GUE, degObj.param_dim, degObj.N, isOrtho );
	H7 = H_GUE( degObj.N );
	HmatFun = (H,xlst)-> Hmat_3comb_ratio_Hoffset!(H,xlst,H_GUElst, ratio, H7; c1 = sqrt(1-alpha) / sqrt(3), cOff = sqrt(alpha) );
	
	if maxDepth == 0
		posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div( degObj, HmatFun ); 
	else 
		posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div_fine( N, HmatFun, param_dim, param_divide_num, param_max, param_min; maxDepth = maxDepth, param_divide_next = param_divide_next );
	end
	
	return posCount, negCount, posLocLst, negLocLst, H_GUElst, BfieldLst, divBlstReal;
end

function locator_div_GOE_ratio( degObj; ratio=[1,1], Hmat_seed = 1234, polarity = +1, maxDepth = 0, param_divide_next = [1,2], isOrtho = false, alpha = 0 )
	isGOE3d = false;
	if degObj.param_dim == 3
		isGOE3d = true;
	end
	H_GOElst = HlstFunc( H_GOE, degObj.param_dim, degObj.N, isOrtho );
	H7 = H_GOE( degObj.N );
	HmatFun = (H,xlst) -> Hmat_3comb_ratio_Hoffset!(H,xlst,H_GOElst, ratio, H7; c1 = sqrt(1-alpha) / sqrt(3), cOff = sqrt(alpha) );
	# @infiltrate
	
	if maxDepth == 0
		posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div( degObj, HmatFun; isGOE3d = isGOE3d ); 
	else 
		posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div_fine( N, HmatFun, param_dim, param_divide_num, param_max, param_min; maxDepth = maxDepth, param_divide_next = param_divide_next );
	end
	
	if alpha == 0
		H_lst_full = H_GOElst;
	else
		H_lst_full = [H_GOElst, H7];
	end
	
	return posCount, negCount, posLocLst, negLocLst, H_lst_full, BfieldLst, divBlstReal;
end

# function locator_div_GOE_ratio( degObj; ratio=[1,1], Hmat_seed = 1234, polarity = +1, maxDepth = 0, param_divide_next = [1,2], isOrtho = false )
	# H_GOElst = HlstFunc( H_GOE, degObj.param_dim, degObj.N, isOrtho );
	# HmatFun = (H,xlst)-> Hmat_3comb_ratio!(H,xlst,H_GOElst, ratio);
	
	# if maxDepth == 0
		# posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div( degObj, HmatFun ); 
	# else 
		# posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div_fine( N, HmatFun, param_dim, param_divide_num, param_max, param_min; maxDepth = maxDepth, param_divide_next = param_divide_next );
	# end
	
	# return posCount, negCount, posLocLst, negLocLst, H_GOElst, BfieldLst, divBlstReal;
# end

function locator_div_sin3( param_divide_num )
	N=2;	
	HmatFun = Hmat_3sin;
	
	param_dim = 3;
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	
	return locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min );
end