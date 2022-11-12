using Utils
using Statistics
# using JLD
# using Infiltrator
using JLD2
using Logging; using LoggingExtras;

function shLocLst( locLst, param_divide_num, dim, shId )
	if shId == 1
		return locLst;
	end
	
	locLstCell = broadcast.( -, locLst, 0.5 );
	if shId == 2
		for locs in locLstCell
			locs[:,dim] = mod.( param_divide_num .- locs[:,dim], param_divide_num );
		end
	elseif shId == 3
		for locs in locLstCell
			locs[:,dim] = mod.( param_divide_num/2 .- locs[:,dim], param_divide_num );
		end
	elseif shId == 4
		for locs in locLstCell
			locs[:,dim] = mod.( param_divide_num/2 .+ locs[:,dim], param_divide_num );
		end
	end
	
	locLst = broadcast.( Int64, broadcast.( +, locLstCell, 0.5 ) );
	return locLst;
end

function divB_profile( mSz, divLst, itNum, seedFed; nDim = 3 )
	if seedFed > 0
		Random.seed!(seedFed);
	end
	minNum = 0;
	maxNum = 2*pi;
	paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, nDim );
	matsFull = matsGridHThreaded( paramsFull, threaded_zeros(ComplexF64,mSz,mSz) );
	degBerrysFull = degBerrysInit( paramsFull, matsFull );
	
	totalNumLst = zeros(Int64,itNum);
	
	for it = 1 : itNum
		print( "\r Iteration: " * string(it) * " / " * string(itNum) );
		Hlst = DegLocatorDiv.HlstFunc(H_GUE,nDim,mSz);
		HmatFun = (H,xLst) -> Hmat_3comb_ratio!( H, xLst, Hlst );
		startNextEigen( matsFull );
		
		@info("Eigen:")
		tFull = @timed eigenAll( matsFull; HmatFun = HmatFun );
		@info( timeMemStr( tFull.time, tFull.bytes ) )

		@info("Link:")
		tFull = @timed DegLocatorDiv.linksCalcAll( degBerrysFull );
		@info( timeMemStr( tFull.time, tFull.bytes ) )
		
		@info("Bfield:")
		tFull = @timed DegLocatorDiv.BfieldCalcAll( degBerrysFull );
		@info( timeMemStr( tFull.time, tFull.bytes ) )
		
		@info("DivB:")
		tFull = @timed DegLocatorDiv.divBCalcAll( degBerrysFull );
		@info( timeMemStr( tFull.time, tFull.bytes ) )
		
		totalNumLst[it] = sum(sum( (x->abs.(x).>1e-9).(degBerrysFull.divBLst) ));
	end
	return totalNumLst;
end

function divB_profile( param_dim = 2, Nlst=[20], num_it_main=2, param_divide=50, seedFedstart = -1, enumSaveMem = memNone; scale=degOptDefaultLst[1], ratioFed=degOptDefaultLst[2], fileNameMod = "", maxDepth = 0, param_divide_next = nothing, isOrtho = false, alpha = degOptDefaultLst[3], thresNM = rtFndDefaultLst[1], thresEDeg = rtFndDefaultLst[2], thresDegCollision = rtFndDefaultLst[2] )
	ratio = ratioFed;
	if ratioFed == degOptDefaultLst[2]
		ratio = ones(Int64,param_dim);
	end
	if maxDepth > 0 && param_divide_next == nothing
		param_divide_next = 2 .* ones(Int64, param_dim);
		param_divide_next[1] = 1;
	end
	
	if param_dim == 2
		locator_func_ratio = locator_div_GOE_ratio;
	else
		locator_func_ratio = locator_div_GUE_ratio;
	end
	param_min = fill(0.0,param_dim);
	param_max = fill(2*pi,param_dim);
	
	# symFact = 4;
	# indCartPosNeg = CartesianIndices( Tuple( symFact * ones(Int64,param_dim) ) );
	seedFed = seedFedstart;
	if seedFed > 0
		Random.seed!(seedFed);
	end
	for itN = 1 : length(Nlst)
		tFullN = @timed begin
		N = Nlst[itN];
		# num_it_sym = symFact ^ param_dim;
		# num_it = num_it_main * num_it_sym;
		num_it = num_it_main;
		numOfPosN = ( enumSaveMem == rootFind ? N-1 : N );
		posNlst = zeros( Int64, num_it, numOfPosN );
		negNlst = zeros( Int64, num_it, numOfPosN );
		if enumSaveMem != rootFind
			locTp = Int64
		else 
			locTp = Float64
		end
		posLocLst = Vector{ Vector{ Array{ locTp,2 } } }(undef, num_it);
		negLocLst = Vector{ Vector{ Array{ locTp,2 } } }(undef, num_it);
		# posTotalLst = zeros( Int64, num_it );
		# negTotalLst = zeros( Int64, num_it );
		
		degObj = DegObj( N, param_dim, param_divide, param_min, param_max; enumSaveMem = enumSaveMem );
		
		H_GUE_lst = [ [ zeros(Complex{Float64},N,N) for iCosSin = 1:2, iDim = 1 : param_dim ] for it = 1:num_it ];
		for it = 1:num_it
			# itStart = itMain * num_it_sym + 1;
			# it = itStart;
			print( "\r Iteration: " * string(it) * " / " * string(num_it) );
			posNlst[it,:], negNlst[it,:], posLocLst[it], negLocLst[it], H_GUE_lst[it] = locator_func_ratio( degObj; ratio = ratio, Hmat_seed = seedFed, maxDepth = maxDepth, param_divide_next = param_divide_next, isOrtho = isOrtho, alpha = alpha );
			# posTotalLst[it] = sum(@view(posNlst[it,:]));
			# negTotalLst[it] = sum(@view(negNlst[it,:]));
			@info( "GC: " );
			tFull = @timed GC.gc();
			@info( timeMemStr( tFull.time, tFull.bytes ) )
		end
		print("\n");
		fMain = "deg";
		attrLst, valLst = fAttrOptLstFunc( N, param_divide, num_it, seedFedstart; dim = param_dim, scale = scale, ratio = ratioFed, alpha = alpha, enumSaveMem = enumSaveMem, thresNM = thresNM, thresEDeg = thresEDeg );
		# attrLst, valLst = fNameAttrLstFunc( N, param_divide, num_it, seedFedstart; dim = param_dim );
		# @infiltrate
		
		varFileName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fileNameMod );
		
		logFileName = fNameFunc( "log_" * fMain, attrLst, valLst, ".log"; fMod = fileNameMod );
		logIO = open(logFileName, "w");
		
		printlnLog = (msg...)->printlnLogFile( logFileName, msg... );
		
		posnegN = [ posNlst, negNlst ];
		# posnegTotal = [ posTotalLst, negTotalLst ];
		itDim = 1;
		eLevDim = 2;
		meanPosNegN = [ funArrDims( mean, Nlst, itDim ) for Nlst in posnegN ];
		stdPosNegN = [ funArrDims( std, Nlst, itDim ) for Nlst in posnegN ];
		posnegTotal = [ funArrDims( sum, Nlst, eLevDim ) for Nlst in posnegN ];
		# meanPosNegN = [ [mean(col) for col in	eachcol(Nlst)] for Nlst in posnegN ];
		# stdPosNegN = [ [std(col) for col in eachcol(Nlst)] for Nlst in posnegN ];
		meanPosNegTotal = [ mean(totalLst) for totalLst in posnegTotal ];
		stdPosNegTotal = [ std(totalLst) for totalLst in posnegTotal ];
		
		printlnLog( "matrix dimension: ", N );
		NnameLst = ["N+", "N-"];
		NEnd = ( enumSaveMem != rootFind ? N : N-1 );
		for it = 1:2
			printlnLog( NnameLst[it], ": " );
			for n = 1:NEnd
				printlnLog( "level ", n, ": ", round( meanPosNegN[it][n]; digits=3 ), " +- ", round( stdPosNegN[it][n]; digits=3 ) );
			end
			printlnLog( "total  : ", round( meanPosNegTotal[it]; digits=3 ), " +- ", round( stdPosNegTotal[it]; digits=3 )  );
		end
		
		meanPosNegNarr = cat( meanPosNegN...; dims = ndims(meanPosNegN[1])+1 );
		@info( size(meanPosNegN) );
		@info(length(size(meanPosNegN[1]))+1);
		@info(ndims(meanPosNegN[1])+1);
		@info( meanPosNegN );
		meanPosNegNarr = permutedims( meanPosNegNarr, [2,1] );
		
		save( varFileName, "N", N, "seed", seedFedstart, "posNlst", posNlst, "negNlst", negNlst, "posLocLst", posLocLst, "negLocLst", negLocLst, "H_GUE_lst", H_GUE_lst, "meanPosNegN", meanPosNegN );
		# @infiltrate
		
		gLstNameMain = "gLst";
		# if degObj.dim == 2
			# gLstNameMain = string( gLstNameMain, "_dim_", dim );
		gNpyName = fNameFunc( gLstNameMain, attrLst, valLst, npyType; fMod = fileNameMod );
		
		npzwrite( gNpyName, meanPosNegNarr );
		end
		tFileName = fNameFunc( "t_" * fMain, attrLst, valLst, jld2Type; fMod = fileNameMod );
		save( tFileName, "tFullN", tFullN );
	end
end

function divB_profile_GOE_3d_test( Nlst=[10], param_divide_num=80, itNum = 1, seedFedstart = -1; scale=1, ratio=nothing, fileNameMod = "", maxDepth = 0, param_divide_next = nothing, isOrtho = false, isSaveMem = false )
	if isSaveMem
		fModSaveMem = "memLayers";
		if fileNameMod == ""
			fileNameMod =fModSaveMem;
		else
			fileNameMod = fileNameMod * "_" * fModSaveMem;
		end
	end
	param_dim = 3;
	if ratio == nothing
		ratio = ones(Int64,param_dim);
	end
	if maxDepth > 0 && param_divide_next == nothing
		param_divide_next = 2 .* ones(Int64, param_dim);
		param_divide_next[1] = 1;
	end
	
	locator_func_ratio = locator_div_GOE_ratio;
	
	param_min = fill(0.0,param_dim);
	param_max = fill(2*pi,param_dim);
	
	param_divide = fill( param_divide_num, param_dim );
	
	seedFed = seedFedstart;
	if seedFed > 0
		Random.seed!(seedFed);
	end
	for itN = 1 : length(Nlst)
		N = Nlst[itN];
		posTotalLst = zeros( Int64, param_dim, param_divide_num );
		negTotalLst = zeros( Int64, param_dim, param_divide_num );
		
		posNlst = Vector{Vector{Vector{Vector{Int64}}}}(undef,itNum);
		negNlst = Vector{Vector{Vector{Vector{Int64}}}}(undef,itNum);
		posLocLst = Vector{Vector{Array{ Array{ Int64,2 }, 2 }}}(undef,itNum);
		negLocLst = Vector{Vector{Array{ Array{ Int64,2 }, 2 }}}(undef,itNum);
		H_GOElst = ( alpha == 0 ? [ [ zeros(Complex{Float64},N,N) for iCosSin = 1:2, iDim = 1 : param_dim ] for it = 1:itNum ] : 
		[ [ [ zeros(Complex{Float64},N,N) for iCosSin = 1:2, iDim = 1 : param_dim ], zeros(Complex{Float64},N,N) ] for it = 1:itNum ] );
		parity2dLst = [ [ [ones(Complex{Float64},N) for id1=1:param_divide_num,id2=1:param_divide_num] for dim = 1 : param_dim ] for it = 1 : itNum];
		
		degObj = DegObj( N, param_dim, param_divide, param_min, param_max; isStrongGC = true, isSaveMem = isSaveMem );
		# @infiltrate
		for it = 1 : itNum
			print( "\r Iteration: " * string(it) );
			colabIOreset();
			# posNlst[it], negNlst[it], posLocLst[it], negLocLst[it], H_GOElst[it] = locator_func_ratio( degObj; ratio = ratio, Hmat_seed = seedFed, maxDepth = maxDepth, param_divide_next = param_divide_next );
			# @infiltrate
			for dim = 1 : param_dim
				# parity2dLst[it][dim] = ones(param_divide_num, param_divide_num);
				for id = 1 : param_divide_num
					print( "\r Iteration: " * string(it) * ", dim = " * string(dim) * ", layer: " * string(id) );
					# @infiltrate
					parity2dLst[it][dim] .= ((a,b)->(a.*b)).( parity2dLst[it][dim], selectdim(degObj.linkLst[dim],dim,id) );
				end
				# parity2dLst[it][dim] .= (x->log.(x)).(parity2dLst[it][dim]) ./ 1im;
			end
			@info( "GC: " );
			tFull = @timed GC.gc();
			@info( timeMemStr( tFull.time, tFull.bytes ) )
			colabIOreset();
		end
		print("\n");
		# num_it = 1;
		fileNameAttr = fileNameAttrFunc( N, param_divide, itNum, seedFedstart );
		fileNameMain = string( "deg_GOE_3d_full", "_dim_", param_dim );
		if scale != 1'
			fileNameAttr = string( fileNameAttr, "_scale_", scale );
		end
		if ratio != ones( Int64, param_dim )
			fileNameAttr = string( fileNameAttr, "_ratio_", ratio );
		end
		if fileNameMod != ""
			fileNameAttr = string( fileNameAttr, "_", fileNameMod );
		end
		oFileName = string( fileNameMain, "_", fileNameAttr );
		varFileName = string( oFileName, jld2Type );
		
		logFileName = string( "log_", oFileName, ".log" );
		logIO = open(logFileName, "w");
		
		printlnLog = (msg...)->printlnLogFile( logFileName, msg... );
		
		# posnegN = [ posNlst, negNlst ];
		# posnegTotal = (lst3->(lst2->(lst1->sum.(lst1)).(lst2)).(lst3)).(posnegN);
		# meanPosNegN = (lst2->(lst->mean.(lst)).(lst2)).(posnegN);
		# stdPosNegN =  (lst2->(lst->std.(lst)).(lst2)).(posnegN);
		# meanPosNegTotal = (lst2->(lst->mean.(lst)).(lst2)).(posnegTotal);
		# stdPosNegTotal = (lst2->(lst->std.(lst)).(lst2)).(posnegTotal);
		# @infiltrate
		
		# printlnLog( "matrix dimension: ", N );
		# NnameLst = ["N+", "N-"];
		# dimRep = 1;
		# it = 1;
		# for ip = 1:2
			# printlnLog( NnameLst[ip], ": " );
			# for n = 1:N
				# printlnLog( "level ", n, ": ", round( meanPosNegN[ip][it][dimRep][n]; digits=3 ), " +- ", round( stdPosNegN[ip][it][dimRep][n]; digits=3 ) );
			# end
			# printlnLog( "total  : ", round( meanPosNegTotal[ip][it][dimRep]; digits=3 ), " +- ", round( stdPosNegTotal[ip][it][dimRep]; digits=3 )  );
		# end
		# @infiltrate
		
		# save( varFileName, "N", N, "seed", seedFedstart, "posNlst", posNlst, "negNlst", negNlst, "posLocLst", posLocLst, "negLocLst", negLocLst, "H_GOElst_All", H_GOElst, "meanPosNegN", meanPosNegN, "parity2dLst", parity2dLst );
		
		gLstNameMain = "gLst_GOE3dfull";
		gNpyNameMain = string( gLstNameMain, "_", fileNameAttr, npyType );
	end
	# @infiltrate
	degObj = nothing;
	# @infiltrate
	@info( "GC: " )
	tFull = @timed GC.gc();
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	# @infiltrate
end

function divB_profile_GOE_3d( Nlst=[10], param_divide_num=80, itNum = 1, seedFedstart = -1; scale=1, ratio=nothing, fileNameMod = "", maxDepth = 0, param_divide_next = nothing, isOrtho = false, enumSaveMem = memNone, alpha = 0 )
	# if isSaveMem
		# fModSaveMem = "memLayers";
		# if fileNameMod == ""
			# fileNameMod =fModSaveMem;
		# else
			# fileNameMod = fileNameMod * "_" * fModSaveMem;
		# end
	# end
	param_dim = 3;
	if ratio == nothing
		ratio = ones(Int64,param_dim);
	end
	if maxDepth > 0 && param_divide_next == nothing
		param_divide_next = 2 .* ones(Int64, param_dim);
		param_divide_next[1] = 1;
	end
	
	locator_func_ratio = locator_div_GOE_ratio;
	
	param_min = fill(0.0,param_dim);
	param_max = fill(2*pi,param_dim);
	
	param_divide = fill( param_divide_num, param_dim );
	
	seedFed = seedFedstart;
	if seedFed > 0
		Random.seed!(seedFed);
	end
	for itN = 1 : length(Nlst)
		N = Nlst[itN];
		posTotalLst = zeros( Int64, param_dim, param_divide_num );
		negTotalLst = zeros( Int64, param_dim, param_divide_num );
		
		posNlst = Vector{Vector{Vector{Vector{Int64}}}}(undef,itNum);
		negNlst = Vector{Vector{Vector{Vector{Int64}}}}(undef,itNum);
		posLocLst = Vector{Vector{Array{ Array{ Int64,2 }, 2 }}}(undef,itNum);
		negLocLst = Vector{Vector{Array{ Array{ Int64,2 }, 2 }}}(undef,itNum);
		H_GOElst = ( alpha == 0 ? [ [ zeros(Complex{Float64},N,N) for iCosSin = 1:2, iDim = 1 : param_dim ] for it = 1:itNum ] : 
		[ [ [ zeros(Complex{Float64},N,N) for iCosSin = 1:2, iDim = 1 : param_dim ], zeros(Complex{Float64},N,N) ] for it = 1:itNum ] );
		parity2dLst = [ [ [ones(Complex{Float64},N) for id1=1:param_divide_num,id2=1:param_divide_num] for dim = 1 : param_dim ] for it = 1 : itNum]
		
		degObj = DegObj( N, param_dim, param_divide, param_min, param_max; isStrongGC = true, enumSaveMem = enumSaveMem, thresNM = thresNM, thresEDeg = thresEDeg, thresDegCollision = thresDegCollision );
		# @infiltrate
		for it = 1 : itNum
			print( "\r Iteration: " * string(it) );
			colabIOreset();
			posNlst[it], negNlst[it], posLocLst[it], negLocLst[it], H_GOElst[it] = locator_func_ratio( degObj; ratio = ratio, Hmat_seed = seedFed, maxDepth = maxDepth, param_divide_next = param_divide_next, alpha = alpha );
			# @infiltrate
			for dim = 1 : param_dim
				# parity2dLst[it][dim] = ones(param_divide_num, param_divide_num);
				for id = 1 : param_divide_num
					print( "\r Iteration: " * string(it) * ", dim = " * string(dim) * ", layer: " * string(id) );
					# @infiltrate
					parity2dLst[it][dim] .= ((a,b)->(a.*b)).( parity2dLst[it][dim], selectdim(degObj.linkLst[dim],dim,id) );
				end
				parity2dLst[it][dim] .= (x->log.(x)).(parity2dLst[it][dim]) ./ 1im;
			end
			@info( "GC: " );
			tFull = @timed GC.gc();
			@info( timeMemStr( tFull.time, tFull.bytes ) )
			colabIOreset();
		end
		print("\n");
		# num_it = 1;
		fMain = "deg_GOE_3d_full";
		attrLst, valLst = fAttrValFunc( N, param_divide, itNum, seedFedstart, param_dim; alpha = alpha, scale = scale );
		fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fileNameMod );
		# fileNameAttr = fileNameAttrFunc( N, param_divide, itNum, seedFedstart );
		# fileNameMain = string( "deg_GOE_3d_full", "_dim_", param_dim );
		# if scale != 1
			# fileNameAttr = string( fileNameAttr, "_scale_", scale );
		# end
		# if ratio != ones( Int64, param_dim )
			# fileNameAttr = string( fileNameAttr, "_ratio_", ratio );
		# end
		# if alpha != 0
			# fileNameAttr = string( fileNameAttr, "_alpha_", alpha );
		# end
		# if fileNameMod != ""
			# fileNameAttr = string( fileNameAttr, "_", fileNameMod );
		# end
		# oFileName = string( fileNameMain, "_", fileNameAttr );
		# varFileName = string( fName, jld2Type );
		
		logFileName = string( "log_", fName, ".log" );
		logIO = open(logFileName, "w");
		
		printlnLog = (msg...)->printlnLogFile( logFileName, msg... );
		
		posnegN = [ posNlst, negNlst ];
		posnegTotal = (lst3->(lst2->(lst1->sum.(lst1)).(lst2)).(lst3)).(posnegN);
		meanPosNegN = (lst2->(lst->mean.(lst)).(lst2)).(posnegN);
		stdPosNegN =  (lst2->(lst->std.(lst)).(lst2)).(posnegN);
		meanPosNegTotal = (lst2->(lst->mean.(lst)).(lst2)).(posnegTotal);
		stdPosNegTotal = (lst2->(lst->std.(lst)).(lst2)).(posnegTotal);
		# @infiltrate
		
		printlnLog( "matrix dimension: ", N );
		NnameLst = ["N+", "N-"];
		dimRep = 1;
		it = 1;
		for ip = 1:2
			printlnLog( NnameLst[ip], ": " );
			for n = 1:N
				printlnLog( "level ", n, ": ", round( meanPosNegN[ip][it][dimRep][n]; digits=3 ), " +- ", round( stdPosNegN[ip][it][dimRep][n]; digits=3 ) );
			end
			printlnLog( "total  : ", round( meanPosNegTotal[ip][it][dimRep]; digits=3 ), " +- ", round( stdPosNegTotal[ip][it][dimRep]; digits=3 )  );
		end
		# @infiltrate
		
		save( fName, "N", N, "seed", seedFedstart, "posNlst", posNlst, "negNlst", negNlst, "posLocLst", posLocLst, "negLocLst", negLocLst, "H_GOElst_All", H_GOElst, "meanPosNegN", meanPosNegN, "parity2dLst", parity2dLst );
		
		gLstNameMain = "gLst_GOE3dfull";
		gNpyNameMain = fNameFunc( gLstNameMain, attrLst, valLst, npyType; fMod = fileNameMod );
	end
	degObj = nothing;
	# @infiltrate
	@info( "GC: " )
	tFull = @timed GC.gc();
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	# @infiltrate
end

function divB_profile_GOE_3d_3itr( Nlst=[10], param_divide_num=80, seedFedstart = -1; scale=1, ratio=nothing, fileNameMod = "", maxDepth = 0, param_divide_next = nothing, isOrtho = false )
	param_dim = 2;
	param_dim_full = 3;
	if ratio == nothing
		ratio = ones(Int64,param_dim);
	end
	if maxDepth > 0 && param_divide_next == nothing
		param_divide_next = 2 .* ones(Int64, param_dim);
		param_divide_next[1] = 1;
	end
	
	locator_func_ratio = locator_div_GOE_ratio;
	
	param_min = fill(0.0,param_dim);
	param_max = fill(2*pi,param_dim);
	
	param_divide = fill( param_divide_num, param_dim );
	
	param_min_3 = 0;
	param_max_3 = 2*pi;
	
	seedFed = seedFedstart;
	if seedFed > 0
		Random.seed!(seedFed);
	end
	for itN = 1 : length(Nlst)
		N = Nlst[itN];
		posNlst = zeros( Int64, param_dim_full, param_divide_num, N );
		negNlst = zeros( Int64, param_dim_full, param_divide_num, N );
		posLocLst = Array{ Vector{ Array{ Int64,2 } },2 }(undef, param_dim_full, param_divide_num);
		negLocLst = Array{ Vector{ Array{ Int64,2 } },2 }(undef, param_dim_full, param_divide_num);;
		posTotalLst = zeros( Int64, param_dim_full, param_divide_num );
		negTotalLst = zeros( Int64, param_dim_full, param_divide_num );
		
		degObj = DegObj( N, param_dim, param_divide, param_min, param_max );
		
		param_grid3 = degObj.param_grids[1];
		
		H_GOElst_All = HlstFunc( H_GOE, param_dim_full, degObj.N, isOrtho );
		dim3 = 1;
		# H_GOElst_3 = HlstFunc( H_GOE, dim3, degObj.N, isOrtho );
		
		# H_lst_full = [ [ zeros(Complex{Float64},N,N) for iCosSin = 1:2, iDim = 1 : param_dim ] for it = 1:param_divide_num ];
		H3 = zeros( Complex{Float64}, N, N );
		for id = 1 : param_dim_full
			H_GOElst_3 = H_GOElst_All[:,id];
			dimLst = collect(1:param_dim_full);
			deleteat!(dimLst,id);
			H_GOElst = H_GOElst_All[:,dimLst];
			for it = 1:param_divide_num
				param3 = param_grid3[it];
				Hmat_3comb_ratio!(H3, [param3],H_GOElst_3,ratio);
				HmatFun = (H,xlst) -> Hmat_3comb_ratio_Hoffset!(H,xlst,H_GOElst, ratio,H3);
				# @infiltrate	
				posNlst[id,it,:], negNlst[id,it,:], posLocLst[id,it], negLocLst[id,it] = locator_div( degObj, HmatFun );
				posTotalLst[id,it] = sum(posNlst[id,it,:]);
				negTotalLst[id,it] = sum(negNlst[id,it,:]);
				# @info( "GC: " );
				# @time GC.gc();
			end
		end
		num_it = 1;
		fileNameAttr = fileNameAttrFunc( N, param_divide, num_it, seedFedstart );
		fileNameMain = string( "deg_GOE_3d_full_itrZ", "_dim_", param_dim_full );
		if scale != 1'
			fileNameAttr = string( fileNameAttr, "_scale_", scale );
		end
		if ratio != ones( Int64, param_dim )
			fileNameAttr = string( fileNameAttr, "_ratio_", ratio );
		end
		if fileNameMod != ""
			fileNameAttr = string( fileNameAttr, "_", fileNameMod );
		end
		oFileName = string( fileNameMain, "_", fileNameAttr );
		varFileName = string( oFileName, jldType );
		
		logFileName = string( "log_", oFileName, ".log" );
		logIO = open(logFileName, "w");
		
		printlnLog = (msg...)->printlnLogFile( logFileName, msg... );
		
		posnegN = [ posNlst, negNlst ];
		posnegTotal = [ posTotalLst, negTotalLst ];
		meanPosNegN = [ dropdims( mean(Nlst; dims=2); dims=2 ) for Nlst in posnegN ];
		stdPosNegN = [ dropdims( std(Nlst; dims=2), dims=2 ) for Nlst in posnegN ];
		meanPosNegTotal = [ mean(totalLst; dims=2) for totalLst in posnegTotal ];
		stdPosNegTotal = [ std(totalLst; dims=2) for totalLst in posnegTotal ];
		
		printlnLog( "matrix dimension: ", N );
		NnameLst = ["N+", "N-"];
		dimRep = 1;
		for it = 1:2
			printlnLog( NnameLst[it], ": " );
			for n = 1:N
				printlnLog( "level ", n, ": ", round( meanPosNegN[it][dimRep,n]; digits=3 ), " +- ", round( stdPosNegN[it][dimRep,n]; digits=3 ) );
			end
			printlnLog( "total  : ", round( meanPosNegTotal[it][dimRep]; digits=3 ), " +- ", round( stdPosNegTotal[it][dimRep]; digits=3 )  );
		end
		# @infiltrate
		
		# meanPosNegNarr = cat( meanPosNegN...; dims = length(size(meanPosNegN[1]))+1 );
		# @info( size(meanPosNegN) );
		# @info(length(size(meanPosNegN[1]))+1);
		# @info( meanPosNegN );
		# meanPosNegNarr = permutedims( meanPosNegNarr, [2,1] )
		
		save( varFileName, "N", N, "seed", seedFedstart, "posNlst", posNlst, "negNlst", negNlst, "posLocLst", posLocLst, "negLocLst", negLocLst, "H_GOElst_All", H_GOElst_All, "meanPosNegN", meanPosNegN );
		
		gLstNameMain = "gLst_GOE3dfull";
		# if degObj.dim == 2
			# gLstNameMain = string( gLstNameMain, "_dim_", dim );
		gNpyNameMain = string( gLstNameMain, "_", fileNameAttr, npyType );
		
		# npzwrite( gNpyNameMain ,meanPosNegNarr );
		# @infiltrate
	end
end

function divB_profile_GOE_3d_oneDir( Nlst=[10], param_divide_num=80, seedFedstart = -1; scale=1, ratio=nothing, fileNameMod = "", maxDepth = 0, param_divide_next = nothing, isOrtho = false )
	param_dim = 2;
	if ratio == nothing
		ratio = ones(Int64,param_dim);
	end
	if maxDepth > 0 && param_divide_next == nothing
		param_divide_next = 2 .* ones(Int64, param_dim);
		param_divide_next[1] = 1;
	end
	
	locator_func_ratio = locator_div_GOE_ratio;
	
	param_min = fill(0.0,param_dim);
	param_max = fill(2*pi,param_dim);
	
	param_divide = fill( param_divide_num, param_dim );
	
	param_min_3 = 0;
	param_max_3 = 2*pi;
	
	seedFed = seedFedstart;
	if seedFed > 0
		Random.seed!(seedFed);
	end
	for itN = 1 : length(Nlst)
		N = Nlst[itN];
		posNlst = zeros( Int64, param_divide_num, N );
		negNlst = zeros( Int64, param_divide_num, N );
		posLocLst = Vector{ Vector{ Array{ Int64,2 } } }(undef, param_divide_num);
		negLocLst = Vector{ Vector{ Array{ Int64,2 } } }(undef, param_divide_num);
		posTotalLst = zeros( Int64, param_divide_num );
		negTotalLst = zeros( Int64, param_divide_num );
		
		degObj = DegObj( N, param_dim, param_divide, param_min, param_max );
		
		param_grid3 = degObj.param_grids[1];
		
		H_GOElst = HlstFunc( H_GOE, degObj.param_dim, degObj.N, isOrtho );
		dim3 = 1;
		H_GOElst_3 = HlstFunc( H_GOE, dim3, degObj.N, isOrtho );
		
		# H_lst_full = [ [ zeros(Complex{Float64},N,N) for iCosSin = 1:2, iDim = 1 : param_dim ] for it = 1:param_divide_num ];
		H3 = zeros( Complex{Float64}, N, N );
		for it = 1:param_divide_num
			param3 = param_grid3[it];
			Hmat_3comb_ratio!(H3, [param3],H_GOElst_3,ratio);
			HmatFun = (H,xlst) -> Hmat_3comb_ratio_Hoffset!(H,xlst,H_GOElst, ratio,H3);
			# @infiltrate	
			posNlst[it,:], negNlst[it,:], posLocLst[it], negLocLst[it] = locator_div( degObj, HmatFun );
			posTotalLst[it] = sum(posNlst[it,:]);
			negTotalLst[it] = sum(negNlst[it,:]);
			# @info( "GC: " );
			# @time GC.gc();
		end	
		num_it = 1;
		fileNameAttr = fileNameAttrFunc( N, param_divide, num_it, seedFedstart );
		fileNameMain = string( "deg_GOE_3d", "_dim_", param_dim );
		if scale != 1'
			fileNameAttr = string( fileNameAttr, "_scale_", scale );
		end
		if ratio != ones( Int64, param_dim )
			fileNameAttr = string( fileNameAttr, "_ratio_", ratio );
		end
		if fileNameMod != ""
			fileNameAttr = string( fileNameAttr, "_", fileNameMod );
		end
		oFileName = string( fileNameMain, "_", fileNameAttr );
		varFileName = string( oFileName, jldType );
		
		logFileName = string( "log_", oFileName, ".log" );
		logIO = open(logFileName, "w");
		
		printlnLog = (msg...)->printlnLogFile( logFileName, msg... );
		
		posnegN = [ posNlst, negNlst ];
		posnegTotal = [ posTotalLst, negTotalLst ];
		meanPosNegN = [ [mean(col) for col in eachcol(Nlst)] for Nlst in posnegN ];
		stdPosNegN = [ [std(col) for col in eachcol(Nlst)] for Nlst in posnegN ];
		meanPosNegTotal = [ mean(totalLst) for totalLst in posnegTotal ];
		stdPosNegTotal = [ std(totalLst) for totalLst in posnegTotal ];
		
		printlnLog( "matrix dimension: ", N );
		NnameLst = ["N+", "N-"];
		for it = 1:2
			printlnLog( NnameLst[it], ": " );
			for n = 1:N
				printlnLog( "level ", n, ": ", round( meanPosNegN[it][n]; digits=3 ), " +- ", round( stdPosNegN[it][n]; digits=3 ) );
			end
			printlnLog( "total  : ", round( meanPosNegTotal[it]; digits=3 ), " +- ", round( stdPosNegTotal[it]; digits=3 )  );
		end
		
		meanPosNegNarr = cat( meanPosNegN...; dims = length(size(meanPosNegN[1]))+1 );
		@info( size(meanPosNegN) );
		@info(length(size(meanPosNegN[1]))+1);
		@info( meanPosNegN );
		meanPosNegNarr = permutedims( meanPosNegNarr, [2,1] )
		
		save( varFileName, "N", N, "seed", seedFedstart, "posNlst", posNlst, "negNlst", negNlst, "posLocLst", posLocLst, "negLocLst", negLocLst, "H_GOElst", H_GOElst, "H_GOElst_3", H_GOElst_3, "meanPosNegN", meanPosNegN );
		
		gLstNameMain = "gLst_GOE3d";
		# if degObj.dim == 2
			# gLstNameMain = string( gLstNameMain, "_dim_", dim );
		gNpyNameMain = string( gLstNameMain, "_", fileNameAttr, npyType );
		
		npzwrite( gNpyNameMain ,meanPosNegNarr );
		# @infiltrate
	end
end
