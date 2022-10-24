using Utils
using Statistics
using JLD
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

function divB_profile( param_dim = 2, Nlst=[20], num_it_main=2, param_divide=50, seedFedstart = -1; scale=1, ratio=nothing, fileNameMod = "", maxDepth = 0, param_divide_next = nothing, isOrtho = false )
	if ratio == nothing
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
		N = Nlst[itN];
		# num_it_sym = symFact ^ param_dim;
		# num_it = num_it_main * num_it_sym;
		num_it = num_it_main;
		posNlst = zeros( Int64, num_it, N );
		negNlst = zeros( Int64, num_it, N );
		# posLocLst = [ [] for it = 1:num_it ];
		# negLocLst = [ [] for it = 1:num_it ];
		posLocLst = Vector{ Vector{ Array{ Int64,2 } } }(undef, num_it);
		negLocLst = Vector{ Vector{ Array{ Int64,2 } } }(undef, num_it);
		posTotalLst = zeros( Int64, num_it );
		negTotalLst = zeros( Int64, num_it );
		
		degObj = DegObj( N, param_dim, param_divide, param_min, param_max );
		
		H_GUE_lst = [ [ zeros(Complex{Float64},N,N) for iCosSin = 1:2, iDim = 1 : param_dim ] for it = 1:num_it ];
		for it = 1:num_it
			# itStart = itMain * num_it_sym + 1;
			# it = itStart;
			posNlst[it,:], negNlst[it,:], posLocLst[it], negLocLst[it], H_GUE_lst[it] = locator_func_ratio( degObj; ratio = ratio, Hmat_seed = seedFed, maxDepth = maxDepth, param_divide_next = param_divide_next, isOrtho );
			posTotalLst[it] = sum(posNlst[it,:]);
			negTotalLst[it] = sum(negNlst[it,:]);
			@info( "GC: " );
			@time GC.gc();
		end		
		fileNameAttr = fileNameAttrFunc( N, param_divide, num_it, seedFedstart );
		fileNameMain = string( "deg", "_dim_", param_dim );
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
		
		save( varFileName, "N", N, "seed", seedFedstart, "posNlst", posNlst, "negNlst", negNlst, "posLocLst", posLocLst, "negLocLst", negLocLst, "H_GUE_lst", H_GUE_lst, "meanPosNegN", meanPosNegN );
		
		gLstNameMain = "gLst";
		# if degObj.dim == 2
			# gLstNameMain = string( gLstNameMain, "_dim_", dim );
		gNpyNameMain = string( "gLst", "_", fileNameAttr, npyType );
		
		npzwrite( gNpyNameMain ,meanPosNegNarr );
		@infiltrate
	end
end

function divB_profile_GOE_3d( Nlst=[10], param_divide_num=80, seedFedstart = -1; scale=1, ratio=nothing, fileNameMod = "", maxDepth = 0, param_divide_next = nothing, isOrtho = false )
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
		posLocLst = Array{ Vector{ Array{ Int64,2 } },2 }(undef,param_dim_full, param_divide_num);
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
				@infiltrate	
				posNlst[id,it,:], negNlst[id,it,:], posLocLst[id,it], negLocLst[id,it] = locator_div( degObj, HmatFun );
				posTotalLst[id,it] = sum(posNlst[id,it,:]);
				negTotalLst[id,it] = sum(negNlst[id,it,:]);
				# @info( "GC: " );
				# @time GC.gc();
			end
		end
		num_it = 1;
		fileNameAttr = fileNameAttrFunc( N, param_divide, num_it, seedFedstart );
		fileNameMain = string( "deg_GOE_3d_full", "_dim_", param_dim_full );
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
		
		gLstNameMain = "gLst_GOE3dfull";
		# if degObj.dim == 2
			# gLstNameMain = string( gLstNameMain, "_dim_", dim );
		gNpyNameMain = string( gLstNameMain, "_", fileNameAttr, npyType );
		
		npzwrite( gNpyNameMain ,meanPosNegNarr );
		@infiltrate
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
			@infiltrate	
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
		@infiltrate
	end
end