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

function divB_profile( param_dim = 2, Nlst=[20], num_it_main=2, param_divide=50, seedFedstart = -1; scale=1, ratio=nothing, fileNameMod = "" )
	# param_dim = 3;
	if ratio == nothing
		ratio = Int64.( ones(param_dim) );
	end
	
	if param_dim == 2
		locator_func_ratio = locator_div_GOE_ratio;
	else
		locator_func_ratio = locator_div_GUE_ratio;
	end
	
	symFact = 4;
	indCartPosNeg = CartesianIndices( Tuple( symFact * ones(Int64,param_dim) ) );
	for itN = 1 : length(Nlst)
		N = Nlst[itN];
		num_it_sym = symFact ^ param_dim;
		num_it = num_it_main * num_it_sym;
		posNlst = zeros( Int64, num_it, N );
		negNlst = zeros( Int64, num_it, N );
		posLocLst = [ [] for it = 1:num_it ];
		negLocLst = [ [] for it = 1:num_it ];
		posTotalLst = zeros( Int64, num_it );
		negTotalLst = zeros( Int64, num_it );
		seedFed = seedFedstart;
		if seedFed > 0
			Random.seed!(seedFed);
		end
		H_GUE_lst = [ [ zeros(Complex{Float64},N,N) for iCosSin = 1:2, iDim = 1 : param_dim ] for it = 1:num_it ];
		# println( typeof( H_GUE_lst ) );
		for itMain = 0:num_it_main-1
			itStart = itMain * num_it_sym + 1;
			it = itStart;
			posNlst[it,:], negNlst[it,:], posLocLst[it], negLocLst[it], H_GUE_lst[it] = locator_func_ratio( N, param_divide; ratio = ratio, Hmat_seed = seedFed ); 
			# display(ratio);
			posTotalLst[it] = sum(posNlst[it,:]);
			negTotalLst[it] = sum(negNlst[it,:]);
			for itSym in indCartPosNeg
				if it == itStart
					it += 1;
					continue;
				end
				posLocLstTmp = deepcopy( posLocLst[itStart] );
				negLocLstTmp = deepcopy( negLocLst[itStart] );
				for dim = 1 : param_dim
						if param_divide isa Array{Int64}
							param_divide_num = param_divide[dim];
						else
							param_divide_num = param_divide;
						end
					itSymDim = itSym[dim];
					posLocLstTmp = shLocLst( posLocLstTmp, param_divide_num, dim, itSymDim );
					negLocLstTmp = shLocLst( negLocLstTmp, param_divide_num, dim, itSymDim );
				end
				posNlst[it,:] = posNlst[itStart,:];
				negNlst[it,:] = negNlst[itStart,:];
				posLocLst[it] = posLocLstTmp;
				negLocLst[it] = negLocLstTmp;
				H_GUE_lst[it] = H_GUE_lst[itStart];
				posTotalLst[it] = sum(posNlst[it,:]);
				negTotalLst[it] = sum(negNlst[it,:]);
				it+= 1;
			end
			# if seedFedstart != -1
				# seedFed += 1;
			# end
			GC.gc();
		end		
		fileNameAttr = fileNameAttrFunc( N, param_divide, num_it, seedFedstart );
		fileNameMain = string( "deg", "_dim_", param_dim );
		if scale != 1'
			fileNameAttr = string( fileNameAttr, "_scale_", scale );
		end
		if fileNameMod != ""
			fileNameAttr = string( fileNameAttr, "_", fileNameMod );
		end
		oFileName = string( fileNameMain, "_", fileNameAttr );
		varFileName = string( oFileName, jldType );
		save( varFileName, "N", N, "seed", seedFedstart, "posNlst", posNlst, "negNlst", negNlst, "posLocLst", posLocLst, "negLocLst", negLocLst, "H_GUE_lst", H_GUE_lst );
		
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
	end
end

