using NPZ
using Utils
using JLD2
using Statistics
using LoggingExtras

function deltaN_avg( avgNum, N, param_divide_num )
	deltaNsqSum = zeros( N, param_divide_num );
	for it in 1:avgNum
		posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div_GUE( N, param_divide_num; Hmat_seed = -1 );	
		
		deltaNsqSum += deltaN( N, param_divide_num, posLocLst, negLocLst );
	end
	deltaNsqAvg = deltaNsqSum ./ avgNum;
	return deltaNsqAvg;
end

function deltaN_avg_fromFile( avgNum, N, param_divide, seed; fileNameMod = "", dim = nothing, pow = 2, opt = nothing )
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed; dim = dim );
	fileNameAttr = string( fileNameAttr, fileNameMod );
	optDistill = "distill";
	optDistill2 = "distill2";
	if opt == optDistill
		fileNameMain = "locDistilled";
		levelLn = N-1;
	elseif opt == optDistill2
		fileNameMain = "deglocDistilled";
		levelLn = N-1;
	else 
		fileNameMain = "deg";
		levelLn = N;
	end
	degFileName = string( fileNameMain, "_", fileNameAttr );
	varFileName = string( degFileName, jld2Type );
	@info(varFileName);
	# @infiltrate
	if opt == optDistill || opt == optDistill2
		locLstPol = load( varFileName, "locDistilledLst" );
	else
		posLocLst = load(varFileName, "posLocLst");
		negLocLst = load(varFileName, "negLocLst");
		locLstPol = [posLocLst, negLocLst];
	end
	
	param_divide = Int64.(param_divide);
	if param_divide isa Array{Int64}
		param_divide_num = param_divide[1];
	else
		param_divide_num = param_divide;
	end
	
	deltaNsqSum = zeros( levelLn, param_divide_num );
	# @infiltrate
	for it in 1:avgNum
		deltaNsqSum += deltaNpower( levelLn, param_divide_num, locLstPol[1][it], locLstPol[2][it], pow );
		# @infiltrate
	end
	deltaNsqAvg = deltaNsqSum ./ avgNum;
	
	return deltaNsqAvg;
end

function deltaNCum_avg_fromFile( avgNum, N, param_divide, seed; fileNameMod = "", dim = nothing )
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed; dim = dim );
	fileNameAttr = string( fileNameAttr, fileNameMod );
	degFileName = string( "deg", "_", fileNameAttr );
	varFileName = string( degFileName, jldType );
	@info(varFileName);
	posLocLst = load(varFileName, "posLocLst");
	negLocLst = load(varFileName, "negLocLst");
	
	param_divide = Int64.(param_divide);
	if param_divide isa Array{Int64}
		param_divide_num = param_divide[1];
	else
		param_divide_num = param_divide;
	end
	
	# deltaN	
	deltaNCumSqSum = zeros( N, param_divide_num );
	deltaNCumOnce = zeros( 1,param_divide_num );
	pow = 1;
	for it in 1:avgNum
		deltaNonce = deltaNpower( N, param_divide_num, posLocLst[it], negLocLst[it], pow );
		deltaNCumOnce = deltaNonce[1,:];
		for n = 1 : N
			deltaNCumSqSum[n,:] += deltaNCumOnce.^2;
			if n<N
				deltaNCumOnce = deltaNCumOnce + deltaNonce[n+1,:];
			end
		end
		# println(sum( deltaNonce; dims=1 ));
	end
	deltaNCumSqAvg = deltaNCumSqSum ./ avgNum;
	
	return deltaNCumSqAvg;
end

function deltaN_var_fromFile( avgNum, N, param_divide, seed; fileNameMod = "", dim = nothing, pow = 2, opt = nothing )
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed; dim = dim );
	fileNameAttr = string( fileNameAttr, fileNameMod );
	optDistill = "distill";
	optDistill2 = "distill2";
	if opt == optDistill
		fileNameMain = "locDistilled";
		levelLn = N-1;
	elseif opt == optDistill2
		fileNameMain = "deglocDistilled";
		levelLn = N-1;
	else 
		fileNameMain = "deg";
		levelLn = N;
	end
	degFileName = string( fileNameMain, "_", fileNameAttr );
	varFileName = string( degFileName, jldType );
	@info(varFileName);
	if opt == optDistill || opt == optDistilled2
		locLstPol = load( varFileName, "locDistilledLst" );
	else
		posLocLst = load(varFileName, "posLocLst");
		negLocLst = load(varFileName, "negLocLst");
		locLstPol = [posLocLst, negLocLst];
	end
	
	param_divide = Int64.(param_divide);
	if param_divide isa Array{Int64}
		param_divide_num = param_divide[1];
	else
		param_divide_num = param_divide;
	end
	
	deltaNsqLst = zeros( levelLn, param_divide_num, avgNum );
	for it in 1:avgNum
		deltaNsqLst[:,:,it] = deltaNpower( levelLn, param_divide_num, locLstPol[1][it], locLstPol[2][it], pow );
	end
	# deltaNsqAvg = deltaNsqSum ./ avgNum;
	deltaNsqVar = std( deltaNsqLst, dims=3 );
	
	return deltaNsqVar;
end

function varG_fromFile( avgNum, N, param_divide, seed; fileNameMod = "", dim = nothing, pow = 2 )
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed; dim = dim );
	fileNameAttr = string( fileNameAttr, fileNameMod );
	degFileName = string( "deg", "_", fileNameAttr );
	varFileName = string( degFileName, jldType );
	@info(varFileName);
	posNLst = load(varFileName, "posNlst");
	negNLst = load(varFileName, "negNlst");
	
	GstdLst = std( 2*posNLst, dims=1 );
	
	return GstdLst;
end

function deltaN_avg_lst( Nlst = [10], avgNum = 100, param_divide_num = 26 )
	filenameMain = "deltaN";

	for N in Nlst
		deltaNsq = deltaN_avg( avgNum, N, param_divide_num );
		# filenameAttr = string( "_N_", N, "_avgNum_", avgNum, "param_divide", param_divide_num );
		filenameAttr = fileNameAttrFunc( N, param_divide, avgNum );
		filenameNpy = string( filenameMain, filenameAttr, npyType );
		npzwrite( filenameNpy, deltaNsq );
	end
end


function deltaN_avg_lst_fromFile( Nlst = [10], avgNum = 100, param_divide = 26, seed = -1; dim = nothing, fileNameMod = "", fileNameOutMod = "", pow = 2, opt=nothing )
	if opt == "distill" || opt == "distill2"
		filenameMain = "deltaN_distill";
	else 
		filenameMain = "deltaN";
	end
	for N in Nlst
		deltaNsq = deltaN_avg_fromFile( avgNum, N, param_divide, seed; fileNameMod = fileNameMod, dim = dim, pow = pow, opt=opt );
		filenameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed );
		filenameAttr = string( filenameAttr, fileNameMod );
		if pow != 2
			filenameAttr = string( "pow_", pow, "_", filenameAttr )
		end
		# filenameAttr = string( "_N_", N, "_avgNum_", avgNum, "param_divide", param_divide_num );
		if fileNameOutMod != ""
			filenameAttr = string( filenameAttr, "_", fileNameOutMod);
		end
		filenameNpy = string( filenameMain, "_",  filenameAttr, npyType );
		npzwrite( filenameNpy, deltaNsq );
	end
end

function deltaNCum_avg_lst_fromFile( Nlst = [10], avgNum = 100, param_divide = 26, seed = -1; dim = nothing, fileNameMod = "", fileNameOutMod = "" )
	filenameMain = "deltaNCum";
	for N in Nlst
		deltaNsq = deltaNCum_avg_fromFile( avgNum, N, param_divide, seed; fileNameMod = fileNameMod, dim = dim );
		filenameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed );
		filenameAttr = string( filenameAttr, fileNameMod );	
		# filenameAttr = string( "_N_", N, "_avgNum_", avgNum, "param_divide", param_divide_num );
		if fileNameOutMod != ""
			filenameAttr = string( filenameAttr, "_", fileNameOutMod);
		end
		filenameNpy = string( filenameMain, "_",  filenameAttr, npyType );
		npzwrite( filenameNpy, deltaNsq );
	end
end

function varG_lst_fromFile( Nlst = [10], avgNum = 100, param_divide = 26, seed = -1; dim = nothing, fileNameMod = "", fileNameOutMod = "", pow = 2 )
	filenameMain = "varG";
	for N in Nlst
		GstdLst = varG_fromFile( avgNum, N, param_divide, seed; fileNameMod = fileNameMod, dim = dim, pow = pow );
		filenameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed );
		filenameAttr = string( filenameAttr, fileNameMod );
		if pow != 2
			filenameAttr = string( "pow_", pow, "_", filenameAttr )
		end
		# filenameAttr = string( "_N_", N, "_avgNum_", avgNum, "param_divide", param_divide_num );
		if fileNameOutMod != ""
			filenameAttr = string( filenameAttr, "_", fileNameOutMod);
		end
		filenameNpy = string( filenameMain, "_",  filenameAttr, npyType );
		npzwrite( filenameNpy, GstdLst );
	end
end

function deltaN_var_lst_fromFile( Nlst = [10], avgNum = 100, param_divide = 26, seed = -1; dim = nothing, fileNameMod = "", fileNameOutMod = "", pow = 2, opt=nothing )
	filenameMain = "deltaVar";
	if opt == "distill" || opt == "distill2"
		filenameMain = string( filenameMain, "_distill" );
	end
	for N in Nlst
		GstdLst = deltaN_var_fromFile( avgNum, N, param_divide, seed; fileNameMod = fileNameMod, dim = dim, pow = pow, opt = opt );
		filenameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed );
		filenameAttr = string( filenameAttr, fileNameMod );
		if pow != 2
			filenameAttr = string( "pow_", pow, "_", filenameAttr )
		end
		# filenameAttr = string( "_N_", N, "_avgNum_", avgNum, "param_divide", param_divide_num );
		if fileNameOutMod != ""
			filenameAttr = string( filenameAttr, "_", fileNameOutMod);
		end
		filenameNpy = string( filenameMain, "_",  filenameAttr, npyType );
		npzwrite( filenameNpy, GstdLst );
	end
end

function deltaN_var_fromFileDirect( N, param_divide, itNum, seed; fMod = "", dim = nothing, pow = 2, opt = nothing, alpha = nothing, fExt = jld2Type )
	attrLst, valLst = fNameAttrLstFunc( N, param_divide, itNum, seed; dim = dim, alpha = alpha );
	optDistill = "distill";
	optDistill2 = "distill2";
	if opt == optDistill
		fileNameMain = "locDistilled";
		levelLn = N-1;
	elseif opt == optDistill2
		fileNameMain = "deglocDistilled";
		levelLn = N-1;
	else 
		fileNameMain = "deg";
		levelLn = N;
	end
	varFileName = fNameFunc( fileNameMain, attrLst, valLst, fExt; fMod = fMod );
	@info(varFileName);
	fGStdMain = "deltaVar";
	if opt == optDistill || opt == optDistill2
		locLstPol = load( varFileName, "locDistilledLst" );
		fGStdMain = fGStdMain * "_distill";
	else
		posLocLst = load(varFileName, "posLocLst");
		negLocLst = load(varFileName, "negLocLst");
		locLstPol = [posLocLst, negLocLst];
	end
	
	param_divide = Int64.(param_divide);
	if param_divide isa Array{Int64}
		param_divide_num = param_divide[1];
	else
		param_divide_num = param_divide;
	end
	
	deltaNsqLst = zeros( levelLn, param_divide_num, itNum );
	for it in 1:itNum
		deltaNsqLst[:,:,it] = deltaNpower( levelLn, param_divide_num, locLstPol[1][it], locLstPol[2][it], pow );
	end
	# deltaNsqAvg = deltaNsqSum ./ avgNum;
	deltaNsqVar = std( deltaNsqLst, dims=3 );
	
	# @infiltrate
	fGStd = fNameFunc( fGStdMain, attrLst, valLst, npyType; fMod = fMod );
	filenameNpy = npzwrite( fGStd, deltaNsqVar );
	
	return deltaNsqVar;
end
