using NPZ
using Utils
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

function deltaN_avg_fromFile( avgNum, N, param_divide, seed; fileNameMod = "", dim = nothing, pow = 2 )
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
	
	deltaNsqSum = zeros( N, param_divide_num );
	for it in 1:avgNum
		deltaNsqSum += deltaNpower( N, param_divide_num, posLocLst[it], negLocLst[it], pow );
	end
	deltaNsqAvg = deltaNsqSum ./ avgNum;
	
	return deltaNsqAvg;
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


function deltaN_avg_lst_fromFile( Nlst = [10], avgNum = 100, param_divide = 26, seed = -1; dim = nothing, fileNameMod = "", fileNameOutMod = "", pow = 2 )
	filenameMain = "deltaN";
	for N in Nlst
		deltaNsq = deltaN_avg_fromFile( avgNum, N, param_divide, seed; fileNameMod = fileNameMod, dim = dim, pow = pow );
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
