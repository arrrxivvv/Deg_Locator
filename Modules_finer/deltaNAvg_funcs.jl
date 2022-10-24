using NPZ
using Utils

function deltaN_avg( avgNum, N, param_divide_num )
	deltaNsqSum = zeros( N, param_divide_num );
	for it in 1:avgNum
		posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div_GUE( N, param_divide_num; Hmat_seed = -1 );	
		
		deltaNsqSum += deltaN( N, param_divide_num, posLocLst, negLocLst );
	end
	deltaNsqAvg = deltaNsqSum ./ avgNum;
	return deltaNsqAvg;
end

function deltaN_avg_fromFile( avgNum, N, param_divide, seed )
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed );
	degFileName = string( "deg", "_", fileNameAttr );
	varFileName = string( degFileName, jldType );
	data = load(varFileName);
	
	param_divide = Int64.(param_divide);
	if param_divide isa Array{Int64}
		param_divide_num = param_divide[1];
	else
		param_divide_num = param_divide;
	end
	
	deltaNsqSum = zeros( N, param_divide_num );
	for it in 1:avgNum
		deltaNsqSum += deltaN( N, param_divide_num, data["posLocLst"][it], data["negLocLst"][it] );
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


function deltaN_avg_lst_fromFile( Nlst = [10], avgNum = 100, param_divide = 26, seed = -1; fileNameMod = "" )
	filenameMain = "deltaN";
	for N in Nlst
		deltaNsq = deltaN_avg_fromFile( avgNum, N, param_divide, seed );
		filenameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed );
		filenameAttr = string( filenameAttr, fileNameMod );
		# filenameAttr = string( "_N_", N, "_avgNum_", avgNum, "param_divide", param_divide_num );
		filenameNpy = string( filenameMain, "_",  filenameAttr, npyType );
		npzwrite( filenameNpy, deltaNsq );
	end
end
