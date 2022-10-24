using Logging; using LoggingExtras;
using NPZ
include("divB.jl")
include("structure_factors.jl")

function struct_abs_avg( avgNum, N, param_divide_num )
	param_dim=3;
	dimLst = Int.( push!( param_divide_num * ones( param_dim ), N ) );
	param_step = 2*pi/param_divide_num;
	structFactAbsTotal = zeros( dimLst... );
	structFactPosAbsTotal = zeros( dimLst... );
	structFactNegAbsTotal = zeros( dimLst... );
	for it = 1:avgNum
		posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div_GUE( N, param_divide_num; Hmat_seed = -1 );
		structFact = structure_factor_true( N, param_dim, param_divide_num, param_step, posLocLst, negLocLst );
		structFactPos =  structure_factor_pol( N, param_dim, param_divide_num, param_step, posLocLst, +1 );
		structFactNeg =  structure_factor_pol( N, param_dim, param_divide_num, param_step, negLocLst, -1 );
		structFactAbsTotal += (abs.( structFact )).^2;
		structFactPosAbsTotal += (abs.( structFactPos )).^2;
		structFactNegAbsTotal += (abs.( structFactNeg )).^2;
		@info( posCount, " ", negCount );
	end
	
	return ( structFactAbsTotal / avgNum, structFactPosAbsTotal / avgNum, structFactNegAbsTotal / avgNum );
end

function struct_abs_avg_lst(  )
	Nlst = [3,5,10,15];
	# Nlst = [5];
	avgNum = 100;
	param_divide_num = 50;
	filenameMain = "structAbsAvg";
	filenameAttrPos = "_pos";
	filenameAttrNeg = "_neg";
	for N in Nlst
		structAbsAvg, structPosAvg, structNegAvg = struct_abs_avg( avgNum, N, param_divide_num );
		filenameAttr = string( "_N_", N, "_avgNum_", avgNum,  npyType );
		filenameNpy = string( filenameMain,filenameAttr );
		filenamePos = string( filenameMain, filenameAttrPos, filenameAttr );
		filenameNeg = string( filenameMain, filenameAttrNeg, filenameAttr );
		npzwrite( filenameNpy, structAbsAvg );
		npzwrite( filenamePos, structPosAvg );
		npzwrite( filenameNeg, structNegAvg );
	end
end
