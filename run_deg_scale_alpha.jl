using DegLocatorDiv

nLst = [10:1:60;];
alphaLst = [0.2:0.2:1;];

nMin = 10;
nMax = 50;
dNres = nMax - nMin;

resMin = 16;
resMax = 32;
dRes = resMax - resMin;

dim = 3;
itNum = 200;
seed = 1000;

param_divide = [80,16,16];
fNameModBase = "alphaTest";

nFinished = 55;
iAstart = 1;
iAend = length(alphaLst);
iAend = 2;
iAnow = 1;

for iAl = iAstart : iAend
	alpha = alphaLst[iAl];
	for n in nLst
		if iAl == iAnow && n <= nFinished
			continue;
		end
		res = 16 + Int64(floor( ( n - nMin ) / dNres * dRes ));
		println(res)
		param_divide .= [80,res,res];
		
		@time divB_profile( 3, [n], itNum, param_divide, seed; alpha = alpha, fileNameMod = fNameModBase );
		fNameMod = fNameModBase;
		if alpha > 0
			fNameMod = "_alpha_" * string(alpha) * "_" * fNameModBase;
		else 
			fNameMod = "_" * fNameModBase;
		end
		
		@time distillLocsFromFile( itNum, n, param_divide, seed; fileNameMod = fNameMod );
		@time deltaN_avg_lst_fromFile( [n], itNum, param_divide, seed; opt = "distill2", dim=3, fileNameMod = fNameMod );
		GC.gc();
		@time locNvarFromFile( n, param_divide, itNum, seed; dim = dim, fMod = fNameModBase, alpha = alpha );
		@time deltaN_var_fromFileDirect( n, param_divide, itNum, seed; dim = dim, fMod = fNameModBase, alpha = alpha, opt = "distill2" );
		
		# @time divB_profile( 3, [n],1000, [80,res,res], 1000 );
		# @time distillLocsFromFile( avgNum, n, param_divide, seed );
		# @time deltaN_avg_lst_fromFile( [n], avgNum, param_divide, seed; dim = 3, pow=2, opt="distill" );
		# @time deltaN_avg_lst_fromFile( [n], 1000, [80,res,res], 1000; dim = 3, fileNameMod = "" );
		# @time deltaNCum_avg_lst_fromFile( [n], 1000, [80,res,res], 1000; dim = 3, fileNameMod = "" );
		# @time fFunc_stat_from_file(n, [80,res,res], 1000, 1000, "", 3, "");
		# @time deltaN_var_lst_fromFile( [n], 1000, [80,res,res], 1000; dim = 3, fileNameMod = "", opt="distill" );
		# @time fFunc_stat_from_file(n, [80,res,res], 1000, 1000, "", 3, isDistilled = true);
		# @time fFunc_stat_across_from_file(n, [80,res,res], 1000, 1000, "", 3, isDistilled = true);
	end
end
