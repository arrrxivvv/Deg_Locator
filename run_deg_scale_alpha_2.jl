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

nFinished = 1;
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
		
		# @time distillLocsFromFile( itNum, n, param_divide, seed; fileNameMod = fNameMod );
		# @time deltaN_avg_lst_fromFile( [n], itNum, param_divide, seed; opt = "distill2", dim=3, fileNameMod = fNameMod );
		# GC.gc();
		@time locNvarFromFile( n, param_divide, itNum, seed; dim = dim, fMod = fNameModBase, alpha = alpha );
		@time deltaN_var_fromFileDirect( n, param_divide, itNum, seed; dim = dim, fMod = fNameModBase, alpha = alpha, opt = "distill2" );
	end
end
