using DegLocatorDiv

nLst = [10];
alphaLst = [0.2:0.2:1;];

nMin = 10;
nMax = 50;
dNres = nMax - nMin;

resMin = 80;
resMax = 80;
dRes = resMax - resMin;

dim = 3;
itNum = 500;
seed = 1000;

param_divide = fill(resMin,dim);
fNameModBase = "";

nFinished = 55;
iAstart = 1;
iAend = length(alphaLst);
# iAend = 2;
iAnow = 1;

for iAl = iAstart : iAend
	alpha = alphaLst[iAl];
	for n in nLst
		if iAl == iAnow && n <= nFinished
			continue;
		end
		res = resMin + Int64(floor( ( n - nMin ) / dNres * dRes ));
		println(res)
		param_divide .= fill(res,dim);
		
		# @time divB_profile_GOE_3d( [n], res, itNum, seed; alpha = alpha, fileNameMod = fNameModBase );
		# fNameMod = fNameModBase;
		@time parityGOE_resave_fromFile( n, param_divide, itNum, seed; fMod = fNameModBase, dim = dim, alpha = alpha );
		@time parityAvg_fromFile( n, param_divide, itNum, seed; fMod = fNameModBase, dim = dim, alpha = alpha );
		
		GC.gc();
	end
end
