using DegLocatorDiv
using Utils

nLst = collect(10:1:60);

nMin = 10;
nMax = 50;
dNres = nMax - nMin;

resMin = 16;
resMax = 32;
dRes = resMax - resMin;

param_divide = [80,16,16];
fMod10 = "degObj";

itNum = 1000;
seed = 1000;

dim = 3;

for n in nLst
	res = 16 + Int64(floor( ( n - nMin ) / dNres * dRes ));
	println(res)
	if n > 10
		fMod = "";
	else
		fMod = fMod10;
	end
	param_divide .= [80,res,res];
	# @time divB_profile( 3, [n],1000, [80,res,res], 1000 );
	# @time distillLocsFromFile( avgNum, n, param_divide, seed );
	# @time deltaN_avg_lst_fromFile( [n], avgNum, param_divide, seed; dim = 3, pow=2, opt="distill" );
	# @time deltaN_avg_lst_fromFile( [n], 1000, [80,res,res], 1000; dim = 3, fileNameMod = "" );
	# @time deltaNCum_avg_lst_fromFile( [n], 1000, [80,res,res], 1000; dim = 3, fileNameMod = "" );
	# @time fFunc_stat_from_file(n, [80,res,res], 1000, 1000, "", 3, "");
	# @time deltaN_var_lst_fromFile( [n], 1000, [80,res,res], 1000; dim = 3, fileNameMod = "", opt="distill" );
	# @time fFunc_stat_from_file(n, [80,res,res], 1000, 1000, "", 3, isDistilled = true);
	# @time fFunc_stat_across_from_file(n, [80,res,res], 1000, 1000, "", 3, isDistilled = true);
	# @time locNvarFromFile( n, param_divide, itNum, seed; dim = dim, fMod = fMod, fExt = jldType, opt = "distill" );
	@time deltaN_var_fromFileDirect( n, param_divide, itNum, seed; dim = dim, fMod = fMod, opt = "distill", fExt = jldType );
	GC.gc();
end
