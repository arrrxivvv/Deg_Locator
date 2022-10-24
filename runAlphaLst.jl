using Utils
alphaLst = [0:0.2:1;];
N = 10;
itNum = 5000;
seed = 1000;
param_divide = [80, 16, 16];
fNameModBase = "alphaTest";

for iAl = 1 : length(alphaLst)
	alpha = alphaLst[iAl];
	println( "alpha = " * string(alpha) );
	
	print("\rdivB\r");
	@info( "divB" )
	tFull = @timed divB_profile( 3, [N], itNum, param_divide, seed; alpha = alpha, fileNameMod = fNameModBase );
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	fNameMod = fNameModBase;
	if alpha > 0
		fNameMod = "_alpha_" * string(alpha) * "_" * fNameModBase;
	else 
		fNameMod = "_" * fNameModBase;
	end
	
	print("\rdistillLocs\r");
	@info( "distillLocs" )
	tFull = @timed distillLocsFromFile( itNum, N, param_divide, seed; fileNameMod = fNameMod );
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	print("\rdeltaN\r");
	@info( "deltaN" )
	tFull = @timed deltaN_avg_lst_fromFile( [N], itNum, param_divide, seed; opt = "distill2", dim=3, fileNameMod = fNameMod );
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	GC.gc();
end
