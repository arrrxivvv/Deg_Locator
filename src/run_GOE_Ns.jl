itNum = 100;
resLow = 64;
resHigh = 64;
dRes = resHigh - resLow;
nDim = 3;

seed = 1000;

nMin = 10;
nMax = 50;
dN = nMax - nMin;

# nLst = [35,40,45,50];
nLst = [10,15,20,25,30];

fMod = "memLayers";
fMod = "";

itNumGc = 2;
nGc = 2;
resGc = 2;

for n in nLst
	println( "N: " * string(n) )
	res = Int64( floor( (n - nMin) * dRes / dN ) ) + resLow;
	# param_divide = fill(res,3);
	divLst = fill(res,nDim);
	# @time divB_profile_GOE_3d( [n], res, itNum, seed; isSaveMem = true );
	# GC.gc();
	# @time divB_profile_GOE_3d( [nGc], resGc, itNumGc, seed );
	# @time parityGOE_resave_fromFile( n, divLst, itNum, seed; dim = nDim, fMod = fMod );
	# @time parity_corr_GOE_arr_from_file( n, divLst, itNum, seed; dim = nDim, fMod = fMod );
	# @time divB_profile_GOE_3d( [4], 2, 1, seed; isSaveMem = true );
	GC.gc();
	# @time parityAvg_fromFile( n, divLst, itNum, seed; dim = 3, fMod = fMod );
	@time parityCorr_GOE_AvgStd_fromFile( n, divLst, itNum, seed; dim = 3 );
end

# @time divB_profile_GOE_3d( [4], 2, 1, seed; isSaveMem = true );
# GC.gc();

# itNum = 500;
# resLow = 80;
# resHigh = 80;
# dRes = resHigh - resLow;
# nDim = 3;

# seed = 1000;

# nMin = 10;
# nMax = 50;
# dN = nMax - nMin;

# # nLst = [35,40,45,50];
# nLst = [10,15,20,25,30];

# itNumGc = 2;
# nGc = 2;
# resGc = 2;

# for n in nLst
	# println( "N: " * string(n) )
	# res = Int64( floor( (n - nMin) * dRes / dN ) ) + resLow;
	# divLst = fill(res,nDim);
	# # @time divB_profile_GOE_3d( [n], res, itNum, seed; isSaveMem = true );
	# # GC.gc();
	# # @time divB_profile_GOE_3d( [nGc], resGc, itNumGc, seed );
	# @time parityGOE_resave_fromFile( n, divLst, itNum, seed; dim = nDim, fMod = fMod );
	# @time parity_corr_GOE_arr_from_file( n, divLst, itNum, seed; dim = nDim, fMod = fMod );
	# # @time divB_profile_GOE_3d( [4], 2, 1, seed; isSaveMem = true );
	# GC.gc();
# end
