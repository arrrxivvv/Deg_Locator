res = 80;
N = 4;
nDim = 3;
divLst = fill(res,nDim);
seed = 1000;
itNum = 2000;
# @time divB_profile_GOE_3d( N, res, itNum, seed )
# implement logging the time to some log file
# ask better way to sync file to ssh server
# @time parityGOE_resave_fromFile( N, divLst, itNum, seed; dim = nDim );
# @time parityAvg_fromFile( N, divLst, itNum, seed; dim = nDim );

N = 30;
itNum = 500;
@time divB_profile_GOE_3d( N, res, itNum, seed; isSaveMem = true );

mSzLst = [35,40,45,50];
for mSz in mSzLst
	@time divB_profile_GOE_3d( mSz, res, itNum, seed; isSaveMem = true );
end
