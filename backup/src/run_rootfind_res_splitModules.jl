using DegLocatorDiv

resLst = [5,10,15,20];
nDim = 3;
# mSz = 5;
mSzLst = [5,10,15,20];
itNum = 100;
thresLst = [1e-9, 1e-7, 1e-8];
divLst = ones(Int64, nDim);
seed = 1000;

for thres in thresLst
	for res in resLst
		divLst .= res;
		for mSz in mSzLst
			# @time divB_profile_rootFind( mSz, divLst, itNum, seed; thresVal = thres, thresSz = thres );
			@time locRootFindRawProfile( mSz, divLst, itNum, seed; thresVal = thres, thresSz = thres );
		end
	end
end
