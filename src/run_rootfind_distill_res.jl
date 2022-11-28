using DegLocatorDiv

# resLst = [5,10,15,20];
resLst = [5]
nDim = 3;
# mSzLst = [5,10,15,20];
mSzLst = [5];
# itNum = 100;
itNum = 2;
thresLst = [1e-9];
thresValRatioLst = [10:10:100;];
thresSzRatioLst = [10:10:100;];
divLst = ones(Int64, nDim);
seed = 1000;

# println(pwd());

for thres in thresLst 
	for res in resLst
		divLst .= res;
		for mSz in mSzLst, thresValRatio in thresValRatioLst, thresSzRatio in thresSzRatioLst
			DegLocatorDiv.pick0gapLocsFromFile( mSz, divLst, itNum, seed; thresVal = thres, thresSz = thres, thresValRatio = thresValRatio );
			DegLocatorDiv.locLstCollisionRemoveFromFile( mSz, divLst, itNum, seed; thresVal = thres, thresSz = thres, thresValRatio = thresValRatio, thresSzRatio = thresSzRatio );
		end
	end
end
