using DegLocatorDiv

nDim = 3;
resLst = [5,10,15,20];
# resLst = [5]
# mSzLst = [5,10,15,20];
mSzLst = [10];
itNum = 100;
# itNum = 2;
thresLst = [1e-9];
thresValRatioLst = [10:10:100;];
# thresSzRatioLst = [10:10:100;];
thresSzRatioLst = [110:10:200;];
divLst = ones(Int64, nDim);
seed = 1000;

# println(pwd());

for thres in thresLst 
	for res in resLst
		divLst .= res;
		for mSz in mSzLst
			for thresValRatio in thresValRatioLst
			DegLocatorDiv.pick0gapLocsFromFile( mSz, divLst, itNum, seed; thresVal = thres, thresSz = thres, thresValRatio = thresValRatio );
				for thresSzRatio in thresSzRatioLst
				print("mSz = $mSz, res = $res, thres = $thres, thresValRatio = $thresValRatio, thresSzRatio = $thresSzRatio", "               \r")
				
				DegLocatorDiv.locLstCollisionRemoveFromFile( mSz, divLst, itNum, seed; thresVal = thres, thresSz = thres, thresValRatio = thresValRatio, thresSzRatio = thresSzRatio );
				end
			end
		end
	end
end
