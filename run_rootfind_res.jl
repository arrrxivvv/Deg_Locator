using DegLocatorDiv

resLst = [5,10,15,20];
nDim = 3;
# mSz = 5;
mSzLst = [5,10,15,20];
itNum = 100;
thresLst = [1e-10, 1e-7, 1e-8];
fModBase = "rootFind";

for thres in thresLst
	for res in resLst
		for mSz in mSzLst
			thresEDeg = 1e2 * thres;
			thresDegCollision = 1e2 * thres;
			@time divB_profile( nDim, [mSz], itNum, fill(res,nDim), 1000, rootFind; fileNameMod = "rootFind", thresNM = thres, thresEDeg = thresEDeg, thresDegCollision = thresDegCollision );
		end
	end
end
