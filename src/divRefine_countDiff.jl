itNum = 10;
nRes = 3;
cntDiffPol = zeros(Int64,2, nRes);
cntPol = zeros(Int64,2, nRes);
thrDivB = 1e-8;

degCnt = 0;

mSz = 10;
for iPol = 1 : 2
	sgn = (-1)^(iPol-1);
	for it = 1 : itNum, n = 1 : mSz
		for iLoc = 1:size(divBLstSeq3[iPol][it][n],1), iRes = 1 : nRes
			global degCnt += 1;
			if sgn * real(divBLstSeq3[iPol][it][n][iLoc,iRes]) < thrDivB
				global cntDiffPol[iPol,iRes] += 1;
			end
			if (sgn * real(divBLstSeq3[iPol][it][n][iLoc,iRes]) - 1) > -thrDivB
				global cntPol[iPol,iRes] += 1;
			end
		end
	end
end
