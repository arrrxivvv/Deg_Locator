using RandomCircle
using FilenameManip
using JLD2
using SharedFNames
using DelimitedFiles


isRenewRand = false;

numRots = 10;
nDim3 = 3;

nSample = 1000;
rCirc = 0.2;

divNum = 128;

if isRenewRand
	ptLst = RandomCircle.genRandPt3dLst( numRots );
	rotMatLst = RandomCircle.genRotMatLst( numRots );
	RandomCircle.scaleRotMatLst!( rotMatLst, rCirc );
		
	axis1Lst = zeros(nDim3, numRots);

	for ii = 1 : numRots, iD = 1 : nDim3
		axis1Lst[iD,ii] = rotMatLst[ii][iD,1];
	end

	pt2dLst = RandomCircle.genRandPt2dFrom3dLst( ptLst );
	rotMat2dLst = RandomCircle.genRotMat2dLst( rotMatLst );
	rotMat2dInvLst = RandomCircle.genRotMat2dInvLst( rotMat2dLst );
end

ellipseSampleLst = RandomCircle.genEllipseSampleLst( rotMat2dLst, pt2dLst, nSample );

widthLst = RandomCircle.genEllipseWidthLst( rotMat2dLst, rotMat2dInvLst );
bndLst = RandomCircle.genEllipseBndFromWidthLst( widthLst, pt2dLst );
bndModLst, bndExtendedLst, idExtendedXYLst, pt2dModXYLst, pt2dExtendedXYLst = RandomCircle.genEllipseBndModExtendedLst( bndLst, pt2dLst );
idSortedBndModLst, idSortedBndExtendedLst = RandomCircle.sortBndLst( bndModLst, bndExtendedLst );

eqCoeffLst = RandomCircle.genEllipseEqCoeffLst( rotMat2dInvLst );

stepTest = 0.001;
xTestLst = [bndLst[1][1][1] + stepTest : stepTest : bndLst[1][1][2] - stepTest; ];
yTestLst = [bndLst[1][2][1] + stepTest : stepTest : bndLst[1][2][2] - stepTest; ];
xIntersectLst = [ RandomCircle.solveOtherXY( eqCoeffLst[1], pt2dLst[1], yVal; xySolved = 'X' ) for yVal in yTestLst ];
yIntersectLst = [ RandomCircle.solveOtherXY( eqCoeffLst[1], pt2dLst[1], xVal; xySolved = 'Y' ) for xVal in xTestLst ];

intersectAtXLst = zeros(2, 2*length(yTestLst));
intersectAtYLst = zeros(2, 2*length(xTestLst));
for ii = 1 : length(yTestLst)
	intersectAtXLst[2,ii] = yTestLst[ii];
	intersectAtXLst[2,ii+length(yTestLst)] = yTestLst[ii];
	intersectAtXLst[1,ii] = xIntersectLst[ii][1];
	intersectAtXLst[1,ii+length(yTestLst)] = xIntersectLst[ii][2];
end
for ii = 1 : length(xTestLst)
	intersectAtYLst[1,ii] = xTestLst[ii];
	intersectAtYLst[1,ii+length(xTestLst)] = xTestLst[ii];
	intersectAtYLst[2,ii] = yIntersectLst[ii][1];
	intersectAtYLst[2,ii+length(xTestLst)] = yIntersectLst[ii][2];
end

zakArr = RandomCircle.genZakArr( divNum, eqCoeffLst, pt2dLst, pt2dModXYLst, pt2dExtendedXYLst, bndModLst, bndExtendedLst, idExtendedXYLst, idSortedBndModLst, idSortedBndExtendedLst );

# fMainPtsRotMatsLst = 

fMainEllipseBndLst = "ellipseBnd";
attrLstBnd = [ "numCircs", "rCirc" ];
valLstBnd = Any[numRots, rCirc];

fNameEllipseBndLst = fNameFunc( fMainEllipseBndLst, attrLstBnd, valLstBnd, jld2Type );
save( fNameEllipseBndLst, "numCircs", numRots, "rCirc", rCirc, "bndLst", bndLst, "bndModLst", bndModLst, "bndExtendedLst", bndExtendedLst, "widthLst", widthLst, "ellipseSampleLst", ellipseSampleLst );

fMainSampleFromSolve = "ellipseSampleFromSolve";
fNameSampleFromSolve = fNameFunc( fMainSampleFromSolve, [], [], jld2Type );
save( fNameSampleFromSolve, "intersectAtXLst", intersectAtXLst, "intersectAtYLst", intersectAtYLst );

fMainAxisLst = "rotAxisLst";
fNameAxisLst = fNameFunc( fMainAxisLst, [],[], jld2Type );
save( fNameAxisLst, "axis1Lst", axis1Lst );

fMainZakArr = "zakArrRandCirc"
fNameZakArr = fNameFunc( fMainZakArr, attrLstBnd, valLstBnd, jld2Type );
save( fNameZakArr, "zakArr", zakArr );

fNameArr = [fNameAxisLst, fNameEllipseBndLst, fNameSampleFromSolve, fNameZakArr];

fLstMain = "fLst" * fMainAxisLst;
fLstAxisLstName = fNameFunc( fLstMain, [],[], jld2Type );

save( fLstAxisLstName, "fNameArr", fNameArr );

open( SharedFNames.dirLog * SharedFNames.fNameTmpNameFileLst, "w" ) do io
	println( io, fLstAxisLstName );
end
