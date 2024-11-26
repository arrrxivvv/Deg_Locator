using RandomCircle
using FilenameManip
using JLD2
using SharedFNames
using DelimitedFiles
using Utils
using LsqFit
using Statistics

using Infiltrator

isRenewRand = false;

isFileNameOnly = true;

nCirc = 10;
nDim3 = 3;

itNum = 10;

itNum1Pass = 10;

nCircStep = 5;
rCircStep = 0.1;
nCircLst = [5:nCircStep:50;];
rCircLst = [0.1:rCircStep:0.5;];
lnNCirc = length( nCircLst );
lnRCirc = length( rCircLst );

rRaw = 1;

nSample = 1000;

divNum1Pass = 128;

fMain = "randCirc";
fMainParams = fMain * "Params";

attrLstBase = ["nCircLst", "rCircLst", "itNum1Pass", "itNum", "divNum1Pass"];
valLstBase = [ FilenameManip.fAttr_arrSummary( nCircLst, nCircStep ), FilenameManip.fAttr_arrSummary( rCircLst, rCircStep ), itNum1Pass, itNum, divNum1Pass ];
fNameRandCirc = fNameFunc( fMain, attrLstBase, valLstBase, jld2Type );
fNameRandCircParams = fNameFunc( fMainParams, attrLstBase, valLstBase, jld2Type );

fNameArr = [fNameRandCirc, fNameRandCircParams];

fMainFNameArr = fMain * "FNameArr";

fNameRandCirc = fNameFunc( fMainFNameArr, attrLstBase, valLstBase, jld2Type );
jldsave( fNameRandCirc; fNameArr );

open( SharedFNames.dirLog * SharedFNames.fNameTmpNameFileLst, "w" ) do io
	println( io, fNameRandCirc );
end



if !isFileNameOnly

	runData = RandomCircle.RunRandCircData( nCircLst, rCircLst, itNum1Pass, itNum, nSample, divNum1Pass );

	RandomCircle.runBaseInfo!( runData, rCircLst );
	RandomCircle.run1Pass!( runData );
	RandomCircle.runFine!( runData );


	corrLenLst, expScaleLst, expShLst, zakMean1dLst, zak1dLst, numCornerMeanLst, numCornerLst, idCornerLstLst = RandomCircle.exportDataDetailed( runData );

	jldsave( fNameRandCirc; corrLenLst, expScaleLst, expShLst, zakMean1dLst, zak1dLst, numCornerMeanLst, numCornerLst, idCornerLstLst );
	jldsave( fNameRandCircParams; nCircLst, rCircLst, divNum1Pass, itNum, itNum1Pass );
end
