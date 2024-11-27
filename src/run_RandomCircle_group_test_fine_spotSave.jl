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

isFileNameOnly = false;
# isFileNameOnly = true;

# isStoreCorrFull = false;
isStoreCorrFull = true;

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
fMainCorrFull = fMain * "CorrFull";

attrLstBase = ["nCircLst", "rCircLst", "itNum1Pass", "itNum", "divNum1Pass"];
valLstBase = [ FilenameManip.fAttr_arrSummary( nCircLst, nCircStep ), FilenameManip.fAttr_arrSummary( rCircLst, rCircStep ), itNum1Pass, itNum, divNum1Pass ];
fNameRandCirc = fNameFunc( fMain, attrLstBase, valLstBase, jld2Type );
fNameRandCircParams = fNameFunc( fMainParams, attrLstBase, valLstBase, jld2Type );
fNameCorrFull = fNameFunc( fMainCorrFull, attrLstBase, valLstBase, jld2Type );



fNameArr = [fNameRandCirc, fNameRandCircParams, fNameCorrFull];

fMainFNameArr = fMain * "FNameArr";
fMainFNameLst = fMain * "FNameLst";

fNameFNameArr = fNameFunc( fMainFNameArr, attrLstBase, valLstBase, jld2Type );
jldsave( fNameFNameArr; fNameArr );

fNameFNameLst = fNameFunc( fMainFNameLst, attrLstBase, valLstBase, txtType );
writedlm( fNameFNameLst, fNameArr );

open( SharedFNames.dirLog * SharedFNames.fNameFileLstJld2Lst, "w" ) do io
	println( io, fNameFNameArr );
end

open( SharedFNames.dirLog * SharedFNames.fNameFileLstLst, "w" ) do io
	println( io, fNameFNameLst );
end



if !isFileNameOnly

	runData = RandomCircle.RunRandCircData( nCircLst, rCircLst, itNum1Pass, itNum, nSample, divNum1Pass; isStoreCorrFull = isStoreCorrFull );

	RandomCircle.runBaseInfo!( runData, rCircLst );
	RandomCircle.run1Pass!( runData );
	RandomCircle.runFine!( runData );


	corrLenLst, expScaleLst, expShLst, zakCorrMean1dLst, zakCorr1dLst, numCornerMeanLst, numCornerLst, idCornerLstLst, zakArrAvgLst, zakArrAvgMeanLst = RandomCircle.exportDataDetailed( runData );
	nCircLst, rCircLst, itNum, divNum1Pass, divNumNxtLst = RandomCircle.exportParams( runData );

	jldsave( fNameRandCirc; corrLenLst, expScaleLst, expShLst, zakCorrMean1dLst, zakCorr1dLst, numCornerMeanLst, numCornerLst, idCornerLstLst, zakArrAvgLst, zakArrAvgMeanLst );
	jldsave( fNameRandCircParams; nCircLst, rCircLst, divNum1Pass, itNum, divNumNxtLst );
	
	if runData.isStoreCorrFull
		zakCorrMeanFullLst = RandomCircle.exportFullCorrMean( runData );
		jldsave( fNameCorrFull; zakCorrMeanFullLst );
	end
end
