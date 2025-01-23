using RandomCircle
using FilenameManip
using JLD2
using SharedFNames
using DelimitedFiles
using Utils
using LsqFit
using Statistics

# using Infiltrator

isRenewRand = false;

isFileNameOnly = false;
# isFileNameOnly = true;

isCornerDetect = false;
# isCornerDetect = true;

isSample = false;
# isSample = true;

# isStoreCorrFull = true;
isStoreCorrFull = false;

nCirc = 10;
nDim3 = 3;

itNum = 10;
# itNum = 100;

itNum1Pass = 10;

# nCircStep = 1;
# rCircStep = 0.025;
nCircStep = 5;
rCircStep = 0.1;
# rCircMax = 0.5;
# rCircMax = 1;
# nMax = 50;
# nMax = 100;
rCircMax = 0.7;
nMax = 70;
nCircLst = [5:nCircStep:nMax;]; 
rCircLst = [0.1:rCircStep:rCircMax;];
lnNCirc = length( nCircLst );
lnRCirc = length( rCircLst );

rRaw = 1;

nSample = 1000;

divNum1Pass = 128;

fMain = "randCirc";
fMainParams = fMain * "Params";
fMainCorrFull = fMain * "CorrFull";
fMainCorner = fMain * "Corner";

fMod = "";
if !isCornerDetect
	fMod = Utils.strAppendWith_( fMod, "noCorner" );
end
if !isStoreCorrFull
	fMod = Utils.strAppendWith_( fMod, "noFull" );
end

attrLstBase = ["nCircLst", "rCircLst", "itNum1Pass", "itNum", "divNum1Pass"];
valLstBase = [ FilenameManip.fAttr_arrSummary( nCircLst, nCircStep ), FilenameManip.fAttr_arrSummary( rCircLst, rCircStep ), itNum1Pass, itNum, divNum1Pass ];
fNameRandCirc = fNameFunc( fMain, attrLstBase, valLstBase, jld2Type; fMod = fMod );
fNameRandCircParams = fNameFunc( fMainParams, attrLstBase, valLstBase, jld2Type; fMod = fMod );
fNameCorrFull = fNameFunc( fMainCorrFull, attrLstBase, valLstBase, jld2Type; fMod = fMod );
fNameCorner = fNameFunc( fMainCorner, attrLstBase, valLstBase, jld2Type; fMod = fMod );



fNameArr = [fNameRandCirc, fNameRandCircParams, fNameCorner, fNameCorrFull];

fMainFNameArr = fMain * "FNameArr";
fMainFNameLst = fMain * "FNameLst";

fNameFNameArr = fNameFunc( fMainFNameArr, attrLstBase, valLstBase, jld2Type; fMod = fMod );
jldsave( fNameFNameArr; fNameArr );

fNameFNameLst = fNameFunc( fMainFNameLst, attrLstBase, valLstBase, txtType; fMod = fMod );
writedlm( fNameFNameLst, fNameArr );

open( SharedFNames.dirLog * SharedFNames.fNameFileLstJld2Lst, "w" ) do io
	println( io, fNameFNameArr );
end

open( SharedFNames.dirLog * SharedFNames.fNameFileLstLst, "w" ) do io
	println( io, fNameFNameLst );
end



if !isFileNameOnly

	runData = RandomCircle.RunRandCircData( nCircLst, rCircLst, itNum1Pass, itNum, nSample, divNum1Pass; isStoreCorrFull = isStoreCorrFull );

	@time RandomCircle.runBaseInfo!( runData, rCircLst; isSample = isSample );
	@time RandomCircle.run1Pass!( runData );
	@time RandomCircle.runFine!( runData; isCornerDetect = isCornerDetect );


	corrLenLst, expScaleLst, expShLst, zakCorrMean1dLst, zakCorr1dLst, zakArrAvgLst, zakArrAvgMeanLst = RandomCircle.exportDataDetailed( runData );
	nCircLst, rCircLst, itNum, divNum1Pass, divNumNxtLst = RandomCircle.exportParams( runData );

	jldsave( fNameRandCirc; corrLenLst, expScaleLst, expShLst, zakCorrMean1dLst, zakArrAvgLst, zakArrAvgMeanLst );
	# , zakCorr1dLst
	jldsave( fNameRandCircParams; nCircLst, rCircLst, divNum1Pass, itNum, divNumNxtLst );
	
	if runData.isStoreCorrFull
		zakCorrMeanFullLst = RandomCircle.exportFullCorrMean( runData );
		jldsave( fNameCorrFull; zakCorrMeanFullLst );
	end
	if isCornerDetect
		numCornerLst, numCornerMeanLst, idCornerLst = RandomCircle.exportDataCorner( runData );
		jldsave( fNameCorner; numCornerLst, numCornerMeanLst, idCornerLst );
	end
end
