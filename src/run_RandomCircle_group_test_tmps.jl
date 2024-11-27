using Statistics

using RandomCircle
using Utils
using SharedFNames
using JLD2
using FilenameManip	
using FFTW

fNameFNameArr = Utils.strReadLastLine( SharedFNames.dirLog * SharedFNames.fNameTmpNameFileLst );

fNameArr = load( fNameFNameArr, "fNameArr" );

# fNameRandCircZak, fNameRandCircCorr, fNameRandCircFine, fNameRandCircCorners = @view fNameArr[1:4];

# fNameRandCircFine = fNameArr[3];
fNameRandCirc = fNameArr[1];
fNameRandCircParams = fNameArr[2];
fNameRandCircCorrFull = fNameArr[3];

# zakCorrLst = load( fNameRandCircFine, "zakCorrLst" );

zakCorrMeanLst = load( fNameRandCircCorrFull, "zakCorrMeanFullLst" );

nCircLst, rCircLst, divNumNxtLst = ( vName -> load( fNameRandCircParams, vName ) ).( ("nCircLst", "rCircLst", "divNumNxtLst") );
lnNCirc, lnRCirc = length.( (nCircLst, rCircLst) );

dIt = 3;
# zakCorrMeanLst = ( corr -> dropdims( mean( corr; dims = dIt ); dims = dIt ) ).( zakCorrLst );

zakCorrFftLst = fft.( zakCorrMeanLst );

freqSqXLst = ( d -> [0:d-1;] ).( divNumNxtLst );
freqSqXLst = ( x -> x .= x.^2 ).(freqSqXLst);
freqSqYLst = transpose.( freqSqXLst );

freqNormLst = ( (fx2,fy2) -> sqrt.( fx2 .+ fy2 ) ).( freqSqXLst, freqSqYLst );

fMainCorrFft = "randCircCorrFft";
fNameCorrFft = fNameFunc( fMainCorrFft, attrLstBase, valLstBase, jld2Type );

jldsave( fNameCorrFft; zakCorrFftLst, freqNormLst );

fNameArr = push!( fNameArr[1:3], fNameCorrFft );

jldsave( fNameFNameArr; fNameArr );
	