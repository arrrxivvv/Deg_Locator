using RandomCircle
using FilenameManip
using JLD2
using SharedFNames
using DelimitedFiles

isRenewRand = false;

nCirc = 10;
nDim3 = 3;

nSample = 1000;
rCirc = 0.2;

divNum = 128;

attrLstBase = ["nCirc", "rCirc"];
valLstBase = Any[nCirc, rCirc];
attrLstSample = push!( copy(attrLstBase), "nSample" );
valLstSample = push!( copy(valLstBase), nSample );

# randCircData = RandomCircle.RandCircData( nCirc, rCirc );

# RandomCircle.refreshPtLst!( randCircData );

# RandomCircle.refreshQuatLst!( randCircData );

# RandomCircle.refreshRotMatFull!( randCircData );

# RandomCircle.refreshEqCoeff!( randCircData );

# RandomCircle.calcSampleLst!( randCircData );

fMainSampleCirc = "sampleCircLst";
fNameSampleCirc = fNameFunc( fMainSampleCirc, attrLstSample, valLstSample, jld2Type );
save( fNameSampleCirc, "sampleCircLst", randCircData.sampleCircLst );

# RandomCircle.refreshBndLst!( randCircData );

# RandomCircle.refreshBndPtsModExtendLst!( randCircData );

# RandomCircle.refreshSortBndLst!( randCircData );

RandomCircle.calcZakArr!( randCircData );

fMainBndCirc = "circBndLst";
fNameBndCirc = fNameFunc( fMainBndCirc, attrLstBase, valLstBase, jld2Type );

save( fNameBndCirc, "bndLst", randCircData.bndLst );



fMainZakCirc = "zakArrCirc"
attrLstZak = push!( deepcopy( attrLstBase ), "divNum" );
valLstZak = push!( deepcopy( valLstBase ), divNum );
fNameZakCirc = fNameFunc( fMainZakCirc, attrLstZak, valLstZak, jld2Type );
save( fNameZakCirc, "zakArr", randCircData.zakArr, "zakXLst", randCircData.zakXLst );

fNameArr = [fNameSampleCirc, fNameBndCirc, fNameZakCirc];
fNameLst = vec(fNameArr);

fMainFNameArr = "fNameArrJLD2";
fNameFNameArr = fNameFunc( fMainFNameArr, attrLstBase, valLstBase, jld2Type );
fMainFNameLst = "fNameLst";
fNameFNameLst = fNameFunc( fMainFNameLst, attrLstBase, valLstBase, txtType );

save( fNameFNameArr, "fNameArr", fNameArr );
writedlm( fNameFNameLst, fNameLst );

open( SharedFNames.dirLog * SharedFNames.fNameTmpNameFileLst, "w" ) do io
	println( io, fNameFNameArr );
end
