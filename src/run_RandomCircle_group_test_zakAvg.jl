using Statistics

using RandomCircle
using Utils
using SharedFNames
using JLD2
using FilenameManip	

fNameFNameArr = Utils.strReadLastLine( SharedFNames.dirLog * SharedFNames.fNameTmpNameFileLst );

fNameArr = load( fNameFNameArr, "fNameArr" );

fNameRandCircZak, fNameRandCircCorr, fNameRandCircFine, fNameRandCircCorners = @view fNameArr[1:4];

# zakArrLst = load( fNameArr[3], "zakArrLst" );
# zakCorrLst = load( fNameArr[3], "zakCorrLst" );

nCircLst = load( fNameArr[1], "nCircLst" );
rCircLst = load( fNameArr[1], "rCircLst" );
fNameRandCircSingleLst = load( fNameArr[3], "fNameRandCircSingleLst" );

lnNCirc, lnRCirc = length.( (nCircLst, rCircLst) );
# itNum = 10;

fMainRandCircShAvgSingle = fMainRandCircZak * "ShAvgSingle";

fNameRandCircShAvgSingleLst = Array{String}(undef, lnRCirc, lnNCirc);

for iR = 1 : lnRCirc, iN = 1 : lnNCirc
	valLstSingle[2] = nCircLst[iN];
	valLstSingle[3] = rCircLst[iR];
	fNameRandCircShAvgSingleLst[iR, iN] = fNameFunc( fMainRandCircShAvgSingle, attrLstSingle, valLstSingle, jld2Type );
end

zakAvgLst = [ zeros(itNum) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
# zakCorrAvgShAvg1dLst = (arr -> arr[:,1]).( zakCorrAvgShAvgLst );
# zakCorrAvgShAvg1dHalfFineLst = [ @view( zakCorrAvgShAvg1dLst[iR,iN][1:divNumNxtHalf[iR,iN]] ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
zakCorrAvgShAvg1dLst = Array{Vector{Float64}}(undef, lnRCirc, lnNCirc);
# zakCorrAvgShAvg1dHalfFineLst = Array{AbstractVector{Float64}}(undef, lnRCirc, lnNCirc);;


for iN = 1 : lnNCirc, iR = 1 : lnRCirc
	zakArrLst = load( fNameRandCircSingleLst[iR,iN], "zakArrLst" );
	zakCorrShAvgLst = load( fNameRandCircSingleLst[iR,iN], "zakCorrLst" );
	GC.gc()
	for it = 1 : itNum
		zakAvgLst[iR,iN][it] = RandomCircle.calcZakAvg( @view( zakArrLst[:,:,it] ) );
	end
	zakCorrShAvgLst .-= reshape( zakAvgLst[iR,iN], 1, 1, itNum ).^2;
	zakCorrAvgShAvgLst = mean( zakCorrShAvgLst; dims = 3 );
	zakCorrAvgShAvg1dLst[iR,iN] = zakCorrAvgShAvgLst[:,1];
	
	jldsave( fNameRandCircShAvgSingleLst[iR, iN]; zakCorrShAvgLst, zakCorrAvgShAvgLst )
end
# zakAvgLst = [ [ RandomCircle.calcZakAvg( @view( zakArrLst[iR,iN][:,:,it] ) ) for it = 1 : itNum ] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

zakAbsAvgLst = mean.( abs, zakAvgLst );

zakCorrAvgShAvg1dHalfFineLst = [ @view( zakCorrAvgShAvg1dLst[iR,iN][1:divNumNxtHalf[iR,iN]] ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

# zakCorrShAvgLst = deepcopy( zakCorrLst );
# ( ( corr, avg )-> corr .= corr .-= reshape(avg, 1, 1, length(avg)).^2 ).( zakCorrShAvgLst, zakAvgLst );

# zakCorrAvgShAvgLst = mean.( zakCorrShAvgLst; dims = 3 );

fittedModelFineAvgShLst = [ curve_fit( modelExp, xLstZakHalfLst[iR,iN], zakCorrAvgShAvg1dHalfFineLst[iR,iN], p0Fit ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
expShFineAvgShLst = ( m -> m.param[3] ).( fittedModelFineAvgShLst );

fMainTmp = fMainRandCircZak * "_tmp";
fNameTmp = fNameFunc( fMainTmp, attrLstBase, valLstBase, jld2Type );

jldsave( fNameTmp; zakAvgLst, zakAbsAvgLst, zakCorrShAvgLst, zakCorrAvgShAvgLst, zakCorrAvgShAvg1dLst, expShFineAvgShLst );

fNameArr = [ fNameRandCircZak, fNameRandCircCorr, fNameRandCircFine, fNameRandCircCorners, fNameTmp ];

jldsave( fNameFNameArr; fNameArr );
