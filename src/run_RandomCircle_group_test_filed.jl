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

nCirc = 10;
nDim3 = 3;

itNum = 100;

nCircStep = 5;
rCircStep = 0.1;
nCircLst = [5:nCircStep:50;];
rCircLst = [0.1:rCircStep:0.5;];
lnNCirc = length( nCircLst );
lnRCirc = length( rCircLst );

rRaw = 1;

nSample = 1000;

divNum1Pass = 128;

zakArrLst1Pass = zeros( Bool, divNum1Pass, divNum1Pass, itNum, lnRCirc, lnNCirc );
zakCorrLst1Pass = similar( zakArrLst1Pass, Float64 );
zakCorrAvg1Pass = zeros( divNum1Pass, divNum1Pass, 1, lnRCirc, lnNCirc );

randCircDataLst = [ RandomCircle.RandCircData( nCircLst[iN], rRaw; divNum = divNum1Pass, nSample = nSample ) for iN = 1 : lnNCirc ];

rotMatBackupLst = [ [ [ zeros(3,3) for iCirc = 1 : nCircLst[iN] ] for it = 1 : itNum ] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
rotMat2dBackupLst = [ [ [ zeros(2,2) for iCirc = 1 : nCircLst[iN] ] for it = 1 : itNum ] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
rotMat2dInvBackupLst = deepcopy( rotMat2dBackupLst );

sampleLstLst = [ [ [ zeros(2, nSample) for iCirc = 1 : nCircLst[iNCirc] ] for it = 1 : itNum ] for iRCirc = 1 : lnRCirc, iNCirc = 1 : lnNCirc ];

for iN = 1 : lnNCirc
	data = randCircDataLst[iN];
	for iR = 1 : lnRCirc
		RandomCircle.setRCircNoRotUpdate!( data, rCircLst[iR] );
		for it = 1 : itNum
			RandomCircle.refreshPtLst!( data );
			RandomCircle.refreshQuatLst!( data );
			RandomCircle.refreshRotMatFull!( data );
			RandomCircle.refreshEqCoeff!( data );
			RandomCircle.refreshBndLstFull!( data );
			RandomCircle.calcZakArr!( data );
			RandomCircle.calcZakCorr!( data );
			
			RandomCircle.backupRotMat!( rotMatBackupLst[iR, iN][it], rotMat2dBackupLst[iR, iN][it], rotMat2dInvBackupLst[iR, iN][it], data );
			RandomCircle.backupZakArrCorr!( @view( zakArrLst1Pass[:, :, it, iR, iN] ), @view( zakCorrLst1Pass[:, :, it, iR, iN] ), data );
			
			RandomCircle.calcSampleLst!( data );
			RandomCircle.backupSampleLst!( sampleLstLst[iR, iN][it], data );
		end
	end
end

zakCorrAvg1Pass = mean( zakCorrLst1Pass; dims = 3 );
zakCorrAvg1d1Pass = @view zakCorrAvg1Pass[:,1,1,:,:];

zakCorrAvg1d1PassSave = zakCorrAvg1Pass[:,1,1,:,:];

xLstZak = [0.0:divNum1Pass-1;];
xLstZak .*= 1 / divNum1Pass;
divNum1PassHalf = div( divNum1Pass,2 );
xLstZakHalf = xLstZak[1:divNum1PassHalf];
zakCorrAvg1dHalf1PassArr = [ @view zakCorrAvg1Pass[1:divNum1PassHalf,1,1,iR,iN] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

modelExp( x, p ) = p[1] .* exp.( - p[2] .* x ) .+ p[3];
p0Fit = [1, 0.5, 0];

fittedModelLst = [ curve_fit( modelExp, xLstZakHalf, zakCorrAvg1dHalf1PassArr[iR, iN], p0Fit ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
# fittedModelLst = [ curve_fit( modelExp, xLstZak, zakCorrAvg1d1PassArr[iR, iN], p0Fit ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

corrLenInvLst = ( x -> x.param[2] ).( fittedModelLst );
corrLenLst = 1 ./ corrLenInvLst;

expScaleLst = ( x -> x.param[1] ).( fittedModelLst );
expShLst = ( x -> x.param[3] ).( fittedModelLst );

divNumNxt = ( x -> x < corrLenLst[1,1] ? Int64( floor( divNum1Pass * corrLenLst[1,1] / x ) ) : divNum1Pass ).(corrLenLst);

# zakArrLst = [ zeros( Bool, divNumNxt[iR,iN], divNumNxt[iR,iN], itNum ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
# zakCorrLst = [ zeros( divNumNxt[iR,iN], divNumNxt[iR,iN], itNum ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
zakArrLst = Array{Bool,3}(undef,0,0,0);
zakCorrLst = Array{Float64,3}(undef,0,0,0);

fMainRandCircZak = "randCircZak";
attrLstBase = [ "nCircLst", "rCircLst", "itNum" ];
valLstBase = Any[FilenameManip.fAttr_arrSummary.( (nCircLst, rCircLst), (nCircStep, rCircStep) )..., itNum];
fNameRandCircZak = fNameFunc( fMainRandCircZak, attrLstBase, valLstBase, jld2Type );
jldsave( fNameRandCircZak; sampleLstLst, zakArrLst, divNumNxt, zakArrLst1Pass, nCircLst, rCircLst, itNum );

fMainRandCircCorr = "randCircCorr";
fNameRandCircCorr = fNameFunc( fMainRandCircCorr, attrLstBase, valLstBase, jld2Type );
jldsave( fNameRandCircCorr; corrLenLst, expScaleLst, expShLst, zakCorrAvg1d = zakCorrAvg1d1PassSave );

fMainRandCircZakSingle = fMainRandCircZak * "single";
fNameRandCircSingleLst = Array{String}(undef, lnRCirc, lnNCirc);
attrLstSingle = ["nCirc", "rCirc", "itNum"];
valLstSingle = Any[nCircLst[1], rCircLst[1], itNum];

for iN = 1 : lnNCirc, iR = 1 : lnRCirc
	valLstSingle[1] = nCircLst[iN];
	valLstSingle[2] = rCircLst[iR];
	fNameRandCircSingleLst[iR, iN] = fNameFunc( fMainRandCircZakSingle, attrLstSingle, valLstSingle, jld2Type );
end

zakCorrAvgFine = Matrix{Float64}(undef, 0, 0);
zakCorrAvg1dFineLst = Matrix{Vector{Float64}}(undef, lnRCirc, lnNCirc);

dAvg = 3;
for iR = 1 : lnRCirc, iN = 1 : lnNCirc
	data = randCircDataLst[iN];
	divNum = divNumNxt[iR, iN];
	RandomCircle.setDivNum!( data, divNumNxt[iR, iN] );
	zakArrLst = zeros( Bool, divNum, divNum, itNum );
	zakCorrLst = similar( zakArrLst, Float64 );
	GC.gc();
	for it = 1 : itNum
		RandomCircle.restoreRotMat!( data, rotMatBackupLst[iR,iN][it], rotMat2dBackupLst[iR,iN][it], rotMat2dInvBackupLst[iR,iN][it] );
		
		RandomCircle.calcZakArr!( data );
		RandomCircle.calcZakCorr!( data );
		
		RandomCircle.backupZakArrCorr!( @view( zakArrLst[:,:,it] ), @view( zakCorrLst[:,:,it] ), data );
	end
	
	zakCorrAvgFine = dropdims( mean( zakCorrLst; dims = dAvg ); dims = dAvg );
	zakCorrAvg1dFineLst[iR, iN] = zakCorrAvgFine[:,1];
	# @infiltrate
	
	jldsave( fNameRandCircSingleLst[iR,iN]; zakArrLst, zakCorrLst, zakCorrAvgFine );
end

dAvg = 3;


xLstZakLst = [ [0:divNumNxt[iR,iN]-1;] ./ divNumNxt[iR,iN] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
divNumNxtHalf = div.( divNumNxt, 2 );
xLstZakHalfLst = ( ( a, i ) -> a[1:i] ).(xLstZakLst, divNumNxtHalf);
zakCorrAvg1dHalfFineLst = [ @view( zakCorrAvg1dFineLst[iR,iN][1:divNumNxtHalf[iR,iN]] ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

fittedModelFineLst = [ curve_fit( modelExp, xLstZakHalfLst[iR,iN], zakCorrAvg1dHalfFineLst[iR,iN], p0Fit ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

expScaleFineLst = ( m -> m.param[1] ).( fittedModelFineLst );
expShFineLst = ( m -> m.param[3] ).( fittedModelFineLst );
corrLenInvFineLst = ( m -> m.param[2] ).( fittedModelFineLst );
corrLenFineLst = 1 ./ corrLenInvFineLst;

fMainRandCircFine = "randCircFine";
fNameRandCircFine = fNameFunc( fMainRandCircFine, attrLstBase, valLstBase, jld2Type );
jldsave( fNameRandCircFine; fNameRandCircSingleLst, expScaleFineLst, expShFineLst, corrLenFineLst, xLstZakLst, zakCorrAvg1dFineLst );





 
fNameArr = [ fNameRandCircZak, fNameRandCircCorr, fNameRandCircFine ];


fMainFNameArr = "fNameArr" * fNameRandCircZak * "_filed";
fNameFNameArr = fNameFunc( fMainFNameArr, attrLstBase, valLstBase, jld2Type );
save( fNameFNameArr, "fNameArr", fNameArr );

open( SharedFNames.dirLog * SharedFNames.fNameTmpNameFileLst, "w" ) do io
	println( io, fNameFNameArr );
end
