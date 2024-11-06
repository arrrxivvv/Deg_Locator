using RandomCircle
using FilenameManip
using JLD2
using SharedFNames
using DelimitedFiles
using Utils
using LsqFit

isRenewRand = false;

nCirc = 10;
nDim3 = 3;

itNum = 10;

nCircLst = [5:5:50;];
rCircLst = [0.1:0.1:0.5;];
lnNCirc = length( nCircLst );
lnRCirc = length( rCircLst );

rRaw = 1;

divNum1Pass = 128;

zakArrLst1Pass = zeros( Bool, divNum1Pass, divNum1Pass, itNum, lnRCirc, lnNCirc );
zakCorrLst1Pass = similar( zakArrLst1Pass, Float64 );
zakCorrAvg1Pass = zeros( divNum1Pass, divNum1Pass, 1, lnRCirc, lnNCirc );

randCircDataLst = [ RandCircData( nCircLst[iN], rRaw; divNum = divNum1Pass ) for iN = 1 : lnNCirc ];

rotMatBackupLst = [ [ [ zeros(3,3) for iCirc = 1 : nCircLst[iN] ] for it = 1 : itNum ] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
rotMat2dBackupLst = [ [ [ zeros(2,2) for iCirc = 1 : nCircLst[iN] ] for it = 1 : itNum ] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
rotMat2dInvBackupLst = deepcopy( rotMat2dBackupLst );

for iN = 1 : lnNCirc
	data = randCircDataLst[iN];
	for iR = 1 : lnRCirc
		setRCircNoRotUpdate!( data, rCircLst[iR] );
		for it = 1 : itNum
			refreshQuatLst!( data );
			refreshRotMatFull!( data );
			refreshEqCoeff!( data );
			refreshBndLstFull!( data );
			calcZakArr!( data );
			calcZakCorr!( data );
			
			backupRotMat!( rotMatBackupLst[iR, iN][it], rotMat2dBackupLst[iR, iN][it], rotMat2dInvBackupLst[iR, iN][it], data );
			backupZakArrCorr!( @view zakArrLst1Pass[:, :, it, iR, iN], @view zakCorrLst1Pass[:, :, it, iR, iN], data );
		end
	end
end

zakCorrAvg1Pass = mean( zakCorrLst1Pass; dims = 3 );
zakCorrAvg1d1Pass = @view zakCorrAvg1Pass[:,1,1,:,:];
zakCorrAvg1d1PassArr = [ @view zakCorrAvg1Pass[:,1,1,iR,iN] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

xLstZak = [0:divNum1Pass-1;];
xLstZak .*= 1 / divNum1Pass;

modelExp( x, p ) = p[1] .* exp.( - p[2] .* x ) .+ p[3];

p0Fit = [1, 0.5, 0];
fittedModelLst = [ curve_fit( modelExp, xLstZak, zakCorrAvg1d1PassArr[iR, iN], p0Fit ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

corrLenLst = ( x -> x.param[2] ).( fittedModelLst );

divNumNxt = Int64.( floor.( divNum1Pass .* corrLenLst[1,1] ./ corrLenLst ) );

for iN = 1 : lNCirc
	setDivNum!( randCircDataLst[iN] );
end
