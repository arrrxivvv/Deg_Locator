using Loops_MC
# using Infiltrator
using DelimitedFiles

@enum RunWhich runFull=1 runFileLst runSaveParam

runChoice = runFull;

isRunSim = false;
isRunFileLst = false;
isRunSaveParam = false;
if runChoice == runFull
	isRunSim = true;
	isRunFileLst = true;
	isRunSaveParam = true;
elseif runChoice == runFileLst
	isRunFileLst = true;
elseif runChoice == runSaveParam
	isRunSaveParam = true;
end

fMod = "";

nDim = 2;

isTestingParam = false;
# isTestingParam = true;

fNameLst = Vector{String}(undef,0);

# isInit0Lst = [false,true];
isInit0Lst = [false];
divNumLst = [8];
itNumLst = [10000];
# itNumLst = [10];
# itNumLstLst = [[10000],[2000]];
# itNumLstLst = [[200],[100]];

# cAreaLst = [-1.0];
# cAreaLst = [-7:1.0:7;];
cAreaLst = [0;];
# cPerimLst = [0:0.2:3;];
cPerimLst = 1 ./ [1:0.1:3.5;];
cFerroLst = [0];
cartesianGrp = Loops_MC.CartesianParamsGroup( cAreaLst, cPerimLst, cFerroLst );

paramsGroupLstFlat = [cartesianGrp];

paramsGroupLstFlat = Loops_MC.ParamsGroup[cartesianGrp];

# fMainCollect = Loops_MC.oFNameLoopsMain;
fMainCollect = Loops_MC.oFNameLoopsNumMain;
# fMainCollect = Loops_MC.oFNameLoopsStartMain;

if isTestingParam
	itNumLst = [100];

	betaLst = [3.0];
	betaBase = 1;
	cRatioLst = exp.( [-2:0.4:0;] );
	# cRatioLst = exp.( [-2:0.4:-2;] );
	cFerroRatioLst = [0];
	sgnAreaLst = [1];
	sgnPerimLst = [1];

	isInit0Lst = [false];
	
	betaGrp = Loops_MC.BetaParamsGroup(betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst);
	
	cRadiusLst = [1];
	cAngleLst = [0.2:0.2:0.8;].*pi/2;
	# cAngleLst = [0.2:0.2:0.2;].*pi/2;
	radiusGrp = Loops_MC.RadiusParamsGroup( cRadiusLst, cAngleLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst );
	
	# cAreaLst = [-1.0];
	cAreaLst = [-7:1.0:7;];
	cPerimLst = [-1.0];
	cFerroLst = [0.0];
	cartesianGrp = Loops_MC.CartesianParamsGroup( cAreaLst, cPerimLst, cFerroLst );
	
	paramsGroupLst = Loops_MC.ParamsGroup[radiusGrp,cartesianGrp];
end

updaterType = Loops_MC.AB2dUpdater;
# updaterType = Loops_MC.SingleUpdater;

itNumLstIn = itNumLst;
paramsGroupLst = paramsGroupLstFlat;

# initType = Loops_MC.BinomialInitializer;
initType = Loops_MC.ConstantInitializer;

if isRunSim
	Loops_MC.runLoopMC_withParamsGroup( updaterType, initType, itNumLstIn, divNumLst, paramsGroupLst; fMod = fMod, nDim = nDim );
end

if isRunFileLst
	fNameNumLst = Loops_MC.genFNameLstLoopMC( updaterType, initType, itNumLstIn, divNumLst, paramsGroupLst; fMain = fMainCollect, fMod = fMod, nDim = nDim );
end

if isRunSaveParam
	fNameParamsSave = Loops_MC.saveParamsLoopMC( updaterType, initType, itNumLstIn, divNumLst, paramsGroupLst; fMod = fMod, nDim = nDim );
end
