using Loops_MC
# using Infiltrator
using DelimitedFiles

@enum RunWhich runFull=1 runFileLst runSaveParam

runChoice = runSaveParam;

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

fModSmartOld = "smart";
fModStag = "stag";
# fMod = fModSmartOld;
fMod = fModStag;
# fMod = "";

# isFModMethod = true;
isFModMethod = false;

isTestingParam = false;
# isTestingParam = true;

itNumLst = [30000];
# itNumLst = [10000];

cAreaLst = [2:1:10;];
cPerimLst = -[2];

# betaLst = [0.2,0.5,1,2,3,10,100];
betaLst = [0.2,1,3,10];
# betaLst = [0.2:0.2:3;];
# betaLst = [5:1:15;];
betaBase = 1;
cRatioLst = exp.( [-2:0.4:2;] );
# cRatioLst = exp.( [-1:1:1;] );
# cRatioLst = exp.( [-1:0.4:1;] );
# cRatioLst = exp.( [-5:1:-1;] );
# cRatioLst = exp.( [-3:1:-3;] );
# cRatioLst = exp.( [-1:1:-1;] );
# cFerroLst = [0:0.1:1;];
# cFerroRatioLst = exp.( [-1:1:5;] );
cFerroRatioLst = [0];
sgnAreaLst = [1,-1];
sgnPerimLst = [1,-1];
# sgnAreaLst = [1];
# sgnPerimLst = [1];

isInit0Lst = [false,true];

divNumLst = [64];

fNameLst = Vector{String}(undef,0);

# fMainCollect = Loops_MC.oFNameLoopsMain;
fMainCollect = Loops_MC.oFNameLoopsNumMain;
# fMainCollect = Loops_MC.oFNameLoopsStartMain;

if isTestingParam
	itNumLst = [100];

	betaLst = [10.0];
	betaBase = 1;
	cRatioLst = exp.( [-1:1:-1;] );
	cFerroRatioLst = [0];
	sgnAreaLst = [1];
	sgnPerimLst = [1];

	isInit0Lst = [false];
end

# updaterType = Loops_MC.ABUpdater;
updaterType = Loops_MC.StaggeredCubeUpdater;

if isRunSim
	Loops_MC.runLoopMC_withParams( updaterType, itNumLst, divNumLst, betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst, isInit0Lst; fMod = fMod );
end

if isRunFileLst
	fNameNumLst = Loops_MC.genFNameLstLoopMC( updaterType, itNumLst, divNumLst, betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst, isInit0Lst; fMain = fMainCollect, fMod = fMod, isFModMethod = isFModMethod );
end

if isRunSaveParam
	fNameParamsSave = Loops_MC.saveParamsLoopMC( updaterType, itNumLst, divNumLst, betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst, isInit0Lst; fMod = fMod, isFModMethod = isFModMethod );
end
