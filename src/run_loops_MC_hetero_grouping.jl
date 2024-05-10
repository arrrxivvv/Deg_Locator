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

fMod = "";

isTestingParam = false;
# isTestingParam = true;

# itNumLst = [30000];
# # itNumLst = [10000];

# cAreaLst = [2:1:10;];
# cPerimLst = -[2];

fNameLst = Vector{String}(undef,0);

# isInit0Lst = [false,true];
isInit0Lst = [false];
divNumLst = [64];
# itNumLst = [10000];
# itNumLst = [2000];
itNumLstLst = [[10000],[2000]];
# itNumLstLst = [[50],[30]];

cFerroRatioLst = [0];
sgnAreaLst = [1];
sgnPerimLst = [1];
# cRadiusLst = [0.2:0.2:3;];
cRadiusLst = [0.6:0.2:3;];
cAngleLst = [0.6:0.2:0.8;].*pi/2;
# cAngleLst = [0.2:0.2:0.8;].*pi/2;
radiusGrp = Loops_MC.RadiusParamsGroup( cRadiusLst, cAngleLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst );

# cAreaLst = [-1.0];
cAreaLst = [-7:1.0:7;];
# cAreaLst = [-7:1.0:-7;];
cPerimLst = [-1.0:-1:-3;];
lnCPerim = length(cPerimLst);
cAreaLstScaledLst = cAreaLst .* [1:lnCPerim;]';
cFerroLst = [0.0];
cartesianGrp = Loops_MC.CartesianParamsGroup( cAreaLst, cPerimLst, cFerroLst );
cartesianGrpLstShort = Loops_MC.ParamsGroup[cartesianGrp];
cartesianGrpLst = Loops_MC.ParamsGroup[Loops_MC.CartesianParamsGroup( cAreaLstScaledLst[:,iGrp], cPerimLst[iGrp:iGrp], copy(cFerroLst) ) for iGrp = 1 : lnCPerim];

paramsGroupLstFlat = [radiusGrp,cartesianGrp];

paramsGroupLst = Loops_MC.ParamsGroup[radiusGrp,cartesianGrp];
radiusGrpLst = Loops_MC.ParamsGroup[radiusGrp];
paramsGroupLstLst = [radiusGrpLst, cartesianGrpLst];

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

# updaterType = Loops_MC.ABUpdater;
updaterType = Loops_MC.StaggeredCubeUpdater;

# paramsGroupLstFlat = Loops_MC.ParamsGroup[Loops_MC.BetaParamsGroup(betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst)];

# itNumLst = itNumLst;
# paramsGroupLst = paramsGroupLstFlat;

itNumLstIn = itNumLstLst;
paramsGroupLst = paramsGroupLstLst;

if isRunSim
	# for paramsGrpLst in paramsGroupLstLst
	Loops_MC.runLoopMC_withParamsGroup( updaterType, itNumLstIn, divNumLst, isInit0Lst, paramsGroupLst; fMod = fMod );
	# end
end

if isRunFileLst
	fNameNumLst = Loops_MC.genFNameLstLoopMC( updaterType, itNumLstIn, divNumLst, isInit0Lst, paramsGroupLst; fMain = fMainCollect, fMod = fMod );
end

if isRunSaveParam
	fNameParamsSave = Loops_MC.saveParamsLoopMC( updaterType, itNumLstIn, divNumLst, isInit0Lst, paramsGroupLst; fMod = fMod );
end
