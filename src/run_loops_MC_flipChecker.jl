using Loops_MC
# using Infiltrator
using DelimitedFiles

@enum RunWhich runFull=1 runFileLst runSaveParam runSaveParamAndFileLst

runChoice = runFull;

isRunSim = false;
isRunFileLst = false;
isRunSaveParam = false;
if runChoice == runFull
	isRunSim = true;
	isRunFileLst = true;
	isRunSaveParam = true;
elseif runChoice == runSaveParamAndFileLst
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
# itNumLst = [10];
itNumLst = [1000];
# itNumLstLst = [[10000],[2000]];
# itNumLstLst = [[200],[100]];

cFerroRatioLst = [0];
sgnAreaLst = [1];
sgnPerimLst = [1];
# cRadiusLst = [0.2:0.2:3;];
# cRadiusLst = [0.6:0.2:3;];
cRadiusLst = [0.4:0.05:1;];
# cAngleLst = [0.8:0.2:0.8;].*pi/2;
cAngleLst = [0.2:0.2:0.8;].*pi/2;
radiusGrp = Loops_MC.RadiusParamsGroup( cRadiusLst, cAngleLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst );

# cAreaLst = [-1.0];
# cAreaLst = [-7:1.0:7;];
# cAreaLst = [-7:1.0:-7;];
cAreaLst = [0:0.1:0.3;];
# cPerimLst = [-1.0:-1:-3;];
cPerimLst = [0];
lnCPerim = length(cPerimLst);
cAreaLstScaledLst = cAreaLst .* [1:lnCPerim;]';
cFerroLst = [0.0];
cartesianGrp = Loops_MC.CartesianParamsGroup( cAreaLst, cPerimLst, cFerroLst );
cartesianGrpLstShort = Loops_MC.ParamsGroup[cartesianGrp];
cartesianGrpLst = Loops_MC.ParamsGroup[Loops_MC.CartesianParamsGroup( cAreaLstScaledLst[:,iGrp], cPerimLst[iGrp:iGrp], copy(cFerroLst) ) for iGrp = 1 : lnCPerim];

# paramsGroupLstFlat = [radiusGrp,cartesianGrp];
# paramsGroupLstFlat = Loops_MC.ParamsGroup[cartesianGrp];
paramsGroupLstFlat = Loops_MC.ParamsGroup[radiusGrp];

paramsGroupLst = Loops_MC.ParamsGroup[radiusGrp,cartesianGrp];
radiusGrpLst = Loops_MC.ParamsGroup[radiusGrp];
paramsGroupLstLst = [radiusGrpLst, cartesianGrpLst];

# fMainCollect = Loops_MC.oFNameLoopsMain;
# fMainCollect = Loops_MC.oFNameLoopsNumMain;
# fMainCollect = Loops_MC.oFNameLoopsStartMain;
# fMainCollect = Loops_MC.oFMainLoopsSample;
fMainCollect = Loops_MC.getAuxDataSummarySampleName( Loops_MC.ZakArrAuxData );
# fMainCollect = Loops_MC.getAuxDataSummaryNumName( Loops_MC.BLinkAuxData );
# fMainCollect = Loops_MC.getAuxDataSummaryItSampleLstName( Loops_MC.ZakArrAuxData );

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
updaterType = Loops_MC.StaggeredCubeUpdaterBase;
# updaterType = Loops_MC.CubeUpdater;
# updaterType = Loops_MC.CubeStaggeredCubeUpdater;
# updaterType = Loops_MC.SwitchingUpdater{Tuple{Loops_MC.StaggeredCubeUpdaterBase,Loops_MC.StaggeredCubeUpdaterBase}};

# paramsGroupLstFlat = Loops_MC.ParamsGroup[Loops_MC.BetaParamsGroup(betaLst, cRatioLst, cFerroRatioLst, sgnAreaLst, sgnPerimLst)];

itNumLstIn = itNumLst;
paramsGroupLst = paramsGroupLstFlat;

# itNumLstIn = itNumLstLst;
# paramsGroupLst = paramsGroupLstLst;

initType = Loops_MC.BinomialInitializer;
# initType = Loops_MC.ConstantInitializer;

itControllerType = Loops_MC.ItNumItController;
# itControllerType = nothing;

# flipProposerType = Loops_MC.OneFlipProposer;
# flipProposerType = Loops_MC.CubeFlipProposer;
flipProposerType = Loops_MC.SwitchingFlipProposer{Tuple{Loops_MC.OneFlipProposer,Loops_MC.CubeFlipProposer}};
# flipProposerType = Loops_MC.SwitchingFlipChecker{Tuple{Loops_MC.OneFlipProposer,Loops_MC.CubeFlipProposer}};
# flipProposerType = nothing;

# flipCheckerType = Loops_MC.NeighborFlipChecker;
# flipCheckerType = Loops_MC.CubeFlipChecker;
# flipCheckerType = Loops_MC.SwitchingFlipChecker{Tuple{Loops_MC.NeighborFlipChecker,Loops_MC.CubeFlipChecker}};
flipCheckerType = Loops_MC.IsingFlipChecker;

if isRunSim
	Loops_MC.runLoopMC_withParamsGroup( updaterType, initType, itNumLstIn, divNumLst, paramsGroupLst; fMod = fMod, flipCheckerType = flipCheckerType, flipProposerType = flipProposerType, itControllerType = itControllerType );
end

if isRunFileLst
	fNameNumLst = Loops_MC.genFNameLstLoopMC( updaterType, initType, itNumLstIn, divNumLst, paramsGroupLst; fMain = fMainCollect, fMod = fMod, flipCheckerType = flipCheckerType, flipProposerType = flipProposerType, itControllerType = itControllerType );
end

if isRunSaveParam
	fNameParamsSave = Loops_MC.saveParamsLoopMC( updaterType, initType, itNumLstIn, divNumLst, paramsGroupLst; fMod = fMod );
end
