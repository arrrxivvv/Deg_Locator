using Loops_MC

itNumLst = [30000];

# cAreaLst = [1];
cAreaLst = [0.5:0.5:1.5;];
cPerimLst = -[0.5:0.5:3;];
# cPerimLst = [0];

# betaLst = [0.2,0.5,1,2,3,10,100];
betaLst = [10];
betaBase = 1;
# cRatioLst = exp.( [-2:0.2:2;] );
cRatioLst = exp.( [-3:1:-3;] );
# cFerroLst = [2];
cFerroRatioLst = exp.( [-1:1:5;] )
sgnArea = 1;
sgnPerim = -1;
# betaLst = 1;

itLst = 1:10

isInit0Lst = [false,true];

divNumLst = [64];

isValLstFloat = false;

fModSmart = "smart";
fModSmartInit0 = "smart" * "_" * "isInit0";

# for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
for itNum in itNumLst, cFerroRatio in cFerroRatioLst, cRatio in cRatioLst, beta in betaLst, divNum in divNumLst, isInit0 in isInit0Lst
	cRatioSq = sqrt(cRatio);
	cArea = sgnArea * beta * cRatioSq;
	cPerim = sgnPerim * beta / cRatioSq;
	cFerro = cPerim * cFerroRatio;
	
	fMod = isInit0 ? fModSmartInit0 : fModSmart ;

	# Loops_MC.zakAvgFromFile( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta, isValLstFloat = isValLstFloat );
	# Loops_MC.zakResaveFromFile( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta, isValLstFloat = isValLstFloat, fMod = fModSmart );
	# Loops_MC.zakCorrFromFile( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta, isValLstFloat = isValLstFloat, fMod = fModSmart );
	# Loops_MC.loopsSamplesResaveFromFile( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta, isValLstFloat = isValLstFloat, fMod = fMod );
	
	Loops_MC.linkBfieldLstGenerateFromFile( divNum, itNum; cPerim = cPerim, cArea = cArea, cFerro = cFerro, beta = betaBase, fMod = fMod, itLst = itLst );
end
