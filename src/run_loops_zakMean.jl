using Loops_MC

itNumLst = [10000];

# cAreaLst = [1];
cAreaLst = [0.5:0.5:1.5;];
cPerimLst = -[0.5:0.5:3;];
# cPerimLst = [0];

fModSmart = "smart";

# betaLst = [0.2,0.5,1,2,3];
betaLst = [1];

divNumLst = [64];

isValLstFloat = false;

fModSmart = "smart";

for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
	# Loops_MC.zakAvgFromFile( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta, isValLstFloat = isValLstFloat );
	# Loops_MC.zakResaveFromFile( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta, isValLstFloat = isValLstFloat, fMod = fModSmart );
	# Loops_MC.zakCorrFromFile( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta, isValLstFloat = isValLstFloat, fMod = fModSmart );
	Loops_MC.loopsSamplesResaveFromFile( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta, isValLstFloat = isValLstFloat, fMod = fModSmart );
end
