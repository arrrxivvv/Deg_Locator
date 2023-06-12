using Loops_MC

itNumLst = [10000];

cAreaLst = [0:0.5:3;];
cPerimLst = [0:0.5:3;];

betaLst = 1;

divNumLst = [16,64];

isValLstFloat = true;

for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
	Loops_MC.zakAvgFromFile( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta, isValLstFloat = isValLstFloat );
end
