using Loops_MC

itNumLst = [10000,20000,50000];

cAreaLst = [0:0.5:3;];
cPerimLst = [0:0.5:3;];

betaLst = 1;

divNumLst = [16,64];

for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
	loops_MC( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta );
end
