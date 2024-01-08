using Loops_MC
using JLD2
# using Infiltrator
using DelimitedFiles
using FilenameManip

# itNumLst = [30000];
itNumLst = [100];

# cAreaLst = [2:1:10;];
# cPerimLst = -[2];

fModSmart = "smart";

betaLst = [0.2,0.5,1,2,3,10,100];
# betaLst = [5:1:15;];
# betaLst = [10];
betaBase = 1;
cRatioLst = exp.( [-2:0.2:2;] );
# cRatioLst = exp.( [-3:0.1:-1;] );
# cRatioLst = exp.( [-3:1:-3;] );
# cFerroRatioLst = exp.( [-1:1:5;] )
cFerroRatioLst = [0];
sgnArea = -1;
sgnPerim = 1;
# betaLst = 1;

cAreaLst = zeros( length(cRatioLst), length(betaLst) );
cPerimLst = similar(cAreaLst);
cFerroLst = zeros( length(cFerroRatioLst), length(cRatioLst), length(betaLst) );

isInit0Lst = [false,true];

divNumLst = [64];

fNameLst = Vector{String}(undef,0);

attrLst = ["divNum","itNum","beta","cRatio","sgnArea","sgnPerim"];
valLst = [divNumLst[1], itNumLst[1], betaLst[[1,end]],cRatioLst[[1,end]],sgnArea,sgnPerim];

fMod = "withInit0";

fMainParam = "loops_paramLst";
fNameParam = fNameFunc( fMainParam, attrLst, valLst, jld2Type; fMod = fMod );

# for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
for itNum in itNumLst, iCRatio = 1 : length(cRatioLst), iCFerro = 1 : length(cFerroRatioLst), iBeta = 1 : length(betaLst), divNum in divNumLst, isInit0 in isInit0Lst
	beta = betaLst[iBeta];
	cRatio = cRatioLst[iCRatio];
	# fName = loops_MC( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta );
	cRatioSq = sqrt(cRatio);
	cArea = sgnArea * beta * cRatioSq;
	cPerimAbs = beta / cRatioSq;
	cPerim = sgnPerim * cPerimAbs;
	cFerro = cPerimAbs * cFerroRatioLst[iCFerro];
	
	cAreaLst[iCRatio, iBeta] = cArea;
	cPerimLst[iCRatio, iBeta] = cPerim;
	cFerroLst[iCFerro, iCRatio, iBeta] = cFerro;
end

cAreaStrLst = string.( cAreaLst );
cPerimStrLst = string.( cPerimLst );
betaStrLst = string.( betaLst );

save( fNameParam, "betaLst", betaLst, "cRatioLst", cRatioLst, "cAreaLst", cAreaLst, "cPerimLst", cPerimLst, "cAreaStrLst", cAreaStrLst, "cPerimStrLst", cPerimStrLst, "cFerroLst", cFerroLst, "betaStrLst", betaStrLst );
