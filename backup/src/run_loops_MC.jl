using Loops_MC
# using Infiltrator
using DelimitedFiles

itNumLst = [30000];
# itNumLst = [100];

# cAreaLst = [0.5:0.5:1.5;];
cAreaLst = [2:1:10;];
# cAreaLst = [1];
cPerimLst = -[2];
# cPerimLst = [0];

fModSmart = "smart";

betaLst = [0.2,0.5,1,2,3,10,100];
# betaLst = [5:1:15;];
# betaLst = [10];
betaBase = 1;
cRatioLst = exp.( [-2:0.2:2;] );
# cRatioLst = exp.( [-5:1:-1;] );
# cRatioLst = exp.( [-3:1:-3;] );
# cFerroLst = [0:0.1:1;];
# cFerroLst = [0];
# cFerroRatioLst = exp.( [-1:1:5;] );
cFerroRatioLst = [0];
# sgnArea = 1;
# sgnPerim = -1;
sgnArea = -1;
sgnPerim = 1;
# betaLst = 1;

isInit0Lst = [false,true];

divNumLst = [64];

fNameLst = Vector{String}(undef,0);

# methodLoops = loops_MC_staggeredCube;
methodLoops = loops_MC_smart;

# for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
@time for itNum in itNumLst, cFerroRatio in cFerroRatioLst, cRatio in cRatioLst, beta in betaLst, divNum in divNumLst, isInit0 in isInit0Lst
	# fName = loops_MC( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta );
	cRatioSq = sqrt(cRatio);
	cArea = sgnArea * beta * cRatioSq;
	cPerimAbs = beta / cRatioSq;
	cPerim = sgnPerim * cPerimAbs;
	cFerro = cPerimAbs * cFerroRatio;
	fNameSmart = methodLoops( divNum, itNum; cPerim = cPerim, cArea = cArea, cFerro = cFerro, beta = betaBase, fMod = fModSmart, isInit0 = isInit0 );
	# @infiltrate
	# push!( fNameLst, fName );	
	push!( fNameLst, fNameSmart );
	GC.gc()
end

writedlm( "fileList.txt", fNameLst );
