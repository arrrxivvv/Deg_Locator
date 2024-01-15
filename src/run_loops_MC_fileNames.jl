using Loops_MC
# using Infiltrator
using DelimitedFiles
using FilenameManip

attrLstLoops = ["divNum","itNum","cArea","cPerim","beta"];
attrLstLoopsFerro = deepcopy(attrLstLoops);
attrLstLoopsFerro = push!( attrLstLoopsFerro, "cFerro" );

itNumLst = [30000];
# itNumLst = [1000];

# cAreaLst = [0.5:0.5:1.5;];
cAreaLst = [2:1:10;];
# cAreaLst = [1];
cPerimLst = -[2];
# cPerimLst = [0];

fModSmart = "smart";
fMod = fModSmart;

# betaLst = [0.2,0.5,1,2,3,10,100];
betaLst = [0.2,1,3,10];
# betaLst = [5:1:15;];
# betaLst = [10];
betaBase = 1;
cRatioLst = exp.( [-2:0.4:2;] );
# cRatioLst = exp.( [-5:1:-1;] );
# cRatioLst = exp.( [-3:1:-3;] );
# cRatioLst = exp.( [-1:1:-1;] );
# cFerroLst = [0:0.1:1;];
# cFerroLst = [0];
# cFerroRatioLst = exp.( [-1:1:5;] );
cFerroRatioLst = [0];
# sgnAreaLst = [1,-1];
sgnAreaLst = [1];
# sgnPerimLst = [1,-1];
sgnPerimLst = [-1];
# sgnArea = 1;
# sgnPerim = -1;
# sgnArea = -1;
# sgnPerim = 1;
# betaLst = 1;

# isInit0Lst = [false,true];
isInit0Lst = [false];

divNumLst = [64];

fNameLst = Vector{String}(undef,0);
fNameLoopsLst = similar(fNameLst);

# methodLoops = loops_MC_staggeredCube;
methodLoops = loops_MC_smart;
# methodLoops = Loops_MC.loops_MC_smart_old;
# methodLoops = loops_MC;

oFNameLoopsMain = "loopsSample";

# for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
@time for itNum in itNumLst, cFerroRatio in cFerroRatioLst, cRatio in cRatioLst, beta in betaLst, divNum in divNumLst, isInit0 in isInit0Lst, sgnArea in sgnAreaLst, sgnPerim in sgnPerimLst
	# fName = loops_MC( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta );
	cRatioSq = sqrt(cRatio);
	cArea = sgnArea * beta * cRatioSq;
	cPerimAbs = beta / cRatioSq;
	cPerim = sgnPerim * cPerimAbs;
	cFerro = cPerimAbs * cFerroRatio;
	# push!( fNameLst, fName );	
	
	attrLst = cFerro == 0 ? attrLstLoops : attrLstLoopsFerro;
	# ["divNum","itNum","cArea","cPerim","beta"];
	valLst = Any[divNum, itNum, cArea, cPerim, betaBase];
	if cFerro != 0
		push!( valLst, cFerro );
	end
	fModOut = fMod;
	if isInit0
		if !isempty(fMod)
			fModOut *= "_";
		end
		fModOut *= "isInit0";
	end
	
	oFNameLoops = fNameFunc( oFNameLoopsMain, attrLst, valLst, jld2Type; fMod = fModOut );
	push!( fNameLoopsLst, oFNameLoops );
	GC.gc()
end

writedlm( "fNameLoopsLst.txt", fNameLoopsLst );
