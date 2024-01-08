using DelimitedFiles
using FilenameManip
using Loops_MC

itNumLst = [30000];

betaLst = [0.2,0.5,1,2,3,10,100];
# betaLst = [5:1:15;];
# betaLst = [10];
betaBase = 1;
cRatioLst = exp.( [-2:0.2:2;] );
# cRatioLst = exp.( [-3:0.1:-1;] );
# cRatioLst = exp.( [-3:1:-3;] );
# cFerroLst = [0:0.1:1;];
# cFerroRatioLst = exp.( [-1:1:3;] );
cFerroRatioLst = [0];
sgnArea = -1;
sgnPerim = 1;
# betaLst = 1;

fModSmart = "smart";
# fModSmart = "smartstaggeredCube";
fModSmartInit0 = fModSmart * "_" * "isInit0";

# betaLst = [0.2,0.5,1,2,3];

isInit0Lst = [false,true];
divNumLst = [64];

fNameLstZak = Vector{String}(undef,0);
fNameLstLoopsSample = Vector{String}(undef,0);
fNameLstLoopsStart = Vector{String}(undef,0);

fMainBase = "loops_MC";
fMainZakSample = fMainBase * "_" * "zakLstSampleMean";
fMainLoopsSample = "loopsSample";
fMainLoopsStart = "loopsStartSample";

# attrLstLoops = ["divNum","itNum","cArea","cPerim","beta"];
# attrLstLoopsFerro = deepcopy(attrLstLoops);

# for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
for itNum in itNumLst, cFerroRatio in cFerroRatioLst, cRatio in cRatioLst, beta in betaLst, divNum in divNumLst, isInit0 in isInit0Lst
	cRatioSq = sqrt(cRatio);
	cArea = sgnArea * beta * cRatioSq;
	cPerimAbs = beta / cRatioSq;
	cPerim = sgnPerim * cPerimAbs;
	cFerro = cPerimAbs * cFerroRatio;
	fMod = isInit0 ? fModSmartInit0 : fModSmart ;
	attrLst = cFerro == 0 ? Loops_MC.attrLstLoops : Loops_MC.attrLstLoopsFerro;
	
	valLst = Any[divNum,itNum, cArea, cPerim,betaBase];
	if cFerro != 0
		push!( valLst, cFerro );
	end
	
	fNameZakSample = fNameFunc( fMainZakSample, attrLst, valLst, jld2Type; fMod = fMod );
	fNameLoopsSample = fNameFunc( fMainLoopsSample, attrLst, valLst, jld2Type; fMod = fMod );
	fNameLoopsStart = fNameFunc( fMainLoopsStart, attrLst, valLst, jld2Type; fMod = fMod );
	push!( fNameLstZak, fNameZakSample );
	push!( fNameLstLoopsSample, fNameLoopsSample );
	push!( fNameLstLoopsStart, fNameLoopsStart );
end

writedlm( "fileListZakSample.txt", fNameLstZak );
writedlm( "fileListLoopsSample.txt", fNameLstLoopsSample );
writedlm( "fileListLoopsStart.txt", fNameLstLoopsStart );

