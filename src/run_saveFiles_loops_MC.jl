using DelimitedFiles
using FilenameManip

itNumLst = [10000];

# cAreaLst = [0.5:0.5:1.5;];
cAreaLst = [1];
cPerimLst = -[0.5:0.5:5;];
# cPerimLst = [0];

fModSmart = "smart";

betaLst = [0.2,0.5,1,2,3];

divNumLst = [64];

fNameLst = Vector{String}(undef,0);

fMainBase = "loops_MC";
fMainSample = fMainBase * "_" * "zakLstSampleMean";

attrLstLoops = ["divNum","itNum","cArea","cPerim","beta"];

for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
	valLst = Any[divNum,itNum, cArea, cPerim,beta];
	fName = fNameFunc( fMainSample, attrLstLoops, valLst, jld2Type; fMod = fModSmart );
	push!( fNameLst, fName );
end

writedlm( "fileListSample.txt", fNameLst );
