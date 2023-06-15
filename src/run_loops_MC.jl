using Loops_MC
# using Infiltrator
using DelimitedFiles

itNumLst = [10000];

cAreaLst = [0.5:0.5:1.5;];
cPerimLst = [0.5:0.5:3;];
# cPerimLst = [0];

fModSmart = "smart";

betaLst = 1;

divNumLst = [64];

fNameLst = Vector{String}(undef,0);

for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
	fName = loops_MC( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta );
	fNameSmart = loops_MC_smart( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta, fMod = fModSmart );
	# @infiltrate
	push!( fNameLst, fName );	
	push!( fNameLst, fNameSmart );
end

writedlm( "fileList.txt", fNameLst );