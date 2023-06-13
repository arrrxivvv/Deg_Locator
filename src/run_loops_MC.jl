using Loops_MC
# using Infiltrator
using DelimitedFiles

itNumLst = [300];

cAreaLst = [0:0.5:8;];
# cPerimLst = [0:0.5:3;];
cPerimLst = [0];

betaLst = 1;

divNumLst = [16];

fNameLst = Vector{String}(undef,0);

for itNum in itNumLst, cArea in cAreaLst, cPerim in cPerimLst, beta in betaLst, divNum in divNumLst
	fName = loops_MC( divNum, itNum; cPerim = cPerim, cArea = cArea, beta = beta );
	# @infiltrate
	push!( fNameLst, fName );	
end

writedlm( "fileList.txt", fNameLst );