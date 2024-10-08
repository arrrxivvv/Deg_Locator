using DelimitedFiles
using FilenameManip
using FlameGraphs
using FileIO
using Profile

# isFileNameOnly = false;
isFileNameOnly = true;

jlProfType = ".jlprof";

fMainDataProf = "dataProf";
fMainFlameProf = "flame";

# fNameDataProf = fMainDataProf * "_" * fAttrDataProf * jld2Type;
# fNameFlameProf = fMainFlameProf * "_" * fAttrDataProf * jld2Type;
fNameDataProf = fMainDataProf * "_" * fAttrDataProf * jlProfType;
fNameFlameProf = fMainFlameProf * "_" * fAttrDataProf * jlProfType;

if !isFileNameOnly
	# dataProf = Profile.fetch();
	Profile.clear();
	@profile Loops_MC.profilingWangLandau( ; divNum = 16, nDim = 3, D_hist = 2, dosIncrMin = 0.7 );
	dataProfRet, liDictRet = Profile.retrieve();
	# save( fNameDataProf, "dataProf", dataProf, "liDict", liDictRet );
	save( fNameDataProf, dataProfRet, liDictRet );
	# flameProf = FlameGraphs.flamegraph( dataProf );
	# save( fNameFlameProf, "flameProf", flameProf );
end

fMainFNameLstDataProf = "fNameLst";

fNameFNameLstDataProf = fMainFNameLstDataProf * "_" * fNameDataProf * txtType;

# fLst = [fNameDataProf, fNameFlameProf];
fLst = [fNameDataProf];
writedlm( fNameFNameLstDataProf, fLst );

open( Loops_MC.dirLog * Loops_MC.fNameFileLstLst, "w" ) do io
	println( io, fNameFNameLstDataProf );
end
