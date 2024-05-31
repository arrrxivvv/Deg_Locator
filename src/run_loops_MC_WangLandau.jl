using Loops_MC
using Utils
# using Infiltrator
using DelimitedFiles


isFileNameOnly = false;


divNum = 8;
itNum = 640000;

fLstName = "Loops_MC_WangLandau.txt";
fLstName2d = "Loops_MC_WL2d.txt";
# fNameFileLstLstWL = "fNameFileLstLstWL.txt";

fMainToLoad = Loops_MC.getAuxDataSummarySampleName( Loops_MC.WangLandauAuxData );

fNameLst = Vector{String}(undef,0);

# @time fName = Loops_MC.loops_MC_methods_WangLandau( divNum, (divNum^3)^2 * 100; cAreaInit = 3, dosIncrInit = 1.0, histDivNum = divNum^2 );
# push!(fNameLst, fName)

@time fName = Loops_MC.loops_MC_methods_WL2d( divNum, itNum; cAreaInit = 3, dosIncrInit = 1.0, isFileNameOnly = isFileNameOnly, fMainOutside = fMainToLoad );
push!(fNameLst, fName)

# @time fName = Loops_MC.loops_MC_methods_WangLandauStaggered( divNum, (divNum^3)^2 * 100; cAreaInit = 3, dosIncrInit = 2.0, histDivNum = divNum^2 );
# push!(fNameLst, fName)

writedlm( fLstName, fNameLst );
# writedlm( fLstName2d, fNameLst );

open(Loops_MC.dirLog * Loops_MC.fNameFileLstWL, "w") do io
	println(io, fName);
end
