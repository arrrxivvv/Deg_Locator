using Loops_MC
using Utils
# using Infiltrator
using DelimitedFiles

divNum = 4;

fLstName = "Loops_MC_WangLandau.txt";

fNameLst = Vector{String}(undef,0);

@time fName = Loops_MC.loops_MC_methods_WangLandau( divNum, (divNum^3)^2 * 10; cAreaInit = 3, dosIncrInit = 2.0, histDivNum = divNum^2 );
push!(fNameLst, fName)

# @time fName = Loops_MC.loops_MC_methods_WangLandauStaggered( divNum, (divNum^3)^2 * 100; cAreaInit = 3, dosIncrInit = 2.0, histDivNum = divNum^2 );
# push!(fNameLst, fName)

writedlm( fLstName, fNameLst );

open(Loops_MC.dirLog * Loops_MC.fNameFileLstLst, "w") do io
	println(io, fLstName);
end
