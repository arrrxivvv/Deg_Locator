using Loops_MC
using Utils
# using Infiltrator
using DelimitedFiles


isFileNameOnly = false;
# isFileNameOnly = true;


divNum = 16;
itNum = 640000;

dosIncrMin = 0.1;
wlResetInterval = 2000;
itExchange = 200;

cAreaInit = 0;
# wlHistDosType = Loops_MC.WLHistDos2DZoned;
# wlHistDosType = Loops_MC.WLHistDosJoint2DFull;
# wlHistDosType = Loops_MC.WLHistDosFull1dDos;
wlHistDosType = Loops_MC.WLHistDosFull2dDos;

# wlHistDosArgs = (-1, 0.2);
wlHistDosArgs = ();
D_hist = 2;
histCutoffThres = 0.5;
numZones = 16;
numWalksEach = 3;
EMinRatio = -2.0;
EMaxRatio = 2.0;

fLstName = "Loops_MC_WangLandau.txt";
fLstName2d = "Loops_MC_WL2d.txt";
# fNameFileLstLstWL = "fNameFileLstLstWL.txt";

# fMainToLoad = Loops_MC.getAuxDataSummarySampleName( Loops_MC.WangLandauAuxData );
# fMainToLoad = Loops_MC.getAuxDataSummarySampleName( Loops_MC.WangLandauAuxData );
# fMainToLoad = Loops_MC.getAuxDataSummaryItSampleLstName( Loops_MC.WangLandauAuxData );
# fMainToLoad = "loops_WL2dZonesGluedThrd";
fMainToLoad = "loops_WL2dZonesGluedReplica";

fNameLst = Vector{String}(undef,0);

# @time fName = Loops_MC.loops_MC_methods_WangLandau( divNum, (divNum^3)^2 * 100; cAreaInit = 3, dosIncrInit = 1.0, histDivNum = divNum^2 );
# push!(fNameLst, fName)

# @time fName = Loops_MC.loops_MC_methods_WL2d( divNum, itNum; cAreaInit = 3, dosIncrInit = 1.0, isFileNameOnly = isFileNameOnly, fMainOutside = fMainToLoad );
# push!(fNameLst, fName)

# @time fNameGotten = Loops_MC.loops_MC_methods_WL2d( divNum; cAreaInit = cAreaInit, dosIncrInit = 1.0, dosIncrMin = dosIncrMin, isFileNameOnly = isFileNameOnly, fMainOutside = fMainToLoad, wlHistDosType = wlHistDosType, wlHistDosArgs = wlHistDosArgs, wlResetInterval = wlResetInterval );
# push!(fNameLst, fNameGotten)

# @time fNameGotten = Loops_MC.loops_MC_methods_WL2dZoned( divNum; dosIncrInit = 1.0, dosIncrMin = dosIncrMin, isFileNameOnly = isFileNameOnly, fMainOutside = fMainToLoad );
# if !isFileNameOnly
	# push!(fNameLst, fNameGotten...)
# else
	# push!(fNameLst, fNameGotten)
# end

@time fNameGotten = Loops_MC.loops_MC_methods_WL2dReplica( divNum; dosIncrInit = 1.0, dosIncrMin = dosIncrMin, isFileNameOnly = isFileNameOnly, fMainOutside = fMainToLoad, itExchange = itExchange, D_hist = D_hist, histCutoffThres = histCutoffThres, numZones = numZones, numWalksEach = numWalksEach, EMinRatio = EMinRatio, EMaxRatio = EMaxRatio );
if !isFileNameOnly
	push!(fNameLst, fNameGotten...)
else
	push!(fNameLst, fNameGotten)
end

# @time fName = Loops_MC.loops_MC_methods_WangLandauStaggered( divNum, (divNum^3)^2 * 100; cAreaInit = 3, dosIncrInit = 2.0, histDivNum = divNum^2 );
# push!(fNameLst, fName)

writedlm( fLstName, fNameLst );
# writedlm( fLstName2d, fNameLst );

open(Loops_MC.dirLog * Loops_MC.fNameFileLstWL, "w") do io
	println(io, fNameGotten);
end
