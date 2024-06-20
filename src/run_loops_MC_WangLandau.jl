using Loops_MC
using Utils
# using Infiltrator
using DelimitedFiles


# isFileNameOnly = false;
isFileNameOnly = true;


divNum = 8;
itNum = 640000;
nDim = 2;
D_hist = 2;

dosIncrMin = 0.1001;
wlResetInterval = 1000;
itExchange = 100;

cAreaInit = 0;
# wlHistDosType = Loops_MC.WLHistDos2DZoned;
# wlHistDosType = Loops_MC.WLHistDosJoint2DFull;
# wlHistDosType = Loops_MC.WLHistDosFull1dDos;
# wlHistDosType = Loops_MC.WLHistDosFull2dDos;
# wlHistDosType = Loops_MC.WL3dHistDosFull1dDos;
wlHistDosType = Loops_MC.WLHistDosFull{nDim,D_hist};

# wlHistDosArgs = (-1, 0.2);
wlHistDosArgs = ();
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

# @time fNameGotten = Loops_MC.loops_MC_methods_WL2d( divNum; cAreaInit = cAreaInit, dosIncrInit = 1.0, dosIncrMin = dosIncrMin, isFileNameOnly = isFileNameOnly, fMainOutside = fMainToLoad, wlHistDosType = wlHistDosType, wlHistDosArgs = wlHistDosArgs, wlResetInterval = wlResetInterval, nDim = nDim );
# push!(fNameLst, fNameGotten)

# @time fNameGotten = Loops_MC.loops_MC_methods_WL2dZoned( divNum; dosIncrInit = 1.0, dosIncrMin = dosIncrMin, isFileNameOnly = isFileNameOnly, fMainOutside = fMainToLoad );
# if !isFileNameOnly
	# push!(fNameLst, fNameGotten...)
# else
	# push!(fNameLst, fNameGotten)
# end

@time fNameGotten = Loops_MC.loops_MC_methods_WL2dReplica( divNum; dosIncrInit = 1.0, dosIncrMin = dosIncrMin, isFileNameOnly = isFileNameOnly, fMainOutside = fMainToLoad, itExchange = itExchange, D_hist = D_hist, histCutoffThres = histCutoffThres, numZones = numZones, numWalksEach = numWalksEach, EMinRatio = EMinRatio, EMaxRatio = EMaxRatio, nDim = nDim );
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

open(Loops_MC.dirLog * Loops_MC.fNameFileLstLst, "w") do io
	println(io, fNameGotten);
end
