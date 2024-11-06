using CornerDetector
using ImgProcessing
using DataStructures
using JLD2
using SharedFNames
using Utils
using FilenameManip
using DelimitedFiles

using Infiltrator

isFileOnly = true;

lnFiltHarris = 2;
lnFiltSteer = 2;
lnConnected = 1;
wdHarrisNearWind = 4;

thresHarris = 0.3;
thresSteer = 0.04;

dirLog = SharedFNames.dirLog;
fNameFileLstJld2Lst = SharedFNames.fNameFileLstJld2Lst;
fNameSaveParamsLst = SharedFNames.fNameSaveParamsLst;

fNameFNameLstJld2 = Utils.strReadLastLine( dirLog * fNameFileLstJld2Lst );
fNameSaveParams = Utils.strReadLastLine( dirLog * fNameSaveParamsLst );

fNameArr = load( fNameFNameLstJld2, "fNameArr" );

nLst = load( fNameSaveParams, "nLst" );
itNum = load( fNameSaveParams, "itNum" );
divLstLst = load( fNameSaveParams, "divLstLst" );

lnNLst = length( nLst );

# itNumTmp = 1
lnNLstTmp = 8;

lnNLst = lnNLstTmp;
# itNum = itNumTmp;

idCornerSteerHarrisMergedLstNLst = [ Matrix{Matrix{Float64}}(undef,nLst[iN],itNum) for iN = 1 : lnNLst ];

fMainOut = "idCornerZakLst";
attrLst = ["mSzLst", "divNumBase"];
valLst = [nLst[[1,end]], divLstLst[1,1]];
fNameIdLstOut = fNameFunc( fMainOut, attrLst, valLst, jld2Type );

if !isFileOnly 

	for iN = 1 : lnNLst
		mSz = nLst[iN];
		
		zakLstLst = nothing;
		GC.gc();
		zakLstLst = load( fNameArr[iN], "zakArr" );
		divNum = divLstLst[1,iN];
		
		steerFiltData = CornerDetector.SteerFiltHelperData( divNum );
		harrisFiltData = CornerDetector.HarrisFiltHelperData( divNum );
		
		for it = 1 : itNum, iM = 1 : mSz
			zakLstCurrent = @view( zakLstLst[iM,:,:,it] );
			
			CornerDetector.genSteerCornerFilt!( steerFiltData, zakLstCurrent, lnFiltSteer );
			
			CornerDetector.genCovMat!( harrisFiltData, zakLstCurrent );
			CornerDetector.genHarrisCornerFiltFromCovMatBox!( harrisFiltData; filtLen = lnFiltHarris );
			
			CornerDetector.nonMaxSuppress!( steerFiltData );
			CornerDetector.nonMaxSuppress!( harrisFiltData );
			
			idCornerSteerMergedLst = CornerDetector.extractMaxIdWithConnectedComp!( steerFiltData, thresSteer, lnConnected );
			idCornerHarrisMergedLst = CornerDetector.extractMaxIdWithConnectedComp!( harrisFiltData, thresHarris, lnConnected );
			
			idCornerSteerHarrisMergedLstCurrent = CornerDetector.genCornerIdSteerHarrisMerged( idCornerSteerMergedLst, CornerDetector.getIsMaxedArr( harrisFiltData ), wdHarrisNearWind );
			
			idCornerSteerHarrisMergedLstNLst[iN][iM,it] = zeros( Float64, 2, length(idCornerSteerHarrisMergedLstCurrent) );
			for ii = 1 : length(idCornerSteerHarrisMergedLstCurrent)
				idCornerSteerHarrisMergedLstNLst[iN][iM,it][:,ii] .= idCornerSteerHarrisMergedLstCurrent[ii];
			end
		end
	end
	
	save( fNameIdLstOut, "idCornerZakLst", idCornerSteerHarrisMergedLstNLst, "mLst", nLst, "itNum", itNum, "divLstLst", divLstLst );
end

fNameLstTxt = fNameFunc( fMainOut, attrLst, valLst, txtType );

fNameLst = Vector{String}(undef,0);
push!( fNameLst, fNameIdLstOut );

writedlm( fNameLstTxt, fNameLst );

open( dirLog * SharedFNames.fNameTmpNameFileLst, "w" ) do io
	println( io, fNameLstTxt );
end
