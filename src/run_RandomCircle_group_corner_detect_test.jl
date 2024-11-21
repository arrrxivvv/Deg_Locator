using RandomCircle
using CornerDetector
using SharedFNames
using Utils
using JLD2
using FilenameManip

fNameFNameArr = Utils.strReadLastLine( SharedFNames.dirLog * SharedFNames.fNameTmpNameFileLst );

fNameArr = load( fNameFNameArr, "fNameArr" );

fNameRandCircFine = fNameArr[3];

zakArrFineLst = load( fNameRandCircFine, "zakArrLst" );

divNumNxtLst = load( fNameArr[1], "divNumNxt" );

nCircLst = load( fNameArr[1], "nCircLst" );
rCircLst = load( fNameArr[1], "rCircLst" );

# nCircLst = [5:5:50;];
# rCircLst = [0.1:0.1:0.5;];
lnNCirc = length( nCircLst );
lnRCirc = length( rCircLst );

lnFiltHarris = 2;
lnFiltSteer = 2;
lnConnected = 1;
thresHarris = 0.3;
thresSteer = 0.04;
wdHarrisNearWind = 4;

itNum = 10;
# lnNCirc = 10;
# lnRCirc = 5;

numCornerLst = zeros(Int64, itNum, lnRCirc, lnNCirc);
idCornerLstLst = [ Vector{Vector{Float64}}(undef,0) for it = 1 : itNum, iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

# for iN = 1 : lnNCirc, iR = 1 : lnRCirc
for iN = 1 : lnNCirc, iR = 1 : lnRCirc
	divNum = divNumNxtLst[iR,iN];
	GC.gc();
	steerHelperTest = CornerDetector.SteerFiltHelperData( divNum );
	harrisHelperTest = CornerDetector.HarrisFiltHelperData( divNum );
	
	# for it = 1 : itNum
	for it = 1 : itNum
		zakArrTest = @view zakArrFineLst[iR,iN][:,:,it];
	
		CornerDetector.genCovMat!( harrisHelperTest, zakArrTest );
		CornerDetector.genHarrisCornerFiltFromCovMatBox!( harrisHelperTest; filtLen = lnFiltHarris );
		CornerDetector.genSteerCornerFilt!( steerHelperTest, zakArrTest, lnFiltSteer );
		
		CornerDetector.nonMaxSuppress!( harrisHelperTest );
		CornerDetector.nonMaxSuppress!( steerHelperTest );
		
		idSteer = CornerDetector.extractMaxIdWithConnectedComp!( steerHelperTest, thresSteer, lnConnected );
		idHarris = CornerDetector.extractMaxIdWithConnectedComp!( harrisHelperTest, thresHarris, lnConnected );
		
		idSteerHarrisMergedLst = CornerDetector.genCornerIdSteerHarrisMerged( idSteer, CornerDetector.getIsMaxedArr( harrisHelperTest ), wdHarrisNearWind );
		
		numCorner = length( idSteerHarrisMergedLst );
		
		numCornerLst[it, iR, iN] = numCorner;
		idCornerLstLst[it, iR, iN] = idSteerHarrisMergedLst;
	end
end

jldsave( "randCircCorners.jld2"; idSteerHarrisMergedLst, zakArrTest, numCorner );

attrLstBase = [ "nCircLst", "rCircLst", "itNum" ];
valLstBase = Any[nCircLst[[1,end]], rCircLst[[1,end]], itNum];

fMainRandCircCorners = "randCircCorners";
fNameRandCircCorners = fNameFunc( fMainRandCircCorners, attrLstBase, valLstBase, jld2Type );

jldsave( fNameRandCircCorners; numCornerLst, idCornerLstLst );

fNameArr = [ fNameRandCircZak, fNameRandCircCorr, fNameRandCircFine, fNameRandCircCorners ];
# fNameArr = [ fNameRandCircCorners ];

# fMainFNameArr = "fNameArr" * "_" * fMainRandCircCorners;
# fNameFNameArr = fNameFunc( fMainFNameArr, attrLstBase, valLstBase, jld2Type );
jldsave( fNameFNameArr; fNameArr );

open( SharedFNames.dirLog * SharedFNames.fNameTmpNameFileLst, "w" ) do io
	println( io, fNameFNameArr );
end
