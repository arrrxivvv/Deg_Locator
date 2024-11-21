using RandomCircle
using CornerDetector
using SharedFNames
using Utils
using JLD2

fNameFNameArr = Utils.strReadLastLine( SharedFNames.dirLog * SharedFNames.fNameTmpNameFileLst );

fNameArr = load( fNameFNameArr, "fNameArr" );

fNameRandCircFine = fNameArr[3];

zakArrFineLst = load( fNameRandCircFine, "zakArrLst" );

iRTest = 3;
iNTest = 5;
it = 1;
zakArrTest = @view zakArrFineLst[iRTest,iNTest][:,:,it];
divNumTest = size( zakArrTest, 1 );

lnFiltHarris = 2;
lnFiltSteer = 2;
lnConnected = 1;
thresHarris = 0.3;
thresSteer = 0.04;
wdHarrisNearWind = 4;

steerHelperTest = CornerDetector.SteerFiltHelperData( divNumTest );

harrisHelperTest = CornerDetector.HarrisFiltHelperData( divNumTest );

# zakArrFloatTest = RandomCircle.boolToIntPosNeg.( zakArrTest );
zakArrFloatTest = zakArrTest;

CornerDetector.genCovMat!( harrisHelperTest, zakArrFloatTest );
CornerDetector.genHarrisCornerFiltFromCovMatBox!( harrisHelperTest; filtLen = lnFiltHarris );

CornerDetector.genSteerCornerFilt!( steerHelperTest, zakArrFloatTest, lnFiltSteer );

CornerDetector.nonMaxSuppress!( harrisHelperTest );
CornerDetector.nonMaxSuppress!( steerHelperTest );

idHarrisLst = CornerDetector.extractMaxIdWithConnectedComp!( harrisHelperTest, thresHarris, lnConnected );
idSteerLst = CornerDetector.extractMaxIdWithConnectedComp!( steerHelperTest, thresSteer, lnConnected );

idSteerHarrisMergedLst = CornerDetector.genCornerIdSteerHarrisMerged( idSteerLst, CornerDetector.getIsMaxedArr( harrisHelperTest ), wdHarrisNearWind );

numCorner = length( idSteerHarrisMergedLst );

jldsave( "randCircCorners.jld2"; idSteerHarrisMergedLst, zakArrTest, numCorner );
