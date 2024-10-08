using CornerDetector
using ImgProcessing
using DataStructures
using JLD2

fName = "deg_GOE3_zakArr_dim_3_N_10_param_divide_[512, 512, 64]_instanceNum_10_seed_1000.jld2";

zakLstLst = load( fName, "zakArr" );

lnFiltHarris = 2;
lnFiltSteer = 2;
lnConnected = 1;
wdHarrisNearWind = 4;

itNum = 10;
mSz = 10;

mLev = 4;
it = 1;
zakLstTest = @view( zakLstLst[mLev,:,:,it] );

szZak = size( zakLstTest, 1 );

steerData = CornerDetector.SteerFiltHelperData( szZak );

thresHarris = 0.3;
thresSteer = 0.04;

diffLstTest = [ similar(zakLstTest) for ii = 1 : 2 ];
covMatTest = [ similar(zakLstTest) for ii = 1 : 2, jj = 1 : 2 ];
covMatFiltedTest = deepcopy( covMatTest );
CornerDetector.genCovMat!( covMatTest, diffLstTest, zakLstTest );

steerFiltTest = similar( zakLstTest, Float64 );
steerFiltedLstTest = [ similar(steerFiltTest) for ii = 1 : 2 ];
harrisFiltTest = similar( steerFiltTest );

CornerDetector.genHarrisCornerFiltFromCovMatBox!( harrisFiltTest, covMatFiltedTest, covMatTest; filtLen = lnFiltHarris );
# CornerDetector.genSteerCornerFilt!( steerFiltTest, steerFiltedLstTest, zakLstTest, lnFiltSteer );
CornerDetector.genSteerCornerFilt!( steerData, zakLstTest, lnFiltSteer );

maxedHarrisFiltTest = similar(harrisFiltTest);
# maxedSteerFiltTest = similar(steerFiltTest);

ImgProcessing.nonMaxSuppress!( maxedHarrisFiltTest, harrisFiltTest );
# ImgProcessing.nonMaxSuppress!( maxedSteerFiltTest, steerFiltTest );
CornerDetector.nonMaxSuppress!( steerData );

isMaxedHarrisArrTest = similar( harrisFiltTest, Bool );
isCornerHarrisArr = isMaxedHarrisArrTest;
# isMaxedSteerArrTest = similar( steerFiltTest, Bool );
isVisitedArr = similar( isMaxedHarrisArrTest );

# idCornerHarrisMergedLst, isCornerHarrisArr = CornerDetector.extractMaxIdWithConnectedComp!( isMaxedHarrisArrTest, isVisitedArr, maxedHarrisFiltTest, thresHarris, lnConnected );

# idCornerSteerMergedLst, isCornerSteerArr = CornerDetector.extractMaxIdWithConnectedComp!( isMaxedSteerArrTest, isVisitedArr, maxedSteerFiltTest, thresSteer, lnConnected );

idCornerHarrisMergedLst = CornerDetector.extractMaxIdWithConnectedComp!( isMaxedHarrisArrTest, isVisitedArr, maxedHarrisFiltTest, thresHarris, lnConnected );

# idCornerSteerMergedLst = CornerDetector.extractMaxIdWithConnectedComp!( isMaxedSteerArrTest, isVisitedArr, maxedSteerFiltTest, thresSteer, lnConnected );
idCornerSteerMergedLst = CornerDetector.extractMaxIdWithConnectedComp!( steerData, thresSteer, lnConnected );

idSteerHarrisMergedLst = CornerDetector.genCornerIdSteerHarrisMerged( idCornerSteerMergedLst, isCornerHarrisArr, wdHarrisNearWind );

idSteerHarrisMergedLstStacked = zeros( Float64, 2, length(idSteerHarrisMergedLst) );
for ii = 1 : length(idSteerHarrisMergedLst)
	idSteerHarrisMergedLstStacked[:,ii] .= idSteerHarrisMergedLst[ii];
end

save( "idSteerHarrisMergedTest.jld2", "idSteerHarrisMergedLstTest", idSteerHarrisMergedLstStacked );
