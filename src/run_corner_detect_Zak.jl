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
harrisFiltData = CornerDetector.HarrisFiltHelperData( szZak);

thresHarris = 0.3;
thresSteer = 0.04;

CornerDetector.genCovMat!( harrisFiltData, zakLstTest );
CornerDetector.genHarrisCornerFiltFromCovMatBox!( harrisFiltData; filtLen = lnFiltHarris );
CornerDetector.genSteerCornerFilt!( steerData, zakLstTest, lnFiltSteer );

CornerDetector.nonMaxSuppress!( harrisFiltData );
CornerDetector.nonMaxSuppress!( steerData );

idCornerHarrisMergedLst = CornerDetector.extractMaxIdWithConnectedComp!( harrisFiltData, thresHarris, lnConnected );

idCornerSteerMergedLst = CornerDetector.extractMaxIdWithConnectedComp!( steerData, thresSteer, lnConnected );

idSteerHarrisMergedLst = CornerDetector.genCornerIdSteerHarrisMerged( idCornerSteerMergedLst, CornerDetector.getIsMaxedArr( harrisFiltData ), wdHarrisNearWind );

idSteerHarrisMergedLstStacked = zeros( Float64, 2, length(idSteerHarrisMergedLst) );
for ii = 1 : length(idSteerHarrisMergedLst)
	idSteerHarrisMergedLstStacked[:,ii] .= idSteerHarrisMergedLst[ii];
end

save( "idSteerHarrisMergedTest.jld2", "idSteerHarrisMergedLstTest", idSteerHarrisMergedLstStacked );
