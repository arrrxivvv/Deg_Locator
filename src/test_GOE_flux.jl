using DegLocatorDiv
using Utils

mSz = 10;
itNum = 2;
seedFed = 1000;
# divLst = 4 .* [128,128,16];
divLst = [128,128,16];
nDim = 3;
	nDimLayer = nDim-1;
	div3 = divLst[nDim];
	HLstFun = DegLocatorDiv.H_GOE;
	divLstLayer = divLst[1:nDim-1];
	
	minNum = 0;
	maxNum = 2*pi;
	
	paramsFull = degParamsPeriodic( mSz, divLst, minNum, maxNum, nDim; isMesh = false );
	paramsLayer = degParamsPeriodic( mSz, divLstLayer, minNum, maxNum, nDim-1 );
	
	degBerrysLayer, non0Lst = DegLocatorDiv.degTmpArrs( paramsLayer, memNone );
	vLstPrev = deepcopy( degBerrysLayer.degMats.vLst );
	vLstFirst = deepcopy( vLstPrev );
	linkLst3 = deepcopy( degBerrysLayer.linkLst[1] );
	linkLstThrough = deepcopy(linkLst3);
	assignArrOfArrs!( linkLstThrough, 1 );
	zakLstLst = [ [ zeros(mSz) for pos in paramsLayer.posLst ] for it = 1 : itNum ];
	
	HLstLst::Array{Matrix{Float64},3} = DegLocatorDiv.HLstRandomGen( mSz, itNum, nDim, seedFed; HRandFun = DegLocatorDiv.H_GOE );
	
	H3sum = zeros( Float64, mSz, mSz );
	x3LstLst = [ [ paramsFull.gridLst[nDim][i3] ] for i3 = 1:div3 ];
	
	NLstPolLayer = [zeros(Int64, mSz) for iPol = 1:2];
	locLstPolLayer = [ Vector{Matrix{Int64}}(undef,0) for iPol = 1:2 ];
	
	NLst = zeros(Int64, mSz, div3, itNum);
	locLst = [ [ zeros(Int64, 0, nDim) for iM = 1:mSz ] for it = 1 : itNum, i3 = 1 : div3 ];
	
	iPol2d = 1;

i3 = 1;
it = 1;

DegLocatorDiv.Hmat_3comb!( H3sum, x3LstLst[i3], @view(HLstLst[:,nDim:nDim,it]) );
HmatFunLayer = (H, xLst2) -> DegLocatorDiv.Hmat_3comb_offset!( H, xLst2, @view(HLstLst[:,1:nDimLayer,it]), H3sum ); 

nDimLayer2::Int64 = nDim-1;
HLstLst2::Array{Matrix{Float64},3} = DegLocatorDiv.HLstRandomGen( mSz, itNum, nDim, seedFed; HRandFun = DegLocatorDiv.H_GOE );
DegLocatorDiv.Hmat_3comb!( H3sum, x3LstLst[i3], @view(HLstLst2[:,nDim:nDim,it]) );
H3sum2::Matrix{Float64} = zeros( Float64, mSz, mSz );
DegLocatorDiv.Hmat_3comb!( H3sum2, x3LstLst[i3], @view(HLstLst2[:,nDim:nDim,it]) );
HmatFunLayer2 = (H, xLst2) -> DegLocatorDiv.Hmat_3comb_offset!( H, xLst2, @view(HLstLst2[:,1:nDimLayer2,it]), H3sum2 ); 
HLstTest = @view(HLstLst2[:,1:nDimLayer2,it]);
HLstTest1 = @view(HLstLst2[:,1:nDimLayer,it]);
HmatFunLayer3 = (H, xLst2) -> DegLocatorDiv.Hmat_3comb_offset!( H, xLst2, HLstTest, H3sum2 ); 

HmatFunLayer12 = (H, xLst2) -> DegLocatorDiv.Hmat_3comb_offset!( H, xLst2, HLstTest, H3sum ); 
HmatFunLayer11 = (H, xLst2) -> DegLocatorDiv.Hmat_3comb_offset!( H, xLst2, HLstTest1, H3sum ); 