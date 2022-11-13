using DegLocatorDiv
using ThreadedArrays
using UMat

mSz = 5;
minNum = 0;
maxNum = 1;
divNum = 2;
nDim = 3;
paramsCube2 = degParamsInit( mSz, fill(divNum,nDim), minNum, maxNum, nDim; isNonPeriodic = true );
matsCube2 = matsGridHThreaded( paramsCube2, threaded_zeros(ComplexF64,mSz,mSz) );
degBerrysCube2 = degBerrysInit( paramsCube2, matsCube2 );
Hlst = DegLocatorDiv.HlstFunc(H_GUE,nDim,mSz);
HmatFun = (H, xLst) -> Hmat_3comb_ratio!( H, xLst, Hlst );

divNum = 4;
paramsCube4 = degParamsInit( mSz, fill(divNum,nDim), minNum, maxNum, nDim; isNonPeriodic = true );
matsCube4 = matsGridHThreaded( paramsCube4, threaded_zeros(ComplexF64,mSz,mSz) );
degBerrysCube4 = degBerrysInit( paramsCube4, matsCube4 );

startNextEigen( matsCube2 );
eigenOnSurface( matsCube2; HmatFun = HmatFun );
DegLocatorDiv.linksCalcSurface( degBerrysCube2 );
DegLocatorDiv.BfieldCalcSurface( degBerrysCube2 );
DegLocatorDiv.divBCalcSurface( degBerrysCube2 );
divBSurfaceOutput( degBerrysCube2, HmatFun );

startNextEigen( matsCube4 );
matsGridTransferSurfaceDouble!( matsCube4, matsCube2 );
eigenOnSurface( matsCube4; HmatFun = HmatFun );
DegLocatorDiv.linksCalcSurface( degBerrysCube4 );
DegLocatorDiv.BfieldCalcSurface( degBerrysCube4 );
DegLocatorDiv.divBCalcSurface( degBerrysCube4 );

# mSz = 10;
# divLst = [80,16,16];
# minNum = 0;
# maxNum = 2*pi;
# paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, nDim );
# matsFull = matsGridHThreaded( paramsFull, threaded_zeros(ComplexF64,mSz,mSz) );
# Hlst = DegLocatorDiv.HlstFunc(H_GUE,nDim,mSz);
# HmatFun = (H,xLst) -> Hmat_3comb_ratio!( H, xLst, Hlst );
