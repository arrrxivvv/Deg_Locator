using DegLocatorDiv
using ThreadedArrays
using UMat

mSz = 10;
divLst = [80,16,16];
minNum = 0;
maxNum = 2*pi;
paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, nDim );
matsFull = matsGridHThreaded( paramsFull, threaded_zeros(ComplexF64,mSz,mSz) );
Hlst = DegLocatorDiv.HlstFunc(H_GUE,nDim,mSz);
HmatFun = (H,xLst) -> Hmat_3comb_ratio!( H, xLst, Hlst );
degBerrysFull = degBerrysInit( paramsFull, matsFull );

startNextEigen( matsFull );
eigenAll( matsFull; HmatFun = HmatFun );

DegLocatorDiv.linksCalcAll( degBerrysFull );
DegLocatorDiv.BfieldCalcAll( degBerrysFull );
DegLocatorDiv.divBCalcAll( degBerrysFull );