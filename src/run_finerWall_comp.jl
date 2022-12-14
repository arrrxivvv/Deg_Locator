using DegLocatorDiv
using Utils
using JLD2
using ShiftedArrays

function makeArrCell( arr, nDim, dimOffset, shDimLst, posLst, itNum, mSz )
	arrShLst = [ 
		ShiftedArrays.circshift( arr, ntuple( x -> x==dimOffset+iDim ? -1 : 0, ndims(arr) ) )
		for iDim = 1 : nDim];
	
	arrCell = [
		iBnd == 1 ? - arr[iDim, n, pos, it] : arrShLst[shDimLst[iDim]][iDim, n, pos, it]
		for iDim = 1:nDim, iBnd = 1:2, n = 1:mSz, pos in posLst, it = 1:itNum];
end

mSz = 10;
divLst = [80, 16, 16];
itNum = 10;
seed = 1000;
nDim = 3;
dimRevLst = [nDim:-1:1;];

fModMethod = "flux";
fModMethodFiner = "fluxCell";
fMain = "deg_detailed";

attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; dim = nDim );

fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModMethod );
fNameFiner = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModMethodFiner );

divBLst1 = load( fName, "divBLstLst" );
divBLstF = load( fNameFiner, "divBLstLst" );

NLstPol1 = load( fName, "NLstPol" );
NLstPolF = load( fNameFiner, "NLstPol" );

locLstPol1 = load( fName, "locLstPol" );

posLst = CartesianIndices( divBLst1[1] );

BfieldLst1 = load( fName, "BfieldLstLst" );
BfieldLstF = load( fNameFiner, "BfieldLstLst" );

dimOffset = 2;

BfieldArr1 = [
	BfieldLst1[it][iDim][pos][n]
	for iDim = 1 : nDim, n = 1:mSz, pos in posLst, it = 1:itNum];
BfieldArrF = [
	BfieldLstF[it][iDim][pos][n]
	for iDim = 1 : nDim, n = 1:mSz, pos in posLst, it = 1:itNum];
	
BfieldLstCell1 = makeArrCell( BfieldArr1, nDim, dimOffset, dimRevLst, posLst, itNum, mSz );
BfieldLstCellF = makeArrCell( BfieldArrF, nDim, dimOffset, dimRevLst, posLst, itNum, mSz );

divBLstDiff = zeros( mSz, divLst..., itNum );
for it = 1 : itNum
	for pos in posLst
		@view( divBLstDiff[:,pos,it] ) .= real.( divBLstF[it][pos] .- divBLst1[it][pos] );
	end
end

BfieldLstDiff = zeros( nDim, mSz, divLst..., itNum );
BfieldLstDiff .= real.( BfieldArrF .- BfieldArr1 );
# for it = 1:itNum, iDim = 1:nDim
	# for pos in posLst
		# @view( BfieldLstDiff[iDim,:, pos, it] ) .= real.( BfieldLstF[it][iDim][pos] .- BfieldLst1[it][iDim][pos] );
	# end
# end

BfieldLstDiffShLst = [ ShiftedArrays.circshift( BfieldLstDiff, ntuple( x-> x==2+iDim ? -1 : 0, 6 ) ) for iDim = 1 : nDim ];

# BfieldLstDiffCell = [
	# iBnd == 1 ? - BfieldLstDiff[iDim,n,pos,it] : BfieldLstDiffShLst[dimRevLst[iDim]][iDim,n,pos,it]
	# for iDim = 1:nDim, iBnd = 1:2, n = 1:mSz, pos in posLst, it =1:itNum];
BfieldLstDiffCell = makeArrCell( BfieldLstDiff, nDim, dimOffset, dimRevLst, posLst, itNum, mSz );

thres = 1e-6;

idSameLst = [
	findall( x->abs(x)<1e-6, @view( divBLstDiff[n, :,:,:, it] ) )
	for n = 1:mSz, it = 1:itNum];
	
BsameLst = [[
	iBnd == 1 ? - BfieldLstDiff[iDim, n, idSameLst[n,it][iLoc], it] : BfieldLstDiffShLst[dimRevLst[iDim]][iDim, n, idSameLst[n,it][iLoc], it]
	for iLoc = 1 : length(idSameLst[n, it]), iDim = 1:nDim, iBnd = 1:2 ]
	for n = 1:mSz, it = 1:itNum];