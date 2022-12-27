using Utils
using Statistics
using JLD2
using Logging; using LoggingExtras;

function HLstRandomGen( mSz::Int64, itNum::Int64, nDim::Int64, seedFed )
	if seedFed > 0
		Random.seed!(seedFed);
	end
	if nDim == 3
		HRandFun = H_GUE;
	elseif nDim == 2
		HRandFun = H_GOE;
	end
	HLstLst = [ HRandFun(mSz) 
		for iCos = 1:2, iDim = 1 : nDim, it = 1:itNum];
	return HLstLst;
end

function divB_profile_new( mSz, divLst, itNum, seedFed; nDim = 3, enumSaveMem = memNone )
	if seedFed > 0
		Random.seed!(seedFed);
	end
	
	minNum = 0;
	maxNum = 2*pi;
	paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, nDim );
	matsFull = matsGridHThreaded( paramsFull, threaded_zeros(ComplexF64,mSz,mSz) );
	matsGridInitAll( matsFull );
	degBerrysFull = degBerrysInit( paramsFull, matsFull; isFullInit = true );
	non0Arr = zeros(divLst...,mSz);
	
	HLstLst = Vector{Array{Array{ComplexF64}}}(undef,itNum);
	locLstPol = [
		Vector{Vector{Array{Int64}}}(undef,itNum)
		for iPol = 1:2];
	NLstPol = [
		zeros(Int64, itNum, mSz)
		for iPol = 1:2];
	
	if nDim == 3
		HRandFun = H_GUE;
	elseif nDim == 2
		HRandFun = H_GOE;
	end
	
	for it = 1 : itNum
		HLstLst[it] = DegLocatorDiv.HlstFunc(HRandFun,paramsFull.nDim,paramsFull.N);
	end
	
	for it = 1 : itNum
		print( "\rIteration: $it / $itNum         " )
		HmatFun = (H,xLst) -> Hmat_3comb!( H, xLst, HLstLst[it] );
		NLstPol[1][it,:], NLstPol[2][it,:], locLstPol[1][it], locLstPol[2][it] = locateDiv( degBerrysFull, non0Arr; HmatFun = HmatFun );
	end
	
	posLstAvg = mean(NLstPol[1]; dims = 1);
	posLstStd = std( NLstPol[1]; dims = 1 );
	posTotalLst = sum(NLstPol[1]; dims = 2);
	posTotalAvg = mean( posTotalLst );
	posTotalStd = std( posTotalLst );
	
	for n = 1 : mSz
		println( "$(n): $(posLstAvg[n]) +/- $(round(posLstStd[n]; digits=3))" );
	end
	println( "Total: $posTotalAvg +/- $(round(posTotalStd; digits=2))" );
		
	@info("GC")
	Utils.@timeInfo GC.gc();
end

function divB_profile_rootFind( mSz, divLst, itNum, seedFed; nDim = 3, thresVal = 1e-9, thresSz = 1e-9, thresRelaxRatio = 30, fMod = "" )
	locFun( tmpArrs...; HmatFun = HmatFun ) = locateRootFind( tmpArrs...; HmatFun = HmatFun, thresVal = thresVal, thresSz = thresSz, thresRelaxRatio = thresRelaxRatio );
	
	tmpArrsFun = tmpArrsDegRootFind;
	
	fModMethod = "rootFind";
	attrMoreLst = ["thresVal", "thresSz", "thresRelaxRatio"];
	valMoreLst = [thresVal, thresSz, thresRelaxRatio];
	
	divB_profile_base( mSz, divLst, itNum, seedFed; nDim = nDim, locFun = locFun, tmpArrsFun = tmpArrsFun, isOnlyBetween = true, locType = Float64, fMod = [fModMethod, fMod], attrMoreLst = attrMoreLst, valMoreLst = valMoreLst );
end

function divB_profile_flux(mSz, divLst, itNum, seedFed; nDim = 3, enumSaveMem = memNone, fMod = "")
	locFun = locateDiv;
	tmpArrsFun( paramsFull ) = degTmpArrs( paramsFull, enumSaveMem );
	fModMethod = "flux";
	
	divB_profile_base( mSz, divLst, itNum, seedFed; nDim = nDim, locFun = locFun, tmpArrsFun = tmpArrsFun, fMod = [fModMethod,fMod] );
end

function divB_profile_flux_cell( mSz, divLst, itNum, seedFed; ratFine = 4, nDim = 3, fMod = "" )
	locFun = locateDivCell;
	tmpArrsFun( paramsFull ) = degTmpArrsCell( paramsFull, ratFine );
	fModMethod = "fluxCell";
	
	divB_profile_base( mSz, divLst, itNum, seedFed; nDim = nDim, locFun = locFun, tmpArrsFun = tmpArrsFun, fMod = [fModMethod, fMod] );
end

function divB_profile_flux_cell_rerun( mSz, divLst, itNum, seedFed; ratFine = 4, nDim = 3, fMod = "" )
	fModMethodOrg = "flux";
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seedFed; dim = nDim );
	fName = fNameFunc( fDeg, attrLst, valLst, jld2Type; fMod = [fModMethodOrg, fMod] );
	HLstLst = load( fName, "HLstLst" );
	
	locFun = locateDivCell;
	tmpArrsFun( paramsFull ) = degTmpArrsCell( paramsFull, ratFine );
	fModMethod = "fluxCell";
	
	attrMoreLst = ["ratFine"];
	valMoreLst = [ratFine];
	
	divB_profile_base( mSz, divLst, itNum, seedFed, HLstLst; nDim = nDim, locFun = locFun, tmpArrsFun = tmpArrsFun, fMod = [fModMethod, fMod], attrMoreLst = attrMoreLst, valMoreLst = valMoreLst );
end

function divB_profile_base( mSz, divLst, itNum, seedFed; nDim = 3, fMod = "", attrMoreLst = [], valMoreLst = [], fExt = jld2Type, locFun, tmpArrsFun, isOnlyBetween = false, locType = Int64 )
	if seedFed > 0
		Random.seed!(seedFed);
	end
	if nDim == 3
		HRandFun = H_GUE;
	elseif nDim == 2
		HRandFun = H_GOE;
	end
	HLstLst = [ HRandFun(mSz) 
		for iCos = 1:2, iDim = 1 : nDim, it = 1:itNum];
	
	divB_profile_base( mSz, divLst, itNum, seedFed, HLstLst; nDim = nDim, fMod = fMod, attrMoreLst = attrMoreLst, valMoreLst = valMoreLst, fExt = fExt, locFun = locFun, tmpArrsFun = tmpArrsFun, isOnlyBetween = isOnlyBetween, locType = locType );
end

function divB_profile_base( mSz, divLst, itNum, seedFed, HLstLst; nDim = 3, fMod = "", attrMoreLst = [], valMoreLst = [], fExt = jld2Type, locFun, tmpArrsFun, isOnlyBetween = false, locType = Int64 )
	minNum = 0;
	maxNum = 2*pi;
	paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, nDim );
	tmpArrs = tmpArrsFun( paramsFull );
	
	nLevels = isOnlyBetween ? mSz-1 : mSz;
	totalNumLst = zeros(Int64,itNum);
	locLstPol = [
		Vector{Vector{Array{locType}}}(undef,itNum)
		for iPol = 1:2];
	NLstPol = [
		zeros(Int64, itNum, nLevels)
		for iPol = 1:2];
	
	for it = 1 : itNum
		print( "\rIteration: $it / $itNum         " )
		HLst = @view(HLstLst[:,:,it]);
		# HLst = DegLocatorDiv.HlstFunc(H_GUE,paramsFull.nDim,paramsFull.N);
		# @infiltrate
		HmatFun = (H,xLst) -> Hmat_3comb!( H, xLst, HLst );
		NLstPol[1][it,:], NLstPol[2][it,:], locLstPol[1][it], locLstPol[2][it] = locFun( tmpArrs...; HmatFun = HmatFun );
	end
	
	posLstAvg = mean(NLstPol[1]; dims = 1);
	posLstStd = std( NLstPol[1]; dims = 1 );
	posTotalLst = sum(NLstPol[1]; dims = 2);
	posTotalAvg = mean( posTotalLst );
	posTotalStd = std( posTotalLst );
	
	for n = 1 : nLevels
		println( "$(n): $(posLstAvg[n]) +/- $(round(posLstStd[n]; digits=3))" );
	end
	println( "Total: $posTotalAvg +/- $(round(posTotalStd; digits=2))" );
	
	fMain = "deg";
	attrLst, valLst = ttrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seedFed; dim = nDim, attrMoreLst = attrMoreLst, valMoreLst = valMoreLst );
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	save( fName, "N", mSz, "seed", seedFed, "NLstPol", NLstPol, "locLstPol", locLstPol, "HLstLst", HLstLst );
		
	@info("GC")
	Utils.@timeInfo GC.gc();
end

function divB_profile_flux_detailedOutput(mSz, divLst, itNum, seedFed; nDim = 3, enumSaveMem = memNone, fMod = "")
	locFun = locateDiv_detailedOutput;
	tmpArrsFun( paramsFull ) = degTmpArrs( paramsFull, enumSaveMem );
	fModMethod = "flux";
	
	divB_profile_base_detailedOutput( mSz, divLst, itNum, seedFed; nDim = nDim, locFun = locFun, tmpArrsFun = tmpArrsFun, fMod = [fModMethod,fMod] );
end

function divB_profile_flux_cell_rerun_detailedOutput( mSz, divLst, itNum, seedFed; ratFine = 4, nDim = 3, fMod = "" )
	fModMethodOrg = "flux";
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seedFed; dim = nDim );
	fName = fNameFunc( fDeg, attrLst, valLst, jld2Type; fMod = [fModMethodOrg, fMod] );
	HLstLst = load( fName, "HLstLst" );
	
	locFun = locateDivCell_detailedOutput;
	tmpArrsFun( paramsFull ) = degTmpArrsCell( paramsFull, ratFine );
	fModMethod = "fluxCell";
	attrMoreLst = ["ratFine"];
	valMoreLst = [ratFine];
	
	divB_profile_base_detailedOutput( mSz, divLst, itNum, seedFed, HLstLst; nDim = nDim, locFun = locFun, tmpArrsFun = tmpArrsFun, fMod = [fModMethod, fMod], attrMoreLst = attrMoreLst, valMoreLst = valMoreLst );
end

function divB_profile_base_detailedOutput( mSz, divLst, itNum, seedFed; nDim = 3, fMod = "", attrMoreLst = [], valMoreLst = [], fExt = jld2Type, locFun, tmpArrsFun, isOnlyBetween = false, locType = Int64 )
	HLstLst = HLstRandomGen( mSz, itNum, nDim, seedFed );
	
	return divB_profile_base_detailedOutput( mSz, divLst, itNum, seedFed, HLstLst; nDim = nDim, fMod = fMod, attrMoreLst = attrMoreLst, valMoreLst = valMoreLst, fExt = fExt, locFun = locFun, tmpArrsFun = tmpArrsFun, isOnlyBetween = isOnlyBetween, locType = locType );
end

function divB_profile_base_detailedOutput( mSz, divLst, itNum, seedFed, HLstLst; nDim = 3, fMod = "", attrMoreLst = [], valMoreLst = [], fExt = jld2Type, locFun, tmpArrsFun, isOnlyBetween = false, locType = Int64 )
	minNum = 0;
	maxNum = 2*pi;
	paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, nDim );
	tmpArrs = tmpArrsFun( paramsFull );
	
	nLevels = isOnlyBetween ? mSz-1 : mSz;
	totalNumLst = zeros(Int64,itNum);
	locLstPol = [
		Vector{Vector{Array{locType}}}(undef,itNum)
		for iPol = 1:2];
	NLstPol = [
		zeros(Int64, itNum, nLevels)
		for iPol = 1:2];
	BfieldLstLst = Vector{Vector{Array{Vector{ComplexF64},nDim}}}(undef,itNum);
	divBLstLst = Vector{Array{Vector{ComplexF64},nDim}}(undef,itNum);
	
	for it = 1 : itNum
		print( "\rIteration: $it / $itNum         " )
		HLst = @view(HLstLst[:,:,it]);
		# HLst = DegLocatorDiv.HlstFunc(H_GUE,paramsFull.nDim,paramsFull.N);
		# @infiltrate
		HmatFun = (H,xLst) -> Hmat_3comb!( H, xLst, HLst );
		NLstPol[1][it,:], NLstPol[2][it,:], locLstPol[1][it], locLstPol[2][it], BfieldLstTmp, divBLstTmp = locFun( tmpArrs...; HmatFun = HmatFun );
		BfieldLstLst[it] = deepcopy(BfieldLstTmp);
		divBLstLst[it] = deepcopy(divBLstTmp);
	end
	
	posLstAvg = mean(NLstPol[1]; dims = 1);
	posLstStd = std( NLstPol[1]; dims = 1 );
	posTotalLst = sum(NLstPol[1]; dims = 2);
	posTotalAvg = mean( posTotalLst );
	posTotalStd = std( posTotalLst );
	
	for n = 1 : nLevels
		println( "$(n): $(posLstAvg[n]) +/- $(round(posLstStd[n]; digits=3))" );
	end
	println( "Total: $posTotalAvg +/- $(round(posTotalStd; digits=2))" );
	
	fMain = fDeg * "_detailed";
	attrLst, valLst = ttrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seedFed; dim = nDim, attrMoreLst = attrMoreLst, valMoreLst = valMoreLst );
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	save( fName, "N", mSz, "seed", seedFed, "NLstPol", NLstPol, "locLstPol", locLstPol, "HLstLst", HLstLst, "BfieldLstLst", BfieldLstLst, "divBLstLst", divBLstLst );
		
	@info("GC")
	Utils.@timeInfo GC.gc();
end

function locRootFindRawProfile( mSz, divLst, itNum, seedFed; nDim = 3, fMod = "", fExt = jld2Type, thresVal = 1e-9, thresSz = 1e-9 )
	if seedFed > 0
		Random.seed!(seedFed);
	end
	
	minNum = 0;
	maxNum = 2*pi;
	paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, nDim );
	degMats, degSmplx, nmArrsThr = tmpArrsDegRootFind( paramsFull );
	
	nLevels = mSz-1;
	
	HLstLst = Vector{Array{Array{ComplexF64}}}(undef,itNum);
	locLstRawLst = zeros( degSmplx.params.nDim, degSmplx.lnSimpAll, degSmplx.params.divLst..., degSmplx.params.N-1, itNum ); 
	gapLstRawLst = zeros( degSmplx.lnSimpAll, degSmplx.params.divLst..., degSmplx.params.N-1, itNum ); 
	dLastLocs = ndims( locLstRawLst );
	dLastGaps = ndims( gapLstRawLst );
	
	for it = 1 : itNum
		HLstLst[it] = DegLocatorDiv.HlstFunc(H_GUE,paramsFull.nDim,paramsFull.N);
	end
	
	for it = 1 : itNum
		print( "\rIteration: $it / $itNum         " )
		HmatFun = (H,xLst) -> Hmat_3comb!( H, xLst, HLstLst[it] );
		
		locLstRaw, gapLstRaw = locateRootFindRaw( degMats, degSmplx, nmArrsThr; HmatFun = HmatFun, thresVal = thresVal, thresSz = thresSz );
		selectdim( locLstRawLst, dLastLocs, it ) .= locLstRaw;
		selectdim( gapLstRawLst, dLastGaps, it ) .= gapLstRaw;
		# @infiltrate
	end
	
	attrMoreLst = ["thresVal", "thresSz"];
	valMoreLst = [thresVal, thresSz];
	fMain = "locsRaw";
	attrLst, valLst = ttrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seedFed; dim = nDim, attrMoreLst = attrMoreLst, valMoreLst = valMoreLst );
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	save( fName, "itNum", itNum, "locLstRawLst", locLstRawLst, "gapLstRawLst", gapLstRawLst );
		
	@info("GC")
	Utils.@timeInfo GC.gc();
end
