using Utils
using Statistics
using JLD2
using Logging; using LoggingExtras;

function divB_profile_new( mSz, divLst, itNum, seedFed; nDim = 3, enumSaveMem = memNone )
	if seedFed > 0
		Random.seed!(seedFed);
	end
	
	minNum = 0;
	maxNum = 2*pi;
	paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, nDim );
	tmpArrs = degTmpArrs( paramsFull, enumSaveMem );
	
	totalNumLst = zeros(Int64,itNum);
	HLstLst = Vector{Array{Array{ComplexF64}}}(undef,itNum);
	locLstPol = [
		Vector{Vector{Array{Int64}}}(undef,itNum)
		for iPol = 1:2];
	NLstPol = [
		zeros(Int64, itNum, mSz)
		for iPol = 1:2];
	
	for it = 1 : itNum
		print( "\rIteration: $it / $itNum         " )
		NLstPol[1][it,:], NLstPol[2][it,:], locLstPol[1][it], locLstPol[2][it], HLstLst[it] = locateDiv( tmpArrs );
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

function divB_profile_rootFind( mSz, divLst, itNum, seedFed; nDim = 3, thresVal = 1e-9, thresSz = 1e-9 )
	locFun( tmpArrs...; HmatFun = HmatFun ) = locateRootFind( tmpArrs...; HmatFun = HmatFun, thresVal = thresVal, thresSz = thresSz );
	
	tmpArrsFun = tmpArrsDegRootFind;
	
	divB_profile_base( mSz, divLst, itNum, seedFed; nDim = nDim, locFun = locFun, tmpArrsFun = tmpArrsFun, isOnlyBetween = true, locType = Float64 );
end

function divB_profile_flux(mSz, divLst, itNum, seedFed; nDim = 3, enumSaveMem = memNone)
	locFun = locateDiv;
	tmpArrsFun( paramsFull ) = degTmpArrs( paramsFull, enumSaveMem );
	
	divB_profile_base( mSz, divLst, itNum, seedFed; nDim = nDim, locFun = locFun, tmpArrsFun = tmpArrsFun );
end

function divB_profile_base( mSz, divLst, itNum, seedFed; nDim = 3, locFun, tmpArrsFun, isOnlyBetween = false, locType = Int64 )
	if seedFed > 0
		Random.seed!(seedFed);
	end
	
	minNum = 0;
	maxNum = 2*pi;
	paramsFull = degParamsInit( mSz, divLst, minNum, maxNum, nDim );
	tmpArrs = tmpArrsFun( paramsFull );
	
	nLevels = isOnlyBetween ? mSz-1 : mSz;
	totalNumLst = zeros(Int64,itNum);
	HLstLst = Vector{Array{Array{ComplexF64}}}(undef,itNum);
	locLstPol = [
		Vector{Vector{Array{locType}}}(undef,itNum)
		for iPol = 1:2];
	NLstPol = [
		zeros(Int64, itNum, nLevels)
		for iPol = 1:2];
	
	for it = 1 : itNum
		HLstLst[it] = DegLocatorDiv.HlstFunc(H_GUE,paramsFull.nDim,paramsFull.N);
	end
	
	for it = 1 : itNum
		print( "\rIteration: $it / $itNum         " )
		HmatFun = (H,xLst) -> Hmat_3comb_ratio!( H, xLst, HLstLst[it] );
		NLstPol[1][it,:], NLstPol[2][it,:], locLstPol[1][it], locLstPol[2][it] = locFun( tmpArrs...; HmatFun = HmatFun );
		# , HLstLst[it]
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
		
	@info("GC")
	Utils.@timeInfo GC.gc();
end
