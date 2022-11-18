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
	if enumSaveMem == memNone
		matsFull = matsGridHThreaded( paramsFull, threaded_zeros(ComplexF64,mSz,mSz) );
		degBerrysFull = degBerrysInit( paramsFull, matsFull; isFullInit = true );
	elseif enumSaveMem >= memEig
		degBerrysFull = degBerrysEigLayered( paramsFull );
	end
	
	totalNumLst = zeros(Int64,itNum);
	HLstLst = Vector{Array{Array{ComplexF64}}}(undef,itNum);
	locLstPol = [
		Vector{Vector{Array{Int64}}}(undef,itNum)
		for iPol = 1:2];
	NLstPol = [
		zeros(Int64, itNum, mSz)
		for iPol = 1:2];
	# non0Arr = [
		# zeros(mSz)
		# for pos in paramsFull.posLst];
	non0Arr = zeros(paramsFull.divLst...,paramsFull.N);
	
	for it = 1 : itNum
		print( "\rIteration: $it / $itNum         " )
		NLstPol[1][it,:], NLstPol[2][it,:], locLstPol[1][it], locLstPol[2][it], HLstLst[it] = locateDiv( non0Arr, degBerrysFull );
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
