using Utils
using Statistics
using JLD2
using Logging; using LoggingExtras;

function HLstRandomGen( mSz::Int64, itNum::Int64, nDim::Int64, seedFed=-1; HRandFun = nothing )
	if seedFed > 0
		Random.seed!(seedFed);
	end
	if isnothing(HRandFun)
		if nDim == 3
			HRandFun = H_GUE;
		elseif nDim == 2
			HRandFun = H_GOE;
		end
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

function divB_profile_flux(mSz, divLst, itNum, seedFed; nDim = 3, enumSaveMem = memNone, fMod = "", HRandFun = nothing, alpha = 1 )
	locFun = locateDiv;
	tmpArrsFun( paramsFull ) = degTmpArrs( paramsFull, enumSaveMem );
	fModMethod = "flux";
	
	if HRandFun == nothing
		if nDim == 3
			HRandFun = H_GUE;
		else
			HRandFun = H_GOE;
		end
	end
	
	attrMoreLst = [];
	valMoreLst = [];
	
	if alpha == 1
		HLstFun = ( (m) -> Hmat_3comb_HLstGen( m, nDim, HRandFun ) );
		HmatFun = Hmat_3comb!;
	else
		HLstFun = ( (m) -> Hmat_3comb_offset_HLstGen( m, nDim, HRandFun ) );
		HmatFun = (H, xLst, HLst) -> Hmat_3comb_offset!( H, xLst, HLst; c1 = sqrt(alpha/3), cOff = sqrt(1-alpha)  );
		push!( attrMoreLst, "alpha" );
		push!( valMoreLst, alpha );
	end
	
	HmatFunLst, HLstLst = HmatFunGen( mSz, itNum, HLstFun, HmatFun; seedFed = seedFed );
	
	divB_profile_base( mSz, divLst, itNum; nDim = nDim, locFun = locFun, tmpArrsFun = tmpArrsFun, fMod = [fModMethod,fMod], HmatFunLst = HmatFunLst, HLstLst = HLstLst, seedFed = seedFed, attrMoreLst = attrMoreLst, valMoreLst = valMoreLst );
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

function divB_profile_base( mSz, divLst, itNum; nDim = 3, fMod = "", attrMoreLst = [], valMoreLst = [], fExt = jld2Type, locFun, tmpArrsFun, HmatFunLst, isOnlyBetween = false, locType = Int64, HLstLst = nothing, seedFed = -1 )
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
		NLstPol[1][it,:], NLstPol[2][it,:], locLstPol[1][it], locLstPol[2][it] = locFun( tmpArrs...; HmatFun = HmatFunLst[it] );
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

function divB_profile_base_detailedOutput( mSz, divLst, itNum, seedFed, HLstLst; nDim = 3, fMod = "", attrMoreLst = [], valMoreLst = [], fExt = jld2Type, locFun, tmpArrsFun, isOnlyBetween = false, locType = Int64, alpha = 1 )
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

function divB_profile_GOE_layered( mSz, divLst, itNum, seedFed; fMod = "", fExt = jld2Type )
	nDim = 3;
	nDimLayer = nDim-1;
	div3 = divLst[nDim];
	HLstFun = H_GOE;
	divLstLayer = divLst[1:nDim-1];
	
	minNum = 0;
	maxNum = 2*pi;
	
	paramsFull = degParamsPeriodic( mSz, divLst, minNum, maxNum, nDim );
	paramsLayer = degParamsPeriodic( mSz, divLstLayer, minNum, maxNum, nDim-1 );
	
	degBerrysLayer, non0Lst = degTmpArrs( paramsLayer, memNone );
	vLstPrev = deepcopy( degBerrysLayer.degMats.vLst );
	vLstFirst = deepcopy( vLstPrev );
	linkLst3 = deepcopy( degBerrysLayer.linkLst[1] );
	linkLstThrough = deepcopy(linkLst3);
	assignArrOfArrs!( linkLstThrough, 1 );
	zakLstLst = [ [ zeros(mSz) for pos in paramsLayer.posLst ] for it = 1 : itNum ];
	
	HLstLst = HLstRandomGen( mSz, itNum, nDim, seedFed; HRandFun = H_GOE );
	
	H3sum = zeros( Float64, mSz, mSz );
	x3LstLst = [ [ paramsFull.gridLst[nDim][i3] ] for i3 = 1:div3 ];
	
	NLstPolLayer = [zeros(Int64, mSz) for iPol = 1:2];
	locLstPolLayer = [ Vector{Matrix{Int64}}(undef,0) for iPol = 1:2 ];
	
	NLst = zeros(Int64, mSz, div3, itNum);
	locLst = [ [ zeros(Int64, 0, nDim) for iM = 1:mSz ] for it = 1 : itNum, i3 = 1 : div3 ];
	
	iPol2d = 1;
	for it = 1 : itNum
		print( "\rIteration: $it / $itNum         " )
		assignArrOfArrs!(linkLstThrough,1);
		for i3 = 1 : div3
			Hmat_3comb!( H3sum, x3LstLst[i3], @view(HLstLst[:,nDim:nDim,it]) );
			HmatFunLayer = (H, xLst2) -> Hmat_3comb_offset!( H, xLst2, @view(HLstLst[:,1:nDimLayer,it]), H3sum );
			NLst[:,i3,it], NLstPolLayer[2], locLstPolLayer[1], locLstPolLayer[2] = locateDiv( degBerrysLayer, non0Lst; HmatFun = HmatFunLayer, yesGC = false );
			for iM = 1 : mSz
				locLst[it, i3][iM] = hcat( locLstPolLayer[iPol2d][iM], fill(i3, NLst[iM,i3,it]) );
			end
			
			if i3 == 1
				assignArrOfArrs!( vLstFirst, degBerrysLayer.degMats.vLst );
			else
				Threads.@threads for pos in paramsLayer.posLst
					dotEachCol!( linkLst3[pos], degBerrysLayer.degMats.vLst[pos], vLstPrev[pos] );
					linkLst3[pos] .= linkLst3[pos] ./ abs.(linkLst3[pos]);
					
					linkLstThrough[pos] .*= linkLst3[pos];
				end
			end
			if i3 == div3
				Threads.@threads for pos in paramsLayer.posLst
					dotEachCol!( linkLst3[pos], vLstFirst[pos], degBerrysLayer.degMats.vLst[pos] );
					linkLst3[pos] .= linkLst3[pos] ./ abs.(linkLst3[pos]);
					
					linkLstThrough[pos] .*= linkLst3[pos];
					zakLstLst[it][pos] .= real.( log.(linkLstThrough[pos]) ./ (pi*1im) );
				end
			else
				assignArrOfArrs!( vLstPrev, degBerrysLayer.degMats.vLst );
			end
		end
		@info("GC: ")
		Utils.@timeInfo GC.gc();
	end
	
	fMain = "deg_GOE3";
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seedFed; dim = nDim );
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	save( fName, "N", mSz, "seed", seedFed, "NLst", NLst, "locLst", locLst, "HLstLst", HLstLst, "zakLstLst", zakLstLst );
end

function divB_profile_GOE_layered_base( HLstLst::Array{Matrix{Float64}}, mSz::Int64, divLst::Vector{Int64} )
	nDim = 3;
	nDimLayer = nDim - 1;
	divLstLayer = divLst[1:nDim-1];
	div3 = divLst[nDim];
	iPol2d = 1;
	
	minNum = 0;
	maxNum = 2*pi;
	
	paramsFull = degParamsPeriodic( mSz, divLst, minNum, maxNum, nDim );
	paramsLayer = degParamsPeriodic( mSz, divLstLayer, minNum, maxNum, nDim-1 );
	
	degBerrysLayer, non0Lst = degTmpArrs( paramsLayer, memNone );
	vLstPrev = deepcopy( degBerrysLayer.degMats.vLst );
	vLstFirst = deepcopy( vLstPrev );
	linkLst3 = deepcopy( degBerrysLayer.linkLst[1] );
	linkLstThrough = deepcopy(linkLst3);
	assignArrOfArrs!( linkLstThrough, 1 );
	zakLstLst = [ zeros(mSz) for pos in paramsLayer.posLst ];
	
	assignArrOfArrs!(linkLstThrough,1);
	
	H3sum = zeros( Float64, mSz, mSz );
	x3LstLst = [ [ paramsFull.gridLst[nDim][i3] ] for i3 = 1:div3 ];
	
	NLstPolLayer = [zeros(Int64, mSz) for iPol = 1:2];
	locLstPolLayer = [ Vector{Matrix{Int64}}(undef,0) for iPol = 1:2 ];
	
	NLst = zeros(Int64, mSz, div3);
	locLst = [ [ zeros(Int64, 0, nDim) for iM = 1:mSz ] for i3 = 1 : div3 ];	
	
	for i3 = 1 : div3
		Hmat_3comb!( H3sum, x3LstLst[i3], @view(HLstLst[:,nDim:nDim]) );
		HmatFunLayer = (H, xLst2) -> Hmat_3comb_offset!( H, xLst2, @view(HLstLst[:,1:nDimLayer]), H3sum );
		NLst[:,i3], NLstPolLayer[2], locLstPolLayer[1], locLstPolLayer[2] = locateDiv( degBerrysLayer, non0Lst; HmatFun = HmatFunLayer, yesGC = false );
		for iM = 1 : mSz
			locLst[i3][iM] = hcat( locLstPolLayer[iPol2d][iM], fill(i3, NLst[iM,i3]) );
		end
		
		if i3 == 1
			assignArrOfArrs!( vLstFirst, degBerrysLayer.degMats.vLst );
		else
			Threads.@threads for pos in paramsLayer.posLst
				dotEachCol!( linkLst3[pos], degBerrysLayer.degMats.vLst[pos], vLstPrev[pos] );
				linkLst3[pos] .= linkLst3[pos] ./ abs.(linkLst3[pos]);
				
				linkLstThrough[pos] .*= linkLst3[pos];
			end
		end
		if i3 == div3
			Threads.@threads for pos in paramsLayer.posLst
				dotEachCol!( linkLst3[pos], vLstFirst[pos], degBerrysLayer.degMats.vLst[pos] );
				linkLst3[pos] .= linkLst3[pos] ./ abs.(linkLst3[pos]);
				
				linkLstThrough[pos] .*= linkLst3[pos];
				zakLstLst[pos] .= real.( log.(linkLstThrough[pos]) ./ (pi*1im) );
			end
		else
			assignArrOfArrs!( vLstPrev, degBerrysLayer.degMats.vLst );
		end
	end	
	
	return zakLstLst;
end

function divB_profile_GOE_layered_div3Lst( mSz::Int64, divNum::Int64, div3Lst::Vector{Int64}, seedFed; fMod = "" )
	itNumDummy = 1;
	nDim = 3;
	HLstLst = HLstRandomGen( mSz, itNumDummy, nDim, seedFed; HRandFun = H_GOE )[:,:,1];
	divLst = fill(divNum, nDim);
	
	zakLstLst = [ [ zeros(mSz) for i1 = 1 : divNum, i2 = 1 : divNum ] for i3 = 1 : length(divLst) ];
	
	for iDiv3 = 1 : length(div3Lst)
		divLst[nDim] = div3Lst[iDiv3];
		
		zakLstLst[iDiv3] = divB_profile_GOE_layered_base( HLstLst, mSz, divLst );
	end
	
	fMain = "deg_GOE3_div3Test";
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst[1:2], itNumDummy, seedFed; dim = nDim );
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fMod );
	
	save( fName, "zakLstLst", zakLstLst, "div3Lst", div3Lst );
end

function deg_GOE3_zak_resave( mSz, divLst, itNum, seedFed; nDim = 3, fMod = "", fExt = jld2Type )
	fMain = "deg_GOE3";
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seedFed; dim = nDim );
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	zakLstLst = load( fName, "zakLstLst" );
	
	fMainOut = "zakLstLst";
	fNameOut = fNameFunc( fMainOut, attrLst, valLst, fExt; fMod = fMod );
	
	save( fNameOut, "zakLstLst", zakLstLst );
end

function deg_GOE3_zak_lstToArr_resave( mSz, divLst, itNum, seedFed; nDim = 3, fMod = "", fExt = jld2Type, fMain = "deg_GOE3", itNumAlt = 0 )
	# fMain = "deg_GOE3";
	if itNumAlt != 0
		itNumTrue = itNumAlt;
	else
		itNumTrue = itNum;
	end
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seedFed; dim = nDim );
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	zakLstLst = load( fName, "zakLstLst" );
	
	zakArr = zeros( mSz, size(zakLstLst[1])..., itNumTrue );
	
	idLst = CartesianIndices( zakLstLst[1] );
	
	for it = 1 : itNumTrue
		for idCart in idLst
			zakArr[:,idCart,it] .= zakLstLst[it][idCart];
		end
	end
	
	fMainOut = "deg_GOE3_zakArr";
	fNameOut = fNameFunc( fMainOut, attrLst, valLst, fExt; fMod = fMod );
	
	save( fNameOut, "zakArr", zakArr );
end

# function divB_profile_GOE_layered( mSz, divLst, itNum, seedFed; fMod = "", fExt = jld2Type )
	# nDim = 3;
	# nDimLayer = nDim-1;
	# div3 = divLst[nDim];
	# HLstFun = H_GOE;
	# divLstLayer = divLst[1:nDim-1];
	
	# minNum = 0;
	# maxNum = 2*pi;
	
	# paramsFull = degParamsPeriodic( mSz, divLst, minNum, maxNum, nDim );
	# paramsLayer = degParamsPeriodic( mSz, divLstLayer, minNum, maxNum, nDim-1 );
	
	# degBerrysLayer, non0Lst = degTmpArrs( paramsLayer, memNone );
	# vLstPrev = deepcopy( degBerrysLayer.degMats.vLst );
	# vLstFirst = deepcopy( vLstPrev );
	# linkLst3 = deepcopy( degBerrysLayer.linkLst[1] );
	# linkLstThrough = deepcopy(linkLst3);
	# assignArrOfArrs!( linkLstThrough, 1 );
	# zakLstLst = [ [ zeros(mSz) for pos in paramsLayer.posLst ] for it = 1 : itNum ];
	# # for it = 1 : itNum
		# # assignArrOfArrs!( zakLstLst[it], 1 );
	# # end
	
	# HLstLst = HLstRandomGen( mSz, itNum, nDim, seedFed; HRandFun = H_GOE );
	
	# H3sum = zeros( Float64, mSz, mSz );
	# x3LstLst = [ [ paramsFull.gridLst[nDim][i3] ] for i3 = 1:div3 ];
	
	# NLstPolLayer = [zeros(Int64, mSz) for iPol = 1:2];
	# locLstPolLayer = [ Vector{Matrix{Int64}}(undef,0) for iPol = 1:2 ];
	
	# NLst = zeros(Int64, mSz, div3, itNum);
	# locLst = [ [ zeros(Int64, 0, nDim) for iM = 1:mSz ] for it = 1 : itNum ];
	
	# iPol2d = 1;
	# for it = 1 : itNum
		# print( "\rIteration: $it / $itNum         " )
		# for i3 = 1 : div3
			# Hmat_3comb!( H3sum, x3LstLst[i3], @view(HLstLst[:,nDim:nDim,it]) );
			# HmatFunLayer = (H, xLst2) -> Hmat_3comb_offset!( H, xLst2, @view(HLstLst[:,1:nDimLayer,it]), H3sum );
			# NLst[:,i3,it], NLstPolLayer[2], locLstPolLayer[1], locLstPolLayer[2] = locateDiv( degBerrysLayer, non0Lst; HmatFun = HmatFunLayer, yesGC = false );
			# for iM = 1 : mSz
				# locLst[it][iM] = vcat( locLst[it][iM], hcat( locLstPolLayer[iPol2d][iM], fill(i3, NLst[iM,i3,it]) ) );
			# end
			# # @infiltrate
			
			# if i3 == 1
				# assignArrOfArrs!( vLstFirst, degBerrysLayer.degMats.vLst );
			# else
				# Threads.@threads for pos in paramsLayer.posLst
					# dotEachCol!( linkLst3[pos], degBerrysLayer.degMats.vLst[pos], vLstPrev[pos] );
					# linkLst3[pos] .= linkLst3[pos] ./ abs.(linkLst3[pos]);
					
					# linkLstThrough[pos] .*= linkLst3[pos];
				# end
			# end
			# if i3 == div3
				# Threads.@threads for pos in paramsLayer.posLst
					# dotEachCol!( linkLst3[pos], vLstFirst[pos], degBerrysLayer.degMats.vLst[pos] );
					# linkLst3[pos] .= linkLst3[pos] ./ abs.(linkLst3[pos]);
					
					# linkLstThrough[pos] .*= linkLst3[pos];
					# zakLstLst[it][pos] .= real.( log.(linkLstThrough[pos]) ./ (2*pi*1im) );
				# end
			# else
				# assignArrOfArrs!( vLstPrev, degBerrysLayer.degMats.vLst );
			# end
		# end
		# # @infiltrate
		# @info("GC: ")
		# Utils.@timeInfo GC.gc();
	# end
	
	# fMain = "deg_GOE3";
	# attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seedFed; dim = nDim );
	# fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	# save( fName, "N", mSz, "seed", seedFed, "NLst", NLst, "locLst", locLst, "HLstLst", HLstLst, "zakLstLst", zakLstLst );
# end

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

function HmatFunGen( mSz::Int64, itNum::Int64, HLstFun, HmatFun; isNewGen = true, HLstLst = nothing, seedFed = -1 )
	if isNewGen
		if seedFed > 0
			Random.seed!(seedFed);
		end
		HLstLst = [ 
			HLstFun( mSz ) 
			for it = 1 : itNum ];
	end
	HmatFunLst = [
		( H, xLst ) -> HmatFun( H, xLst, HLstLst[it] )
		for it = 1 : itNum ];
	return HmatFunLst, HLstLst;
end
