# Loops_MC module

abstract type AbstractWangLandauFlipChecker <: FlipChecker end

fNameFileLstWL = "fNameFileLstWL.txt";

function flipCheck( flipChecker::AbstractWangLandauFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	isFlip = flipCheckHistDos( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	
	if wangLandauHistResetCheck( flipChecker )
		wangLandauCheckDosIncr( flipChecker );
	end
	
	return isFlip;
end

function wangLandauHistResetCheck( flipChecker::AbstractWangLandauFlipChecker )
	error("WangLandau Flipper not defined")
end

function wangLandauCheckDosIncr( flipChecker::AbstractWangLandauFlipChecker )
	flipChecker.histArr .= 0
	flipChecker.dosIncrRef[] /= 2;
end

function flipCheck( flipChecker::AbstractWangLandauFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	dS = boolToFlipChange( BfieldLst[dim][pos] );
	dL = 0;
	for iLnkDim in 1 : params.nDimLayer
		dimLink = params.linkDimLst[dim][iLnkDim];
		dimLinkSh = params.linkDimShLst[dim][iLnkDim];
		dL += boolToFlipChange( linkLst[dimLink][pos] );
		dL += boolToFlipChange( linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] );
	end
	
	sTotal = sum( sum.( BfieldLst ) );
	lTotal = sum( sum.( linkLst ) );
	
	sTotalNxt = sTotal + dS;
	lTotalNxt = lTotal + dL;
	
	sRatio = sTotal / params.nDim / params.grdNum;
	lRatio = lTotal / params.nDim / params.grdNum;
	
	sRatioNxt = sTotalNxt / params.nDim / params.grdNum;
	lRatioNxt = lTotalNxt / params.nDim / params.grdNum;
	
	sId = ratioToBinId( sRatio, flipChecker.histDivNum );
	lId = ratioToBinId( lRatio, flipChecker.histDivNum );
	sIdNxt = ratioToBinId( sRatioNxt, flipChecker.histDivNum );
	lIdNxt = ratioToBinId( lRatioNxt, flipChecker.histDivNum );
	
	p = exp( flipChecker.dosArr[sId,lId] - flipChecker.dosArr[sIdNxt,lIdNxt] );
	isFlip = rand() < p;
	if isFlip
		flipChecker.histArr[sIdNxt,lIdNxt] += 1;
		flipChecker.dosArr[sIdNxt,lIdNxt] += flipChecker.dosIncrRef[];
	else
		flipChecker.histArr[sId,lId] += 1;
		flipChecker.dosArr[sId,lId] += flipChecker.dosIncrRef[];
	end
	
	if wangLandauHistResetCheck( flipChecker )
		wangLandauCheckDosIncr( flipChecker );
	end
	
	return isFlip;
end

struct WangLandauNoResetFlipChecker <: AbstractWangLandauFlipChecker
	histDivNum::Int64;
	histArr::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	
	function WangLandauNoResetFlipChecker( histDivNum::Int64; dosIncrVal::Float64 = 1.0 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		dosArr = zeros(Int64, histDivNum, histDivNum);
		
		new( histDivNum, histArr, dosArr, Ref(dosIncrVal) );
	end
end

function getFlipCheckerName( flipType::Type{WangLandauNoResetFlipChecker} )
	return "WangLandauNoResetFlip";
end

function getFlipCheckerAttrLst( flipChecker::WangLandauNoResetFlipChecker )
	return ["histDiv"];
end

function getFlipCheckerAttrValLst( flipChecker::WangLandauNoResetFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst(flipChecker);
	valLst = roundKeepInt.( [flipChecker.histDivNum]; digits = rndDigs );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WangLandauNoResetFlipChecker )
	return false;
end

struct WangLandauFlipChecker <: AbstractWangLandauFlipChecker
	histDivNum::Int64;
	histArr::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	histCoverThreshold::Float64;
	histLowRatioThreshold::Float64;
	
	function WangLandauFlipChecker( histDivNum::Int64; histLowRatioThreshold::Float64 = 0.8, histCoverThreshold::Float64 = 0.8, dosIncrVal::Float64 = 1.0 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		dosArr = zeros(Int64, histDivNum, histDivNum);
		
		new( histDivNum, histArr, dosArr, Ref(dosIncrVal), histCoverThreshold, histLowRatioThreshold );
	end
end

function getFlipCheckerName( flipType::Type{WangLandauFlipChecker} )
	return "WangLandauFlip";
end

function getFlipCheckerAttrLst( flipType::WangLandauFlipChecker )
	attrLst = ["histCoverThres", "histLowRatioThres"];
end

function getFlipCheckerAttrValLst( flipChecker::WangLandauFlipChecker )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = Any[ flipChecker.histCoverThreshold, flipChecker.histLowRatioThreshold ];
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WangLandauFlipChecker )
	isHistRest = false;
	coverNum = sum( x-> x>0, flipChecker.histArr );
	coverRatio = coverNum / flipChecker.histDivNum^2;
	if coverRatio > flipChecker.histCoverThreshold
		histAvg = sum( flipChecker.histArr ) / coverNum;
		histMin = minimum( x -> x == 0 ? Inf : x, flipChecker.histArr );
		isHistRest = histMin > flipChecker.histLowRatioThreshold * histAvg
	end
	
	return isHistRest;
end

function wangLandau_dSdL( flipChecker::WangLandauFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	dS = boolToFlipChange( BfieldLst[dim][pos] );
	dL = 0;
	for iLnkDim in 1 : params.nDimLayer
		dimLink = params.linkDimLst[dim][iLnkDim];
		dimLinkSh = params.linkDimShLst[dim][iLnkDim];
		dL += boolToFlipChange( linkLst[dimLink][pos] );
		dL += boolToFlipChange( linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] );
	end
	
	return dS, dL;
end

function wangLandau_flipCheck_dSdL( dS::Int64, dL::Int64, flipChecker::WangLandauFlipChecker, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	sTotal = sum( sum.( BfieldLst ) );
	lTotal = sum( sum.( linkLst ) );
	
	sTotalNxt = sTotal + dS;
	lTotalNxt = lTotal + dL;
	
	sRatio = sTotal / params.nDim / params.grdNum;
	lRatio = lTotal / params.nDim / params.grdNum;
	
	sRatioNxt = sTotalNxt / params.nDim / params.grdNum;
	lRatioNxt = lTotalNxt / params.nDim / params.grdNum;
	
	sId = ratioToBinId( sRatio, flipChecker.histDivNum );
	lId = ratioToBinId( lRatio, flipChecker.histDivNum );
	sIdNxt = ratioToBinId( sRatioNxt, flipChecker.histDivNum );
	lIdNxt = ratioToBinId( lRatioNxt, flipChecker.histDivNum );
	
	p = exp( flipChecker.dosArr[sId,lId] - flipChecker.dosArr[sIdNxt,lIdNxt] );
	isFlip = rand() < p;
	if isFlip
		flipChecker.histArr[sIdNxt,lIdNxt] += 1;
		flipChecker.dosArr[sIdNxt,lIdNxt] += flipChecker.dosIncrRef[];
	else
		flipChecker.histArr[sId,lId] += 1;
		flipChecker.dosArr[sId,lId] += flipChecker.dosIncrRef[];
	end
	
	return isFlip;
end

struct WangLandauStdFlipChecker <: AbstractWangLandauFlipChecker
	histDivNum::Int64;
	histArr::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	
	histMaxCntThres::Int64;
	histCntRatioCutoff::Float64;
	histStdRatioThres::Float64;
	histCoverCntThres::Int64;
	
	histIsOverThresArr::Array{Bool};
	histMaskedOverArr::Array{Int64};
	
	function WangLandauStdFlipChecker( histDivNum::Int64; dosIncrVal::Float64 = 1.0, histMaxCntThres::Int64 = 20, histCntRatioCutoff::Float64 = 0.5, histStdRatioThres::Float64 = 0.1, histCoverCntThres::Int64 = 5 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		dosArr = zeros(Int64, histDivNum, histDivNum);
		
		histIsOverThresArr = similar( histArr, Bool );
		histMaskedOverArr = similar( histArr );
		
		new( histDivNum, histArr, dosArr, Ref(dosIncrVal), histMaxCntThres, histCntRatioCutoff, histStdRatioThres, histCoverCntThres, histIsOverThresArr, histMaskedOverArr );
	end
end

function getFlipCheckerName( flipCheckerType::Type{WangLandauStdFlipChecker} )
	return "WangLandauStdFlip";
end

function getFlipCheckerAttrLst( flipChecker::WangLandauStdFlipChecker )
	return ["histDiv", "histMaxCntThres","histCntRatioCut","histStdRatio","histCoverCntThres"];
end

function getFlipCheckerAttrValLst( flipChecker::WangLandauStdFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	
	valLst = roundKeepInt.( [ flipChecker.histDivNum, flipChecker.histMaxCntThres, flipChecker.histCntRatioCutoff, flipChecker.histStdRatioThres, flipChecker.histCoverCntThres ]; digits = rndDigs );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WangLandauStdFlipChecker )
	isReset = false;
	
	maxCnt, idMax = findmax( flipChecker.histArr );
	# maxCnt = flipChecker.histArr[idMax];
	maxDos = flipChecker.dosArr[idMax];
	if maxCnt > flipChecker.histMaxCntThres
		cntThres = maxCnt * flipChecker.histCntRatioCutoff;
		dosThres = maxDos * flipChecker.histCntRatioCutoff;
		flipChecker.histIsOverThresArr .= flipChecker.histArr .> cntThres .|| flipChecker.dosArr .> dosThres;
		
		cntOverThres = sum( flipChecker.histIsOverThresArr );
		if cntOverThres > flipChecker.histCoverCntThres
			flipChecker.histMaskedOverArr .= flipChecker.histArr .* flipChecker.histIsOverThresArr;
			histMean = sum( flipChecker.histMaskedOverArr ) / cntOverThres;
			histStd = sqrt.( sum( flipChecker.histIsOverThresArr.* ( flipChecker.histMaskedOverArr .- histMean ).^2 ) / (cntOverThres - 1) );
			@infiltrate flipChecker.dosIncrRef[] < 1
			if histStd < histMean * flipChecker.histStdRatioThres
				isReset = true;
				@infiltrate
			end
		end
	end
	
	
	return isReset;
end

struct StaggeredCubeUpdaterWangLandau{N,Nplus1} <: LoopsUpdater
	posLstSh0::CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}};
	posLstAdvOrNot::Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}} }};
	posShOrNotLst::Vector{Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}}}};
	posStagCubeLst::Vector{Array{CartesianIndex{N},Nplus1}};
	idStagLst::CartesianIndices{Nplus1,NTuple{Nplus1,Base.OneTo{Int64}}};
	randIsFlipStagLst::Array{Bool};
	
	iDimLst::UnitRange{Int64};
	iIsShLst::UnitRange{Int64};
	randIDimLst::Array{Int64,Nplus1};
	randIShLst::Array{Int64,Nplus1};
	
	function StaggeredCubeUpdaterWangLandau( params::ParamsLoops; )
		posLstSh0 = ShiftedArrays.circshift( params.posLst, ntuple(x->0,params.nDim) );
		posLstAdvanced = [ ShiftedArrays.circshift( params.posLst, ntuple( ( i -> ( i == dim ? -1 : 0 ) ), params.nDim ) ) for dim = 1 : params.nDim ];
		posLstAdvOrNot = pushfirst!( posLstAdvanced, ShiftedArrays.circshift( params.posLst, ntuple(x->0,params.nDim) ) );
		posShOrNotLst = [ [ posLstSh0, posLstAdvanced[dim] ] for dim = 1 : params.nDim ];
		posShOrNotArr = [ (iSh == 1 ? posLstSh0 : posLstAdvanced[dim]) for dim = 1 : params.nDim, iSh = 1:2 ];
		
		coordsStagA = ntuple( iDim -> 1:2:params.divLst[iDim], params.nDim );
		coordsStagB = ntuple( iDim -> 2:2:params.divLst[iDim], params.nDim );
		posStagCubeLstA = @view( params.posLst[coordsStagA...] );
		posStagCubeLstB = @view( params.posLst[coordsStagB...] );
		posStagCubeLst = [ cat( @view( posLstAdvOrNot[iAdv][coordsStagA...]), @view(posLstAdvOrNot[iAdv][coordsStagB...]); dims = params.nDim+1 ) for iAdv = 1 : params.nDim+1 ];
		idStagLst = CartesianIndices( posStagCubeLst[1] );
		randIsFlipStagLst = similar( idStagLst, Bool );
		
		iDimLst = 1:params.nDim;
		iIsShLst = 1:2;
		randIDimLst = similar( posStagCubeLst[1], Int64 );
		randIShLst = similar(randIDimLst);
		
		new{params.nDim,params.nDim+1}( posLstSh0, posLstAdvOrNot, posShOrNotLst, posStagCubeLst, idStagLst, randIsFlipStagLst, iDimLst, iIsShLst, randIDimLst, randIShLst );
	end
end

function getUpdaterFMod( updaterType::Type{StaggeredCubeUpdaterWangLandau} )
	return "upStagCube";
end

function updateLoops( updater::StaggeredCubeUpdaterWangLandau, flipChecker::WangLandauFlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for iAdv = 1 : params.nDim+1
		rand!( updater.randIDimLst, updater.iDimLst );
		rand!( updater.randIShLst, updater.iIsShLst );
		rand!( updater.randIsFlipStagLst );
		dSSumRef = Threads.Atomic{Int}(0);
		dLSumRef = Threads.Atomic{Int}(0);
		Threads.@threads for idStag in updater.idStagLst
			posCube = updater.posStagCubeLst[iAdv][idStag];
			pos = updater.posShOrNotLst[updater.randIDimLst[idStag]][updater.randIShLst[idStag]][posCube];
			
			if updater.randIsFlipStagLst[idStag]
				dS, dL = wangLandau_dSdL( flipChecker, params, updater.randIDimLst[idStag], pos, BfieldLst, linkLst, linkFerroLst );
				Threads.atomic_add!( dSSumRef, dS );
				Threads.atomic_add!( dLSumRef, dL );
			end
		end
		if wangLandau_flipCheck_dSdL( dSSumRef[], dLSumRef[], flipChecker, params, BfieldLst, linkLst, linkFerroLst )
			Threads.@threads for idStag in updater.idStagLst
				if updater.randIsFlipStagLst[idStag]
					posCube = updater.posStagCubeLst[iAdv][idStag];
					pos = updater.posShOrNotLst[updater.randIDimLst[idStag]][updater.randIShLst[idStag]][posCube];
					
					flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = updater.randIDimLst[idStag] );
				end
			end
		end
	end
end


function loops_MC_methods_WangLandauStaggered( divNum = 64, itNum = 10000; histDivNum = 64, itNumSample = 100, itStartSample = 50, isInit0 = false, dosIncrInit = 1, cAreaInit = 0 )
	flipChecker = WangLandauFlipChecker( histDivNum; dosIncrVal = dosIncrInit );
	updaterType = StaggeredCubeUpdaterWangLandau;
	
	numBfieldLst, numLinkLst, BfieldSampleLst, linkSampleLst, BfieldStartSampleLst, linkStartSampleLst, zakMeanLst, zakLstSampleLst = loops_MC_methods_inJulia( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, itNumSample = itNumSample, itStartSample = itStartSample, isInit0 = isInit0, cAreaInit = cAreaInit );
	
	fMain = "loops_WL_staggered";
	attrLst = ["divNum", "itNum", "histDivNum", "cAreaInit", "dosIncrInit"];
	valLst = [divNum, itNum, histDivNum, cAreaInit, dosIncrInit];
	fNameWL = fNameFunc( fMain, attrLst, valLst, jld2Type );
	
	save( fNameWL, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	
	return fNameWL;
end

function loops_MC_methods_WangLandau( divNum = 64, itNum = 10000; histDivNum = 64, itNumSample = 100, itStartSample = 50, isInit0 = false, cAreaInit = 0, dosIncrInit = 1, nDim = 3 )
	# flipChecker = WangLandauFlipChecker( histDivNum; dosIncrVal = dosIncrInit );
	# flipChecker = WangLandauNoResetFlipChecker( histDivNum );
	flipChecker = WangLandauStdFlipChecker( histDivNum; histMaxCntThres = 10, histStdRatioThres = 0.15 );
	initializer = genMeanFieldInitializer( cAreaInit );
	updaterType = SingleUpdater;
	
	# numBfieldLst, numLinkLst, BfieldSampleLst, linkSampleLst, BfieldStartSampleLst, linkStartSampleLst, zakMeanLst, zakLstSampleLst = loops_MC_methods_inJulia( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, itNumSample = itNumSample, itStartSample = itStartSample, isInit0 = isInit0, cAreaInit = cAreaInit );
	fName = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, initializer = initializer, itNumSample = itNumSample, itStartSample = itStartSample );
	
	fMain = "loops_WL_single";
	attrLst, valLst = genAttrLstLttcFlipInit( divNum, itNum, nDim, flipChecker, initializer );
	# attrLst = ["divNum", "itNum", "histDivNum", "cAreaInit", "dosIncrInit"];
	# valLst = Any[divNum, itNum, histDivNum, cAreaInit, dosIncrInit];
	fNameWL = fNameFunc( fMain, attrLst, valLst, jld2Type );
	
	save( fNameWL, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	
	open( dirLog * fNameFileLstWL, "w" ) do io
		println( io, fNameWL );
	end
	
	return fNameWL;
end
