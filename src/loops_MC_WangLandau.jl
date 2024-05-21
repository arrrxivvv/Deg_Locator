# Loops_MC module

function calcDL( params::ParamsLoops, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D}
	dL = 0;
	for iLnkDim in 1 : params.nDimLayer
		dimLink = params.linkDimLst[dim][iLnkDim];
		dimLinkSh = params.linkDimShLst[dim][iLnkDim];
		dL += boolToFlipChange( linkLst[dimLink][pos] );
		dL += boolToFlipChange( linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] );
	end
	
	return dL;
end

abstract type AbstractWangLandauFlipChecker <: FlipChecker end

fNameFileLstWL = "fNameFileLstWL.txt";

function flipCheck( flipChecker::AbstractWangLandauFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	isFlip = wangLandauUpdateHistDos( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	
	if wangLandauHistResetCheck( flipChecker )
		wangLandauUpdateDosIncr( flipChecker );
	end
	
	return isFlip;
end

function wangLandauUpdateHistDos( flipChecker::AbstractWangLandauFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} )where {D}
	error("WangLandau Flipper not defined")
end

function wangLandauHistResetCheck( flipChecker::AbstractWangLandauFlipChecker )
	error("WangLandau Flipper not defined")
end

function wangLandauUpdateDosIncr( flipChecker::AbstractWangLandauFlipChecker )
	flipChecker.histArr .= 0
	flipChecker.dosIncrRef[] /= 2;
	# @infiltrate
end

abstract type AbstractWangLandau2dLinkOnlyFlipChecker <: AbstractWangLandauFlipChecker end

function getLinkHistId( flipChecker::AbstractWangLandau2dLinkOnlyFlipChecker, lnkVal::Int64 )
	# lnkId = Int64( ( lnkVal + 2*flipChecker.grdNum ) / 4 );
	lnkId = Int64( lnkVal / 2 );
	if lnkId == flipChecker.grdNum
		lnkId -= 2;
	elseif lnkId > 0
		lnkId -= 1;
	end
	
	lnkId += 1;
	
	return lnkId;
end

function wangLandauUpdateHistDos( flipChecker::AbstractWangLandau2dLinkOnlyFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} )where {D}
	dL = calcDL( params, linkLst, dim, pos );
	
	# lTotal = boolToOnePN( sum( sum.(linkLst) ); cnt = 2*flipChecker.grdNum );
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
	# lTotal = boolToOnePN( lTotalCnt; cnt = 2*flipChecker.grdNum );
	# lNxt = boolToOnePN( lNxtCnt; cnt = 2*flipChecker.grdNum );
	
	lId = getLinkHistId( flipChecker, lTotal );
	lIdNxt = getLinkHistId( flipChecker, lNxt );
	
	p = exp( flipChecker.dosArr[lId] - flipChecker.dosArr[lIdNxt] );
	isFlip = rand() < p;
	
	if isFlip
		flipChecker.histArr[lIdNxt] += 1;
		flipChecker.dosArr[lIdNxt] += flipChecker.dosIncrRef[];
	else
		flipChecker.histArr[lId] += 1;
		flipChecker.dosArr[lId] += flipChecker.dosIncrRef[];
	end
	
	return isFlip;
end

struct WL2dLinkFlatFlipChecker <: AbstractWangLandau2dLinkOnlyFlipChecker
	divNum::Int64;
	grdNum::Int64;
	
	histArr::Vector{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	
	histMinRatioThres::Float64;
	
	wlCounter::Ref{Int64};
	wlResetInterval::Int64;
	
	function WL2dLinkFlatFlipChecker( divNum::Int64, nDim::Int64; dosIncrVal::Float64 = 1.0, histMinRatioThres = 0.8, wlResetInterval::Int64 = 100 )
		grdNum = divNum^nDim;
		histArr = zeros( Int64, grdNum+1-2 );
		dosArr = similar(histArr);
		dosArr .= 0;
		wlCounterVal = 1;
		
		new( divNum, grdNum, histArr, dosArr, Ref(dosIncrVal), histMinRatioThres, Ref(wlCounterVal), wlResetInterval );
	end
end

function getFlipCheckerName( flipCheckerType::Type{WL2dLinkFlatFlipChecker} )
	return "WangLandau2dFlatHist";
end

function getFlipCheckerAttrLst( flipChecker::WL2dLinkFlatFlipChecker )
	return ["histMinRatio"];
end

function getFlipCheckerAttrValLst( flipChecker::WL2dLinkFlatFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.([flipChecker.histMinRatioThres]; digits = rndDigs);
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WL2dLinkFlatFlipChecker )
	isHistRest = false;
	flipChecker.wlCounter[] += 1;
	if mod( flipChecker.wlCounter[], flipChecker.wlResetInterval ) == 0
		# maxHist = maximum(flipChecker.histArr);
		meanHist = mean(flipChecker.histArr);
		minHist = minimum(flipChecker.histArr);
		
		if minHist > meanHist * flipChecker.histMinRatioThres
			isHistRest = true;
			# @infiltrate
		end
		# @infiltrate
	end
	
	return isHistRest;
end

struct WL2dLinkPartFlatFlipChecker <: AbstractWangLandau2dLinkOnlyFlipChecker
	divNum::Int64;
	grdNum::Int64;
	
	histArr::Vector{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	
	histCutoffThres::Float64;
	histStdRatioThres::Float64;	
	
	histIsOverArr::Vector{Bool};
	histIsOverArrPrev::Vector{Bool};
	histMaskedArr::Vector{Int64};
	
	wlCounter::Ref{Int64};
	wlResetInterval::Int64;
	
	function WL2dLinkPartFlatFlipChecker( divNum::Int64, nDim::Int64; dosIncrVal::Float64 = 1.0, histCutoffThres = 0.5, histStdRatioThres::Float64 = 0.15, wlResetInterval::Int64 = 100 )
		grdNum = divNum^nDim;
		histArr = zeros( Int64, grdNum+1-2 );
		dosArr = similar(histArr);
		dosArr .= 0;
		histIsOverArr = similar( histArr, Bool );
		histIsOverArrPrev = similar( histIsOverArr );
		histIsOverArr .= false;
		histIsOverArrPrev .= false;
		histMaskedArr = similar( histArr );
		wlCounterVal = 1;
		
		new( divNum, grdNum, histArr, dosArr, Ref(dosIncrVal), histCutoffThres, histStdRatioThres, histIsOverArr, histIsOverArrPrev, histMaskedArr, Ref(wlCounterVal), wlResetInterval );
	end
end

function getFlipCheckerName( flipCheckerType::Type{WL2dLinkPartFlatFlipChecker} ) 
	return "WL2dExploredFlatFlip"
end

function getFlipCheckerAttrLst( flipChecker::WL2dLinkPartFlatFlipChecker )
	return ["histStdRatioThres", "histCutoffThres", "wlResetIntvl"];
end

function getFlipCheckerAttrValLst( flipChecker::WL2dLinkPartFlatFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.( [flipChecker.histStdRatioThres, flipChecker.histCutoffThres, flipChecker.wlResetInterval]; digits = rndDigs );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WL2dLinkPartFlatFlipChecker )
	isHistReset = false;
	flipChecker.wlCounter[] += 1;
	if mod( flipChecker.wlCounter[], flipChecker.wlResetInterval ) == 0
		maxCnt, cntId = findmax( flipChecker.histArr );
		# maxDos = flipChecker.dosArr[ cntId ];
		cntCutoff = maxCnt * flipChecker.histCutoffThres;
		# dosCutoff = maxDos * flipChecker.histCutoffThres;
		# flipChecker.histIsOverArr .= flipChecker.histArr .>= cntCutoff .|| flipChecker.dosArr .>= dosCutoff;
		flipChecker.histIsOverArr .= flipChecker.histArr .>= cntCutoff .|| flipChecker.histIsOverArrPrev;
		flipChecker.histMaskedArr .= flipChecker.histArr .* flipChecker.histIsOverArr;
		
		cntOver = sum(flipChecker.histIsOverArr);
		meanHistOver = sum( flipChecker.histMaskedArr ) / cntOver;
		stdHistOver = sqrt.( sum( ( flipChecker.histMaskedArr .- meanHistOver ).^2 .* flipChecker.histIsOverArr ) / (cntOver - 1) );
		
		if stdHistOver < meanHistOver .* flipChecker.histStdRatioThres
			isHistReset = true;
			flipChecker.histIsOverArrPrev .= flipChecker.histIsOverArr;
			@infiltrate
		end
		# @infiltrate
	end
	
	return isHistReset
end

abstract type AbstractWangLandau3dFlipChecker <: AbstractWangLandauFlipChecker end

function wangLandauUpdateHistDos( flipChecker::AbstractWangLandau3dFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	dS = boolToFlipChange( BfieldLst[dim][pos] );
	# dL = 0;
	# for iLnkDim in 1 : params.nDimLayer
		# dimLink = params.linkDimLst[dim][iLnkDim];
		# dimLinkSh = params.linkDimShLst[dim][iLnkDim];
		# dL += boolToFlipChange( linkLst[dimLink][pos] );
		# dL += boolToFlipChange( linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] );
	# end
	dL = calcDL( params, linkLst, dim, pos );
	
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

struct WangLandauNoResetFlipChecker <: AbstractWangLandau3dFlipChecker
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

struct WangLandauFlipChecker <: AbstractWangLandau3dFlipChecker
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

struct WangLandauStdFlipChecker <: AbstractWangLandau3dFlipChecker
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
			flipChecker.dosIncrRef[] < 1
			if histStd < histMean * flipChecker.histStdRatioThres
				isReset = true;
				# @infiltrate
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
	
	println( "f = ", flipChecker.dosIncrRef[] );
	
	save( fNameWL, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	
	open( dirLog * fNameFileLstWL, "w" ) do io
		println( io, fNameWL );
	end
	
	return fNameWL;
end

function loops_MC_methods_WL2d( divNum = 64, itNum = 10000; itNumSample = 100, itStartSample = 50, isInit0 = false, cAreaInit = 0, dosIncrInit = 1, nDim = 2 )
	# flipChecker = WangLandauFlipChecker( histDivNum; dosIncrVal = dosIncrInit );
	# flipChecker = WangLandauNoResetFlipChecker( histDivNum );
	flipChecker = WL2dLinkFlatFlipChecker( divNum, nDim; histMinRatioThres = 0.8 );
	# flipChecker = WL2dLinkPartFlatFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = 1000 );
	initializer = genMeanFieldInitializer( cAreaInit );
	updaterType = SingleUpdater;
	
	fName = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, initializer = initializer, itNumSample = itNumSample, itStartSample = itStartSample, nDim = nDim );
	
	fMain = "loops_WL2d";
	attrLst, valLst = genAttrLstLttcFlipInit( divNum, itNum, nDim, flipChecker, initializer );
	fNameWL = fNameFunc( fMain, attrLst, valLst, jld2Type );
	
	println( "f = ", flipChecker.dosIncrRef[] );
	
	save( fNameWL, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	
	open( dirLog * fNameFileLstWL, "w" ) do io
		println( io, fNameWL );
	end
	
	# @infiltrate
	
	return fNameWL;
end
