# Loops_MC module

fNameFileLstWL = "fNameFileLstWL.txt";

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




abstract type AbstractWLHistDos end

throwWLHistDosUndefined() = error( "Lops_MC: WLHistDos undefined" );

function getHistArr( histDos::AbstractWLHistDos )
	if getIsHistOutBackup( histDos )
		return histDos.histArrFlatBackup;
	else
		return histDos.histArr;
	end
end

function getDosArr( histDos::AbstractWLHistDos )
	return histDos.dosArr;
end

getIsHistOutBackup( histDos::AbstractWLHistDos ) = histDos.isHistOutBackupRef[];

function setHistOutBackup!( histDos::AbstractWLHistDos )
	histDos.isHistOutBackupRef[] = true;
end

function unsetHistOutBackup!( histDos::AbstractWLHistDos )
	histDos.isHistOutBackupRef[] = false;
end

getHistUpdateId( wlHistDos::AbstractWLHistDos, params::ParamsLoops, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D} = throwWLHistDosUndefined();

getLinkHistId( wlHistDos::AbstractWLHistDos, lnkVal::Int64 ) = throwWLHistDosUndefined();

function histUpdate( histDos::AbstractWLHistDos, params::ParamsLoops, flipProposer::FlipProposer, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D}, dosIncr::Float64 ) where {D}
	id, idNxt = getHistUpdateId( histDos, params, linkLst, dim, pos );
	isFlip = histUpdateWithId( histDos, dosIncr, id, idNxt );
	unsetHistOutBackup!( histDos );
	
	return isFlip;
end

function histUpdateWithId( histDos::AbstractWLHistDos, dosIncr::Float64, lId::CartesianIndex{histD}, lIdNxt::CartesianIndex{histD} ) where {histD}
	p = exp( histDos.dosArr[lId] - histDos.dosArr[lIdNxt] );
	isFlip = rand() < p;
	
	if isFlip
		histDos.histArr[lIdNxt] += 1;
		histDos.dosArr[lIdNxt] += dosIncr;
	else
		histDos.histArr[lId] += 1;
		histDos.dosArr[lId] += dosIncr;
	end
	
	return isFlip;
end

function histResetZero!( histDos::AbstractWLHistDos )
	histDos.histArrFlatBackup .= histDos.histArr;
	histDos.histArr .= 0;
	
	setHistOutBackup!( histDos );
end




struct WLHistDos2DFull <: AbstractWLHistDos
	histArr::Vector{Int64};
	histArrFlatBackup::Vector{Int64};
	dosArr::Vector{Float64};
	
	isHistOutBackupRef::Ref{Bool};
	
	grdNum::Int64;
	
	posLst::CartesianIndices{1,Tuple{Vararg{Base.OneTo{Int64},1}}};
	
	function WLHistDos2DFull( divNum )
		nDim = 2;
		grdNum = divNum^nDim;
		histArr = zeros( Int64, grdNum+1-2 );
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr, Float64);
		dosArr .= 0;
		
		isHistOutBackupVal = false;
		
		posLst = CartesianIndices(histArr);
		
		new( histArr, histArrFlatBackup, dosArr, isHistOutBackupVal, grdNum, posLst );
	end
end

function getHistUpdateId( wlHistDos::WLHistDos2DFull, params::ParamsLoops, linkLst::Vector{Array{Bool,2}}, dim::Int64, pos::CartesianIndex{2} )
	dL = calcDL( params, linkLst, dim, pos );
	
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
	lId = getLinkHistId( wlHistDos, lTotal );
	lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	return wlHistDos.posLst[lId], wlHistDos.posLst[lIdNxt];
end

function getLinkHistId( wlHistDos::WLHistDos2DFull, lnkVal::Int64 )
	lnkId = Int64( lnkVal / 2 );
	if lnkId == wlHistDos.grdNum
		lnkId -= 2;
	elseif lnkId > 0
		lnkId -= 1;
	end
	
	lnkId += 1;
	
	return lnkId;
end



struct WLHistDos2DHalf <: AbstractWLHistDos
	histArr::Vector{Int64};
	histArrFlatBackup::Vector{Int64};
	dosArr::Vector{Float64};
	
	isHistOutBackupRef::Ref{Bool};
	
	grdNum::Int64;
	grdNumHalf::Int64;
	
	posLst::CartesianIndices{1,Tuple{Vararg{Base.OneTo{Int64},1}}};
	
	function WLHistDos2DHalf( divNum )
		nDim = 2;
		grdNum = divNum^nDim;
		grdNumHalf = Int64(grdNum / 2);
		histArr = zeros( Int64, Int64( grdNum / 2) );
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr, Float64);
		dosArr .= 0;
		
		isHistOutBackupVal = false;
		
		posLst = CartesianIndices(histArr);
		
		new( histArr, histArrFlatBackup, dosArr, isHistOutBackupVal, grdNum, grdNumHalf, posLst );
	end
end

function getHistUpdateId( wlHistDos::WLHistDos2DHalf, params::ParamsLoops, linkLst::Vector{Array{Bool,2}}, dim::Int64, pos::CartesianIndex{2} )
	dL = calcDL( params, linkLst, dim, pos );
	
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
	lId = getLinkHistId( wlHistDos, lTotal );
	lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	return wlHistDos.posLst[lId], wlHistDos.posLst[lIdNxt];
end

function getLinkHistId( wlHistDos::WLHistDos2DHalf, lnkVal::Int64 )
	lnkId = Int64( lnkVal / 2 );
	if lnkId > wlHistDos.grdNumHalf
		lnkId = wlHistDos.grdNum zZzZS1111111111111111111111   - lnkId;
	end
	if lnkId > 0
		lnkId -= 1;
	end
	
	lnkId += 1;
	
	return lnkId;
end




abstract type AbstractFlipCheckerWithProposer <: FlipChecker end

function flipCheckDoIt( flipChecker::AbstractFlipCheckerWithProposer, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	if flipCheck( flipChecker, flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
		flipDoIt( flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	end
end

function flipCheck( flipChecker::AbstractFlipCheckerWithProposer, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	error("Loops_MC: flipChecker or flipProposer not defined yet");
end




abstract type AbstractWangLandauFlipChecker <: AbstractFlipCheckerWithProposer end


struct WangLandauItController <: ItController
	dosIncrMin::Float64;
	flipChecker::AbstractWangLandauFlipChecker;
	
	itRef::Ref{Int64};
	
	function WangLandauItController( dosIncrMin::Float64, flipChecker::AbstractWangLandauFlipChecker )
		itVal = 1;
		
		new( dosIncrMin, flipChecker, Ref(itVal) );
	end
end

function getAttrLstItController( itControllerType::Type{WangLandauItController} )
	return ["dosIncrMin"];
end

function getValLstItController( itController::WangLandauItController )
	return [itController.dosIncrMin];
end

function getItNumLst( itController::WangLandauItController )
	return zeros(Int64, 3);
end

function testItNotDone( wlItController::WangLandauItController )
	return getDosIncr( wlItController.flipChecker ) > wlItController.dosIncrMin;
end

function testItDoSample( wlItController::WangLandauItController )
	return getIsHistFlat( wlItController.flipChecker );
end

testItDoStartSample( wlItController::WangLandauItController ) = false;




function flipCheck( flipChecker::AbstractWangLandauFlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	isFlip = wangLandauUpdateHistDos( flipChecker, flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	
	if wangLandauHistResetCheck( flipChecker )
		setIsHistFlat!( flipChecker );
		wangLandauUpdateDosIncr( flipChecker );
	else
		unsetIsHistFlat!(flipChecker);
	end
	
	return isFlip;
end

function getDosIncr( flipChecker::AbstractWangLandauFlipChecker )
	return flipChecker.dosIncrRef[];
end

function getIsHistFlat( flipChecker::AbstractWangLandauFlipChecker )
	return flipChecker.isHistFlatRef[];
end

function setIsHistFlat!( flipChecker::AbstractWangLandauFlipChecker )
	flipChecker.isHistFlatRef[] = true;
	flipChecker.histArrFlatBackup .= flipChecker.histArr;
end

function unsetIsHistFlat!( flipChecker::AbstractWangLandauFlipChecker )
	flipChecker.isHistFlatRef[] = false;
end

function getHistArr( flipChecker::AbstractWangLandauFlipChecker )
	if getIsHistFlat( flipChecker )
		return flipChecker.histArrFlatBackup;
	else
		return flipChecker.histArr;
	end
end

function getDosArr( flipChecker::AbstractWangLandauFlipChecker )
	return flipChecker.dosArr;
end

function getHistDosArr( flipChecker::AbstractWangLandauFlipChecker )
	return [getHistArr( flipChecker ), getDosArr( flipChecker )];
end

function wangLandauUpdateHistDos( flipChecker::AbstractWangLandauFlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} )where {D}
	error("WangLandau Flipper or FlipProser not defined")
end

function wangLandauHistResetCheck( flipChecker::AbstractWangLandauFlipChecker )
	error("WangLandau Flipper not defined")
end

function wangLandauUpdateDosIncr( flipChecker::AbstractWangLandauFlipChecker )
	flipChecker.histArr .= 0
	flipChecker.dosIncrRef[] /= 2;
end



abstract type AbstractWLHistStructFlipChecker <: AbstractWangLandauFlipChecker end

function flipCheck( flipChecker::AbstractWLHistStructFlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	# wangLandauUpdateHistDos( flipChecker, flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	isFlip = histUpdate( flipChecker.histDos, params, flipProposer, linkLst, dim, pos, flipChecker.dosIncrRef[] );
	
	if wangLandauHistResetCheck( flipChecker )
		wangLandauUpdateDosIncr( flipChecker );
		histResetZero!( flipChecker.histDos );
	end
	
	return isFlip;
end

function getDosIncr( flipChecker::AbstractWLHistStructFlipChecker )
	return flipChecker.dosIncrRef[];
end

function getIsHistFlat( flipChecker::AbstractWLHistStructFlipChecker )
	return getIsHistOutBackup( flipChecker.histDos );
end

function setIsHistFlat!( flipChecker::AbstractWLHistStructFlipChecker )
	nothing;
end

function unsetIsHistFlat!( flipChecker::AbstractWLHistStructFlipChecker )
	nothing;
end

function getHistArr( flipChecker::AbstractWLHistStructFlipChecker )
	return getHistArr( flipChecker.histDos );
end

function getDosArr( flipChecker::AbstractWLHistStructFlipChecker )
	return getDosArr( flipChecker.histDos );
end

function getHistDosArr( flipChecker::AbstractWLHistStructFlipChecker )
	return [getHistArr( flipChecker ), getDosArr( flipChecker )];
end

# function wangLandauUpdateHistDos( flipChecker::AbstractWLHistStructFlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} )where {D}
	# error("WangLandau Flipper or FlipProser not defined")
# end

function wangLandauHistResetCheck( flipChecker::AbstractWLHistStructFlipChecker )
	error("WangLandau Flipper not defined")
end

function wangLandauUpdateDosIncr( flipChecker::AbstractWLHistStructFlipChecker )
	flipChecker.dosIncrRef[] /= 2;
end



struct WL2dHistStructFlipChecker{T_histDos} <: AbstractWLHistStructFlipChecker
	divNum::Int64;
	grdNum::Int64;
	
	histDos::T_histDos;
	dosIncrRef::Ref{Float64};
	
	histCutoffThres::Float64;
	histStdRatioThres::Float64;	
	
	histIsOverArr::Vector{Bool};
	histIsOverArrPrev::Vector{Bool};
	histMaskedArr::Vector{Int64};
	
	wlCounter::Ref{Int64};
	wlResetInterval::Int64;
	
	function WL2dHistStructFlipChecker( divNum::Int64, nDim::Int64; dosIncrVal::Float64 = 1.0, histCutoffThres = 0.5, histStdRatioThres::Float64 = 0.15, wlResetInterval::Int64 = 100, wlHistDosType = WLHistDos2DFull )
		grdNum = divNum^nDim;
		histDos = wlHistDosType( divNum );
		histIsOverArr = similar( getHistArr(histDos), Bool );
		histIsOverArrPrev = similar( histIsOverArr );
		histIsOverArr .= false;
		histIsOverArrPrev .= false;
		histMaskedArr = similar( getHistArr(histDos) );
		wlCounterVal = 1;
		
		new{wlHistDosType}( divNum, grdNum, histDos, Ref(dosIncrVal), histCutoffThres, histStdRatioThres, histIsOverArr, histIsOverArrPrev, histMaskedArr, Ref(wlCounterVal), wlResetInterval );
	end
end

function getFlipCheckerName( flipCheckerType::Type{<:WL2dHistStructFlipChecker} ) 
	return "WL2dExploredFlatFlip";
end

function getFlipCheckerAttrLst( flipChecker::WL2dHistStructFlipChecker )
	return ["histStdRatioThres", "histCutoffThres", "wlResetIntvl"];
end

function getFlipCheckerAttrValLst( flipChecker::WL2dHistStructFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.( [flipChecker.histStdRatioThres, flipChecker.histCutoffThres, flipChecker.wlResetInterval]; digits = rndDigs );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WL2dHistStructFlipChecker )
	isHistReset = false;
	flipChecker.wlCounter[] += 1;
	if mod( flipChecker.wlCounter[], flipChecker.wlResetInterval ) == 0
		maxCnt, cntId = findmax( getHistArr( flipChecker.histDos ) );
		cntCutoff = maxCnt * flipChecker.histCutoffThres;
		flipChecker.histIsOverArr .= getHistArr( flipChecker.histDos ) .>= cntCutoff .|| flipChecker.histIsOverArrPrev;
		flipChecker.histMaskedArr .= getHistArr( flipChecker.histDos ) .* flipChecker.histIsOverArr;
		
		cntOver = sum(flipChecker.histIsOverArr);
		meanHistOver = sum( flipChecker.histMaskedArr ) / cntOver;
		stdHistOver = sqrt.( sum( ( flipChecker.histMaskedArr .- meanHistOver ).^2 .* flipChecker.histIsOverArr ) / (cntOver - 1) );
		
		if stdHistOver < meanHistOver * flipChecker.histStdRatioThres
			isHistReset = true;
			flipChecker.histIsOverArrPrev .= flipChecker.histIsOverArr;
		end
	end
	
	return isHistReset
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

function wangLandauUpdateHistDos( flipChecker::AbstractWangLandau2dLinkOnlyFlipChecker, flipProposer::OneFlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} )where {D}
	dL = calcDL( params, linkLst, dim, pos );
	
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
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
	histArrFlatBackup::Vector{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	isHistFlatRef::Ref{Bool};
	
	histMinRatioThres::Float64;
	
	wlCounter::Ref{Int64};
	wlResetInterval::Int64;
	
	function WL2dLinkFlatFlipChecker( divNum::Int64, nDim::Int64; dosIncrVal::Float64 = 1.0, histMinRatioThres = 0.8, wlResetInterval::Int64 = 100 )
		grdNum = divNum^nDim;
		histArr = zeros( Int64, grdNum+1-2 );
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr);
		dosArr .= 0;
		wlCounterVal = 1;
		isHistFlatVal = false;
		
		new( divNum, grdNum, histArr, histArrFlatBackup, dosArr, Ref(dosIncrVal), Ref(isHistFlatVal), histMinRatioThres, Ref(wlCounterVal), wlResetInterval );
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
		meanHist = mean(flipChecker.histArr);
		minHist = minimum(flipChecker.histArr);
		
		if minHist > meanHist * flipChecker.histMinRatioThres
			isHistRest = true;
		end
	end
	
	return isHistRest;
end

struct WL2dLinkPartFlatFlipChecker <: AbstractWangLandau2dLinkOnlyFlipChecker
	divNum::Int64;
	grdNum::Int64;
	
	histArr::Vector{Int64};
	histArrFlatBackup::Vector{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	isHistFlatRef::Ref{Bool};
	
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
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr);
		dosArr .= 0;
		histIsOverArr = similar( histArr, Bool );
		histIsOverArrPrev = similar( histIsOverArr );
		histIsOverArr .= false;
		histIsOverArrPrev .= false;
		histMaskedArr = similar( histArr );
		wlCounterVal = 1;
		isHistFlatVal = false;
		
		new( divNum, grdNum, histArr, histArrFlatBackup, dosArr, Ref(dosIncrVal), Ref(isHistFlatVal), histCutoffThres, histStdRatioThres, histIsOverArr, histIsOverArrPrev, histMaskedArr, Ref(wlCounterVal), wlResetInterval );
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
			# @infiltrate
		end
		# @infiltrate
	end
	
	return isHistReset
end

abstract type AbstractWangLandau3dFlipChecker <: AbstractWangLandauFlipChecker end

function wangLandauUpdateHistDos( flipChecker::AbstractWangLandau3dFlipChecker, flipProposer::OneFlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	dS = boolToFlipChange( BfieldLst[dim][pos] );
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
	histArrFlatBackup::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	isHistFlatRef::Ref{Bool};
	
	function WangLandauNoResetFlipChecker( histDivNum::Int64; dosIncrVal::Float64 = 1.0 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		histArrFlatBackup = similar( histArr );
		dosArr = zeros(Int64, histDivNum, histDivNum);
		isHistFlatVal = false;
		
		new( histDivNum, histArr, histArrFlatBackup, dosArr, Ref(dosIncrVal), Ref(isHistFlatVal) );
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
	histArrFlatBackup::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	histCoverThreshold::Float64;
	histLowRatioThreshold::Float64;
	isHistFlatRef::Ref{Bool};
	
	function WangLandauFlipChecker( histDivNum::Int64; histLowRatioThreshold::Float64 = 0.8, histCoverThreshold::Float64 = 0.8, dosIncrVal::Float64 = 1.0 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		histArrFlatBackup = similar(histArr);
		dosArr = zeros(Int64, histDivNum, histDivNum);
		
		isHistFlatVal = false;
		
		new( histDivNum, histArr, histArrFlatBackup, dosArr, Ref(dosIncrVal), Ref(isHistFlatRef), histCoverThreshold, histLowRatioThreshold );
	end
end

function getFlipCheckerName( flipType::Type{WangLandauFlipChecker} )
	return "WangLandauFlip";
end

function getFlipCheckerAttrLst( flipChecker::WangLandauFlipChecker )
	attrLst = ["histCoverThres", "histLowRatioThres"];
end

function getFlipCheckerAttrValLst( flipChecker::WangLandauFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.( [ flipChecker.histCoverThreshold, flipChecker.histLowRatioThreshold ]; digits = rndDigs );
	
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
	histArrFlatBackup::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	isHistFlatRef::Ref{Bool};
	
	histMaxCntThres::Int64;
	histCntRatioCutoff::Float64;
	histStdRatioThres::Float64;
	histCoverCntThres::Int64;
	
	histIsOverThresArr::Array{Bool};
	histMaskedOverArr::Array{Int64};
	
	function WangLandauStdFlipChecker( histDivNum::Int64; dosIncrVal::Float64 = 1.0, histMaxCntThres::Int64 = 20, histCntRatioCutoff::Float64 = 0.5, histStdRatioThres::Float64 = 0.1, histCoverCntThres::Int64 = 5 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		histArrFlatBackup = similar(histArr);
		dosArr = zeros(Float64, histDivNum, histDivNum);
		
		histIsOverThresArr = similar( histArr, Bool );
		histMaskedOverArr = similar( histArr );
		
		isHistFlatVal = false;
		
		new( histDivNum, histArr, dosArr, Ref(dosIncrVal), Ref(isHistFlatVal), histMaxCntThres, histCntRatioCutoff, histStdRatioThres, histCoverCntThres, histIsOverThresArr, histMaskedOverArr );
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

struct WL3dPartFlatStdFlipChecker <: AbstractWangLandau3dFlipChecker
	histDivNum::Int64
	histArr::Array{Int64};
	histArrFlatBackup::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	isHistFlatRef::Ref{Bool};
	
	wlCounterRef::Ref{Int64};
	wlResetInterval::Int64;
	
	histStdRatioThres::Float64;
	histCutoffThres::Float64;
	
	histIsOverArr::Array{Bool};
	histIsOverArrPrev::Array{Bool};
	histMaskedArr::Array{Int64};
	
	function WL3dPartFlatStdFlipChecker( histDivNum::Int64; dosIncrVal::Float64 = 1.0, wlResetInterval = 1000, histStdRatioThres = 0.15, histCutoffThres = 0.5 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		histArrFlatBackup = similar( histArr );
		dosArr = similar(histArr, Float64);
		dosArr .= 0;
		
		wlCounterVal = 1;
		
		histIsOverArr = similar(histArr, Bool);
		histIsOverArrPrev = similar(histArr, Bool);
		histIsOverArr .= false;
		histIsOverArrPrev .= false;
		histMaskedArr = similar(histArr);
		
		new( histDivNum, histArr, dosArr, Ref(dosIncrVal), Ref(isHistFlatVal), Ref(wlCounterVal), wlResetInterval, histStdRatioThres, histCutoffThres, histIsOverArr, histIsOverArrPrev, histMaskedArr );
	end
end

function getFlipCheckerName( flipCheckerType::Type{WL3dPartFlatStdFlipChecker} )
	return "WL3dPartFlatStdFlip";
end

function getFlipCheckerAttrLst( flipChecker::WL3dPartFlatStdFlipChecker )
	return ["histStdRatioThres", "histCutoff", "wlResetIntrvl"];
end

function getFlipCheckerAttrValLst( flipChecker::WL3dPartFlatStdFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.( [flipChecker.histStdRatioThres, flipChecker.histCutoffThres, flipChecker.wlResetInterval]; digits = rndDigsLpsMC );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WL3dPartFlatStdFlipChecker )
	isHistReset = false;
	flipChecker.wlCounterRef[] += 1;
	if mod( flipChecker.wlCounterRef[], flipChecker.wlResetInterval ) == 0
		maxCnt, idMaxCnt = findmax( flipChecker.histArr );
		histCutoff = maxCnt * flipChecker.histCutoffThres;
		flipChecker.histIsOverArr .= flipChecker.histArr .> histCutoff .|| flipChecker.histIsOverArrPrev;
		flipChecker.histMaskedArr .= flipChecker.histArr .* flipChecker.histIsOverArr;
		cntOver = sum(flipChecker.histIsOverArr);
		
		meanHistOver = sum( flipChecker.histMaskedArr ) / cntOver;
		stdHistOver = sqrt( sum( ( flipChecker.histMaskedArr .- meanHistOver ).^2 .* flipChecker.histIsOverArr ) ) / (cntOver - 1);
		
		if stdHistOver < meanHistOver * flipChecker.histStdRatioThres
			isHistReset = true;
			flipChecker.histIsOverArrPrev .= flipChecker.histIsOverArr;
		end
	end
	
	return isHistReset;
end





struct WangLandauAuxData <: AuxData
	itSampleLst::Vector{Int64};
	
	flipChecker::FlipChecker;
	
	dataLst::Vector{Array};
	dataSampleLst::Vector{<:Vector{<:Array}};
	dataStartSampleLst::Vector{<:Vector{<:Array}};
	dataNumLst::Vector{Vector};
	
	dataSampleOutLst::Vector{Array};
	dataStartSampleOutLst::Vector{Array};
	dataNumOutLst::Vector{Array};
	
	jldVarSampleLst::Vector{Any};
	jldVarStartSampleLst::Vector{Any};
	jldVarNumLst::Vector{Any};
	jldVarItSampleLst::Vector{Any};
	
	function WangLandauAuxData( flipChecker::AbstractWangLandauFlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 )
		itSampleLst = zeros(Int64, itNumSample);
		
		dataLst = Array[getHistArr(flipChecker), getDosArr(flipChecker)];
		dataSampleLst = [ Vector[ similar( dataLst[ii] ) for itSample = 1 : itNumSample ] for ii = 1 : length(dataLst) ];
		# dataSampleLst = Vector{Array}[histArrSampleLst, dosArrSampleLst];
		# dataStartSampleLst = Vector{Array}[histArrStartSampleLst, dosArrStartSampleLst];
		dataStartSampleLst = [ Vector[ similar( dataLst[ii] ) for itSample = 1 : itNumStartSample ] for ii = 1 : length(dataLst) ];
		dataNumLst = Vector{Vector}[];
		
		dataSampleOutLst = Vector{Array}(undef, length(dataSampleLst));
		dataStartSampleOutLst = Vector{Array}(undef, length(dataStartSampleLst));
		dataNumOutLst = copy(dataNumLst);
		
		jldVarSampleLst = Vector{Any}(undef,0);
		jldVarStartSampleLst = similar(jldVarSampleLst);
		jldVarNumLst = similar(jldVarSampleLst);
		jldVarItSampleLst = similar(jldVarSampleLst);
		
		new( itSampleLst, flipChecker, dataLst, dataSampleLst, dataStartSampleLst, dataNumLst, dataSampleOutLst, dataStartSampleOutLst, dataNumOutLst, jldVarSampleLst, jldVarStartSampleLst, jldVarNumLst, jldVarItSampleLst );
	end
end

WangLandauAuxData( flipChecker::AbstractWangLandauFlipChecker ) = WangLandauAuxData( flipChecker, 0, 0, 0 );
WangLandauAuxData( params::ParamsLoops, flipChecker::AbstractWangLandauFlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) = WangLandauAuxData( flipChecker, itNum, itNumSample, itNumStartSample );
WangLandauAuxData( params::ParamsLoops, flipChecker::AbstractWangLandauFlipChecker ) = WangLandauAuxData( flipChecker );
genAuxData( auxDataType::Type{WangLandauAuxData}, params::ParamsLoops, flipChecker::AbstractWangLandauFlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) = auxDataType( params, flipChecker, itNum, itNumSample, itNumStartSample );

function getAuxDataSummaryName( wlAuxDataType::Type{WangLandauAuxData} )
	return "WLHistDos";
end

function getAuxDataNameLst( wlAuxDataType::Type{WangLandauAuxData} )
	return ["histArr", "dosArr"];
end

getAuxDataNumNameLst( wlAuxDataType::Type{WangLandauAuxData} ) = String[];

function storeAuxDataSampleNoBndCheck( wlAuxData::WangLandauAuxData, itSample::Int64 )
	for ii = 1 : length(wlAuxData.dataLst)
		wlAuxData.dataSampleLst[ii][itSample] .= wlAuxData.dataLst[ii];
	end
end

function storeAuxDataStartSampleNoBndCheck( wlAuxData::WangLandauAuxData, itStartSample::Int64 )
	for ii = 1 : length(wlAuxData.dataLst)
		wlAuxData.dataSampleLst[ii][itStartSample] .= wlAuxData.dataLst[ii];
	end
end

storeAuxDataNumNoBndCheck( wlAuxData::WangLandauAuxData, it::Int64 ) = nothing;

function flipAuxData!( wlAuxData::WangLandauAuxData, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	nothing;
end

function calcAuxData!( wlAuxData::WangLandauAuxData, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	# nothing;
	wlAuxData.dataLst[1] = getHistArr( wlAuxData.flipChecker );
	wlAuxData.dataLst[2] = getDosArr( wlAuxData.flipChecker );
end

function renewAuxDataOutLst!( wlAuxData::WangLandauAuxData )
	funStackArrs = ( arr -> cat( arr...; dims = ndims(wlAuxData.dataLst[1])+1 ) );
	wlAuxData.dataSampleOutLst .= funStackArrs.(wlAuxData.dataSampleLst);
	wlAuxData.dataStartSampleOutLst .= funStackArrs.(wlAuxData.dataStartSampleLst);
	
	GC.gc();
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
	# flipChecker = WangLandauStdFlipChecker( histDivNum; histMaxCntThres = 10, histStdRatioThres = 0.15 );
	flipChecker = WL3dPartFlatStdFlipChecker( histDivNum; wlResetInterval = 1000, histStdRatioThres = 0.15 );
	flipProposer = OneFlipProposer();
	initializer = genMeanFieldInitializer( cAreaInit );
	updaterType = SingleUpdater;
	
	fName = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, itNumSample = itNumSample, itStartSample = itStartSample );
	
	fMain = "loops_WL_single";
	attrLst, valLst = genAttrLstLttcFlipInit( divNum, itNum, nDim, flipChecker, initializer );
	fNameWL = fNameFunc( fMain, attrLst, valLst, jld2Type );
	
	println( "f = ", flipChecker.dosIncrRef[] );
	
	save( fNameWL, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	
	open( dirLog * fNameFileLstWL, "w" ) do io
		println( io, fNameWL );
	end
	
	return fNameWL;
end

function loops_MC_methods_WL2d( divNum, itNum; itNumSample = 100, itStartSample = 50, isInit0 = false, cAreaInit = 0, dosIncrInit = 1, nDim = 2, isFileNameOnly = false, fMainOutside = "" )
	# flipChecker = WangLandauFlipChecker( histDivNum; dosIncrVal = dosIncrInit );
	# flipChecker = WangLandauNoResetFlipChecker( histDivNum );
	# flipChecker = WL2dLinkFlatFlipChecker( divNum, nDim; histMinRatioThres = 0.8 );
	flipChecker = WL2dLinkPartFlatFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = 1000 );
	flipProposer = OneFlipProposer();
	initializer = genMeanFieldInitializer( cAreaInit );
	updaterType = SingleUpdater;
	auxDataType = WangLandauAuxData;
	
	itController = ItNumItController( itNum, itNumSample, itStartSample );
	
	# fName = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, auxDataType = auxDataType, itNumSample = itNumSample, itStartSample = itStartSample, nDim = nDim, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	
	fName = loops_MC_methods_Base( divNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, auxDataType = auxDataType, itController = itController, nDim = nDim, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	
	fMain = "loops_WL2d";
	attrLst, valLst = genAttrLstLttcFlipInit( divNum, itNum, nDim, flipChecker, initializer );
	fNameWL = fNameFunc( fMain, attrLst, valLst, jld2Type );
	
	println( "f = ", flipChecker.dosIncrRef[] );
	
	save( fNameWL, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	
	open( dirLog * fNameFileLstWL, "w" ) do io
		println( io, fNameWL );
	end
	
	# @infiltrate
	
	return fName;
end

function loops_MC_methods_WL2d( divNum = 64; dosIncrInit = 1, dosIncrMin = 0.001, cAreaInit = 0, nDim = 2, isFileNameOnly = false, fMainOutside = "" )
	# flipChecker = WangLandauFlipChecker( histDivNum; dosIncrVal = dosIncrInit );
	# flipChecker = WangLandauNoResetFlipChecker( histDivNum );
	# flipChecker = WL2dLinkFlatFlipChecker( divNum, nDim; histMinRatioThres = 0.8 );
	# flipChecker = WL2dLinkPartFlatFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = 1000 );
	flipChecker = WL2dHistStructFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = 1000, wlHistDosType = WLHistDos2DHalf );
	flipProposer = OneFlipProposer();
	initializer = genMeanFieldInitializer( cAreaInit );
	updaterType = SingleUpdater;
	auxDataType = WangLandauAuxData;
	
	itController = WangLandauItController( dosIncrMin, flipChecker );
	
	# fName = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, auxDataType = auxDataType, itNumSample = itNumSample, itStartSample = itStartSample, nDim = nDim, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	
	fName = loops_MC_methods_Base( divNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, auxDataType = auxDataType, itController = itController, nDim = nDim, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	
	itNum = itController.itRef[];
	
	fMain = "loops_WL2d";
	attrLst, valLst = genAttrLstLttcFlipInit( divNum, itNum, nDim, flipChecker, initializer );
	fNameWL = fNameFunc( fMain, attrLst, valLst, jld2Type );
	
	println( "f = ", flipChecker.dosIncrRef[] );
	
	save( fNameWL, "histArr", getHistArr( flipChecker ), "dosArr", getDosArr( flipChecker ) );
	
	open( dirLog * fNameFileLstWL, "w" ) do io
		println( io, fNameWL );
	end
	
	# @infiltrate
	
	return fName;
end
