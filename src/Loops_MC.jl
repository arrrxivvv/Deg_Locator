module Loops_MC

using ShiftedArrays
using FilenameManip
using Random
using JLD2
using Statistics
using Distributions
using Utils
 
# using Infiltrator

export loops_MC_methods, loops_MC, loops_MC_smart, loops_MC_staggeredCube

oFNameLoopsMain = "loopsSample";
oFNameLoopsStartMain = "loopsStartSample";
oFNameLoopsNumMain = "loopsNum";

dirLog = "./log/";

fMainLoopsMC = "loops_MC";
attrLstLoops = ["divNum","itNum","cArea","cPerim"];
attrLstLoopsBeta = push!( deepcopy(attrLstLoops), "beta" );
attrFerro = "cFerro";
attrLstLoopsFerro = push!( deepcopy(attrLstLoops), attrFerro );
attrLstLoopsBetaFerro = push!( deepcopy( attrLstLoopsBeta ), attrFerro );

function getAttrValLstLoopsMC( divNum, itNum, cArea, cPerim; cFerro = 0, beta = nothing )
	valLst = Any[divNum, itNum, cArea, cPerim];
	if isnothing(beta)
		attrLstNoFerro = attrLstLoops;
		attrLstFerro = attrLstLoopsFerro;
	else
		attrLstNoFerro = attrLstLoopsBeta;
		attrLstFerro = attrLstLoopsBetaFerro;
		push!(valLst, beta);
	end
	attrLst = cFerro == 0 ? attrLstNoFerro : attrLstFerro;
	if cFerro != 0
		push!( valLst, cFerro );
	end
	
	return attrLst, valLst;
end

struct ParamsLoops{N}
	nDim::Int64;
	nDimLayer::Int64;
	divNum::Int64;
	grdNum::Int64;
	divLst::Vector{Int64};
	posLst::CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}};
	posLstShLst:: Matrix{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}}}; 
	linkDimLst::Vector{Vector{Int64}};
	linkDimShLst::Vector{Vector{Int64}};
	
	function ParamsLoops( divNum::Int64, nDim::Int64 )
		nDimLayer = nDim - 1;
		divLst = fill(divNum, nDim);
		grdNum = divNum^nDim;
		posLst = CartesianIndices(ntuple(x->divNum,nDim));
		posLstShLst = Utils.arrShAdvRetLstFunc( posLst, nDim );
		linkDimLst = [ zeros(Int64, nDim-1) for dim = 1 : nDim ];
		for dim = 1:nDim
			iDLink = 1;
			for dimLink = 1 : nDim
				if dimLink != dim
					linkDimLst[dim][iDLink] = dimLink;
					iDLink += 1;
				end
			end
		end
		linkDimShLst = [ circshift( linkDimLst[dim], 1 ) for dim = 1 : nDim ];
		
		new{nDim}( nDim, nDimLayer, divNum, grdNum, divLst, posLst, posLstShLst, linkDimLst, linkDimShLst );
	end
end

abstract type FlipChecker end

function flipCheck( flipChecker::Type{<:FlipChecker}, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	error( "Loops_MC: FlipChecker not defined yet" );
end

struct NeighborFlipChecker <: FlipChecker
	pFlipLst::Array{Float64};
	
	function NeighborFlipChecker( cArea, cPerim, cFerroSigned )
		pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
		
		new(pFlipLst);
	end
end

function flipCheck( flipChecker::NeighborFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	iArea = BfieldLst[dim][pos] + 1;
	iL = 1;
	for iLnkDim in 1 : params.nDimLayer
		dimLink = params.linkDimLst[dim][iLnkDim];
		dimLinkSh = params.linkDimShLst[dim][iLnkDim];
		iL += linkLst[dimLink][pos];
		iL += linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
	end
	iLFerro = 1;
	for iLnkDim = 1 : params.nDimLayer
		dimLink = params.linkDimLst[dim][iLnkDim];
		dimLinkSh = params.linkDimShLst[dim][iLnkDim];
		iLFerro += linkFerroLst[iLnkDim,dim][pos];
		iLFerro += linkFerroLst[iLnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]];
	end
	
	pSwitchRand = rand();
	return pSwitchRand < flipChecker.pFlipLst[iLFerro, iL, iArea];
end

struct WangLandauFlipChecker <: FlipChecker
	histDivNum::Int64;
	histArr::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	histCoverThreshold::Float64;
	histLowRatioThreshold::Float64;
	
	function WangLandauFlipChecker( histDivNum::Int64; histLowRatioThreshold::Float64 = 0.8, histCoverThreshold::Float64 = 0.5, dosIncrVal::Float64 = 1.0 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		dosArr = zeros(Int64, histDivNum, histDivNum);
		
		new( histDivNum, histArr, dosArr, Ref(dosIncrVal), histCoverThreshold, histLowRatioThreshold );
	end
end

function flipCheck( flipChecker::WangLandauFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
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
	# @infiltrate
	
	return isFlip;
end

function wangLandauCheckDosIncr( flipChecker::WangLandauFlipChecker )
	coverNum = sum( x-> x>0, flipChecker.histArr );
	coverRatio = coverNum / flipChecker.histDivNum^2;
	if coverRatio > flipChecker.histCoverThreshold
		histAvg = sum( flipChecker.histArr ) / coverNum;
		histMin = minimum( x -> x == 0 ? Inf : x, flipChecker.histArr );
		if histMin > histAvg
			flipChecker.histArr .= 0
			flipChecker.dosIncrRef[] /= 2;
		end
	end
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

abstract type LoopsUpdater end

function updateLoops( updater::LoopsUpdater, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	error( "Loops_MC: updater not defined yet" );
end

function updateLoops( updater::LoopsUpdater, flipChecker::FlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	error( "Loops_MC: updater or flipChecker not defined yet" );
end

function getUpdaterFMod( updaterType::(Type{T} where T <: LoopsUpdater) )
	error( "Loops_MC: updater not defined yet" );
end

struct SingleUpdater <: LoopsUpdater
	;
end

SingleUpdater( params::ParamsLoops ) = SingleUpdater();

function updateLoops( updater::SingleUpdater, flipChecker::FlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	pos = rand(params.posLst);
	dim = rand(1:params.nDim);
	
	if flipCheck( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
		flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = dim );
	end
end

struct ABUpdater{N,N_1} <: LoopsUpdater
	posLstDimLst;
	posLayerLst::CartesianIndices{N_1,Tuple{Vararg{Base.OneTo{Int64},N_1}}};
	posABLst::Matrix{Vector{CartesianIndex{N}}};
	
	pFlipLst::Array{Float64};
	
	function ABUpdater( params::ParamsLoops; cArea::Float64, cPerim::Float64, cFerroSigned::Float64 )
		posLstDimLst = [ selectdim( params.posLst, dim, it ) for dim = 1 : params.nDim, it = 1 : params.divNum ];
		posLayerLst = CartesianIndices( ntuple(x->params.divNum,params.nDim-1) );
		
		posABLst = [ Vector{CartesianIndex{params.nDim}}(undef,Int64(params.divNum^(params.nDim)/2)) for iAB = 1 : 2, iDim = 1 : params.nDim ];
		
		pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
		
		for iDim = 1 : params.nDim
			iA = 1; iB = 1;
			for pos in params.posLst
				sumId = 0;
				for dim = 1 : params.nDim-1
					sumId += pos[params.linkDimLst[iDim][dim]];
				end
				if sumId % 2 == 0
					posABLst[1,iDim][iA] = pos;
					iA += 1;
				else
					posABLst[2,iDim][iB] = pos;
					iB += 1;
				end
			end
		end
		
		new{params.nDim,params.nDim-1}( posLstDimLst, posLayerLst, posABLst, pFlipLst );
	end
end

function updateLoops( updater::ABUpdater, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for dim = 1 : params.nDim
		for iAB = 1 : 2
			Threads.@threads for pos in updater.posABLst[iAB,dim]
				# @time begin
				iArea = BfieldLst[dim][pos] + 1;
				iL = 1;
				for iLnkDim in 1 : params.nDimLayer
					dimLink = params.linkDimLst[dim][iLnkDim];
					dimLinkSh = params.linkDimShLst[dim][iLnkDim];
					iL += linkLst[dimLink][pos];
					iL += linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
				end
				iLFerro = 1;
				for iLnkDim = 1 : params.nDimLayer
					dimLink = params.linkDimLst[dim][iLnkDim];
					dimLinkSh = params.linkDimShLst[dim][iLnkDim];
					iLFerro += linkFerroLst[iLnkDim,dim][pos];
					iLFerro += linkFerroLst[iLnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]];
				end
				# end
				# @time begin
				pSwitchRand = rand();
				if pSwitchRand < updater.pFlipLst[iLFerro, iL, iArea];
					flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = dim );
				end
				# end
				# @infiltrate
			end
		end
	end
end

function getUpdaterFMod( updaterType::Type{ABUpdater} )
	return "upAB";
end

struct StaggeredCubeUpdater{N,Nplus1} <: LoopsUpdater
	posLstSh0::CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}};
	posLstAdvOrNot::Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}} }};
	posShOrNotLst::Vector{Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}}}};
	posStagCubeLst::Vector{Array{CartesianIndex{N},Nplus1}};
	idStagLst::CartesianIndices{Nplus1,NTuple{Nplus1,Base.OneTo{Int64}}};
	
	pFlipLst::Array{Float64};
	
	iDimLst::UnitRange{Int64};
	iIsShLst::UnitRange{Int64};
	randIDimLst::Array{Int64,Nplus1};
	randIShLst::Array{Int64,Nplus1};
	
	function StaggeredCubeUpdater( params::ParamsLoops; cArea::Float64, cPerim::Float64, cFerroSigned::Float64 )
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
		
		pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
		
		iDimLst = 1:params.nDim;
		iIsShLst = 1:2;
		randIDimLst = similar( posStagCubeLst[1], Int64 );
		randIShLst = similar(randIDimLst);
		
		new{params.nDim,params.nDim+1}( posLstSh0, posLstAdvOrNot, posShOrNotLst, posStagCubeLst, idStagLst, pFlipLst, iDimLst, iIsShLst, randIDimLst, randIShLst );
	end
end

function updateLoops( updater::StaggeredCubeUpdater, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for iAdv = 1 : params.nDim+1
		rand!( updater.randIDimLst, updater.iDimLst );
		rand!( updater.randIShLst, updater.iIsShLst );
		Threads.@threads for idStag in updater.idStagLst
			posCube = updater.posStagCubeLst[iAdv][idStag];
			pos = updater.posShOrNotLst[updater.randIDimLst[idStag]][updater.randIShLst[idStag]][posCube];
			iArea = BfieldLst[updater.randIDimLst[idStag]][pos] + 1;
			iL = 1;
			for iLnkDim in 1 : params.nDimLayer
				dimLink = params.linkDimLst[updater.randIDimLst[idStag]][iLnkDim];
				dimLinkSh = params.linkDimShLst[updater.randIDimLst[idStag]][iLnkDim];
				iL += linkLst[dimLink][pos];
				iL += linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
			end
			iLFerro = 1;
			for iLnkDim = 1 : params.nDimLayer
				dimLink = params.linkDimLst[updater.randIDimLst[idStag]][iLnkDim];
				dimLinkSh = params.linkDimShLst[updater.randIDimLst[idStag]][iLnkDim];
				iLFerro += linkFerroLst[iLnkDim,updater.randIDimLst[idStag]][pos];
				iLFerro += linkFerroLst[iLnkDim,updater.randIDimLst[idStag]][params.posLstShLst[dimLinkSh,1][pos]];
			end
			pSwitch = rand();
			pFlip = updater.pFlipLst[iLFerro, iL, iArea];
			if pSwitch < pFlip
				flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = updater.randIDimLst[idStag] );
			end
			# @infiltrate
		end
		# @infiltrate
	end
end

function getUpdaterFMod( updaterType::Type{StaggeredCubeUpdater} )
	return "upStagCube";
end

struct StaggeredCubeUpdaterBase{N,Nplus1} <: LoopsUpdater
	posLstSh0::CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}};
	posLstAdvOrNot::Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}} }};
	posShOrNotLst::Vector{Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}}}};
	posStagCubeLst::Vector{Array{CartesianIndex{N},Nplus1}};
	idStagLst::CartesianIndices{Nplus1,NTuple{Nplus1,Base.OneTo{Int64}}};
	
	iDimLst::UnitRange{Int64};
	iIsShLst::UnitRange{Int64};
	randIDimLst::Array{Int64,Nplus1};
	randIShLst::Array{Int64,Nplus1};
	
	function StaggeredCubeUpdater( params::ParamsLoops; cArea::Float64, cPerim::Float64, cFerroSigned::Float64 )
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
		
		iDimLst = 1:params.nDim;
		iIsShLst = 1:2;
		randIDimLst = similar( posStagCubeLst[1], Int64 );
		randIShLst = similar(randIDimLst);
		
		new{params.nDim,params.nDim+1}( posLstSh0, posLstAdvOrNot, posShOrNotLst, posStagCubeLst, idStagLst, iDimLst, iIsShLst, randIDimLst, randIShLst );
	end
end

function updateLoops( updater::StaggeredCubeUpdaterBase, flipChecker::NeighborFlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for iAdv = 1 : params.nDim+1
		rand!( updater.randIDimLst, updater.iDimLst );
		rand!( updater.randIShLst, updater.iIsShLst );
		Threads.@threads for idStag in updater.idStagLst
			posCube = updater.posStagCubeLst[iAdv][idStag];
			pos = updater.posShOrNotLst[updater.randIDimLst[idStag]][updater.randIShLst[idStag]][posCube];
			
			if flipCheck( flipChecker, params, updater.randIDimLst[idStag], pos, BfieldLst, linkLst, linkFerroLst )
				flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = updater.randIDimLst[idStag] );
			end
		end
	end
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
			
			# if flipCheck( flipChecker, params, updater.randIDimLst[idStag], pos, BfieldLst, linkLst, linkFerroLst )
				# flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = updater.randIDimLst[idStag] );
			# end
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

function getFModLoopsMC( fMod::String, updaterType::Type{<:LoopsUpdater}, isInit0::Bool=false; isFModMethod = true )
	fModOut = fMod;
	if isFModMethod
		fModOut = Utils.strAppendWith_( fModOut, getUpdaterFMod(updaterType) );
	end
	if isInit0
		fModOut = Utils.strAppendWith_( fModOut, "isInit0" );
	end
	
	return fModOut;
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

function loops_MC_methods_WangLandau( divNum = 64, itNum = 10000; histDivNum = 64, itNumSample = 100, itStartSample = 50, isInit0 = false, cAreaInit = 0, dosIncrInit = 1 )
	flipChecker = WangLandauFlipChecker( histDivNum; dosIncrVal = dosIncrInit );
	updaterType = SingleUpdater;
	
	numBfieldLst, numLinkLst, BfieldSampleLst, linkSampleLst, BfieldStartSampleLst, linkStartSampleLst, zakMeanLst, zakLstSampleLst = loops_MC_methods_inJulia( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, itNumSample = itNumSample, itStartSample = itStartSample, isInit0 = isInit0, cAreaInit = cAreaInit );
	
	fMain = "loops_WL_single";
	attrLst = ["divNum", "itNum", "histDivNum", "cAreaInit", "dosIncrInit"];
	valLst = Any[divNum, itNum, histDivNum, cAreaInit, dosIncrInit];
	fNameWL = fNameFunc( fMain, attrLst, valLst, jld2Type );
	
	save( fNameWL, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	
	return fNameWL;
end

function loops_MC_methods_cALF( divNum = 64, itNum = 10000; updaterType::Type{<:LoopsUpdater}, fMod = "", cArea = 1, cPerim = 1, cFerro = 0, isBeta = false, beta = nothing, itNumSample = 100, itStartSample = 50, isInit0 = false )
	flipChecker = NeighborFlipChecker( cArea, cPerim, cFerro );
	
	loops_MC_methods( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, fMod = fMod, isBeta = isBeta, beta = beta, itNumSampe = itNumSample, itStartSample = itStartSample, isInit0 = isInit0 );
end

function loops_MC_methods_inJulia( divNum = 64, itNum = 10000; updaterType::Type{<:LoopsUpdater}, flipChecker::FlipChecker, itNumSample = 100, itStartSample = 50, isInit0 = false, cAreaInit = 0 )
	nDim = 3;
	
	params = ParamsLoops( divNum, nDim );
	nDimLayer = nDim-1;
	divLst = params.divLst;
	divTup = Tuple(divLst);
	posLst = params.posLst;
	posLstShLst = params.posLstShLst;
	linkDimLst = params.linkDimLst;
	linkDimShLst = params.linkDimShLst;
	
	itStep = max( Int64( floor(itNum / itNumSample) ), 1 );
	lnSample = Int64( floor( itNum / itStep ) );
	itStartSample = min( itStartSample, itNum );
	
	zakLstLst = zeros( Bool, divNum, divNum, nDim, itNum );

	BfieldLst = [ zeros( Bool, divTup ) for dim = 1 : nDim ];
	linkLst = [ zeros( Bool, divTup ) for dim = 1 : nDim ];
	linkFerroLst = [ zeros( Bool, divTup ) for lnkDim = 1 : nDim-1, BDim = 1 : nDim ];
	
	numBfieldLst = zeros(Int64, itNum, nDim);
	numLinkLst = zeros(Int64, itNum, nDim);
	
	probInit = exp( -cAreaInit ) / ( 1 + exp( -cAreaInit ) );
	distInit = Binomial( 1, probInit );
	if !isInit0
		for dim = 1 : nDim
			rand!( distInit, BfieldLst[dim] );
		end
	end
	updateLinkFrom0ByB( BfieldLst, linkLst, linkFerroLst, params );
	
	BfieldSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	linkSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	BfieldStartSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	linkStartSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];

	updater = updaterType( params );
	
	itSample = 1;
	for it = 1 : itNum
		print( "it = ", it, "         \r" )
		for dim = 1 : nDim
			@view(zakLstLst[:,:,dim,it]) .= dropdims( reduce( xor, BfieldLst[dim]; dims = dim, init = false ); dims = dim );
		end
		
		if mod( it, itStep ) == 0
			Threads.@threads for dim = 1 : nDim
				linkSampleLst[itSample][dim] .= linkLst[dim];
				BfieldSampleLst[itSample][dim] .= BfieldLst[dim];
			end
			itSample += 1;
		end
		
		if it <= itStartSample
			Threads.@threads for dim = 1 : nDim
				linkStartSampleLst[it][dim] .= linkLst[dim];
				BfieldStartSampleLst[it][dim] .= BfieldLst[dim];
			end
		end
		
		updateLoops( updater, flipChecker, BfieldLst, linkLst, linkFerroLst, params );
		
		for dim = 1 : nDim
			numBfieldLst[it,dim] = sum(BfieldLst[dim]);
			numLinkLst[it,dim] = sum( linkLst[dim] );
		end
		# @infiltrate
	end
	
	xyDims = (1,2);
	zakMeanLst = dropdims( mean( zakLstLst; dims = xyDims ); dims = xyDims );
	
	zakLstSampleLst = zakLstLst[:,:,:,[itStep:itStep:itNum;]];
	
	return numBfieldLst, numLinkLst, BfieldSampleLst, linkSampleLst, BfieldStartSampleLst, linkStartSampleLst, zakMeanLst, zakLstSampleLst;
end

function loops_MC_methods( divNum = 64, itNum = 10000; updaterType::Type{<:LoopsUpdater}, flipChecker::FlipChecker, fMod = "", isBeta = false, beta = nothing, itNumSample = 100, itStartSample = 50, isInit0 = false )
	cFerroSigned = cFerro;
	nDim = 3;
	params = ParamsLoops( divNum, nDim );
	nDimLayer = nDim-1;
	divLst = params.divLst;
	divTup = Tuple(divLst);
	posLst = params.posLst;
	posLstShLst = params.posLstShLst;
	linkDimLst = params.linkDimLst;
	linkDimShLst = params.linkDimShLst;
	
	itStep = max( Int64( floor(itNum / itNumSample) ), 1 );
	lnSample = Int64( floor( itNum / itStep ) );
	itStartSample = min( itStartSample, itNum );
	
	zakLstLst = zeros( Bool, divNum, divNum, nDim, itNum );

	BfieldLst = [ zeros( Bool, divTup ) for dim = 1 : nDim ];
	linkLst = [ zeros( Bool, divTup ) for dim = 1 : nDim ];
	linkFerroLst = [ zeros( Bool, divTup ) for lnkDim = 1 : nDim-1, BDim = 1 : nDim ];
	
	numBfieldLst = zeros(Int64, itNum, nDim);
	numLinkLst = zeros(Int64, itNum, nDim);
	
	probInit = exp( -cArea ) / ( 1 + exp( -cArea ) );
	distInit = Binomial( 1, probInit );
	if !isInit0
		for dim = 1 : nDim
			rand!( distInit, BfieldLst[dim] );
		end
	end
	updateLinkFrom0ByB( BfieldLst, linkLst, linkFerroLst, params );
	
	BfieldSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	linkSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	BfieldStartSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	linkStartSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	
	pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );

	updater = updaterType( params; cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
	
	itSample = 1;
	for it = 1 : itNum
		print( "it = ", it, "         \r" )
		for dim = 1 : nDim
			@view(zakLstLst[:,:,dim,it]) .= dropdims( reduce( xor, BfieldLst[dim]; dims = dim, init = false ); dims = dim );
		end
		
		if mod( it, itStep ) == 0
			Threads.@threads for dim = 1 : nDim
				linkSampleLst[itSample][dim] .= linkLst[dim];
				BfieldSampleLst[itSample][dim] .= BfieldLst[dim];
			end
			itSample += 1;
		end
		
		if it <= itStartSample
			Threads.@threads for dim = 1 : nDim
				linkStartSampleLst[it][dim] .= linkLst[dim];
				BfieldStartSampleLst[it][dim] .= BfieldLst[dim];
			end
		end
		
		updateLoops( updater, flipChecker, BfieldLst, linkLst, linkFerroLst, params );
		
		for dim = 1 : nDim
			numBfieldLst[it,dim] = sum(BfieldLst[dim]);
			numLinkLst[it,dim] = sum( linkLst[dim] );
		end
	end
	
	xyDims = (1,2);
	zakMeanLst = dropdims( mean( zakLstLst; dims = xyDims ); dims = xyDims );
	
	zakLstSampleLst = zakLstLst[:,:,:,[itStep:itStep:itNum;]];
	
	fModOut = getFModLoopsMC( fMod, updaterType, isInit0 );
	
	fMain = fMainLoopsMC;
	valLst = Any[divNum, itNum, cArea, cPerim];
	attrLst, valLst = getAttrValLstLoopsMC( divNum, itNum, cArea, cPerim; beta = beta, cFerro = cFerro );
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( fName, "divNum", divNum, "itNum", itNum, "cArea", cArea, "cPerim", cPerim, "beta", beta );
	
	fMainZakLst = fMainLoopsMC * "_" * "zakLstAllMean";
	fMainZakSample = fMainLoopsMC * "_" * "zakLstSampleMean";
	
	fNameZakLst = fNameFunc( fMainZakLst, attrLst, valLst, jld2Type; fMod = fModOut );
	fNameZakSample = fNameFunc( fMainZakSample, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( fNameZakLst, "zakLstLst", zakLstLst, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	save( fNameZakSample, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	
	oFNameLoops = fNameFunc( oFNameLoopsMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoops, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
	
	oFNameLoopsNum = fNameFunc( oFNameLoopsNumMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsNum, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst );
	
	oFNameLoopsStart = fNameFunc( oFNameLoopsStartMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsStart, "linkStartSampleLst", linkStartSampleLst, "BfieldStartSampleLst", BfieldStartSampleLst );
	
	return fName;
end

function loops_MC_methods( divNum = 64, itNum = 10000; updaterType::(Type{T} where T <: LoopsUpdater), fMod = "", cArea = 1, cPerim = 1, cFerro = 0, isBeta = false, beta = nothing, itNumSample = 100, itStartSample = 50, isInit0 = false )
	cFerroSigned = cFerro;
	nDim = 3;
	params = ParamsLoops( divNum, nDim );
	nDimLayer = nDim-1;
	divLst = params.divLst;
	divTup = Tuple(divLst);
	posLst = params.posLst;
	posLstShLst = params.posLstShLst;
	linkDimLst = params.linkDimLst;
	linkDimShLst = params.linkDimShLst;
	
	itStep = max( Int64( floor(itNum / itNumSample) ), 1 );
	lnSample = Int64( floor( itNum / itStep ) );
	itStartSample = min( itStartSample, itNum );
	
	zakLstLst = zeros( Bool, divNum, divNum, nDim, itNum );

	BfieldLst = [ zeros( Bool, divTup ) for dim = 1 : nDim ];
	linkLst = [ zeros( Bool, divTup ) for dim = 1 : nDim ];
	linkFerroLst = [ zeros( Bool, divTup ) for lnkDim = 1 : nDim-1, BDim = 1 : nDim ];
	
	numBfieldLst = zeros(Int64, itNum, nDim);
	numLinkLst = zeros(Int64, itNum, nDim);
	
	probInit = exp( -cArea ) / ( 1 + exp( -cArea ) );
	distInit = Binomial( 1, probInit );
	if !isInit0
		for dim = 1 : nDim
			rand!( distInit, BfieldLst[dim] );
		end
	end
	updateLinkFrom0ByB( BfieldLst, linkLst, linkFerroLst, params );
	
	BfieldSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	linkSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : lnSample ];
	BfieldStartSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	linkStartSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : nDim ] for itSample = 1 : itStartSample ];
	
	pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );

	updater = updaterType( params; cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
	
	itSample = 1;
	for it = 1 : itNum
		print( "it = ", it, "         \r" )
		for dim = 1 : nDim
			@view(zakLstLst[:,:,dim,it]) .= dropdims( reduce( xor, BfieldLst[dim]; dims = dim, init = false ); dims = dim );
		end
		
		if mod( it, itStep ) == 0
			Threads.@threads for dim = 1 : nDim
				linkSampleLst[itSample][dim] .= linkLst[dim];
				BfieldSampleLst[itSample][dim] .= BfieldLst[dim];
			end
			itSample += 1;
		end
		
		if it <= itStartSample
			Threads.@threads for dim = 1 : nDim
				linkStartSampleLst[it][dim] .= linkLst[dim];
				BfieldStartSampleLst[it][dim] .= BfieldLst[dim];
			end
		end
		
		updateLoops( updater, BfieldLst, linkLst, linkFerroLst, params );
		
		for dim = 1 : nDim
			numBfieldLst[it,dim] = sum(BfieldLst[dim]);
			numLinkLst[it,dim] = sum( linkLst[dim] );
		end
	end
	
	xyDims = (1,2);
	zakMeanLst = dropdims( mean( zakLstLst; dims = xyDims ); dims = xyDims );
	
	zakLstSampleLst = zakLstLst[:,:,:,[itStep:itStep:itNum;]];
	
	fModOut = getFModLoopsMC( fMod, updaterType, isInit0 );
	
	fMain = fMainLoopsMC;
	valLst = Any[divNum, itNum, cArea, cPerim];
	attrLst, valLst = getAttrValLstLoopsMC( divNum, itNum, cArea, cPerim; beta = beta, cFerro = cFerro );
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( fName, "divNum", divNum, "itNum", itNum, "cArea", cArea, "cPerim", cPerim, "beta", beta );
	
	fMainZakLst = fMainLoopsMC * "_" * "zakLstAllMean";
	fMainZakSample = fMainLoopsMC * "_" * "zakLstSampleMean";
	
	fNameZakLst = fNameFunc( fMainZakLst, attrLst, valLst, jld2Type; fMod = fModOut );
	fNameZakSample = fNameFunc( fMainZakSample, attrLst, valLst, jld2Type; fMod = fModOut );
	
	save( fNameZakLst, "zakLstLst", zakLstLst, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	save( fNameZakSample, "zakMeanLst", zakMeanLst, "zakLstSampleLst", zakLstSampleLst );
	
	oFNameLoops = fNameFunc( oFNameLoopsMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoops, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
	
	oFNameLoopsNum = fNameFunc( oFNameLoopsNumMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsNum, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst );
	
	oFNameLoopsStart = fNameFunc( oFNameLoopsStartMain, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsStart, "linkStartSampleLst", linkStartSampleLst, "BfieldStartSampleLst", BfieldStartSampleLst );
	
	return fName;
end

function updateLinkFrom0ByB( BfieldLst, linkLst, linkFerroLst, params::ParamsLoops )
	for dim = 1 : params.nDim, pos in params.posLst
		if BfieldLst[dim][pos]
			for dimLink in params.linkDimLst[dim]
				linkLst[dimLink][pos] = !linkLst[dimLink][pos];
				for dimLinkSh in params.linkDimLst[dim]
					if dimLinkSh != dimLink
						linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] = !linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
					end
				end
			end
			for lnkDim = 1 : params.nDim-1
				linkFerroLst[lnkDim,dim][pos] = !linkFerroLst[lnkDim,dim][pos];
				dimLink = params.linkDimLst[dim][lnkDim];
				for dimLinkSh in params.linkDimLst[dim]
					if dimLinkSh != dimLink
						linkFerroLst[lnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]] = !linkFerroLst[lnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]];
					end
				end
			end
			# @infiltrate
		end
	end
end

function genPFlipLst( ; cArea, cPerim, cFerroSigned = 0 )
	numFlipTypes = 2;
	numPerimTypes = 2*numFlipTypes+1;
	flipStep = 2;
	pFlipLst = zeros( numPerimTypes, numPerimTypes, numFlipTypes );
	Bval = 1;
	for iB = 1 : numFlipTypes
		lnkVal = numFlipTypes * flipStep;
		for iL = 1 : numPerimTypes
			lnkFerroVal = numFlipTypes * flipStep;
			for iLFerro = 1 : numPerimTypes
				Eval = - ( cArea * Bval + cPerim * lnkVal + cFerroSigned * lnkFerroVal );
				dE = -2*Eval;
				pFlipLst[iLFerro, iL, iB] = 1 / ( 1 + exp( dE ) );
				lnkFerroVal -= flipStep;
			end
			lnkVal -= flipStep;
		end
		Bval *= -1;
	end
	
	return pFlipLst;
end

function flipBLinkAtPos( params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}; pos::CartesianIndex{D}, dim::Int64 ) where {D}
	BfieldLst[dim][pos] = !BfieldLst[dim][pos];
	for iLnkDim in 1 : params.nDimLayer
		dimLink = params.linkDimLst[dim][iLnkDim];
		dimLinkSh = params.linkDimShLst[dim][iLnkDim];
		linkLst[dimLink][pos] = !linkLst[dimLink][pos];
		linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] = !linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
	end
	for iLnkDim = 1 : params.nDimLayer
		linkFerroLst[iLnkDim,dim][pos] = !linkFerroLst[iLnkDim,dim][pos];
		dimLink = params.linkDimLst[dim][iLnkDim];
		dimLinkSh = params.linkDimShLst[dim][iLnkDim];
		linkFerroLst[iLnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]] = !linkFerroLst[iLnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]];
	end
end

function boolToOnePN( varBool::Bool )
	# return -(-1).^varBool;
	return varBool ? 1 : -1;
end

function boolToFlipChange( sBool::Bool )
	return - boolToOnePN( sBool )
end

function ratioToBinId( ratio::Number, divNum::Int64 )
	return Int64( min( floor( ratio * divNum ) + 1, divNum  ) );
end

include("loops_MC_run_funcs.jl")
export runLoopMC_withParams

include("loops_MC_methods_seperate_funcs.jl")
export loops_MC, loops_MC_smart, loops_MC_staggeredCube

include("loops_MC_resave_funcs.jl")

include("loops_MC_old_funcs.jl")

include("loops_MC_attr_printing.jl");

end #endmodule
