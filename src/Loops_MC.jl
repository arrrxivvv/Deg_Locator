module Loops_MC

using ShiftedArrays
using FilenameManip
using Random
using JLD2
using Statistics
using Distributions
using Utils

const rndDigsLpsMC = 3
 
# using Infiltrator

export loops_MC_methods, loops_MC, loops_MC_smart, loops_MC_staggeredCube

oFNameLoopsMain = "loopsSample";
oFNameLoopsStartMain = "loopsStartSample";
oFNameLoopsNumMain = "loopsNum";

dirLog = "./log/";

fMainLoopsMC = "loops_MC";
oFMainLoopsSample = "loopsSample";
oFMainLoopsStart = "loopsStartSample";
oFMainLoopsNum = "loopsNum";
attrLstLttcBase = ["div", "it", "nDim"];
getAttrLstLttcBase() = deepcopy(attrLstLttcBase);

function genAttrValLstLttc( divNum::Int64, itNum::Int64, nDim::Int64 )
	attrLst = getAttrLstLttcBase();
	valLst = Any[divNum,itNum,nDim];
	
	return attrLst, valLst;
end

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
	nDimB::Int64;
	divNum::Int64;
	grdNum::Int64;
	divLst::Vector{Int64};
	posLst::CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}};
	posLstShLst:: Matrix{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}}}; 
	linkDimLst::Vector{Vector{Int64}};
	linkDimShLst::Vector{Vector{Int64}};
	
	function ParamsLoops( divNum::Int64, nDim::Int64 )
		# nDimLayer = nDim - 1;
		nDimLayer = 2;
		nDimB = Int64( nDim * (nDim - 1) / 2 );
		divLst = fill(divNum, nDim);
		grdNum = divNum^nDim;
		posLst = CartesianIndices(ntuple(x->divNum,nDim));
		posLstShLst = Utils.arrShAdvRetLstFunc( posLst, nDim );
		nDimLink = 2;
		linkDimLst = [ zeros(Int64, nDimLink) for dim = 1 : nDimB ];
		if nDim == 2
			dimBLst = [3];
		elseif nDim == 3
			dimBLst = [1:3;];
		else
			dimBLst = zeros(Int64,0);
		end
		for dim = 1:nDimB
			dimB = dimBLst[dim];
			iDLink = 1;
			for dimLink = 1 : nDim
				if dimLink != dimB
					linkDimLst[dim][iDLink] = dimLink;
					iDLink += 1;
				end
			end
		end
		linkDimShLst = [ circshift( linkDimLst[dim], 1 ) for dim = 1 : nDimB ];
		
		new{nDim}( nDim, nDimLayer, nDimB, divNum, grdNum, divLst, posLst, posLstShLst, linkDimLst, linkDimShLst );
	end
end

function genBfieldLinkArr( params::ParamsLoops )
	divTup = Tuple( params.divLst );
	
	BfieldLst = [ zeros( Bool, divTup ) for dim = 1 : params.nDimB ];
	linkLst = [ zeros( Bool, divTup ) for dim = 1 : params.nDim ];
	linkFerroLst = [ zeros( Bool, divTup ) for lnkDim = 1 : params.nDimLayer, BDim = 1 : params.nDimB ];
	
	return BfieldLst, linkLst, linkFerroLst;
end

function genBfieldLinkNumLst( params::ParamsLoops, itNum::Int64 )
	numBfieldLst = zeros(Int64, itNum, params.nDimB);
	numLinkLst = zeros(Int64, itNum, params.nDim);
	
	return numBfieldLst, numLinkLst
end

function genBfieldLinkArrSample( params::ParamsLoops, lnSample )
	divTup = Tuple(params.divLst);
	
	BfieldSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : params.nDimB ] for itSample = 1 : lnSample ];
	linkSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : params.nDim ] for itSample = 1 : lnSample ];
	
	return BfieldSampleLst, linkSampleLst
end

abstract type BLinkInitializer end

function getInitializerName( initType::Type{<:BLinkInitializer} )
	error("initializer undefined");
end

function getInitializerName( initializer::BLinkInitializer )
	return getInitializerName( typeof(initializer) );
end

function getAttrValInitializer( initializer::BLinkInitializer; rndDigs = rndDigsLpsMC )
	error("initializer undefined");
end

function initializeBL( initializer::BLinkInitializer, BfieldLst, params::ParamsLoops )
	error("initializer undefined");
end

struct BinomialInitializer <: BLinkInitializer
	prob::Float64;
	dist::Distribution;
end

function genMeanFieldInitializer( cArea::Real )
	prob = exp( -cArea ) / ( 1 + exp( -cArea ) );
	
	return genProbInitializer( prob );
end

function genProbInitializer( prob::Real )
	dist = Binomial(1, prob);
	
	return BinomialInitializer( prob, dist );
end

function getInitializerName( initType::Type{BinomialInitializer} )
	return "BinomialInit";
end

function getAttrValInitializer( initializer::BinomialInitializer; rndDigs = rndDigsLpsMC )
	attrLst = [getInitializerName( initializer )];
	valLst = roundKeepInt.( [initializer.prob]; digits = rndDigs );
	
	return attrLst, valLst;
end

function initializeBL( initializer::BinomialInitializer, BfieldLst, params::ParamsLoops )
	for dim = 1 : params.nDimB
		rand!( initializer.dist, BfieldLst[dim] );
	end
end

struct ConstantInitializer <: BLinkInitializer
	initVal::Bool;
end

ConstantInitializer() = ConstantInitializer(false);

function initializeBL( initializer::ConstantInitializer, BfieldLst, params::ParamsLoops )
	for dim = 1 : params.nDimB
		BfieldLst[dim] .= initializer.initVal;
	end
end

function getInitializerName( initType::Type{ConstantInitializer} )
	return "ConstInit";
end

function getAttrValInitializer( initializer::ConstantInitializer; rndDigs = rndDigsLpsMC )
	attrLst = [getInitializerName( initializer )];
	valLst = roundKeepInt.( [initializer.initVal]; digits = rndDigs );
	
	return attrLst, valLst;
end

abstract type FlipChecker end

function flipCheck( flipChecker::Type{<:FlipChecker}, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	error( "Loops_MC: FlipChecker not defined yet" );
end

function getFlipCheckerName( flipCheckerType::Type{<:FlipChecker} )
	error( "Loops_MC: FlipChecker not defined yet" );
end

function getFlipCheckerName( flipChecker::FlipChecker )
	return getFlipCheckerName( typeof(flipChecker) );
end

function getFlipCheckerAttrLst( flipChecker::FlipChecker )
	error( "Loops_MC: FlipChecker not defined yet" );
end

# function getFlipCheckerAttrLst( flipChecker::FlipChecker )
	# return getFlipCheckerAttrLst( typeof(flipChecker) );
# end

function getFlipCheckerAttrValLst( flipChecker::FlipChecker; rndDigs = rndDigsLpsMC )
	error( "Loops_MC: FlipChecker not defined yet" );
end

function genAttrLstLttcFlipChecker( divNum::Int64, itNum::Int64, nDim::Int64, flipChecker::FlipChecker; rndDigs = rndDigsLpsMC )
	attrLst, valLst = genAttrValLstLttc( divNum, itNum, nDim );
	
	attrLstFlip, valLstFlip = getFlipCheckerAttrValLst( flipChecker; rndDigs = rndDigs );
	
	append!( attrLst, attrLstFlip );
	append!( valLst, valLstFlip );
	
	return attrLst, valLst;
end

function genAttrLstLttcFlipInit( divNum::Int64, itNum::Int64, nDim::Int64, flipChecker::FlipChecker, initializer::BLinkInitializer, rndDigs = rndDigsLpsMC )
	attrLst, valLst = genAttrLstLttcFlipChecker( divNum, itNum, nDim, flipChecker; rndDigs = rndDigs );
	
	attrLstInit, valLstInit = getAttrValInitializer( initializer );
	append!( attrLst, attrLstInit );
	append!( valLst, valLstInit );
	
	return attrLst, valLst;
end

abstract type AbstractNeighborFlipChecker <: FlipChecker end

function getFlipCheckerAttrLst( flipChecker::AbstractNeighborFlipChecker )
	return ["cArea", "cPerim", "cFerro"];
end

struct NeighborFlipChecker <: AbstractNeighborFlipChecker
	pFlipLst::Array{Float64};
	cParamLst::Array{Real};
	
	function NeighborFlipChecker( cArea, cPerim, cFerroSigned = 0 )
		pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
		cParamLst = Real[cArea, cPerim, cFerroSigned];
		
		new(pFlipLst,cParamLst);
	end
end

function getFlipCheckerName( flipCheckerType::Type{NeighborFlipChecker} )
	return "cAcLcFFlip";
end

# function getFlipCheckerAttrLst( flipChecker::NeighborFlipChecker )
	# return ["cArea", "cPerim", "cFerro"];
# end

function getFlipCheckerAttrValLst( flipChecker::AbstractNeighborFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.( deepcopy( flipChecker.cParamLst ); digits = rndDigs );
	if flipChecker.cParamLst[3] == 0
		attrLst = attrLst[1:2];
		valLst = valLst[1:2];
	end
	
	return attrLst, valLst;
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

struct CubeFlipChecker <: FlipChecker
	pFlipLst::Array{Float64};
	cArea::Real;
	
	function CubeFlipChecker( cArea, cPerim, cFerroSigned = 0 )
		pFlipLst = genPFlipLstCube( cArea = cArea );
		
		new(pFlipLst,cArea);
	end
end

function getFlipCheckerName( flipCheckerType::Type{CubeFlipChecker} )
	return "cAcLcFCubeFlip";
end

function getFlipCheckerAttrLst( flipChecker::CubeFlipChecker )
	return ["cArea"];
end

function getFlipCheckerAttrValLst( flipChecker::CubeFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.( [flipChecker.cArea]; digits = rndDigs );
	
	return attrLst, valLst;
end

function flipCheck( flipChecker::CubeFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	iArea = 1;
	for dimB = 1 : params.nDimB
		iArea += BfieldLst[dimB][pos];
		iArea += BfieldLst[dimB][params.posLstShLst[dimB,1][pos]];
	end
	
	pSwitchRand = rand();
	return pSwitchRand < flipChecker.pFlipLst[iArea];
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

function getUpdaterFMod( updaterType::Type{SingleUpdater} )
	return "upSingle";
end

function updateLoops( updater::SingleUpdater, flipChecker::FlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	pos = rand(params.posLst);
	dim = rand(1:params.nDimB);
	
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

struct AB2dUpdater{D} <: LoopsUpdater
	posABLst::Vector{Vector{CartesianIndex{D}}};
	
	function AB2dUpdater( params::ParamsLoops )
		if params.nDim != 2
			error( "nDim != 2" );
		end
		if iseven(params.divNum)
			szAB = Int64(params.divNum.^params.nDim/2);
			posABLst = [ Vector{CartesianIndex{params.nDim}}(undef,szAB) for iAB = 1 : 2 ];
			iA = 1; iB = 1;
			for ii = 1 : params.divNum, jj = 1 : params.divNum
				if (ii+jj)%2 != 0
					posABLst[1][iA] = params.posLst[ii,jj];
					iA+=1;
				else
					posABLst[2][iB] = params.posLst[ii,jj];
					iB+=1;
				end
			end
		else
			posABLst = [ params.posLst[iAB:2:end] for iAB = 1 : 2 ];
		end
		
		new{params.nDim}( posABLst );
	end
end

function getUpdaterFMod( updaterType::Type{AB2dUpdater} )
	return "upAB2d";
end

function updateLoops( updater::AB2dUpdater, flipChecker::FlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	if params.nDim != 2
		error( "nDim != 2" );
	end
	
	dim = 1;
	for iAB = 1 : 2
		for pos in updater.posABLst[iAB]
			if flipCheck( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
				flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = dim );
			end
		end
	end
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
	
	function StaggeredCubeUpdaterBase( params::ParamsLoops )
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

function getUpdaterFMod( updaterType::Type{StaggeredCubeUpdaterBase} )
	return "upStagCube";
end

function updateLoops( updater::StaggeredCubeUpdaterBase, flipChecker::FlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
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

struct CubeUpdater <: LoopsUpdater
	
end

CubeUpdater( params::ParamsLoops ) = CubeUpdater();

function getUpdaterFMod( updaterType::Type{CubeUpdater} )
	return "cubeOnly";
end

function updateLoops( updater::CubeUpdater, flipChecker::FlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	pos = rand( params.posLst );
	dim = 1;
	
	if flipCheck( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
		for dimB = 1 : params.nDimB
			BfieldLst[dimB][pos] = !BfieldLst[dimB][pos];
			BfieldLst[dimB][params.posLstShLst[dimB,1][pos]] = !BfieldLst[dimB][params.posLstShLst[dimB,1][pos]];
		end
	end
end

function getFModLoopsMC( fMod::String, updaterType::Type{<:LoopsUpdater}, isInit0::Bool=false; isFModMethod = true, flipChecker::Union{FlipChecker,Nothing} = nothing )
	fModOut = fMod;
	if isFModMethod
		fModOut = Utils.strAppendWith_( fModOut, getUpdaterFMod(updaterType) );
	end
	if !isnothing( flipChecker )
		fModOut = Utils.strAppendWith_( fModOut, getFlipCheckerName(flipChecker) );
	end
	if isInit0
		fModOut = Utils.strAppendWith_( fModOut, "isInit0" );
	end
	
	return fModOut;
end

function getFModLoopsMC( fMod::String, updaterType::Union{Type{<:LoopsUpdater},LoopsUpdater}; isInit0::Bool=false, isFModMethod = true, flipChecker::Union{FlipChecker,Nothing} = nothing )
	fModOut = fMod;
	if isFModMethod
		fModOut = Utils.strAppendWith_( fModOut, getUpdaterFMod(updaterType) );
	end
	if !isnothing( flipChecker )
		fModOut = Utils.strAppendWith_( fModOut, getFlipCheckerName(flipChecker) );
	end
	if isInit0
		fModOut = Utils.strAppendWith_( fModOut, "isInit0" );
	end
	
	return fModOut;
end


function loops_MC_methods_cALFCube( divNum = 64, itNum = 10000; updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, fMod = "", cArea = 1, cPerim = 1, cFerro = 0, itNumSample = 100, itStartSample = 50, nDim = 3, probInit = nothing, isFileNameOnly::Bool = false, fMainOutside::String = "" )
	flipChecker = CubeFlipChecker( cArea, cPerim, cFerro );
	if initializerType == BinomialInitializer
		if isnothing(probInit)
			initializer = genMeanFieldInitializer( cArea );
		else
			initializer = genProbInitializer( probInit );
		end
	elseif initializerType == ConstantInitializer
		initializer = ConstantInitializer();
	end
	
	fNameOutside = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, initializer = initializer, fMod = fMod, itNumSample = itNumSample, itStartSample = itStartSample, isFileNameOnly, fMainOutside );
	
	return fNameOutside;
end

function loops_MC_methods_cALF( divNum = 64, itNum = 10000; updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, fMod = "", cArea = 1, cPerim = 1, cFerro = 0, itNumSample = 100, itStartSample = 50, nDim = 3, probInit = nothing, isFileNameOnly::Bool = false, fMainOutside::String = "" )
	flipChecker = NeighborFlipChecker( cArea, cPerim, cFerro );
	if initializerType == BinomialInitializer
		if isnothing(probInit)
			initializer = genMeanFieldInitializer( cArea );
		else
			initializer = genProbInitializer( probInit );
		end
	elseif initializerType == ConstantInitializer
		initializer = ConstantInitializer();
	end
	
	fNameOutside = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, initializer = initializer, fMod = fMod, itNumSample = itNumSample, itStartSample = itStartSample, isFileNameOnly, fMainOutside = fMainOutside, nDim = nDim );
	
	return fNameOutside;
end

function loops_MC_methods_Base( divNum = 64, itNum = 10000; updaterType::Type{<:LoopsUpdater}, flipChecker::FlipChecker, initializer::BLinkInitializer, fMod = "", itNumSample = 100, itStartSample = 50, nDim = 3, isFileNameOnly::Bool = false, fMainOutside::Union{String, Nothing}= "" )
	fModOut = getFModLoopsMC( fMod, updaterType; flipChecker = flipChecker );
	
	fMain = fMainLoopsMC;
	attrLst, valLst = genAttrLstLttcFlipInit( divNum, itNum, nDim, flipChecker, initializer );
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
	
	if isFileNameOnly
		fNameOutside = fNameFunc( fMainOutside, attrLst, valLst, jld2Type; fMod = fModOut );
		return fNameOutside;
	end
	
	params = ParamsLoops( divNum, nDim );
	
	BfieldLst, linkLst, linkFerroLst = genBfieldLinkArr( params );
	
	itStep = max( Int64( floor(itNum / itNumSample) ), 1 );
	lnSample = Int64( floor( itNum / itStep ) );
	itStartSample = min( itStartSample, itNum );
	
	numBfieldLst, numLinkLst = genBfieldLinkNumLst( params, itNum );
	BfieldSampleLst, linkSampleLst = genBfieldLinkArrSample( params, lnSample );
	BfieldStartSampleLst, linkStartSampleLst = genBfieldLinkArrSample( params, itStartSample );

	initializeBL( initializer, BfieldLst, params );
	updateLinkFrom0ByBAllDims( BfieldLst, linkLst, linkFerroLst, params );
	updater = updaterType( params );
	
	itSample = 1;
	for it = 1 : itNum
		print( "it = ", it, "         \r" )
		
		if mod( it, itStep ) == 0
			Threads.@threads for dim = 1 : params.nDim
				linkSampleLst[itSample][dim] .= linkLst[dim];
			end
			Threads.@threads for dim = 1 : params.nDimB
				BfieldSampleLst[itSample][dim] .= BfieldLst[dim];
			end
			itSample += 1;
		end
		
		if it <= itStartSample
			Threads.@threads for dim = 1 : params.nDim
				linkStartSampleLst[it][dim] .= linkLst[dim];
			end
			Threads.@threads for dim = 1 : params.nDimB
				BfieldStartSampleLst[it][dim] .= BfieldLst[dim];
			end
		end
		
		updateLoops( updater, flipChecker, BfieldLst, linkLst, linkFerroLst, params );
		
		for dim = 1 : params.nDim
			numLinkLst[it,dim] = sum( linkLst[dim] );
		end
		for dim = 1 : params.nDimB
			numBfieldLst[it,dim] = sum(BfieldLst[dim]);
		end
	end
	
	save( fName, "divNum", divNum, "itNum", itNum );
	
	oFNameLoops = fNameFunc( oFMainLoopsSample, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoops, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst, "linkSampleLst", linkSampleLst, "BfieldSampleLst", BfieldSampleLst );
	
	oFNameLoopsNum = fNameFunc( oFMainLoopsNum, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsNum, "numBfieldLst", numBfieldLst, "numLinkLst", numLinkLst );
	
	oFNameLoopsStart = fNameFunc( oFMainLoopsStart, attrLst, valLst, jld2Type; fMod = fModOut );
	save( oFNameLoopsStart, "linkStartSampleLst", linkStartSampleLst, "BfieldStartSampleLst", BfieldStartSampleLst );
	
	if isFileNameOnly
		return fNameOutside;
	end
	
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

function updateLinkFrom0ByBAllDims( BfieldLst, linkLst, params::ParamsLoops )
	for dim = 1 : params.nDimB, pos in params.posLst
		if BfieldLst[dim][pos]
			for iLnkDim in 1 : params.nDimLayer
				dimLink = params.linkDimLst[dim][iLnkDim];
				dimLinkSh = params.linkDimShLst[dim][iLnkDim];
				linkLst[dimLink][pos] = !linkLst[dimLink][pos];
				linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] = !linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]];
			end
		end
	end
end

function updateLinkFrom0ByBAllDims( BfieldLst, linkLst, linkFerroLst, params::ParamsLoops )
	updateLinkFrom0ByBAllDims( BfieldLst, linkLst, params );
	for dim = 1 : params.nDimB, pos in params.posLst
		if BfieldLst[dim][pos]
			for iLnkDim in 1 : params.nDimLayer
				dimLink = params.linkDimLst[dim][iLnkDim];
				dimLinkSh = params.linkDimShLst[dim][iLnkDim];
				linkFerroLst[iLnkDim,dim][pos] = !linkFerroLst[iLnkDim,dim][pos];
				linkFerroLst[iLnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]] = !linkFerroLst[iLnkDim,dim][params.posLstShLst[dimLinkSh,1][pos]];
			end
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

function genPFlipLstCube( ; cArea )
	nDimB = 3;
	numFlipTypes = 2*nDimB+1;
	pFlipLst = 1 ./ (1 .+ Float64[2*nDimB : -2 : -2*nDimB;] .* cArea );
	
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

include("loops_MC_WangLandau.jl")

include("loops_MC_run_funcs.jl")
export runLoopMC_withParams

include("loops_MC_run_funcs_flipChecker.jl")

include("loops_MC_methods_seperate_funcs.jl")
export loops_MC, loops_MC_smart, loops_MC_staggeredCube

include("loops_MC_resave_funcs.jl")

include("loops_MC_old_funcs.jl")

include("loops_MC_attr_printing.jl");

include("loops_MC_2d.jl");

end #endmodule
