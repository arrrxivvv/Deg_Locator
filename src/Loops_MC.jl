module Loops_MC

using ShiftedArrays
using FilenameManip
using Random
using JLD2
using Statistics
using Distributions
using Utils

export loops_MC_methods, loops_MC, loops_MC_smart, loops_MC_staggeredCube

const rndDigsLpsMC = 3

const dirLog = "./log/";
 
using Infiltrator

const oFNameLoopsMain = "loopsSample";
const oFNameLoopsStartMain = "loopsStartSample";
const oFNameLoopsNumMain = "loopsNum";

const fMainLoopsMC = "loops_MC";
const oFMainLoopsSample = "loopsSample";
const oFMainLoopsStart = "loopsStartSample";
const oFMainLoopsNum = "loopsNum";
const attrLstLttcBase = ["div", "it", "nDim"];
const attrLstLttcNoItBase = ["div", "nDim"];
getAttrLstLttcBase() = deepcopy(attrLstLttcBase);
getAttrLstLttcNoItBase() = deepcopy(attrLstLttcNoItBase);

function genAttrValLstLttc( divNum::Int64, itNum::Int64, nDim::Int64 )
	attrLst = getAttrLstLttcBase();
	valLst = Any[divNum,itNum,nDim];
	
	return attrLst, valLst;
end

function genAttrValLstLttc( divNum::Int64, nDim::Int64 )
	attrLst = getAttrLstLttcNoItBase();
	valLst = Any[divNum,nDim];
	
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
	divTup::Tuple{Vararg{Int64,N}};
	posLst::CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}};
	posSlcLst::Vector{Vector{SubArray{CartesianIndex{N}}}};
	posLstShLst:: Matrix{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}}}; 
	linkDimLst::Vector{Vector{Int64}};
	linkDimShLst::Vector{Vector{Int64}};
	
	function ParamsLoops( divNum::Int64, nDim::Int64 )
		# nDimLayer = nDim - 1;
		nDimLayer = 2;
		nDimB = Int64( nDim * (nDim - 1) / 2 );
		divLst = fill(divNum, nDim);
		divTup = Tuple(divLst);
		grdNum = divNum^nDim;
		posLst = CartesianIndices(ntuple(x->divNum,nDim));
		posLstShLst = Utils.arrShAdvRetLstFunc( posLst, nDim );
		
		posSlcLst = [ [ selectdim( posLst, dim, iD ) for iD = 1 : divLst[dim] ] for dim = 1 : nDim ];
		
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
		
		new{nDim}( nDim, nDimLayer, nDimB, divNum, grdNum, divLst, divTup, posLst, posSlcLst, posLstShLst, linkDimLst, linkDimShLst );
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

function genBfieldLinkFerroNumLst( params::ParamsLoops, itNum::Int64 )
	numBfieldLst = zeros(Int64, itNum, params.nDimB);
	numLinkLst = zeros(Int64, itNum, params.nDim);
	numLinkFerroLst = zeros(Int64, itNum, params.nDimLayer, params.nDimB);
	
	return numBfieldLst, numLinkLst, numLinkFerroLst;
end

# function genBfieldLinkNumLst( params::ParamsLoops, itNum::Int64 )
	# numBfieldLst = zeros(Int64, itNum, params.nDimB);
	# numLinkLst = zeros(Int64, itNum, params.nDim);
	
	# return numBfieldLst, numLinkLst
# end

function genBfieldLinkArrSample( params::ParamsLoops, lnSample )
	divTup = Tuple(params.divLst);
	
	BfieldSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : params.nDimB ] for itSample = 1 : lnSample ];
	linkSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : params.nDim ] for itSample = 1 : lnSample ];
	
	return BfieldSampleLst, linkSampleLst
end

function genBfieldLinkFerroArrSample( params::ParamsLoops, lnSample )
	divTup = Tuple(params.divLst);
	
	BfieldSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : params.nDimB ] for itSample = 1 : lnSample ];
	linkSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : params.nDim ] for itSample = 1 : lnSample ];
	linkFerroSampleLst = [[ zeros( Bool, divTup ) for dim = 1 : params.nDimLayer, BDim = 1 : params.nDimB ] for itSample = 1 : lnSample ];
	
	return BfieldSampleLst, linkSampleLst, linkFerroSampleLst;
end




abstract type ItController end

throwItControllerUndefined() = error("Loops_MC: itController undefined");

getItControllerName( itControllerType::Type{<:ItController} ) = throwItControllerUndefined();
getItControllerName( itController::ItController ) = getItControllerName( typeof(itController) );

getAttrLstItController( itControllerType::Type{<:ItController} ) = throwItControllerUndefined();
getAttrLstItController( itController::ItController ) = getAttrLstItController( typeof(itController) );

getValLstItController( itController::ItController ) = throwItControllerUndefined();

function getAttrValLstItController( itController::ItController )
	attrLst = getAttrLstItController( itController );
	valLst = getValLstItController( itController );
	
	return attrLst, valLst;
end

function advanceItControl( itController::ItController )
	itController.itRef[] += 1;
end

function resetItControl( itController::ItController )
	itController.itRef[] = 1;
end

testItNotDone( itController::ItController ) = throwItControllerUndefined();

testItDoSample( itController::ItController ) = throwItControllerUndefined();

testItDoStartSample( itController::ItController ) = throwItControllerUndefined();

getItNumLst( itController::ItController ) = throwItControllerUndefined();



struct ItNumItController <: ItController
	itNum::Int64;
	itNumSample::Int64;
	itNumStartSample::Int64;
	
	itSampleStep::Int64;
	
	itRef::Ref{Int64}
	
	function ItNumItController( itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 )
		itSampleStep = max( 1, Int64( floor( itNum / itNumSample ) ) );
		itVal = 1;
		
		new( itNum, itNumSample, itNumStartSample, itSampleStep, Ref(itVal) );
	end
end




abstract type BLinkInitializer end

function throwInitializerUndefined()
	error("initializer undefined");
end

function getInitializerName( initType::Type{<:BLinkInitializer} )
	throwInitializerUndefined();
end

function getInitializerName( initializer::BLinkInitializer )
	return getInitializerName( typeof(initializer) );
end

function getAttrValInitializer( initializer::BLinkInitializer; rndDigs = rndDigsLpsMC )
	throwInitializerUndefined();
end

function initializeBL( initializer::BLinkInitializer, BfieldLst, params::ParamsLoops )
	throwInitializerUndefined();
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

struct ConstantInitializer <: BLinkInitializer
	initVal::Bool;
end

ConstantInitializer() = ConstantInitializer(false);


struct PresetInitializer{D} <: BLinkInitializer
	BfieldLst::Vector{Array{Bool,D}};
end




abstract type FlipProposer end

function throwFlipProposerUndefined()
	error("Loops_MC: FlipProposer not defined")
end

function flipDoIt( flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	throwFlipProposerUndefined();
end

getFlipProposerName( flipProposerType::Type{<:FlipProposer} ) = throwFlipProposerUndefined();
getFlipProposerName( flipProposer::FlipProposer ) = getFlipProposerName( typeof(flipProposer) );

struct SwitchingFlipProposer{T_tuple} <: FlipProposer
	flipProposerTup::Tuple{Vararg{FlipProposer}};
	counterRef::Ref{Int64};
	
	function SwitchingFlipProposer( flipProposerTup::Tuple{Vararg{FlipProposer}} )
		counterVal = 1;
		
		new{typeof(flipProposerTup)}( flipProposerTup, Ref(counterVal) );
	end
end

SwitchingFlipProposer( flipProposers::FlipProposer... ) = SwitchingFlipProposer( flipProposers );

function SwitchingFlipProposer{T_tuple}() where {T_tuple}
	flipProposerTup = ntuple( ii -> T_tuple.parameters[ii](), length(T_tuple.parameters) );
	
	return SwitchingFlipProposer( flipProposerTup );
end

function getSwitchLstLen( flipProposer::SwitchingFlipProposer )
	return length( flipProposer.flipProposerTup );
end

abstract type AbstractConcreteFlipProposer <: FlipProposer end

struct OneFlipProposer <: AbstractConcreteFlipProposer
	
end

struct CubeFlipProposer <: AbstractConcreteFlipProposer

end




abstract type FlipChecker end

abstract type AuxData end

function throwFlipCheckerUndefined()
	error( "Loops_MC: FlipChecker not defined yet" );
end

function flipCheck( flipChecker::FlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	throwFlipCheckerUndefined();
end

function flipCheck( flipChecker::FlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	throwFlipCheckerUndefined();
end

function flipDoIt( flipChecker::FlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	throwFlipCheckerUndefined();
end

function flipCheckDoIt( flipChecker::FlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	if flipCheck( flipChecker, flipProposer::FlipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
		flipDoIt( flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	end
end

function flipCheckDoIt( flipChecker::FlipChecker, flipProposer::FlipProposer, auxData::AuxData, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	if flipCheck( flipChecker, flipProposer::FlipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
		flipDoIt( flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
		flipAuxData!( auxData, flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	end
	# calcAuxData( auxData, params, BfieldLst, linkLst, linkFerroLst );
end

function flipCheckDoIt( flipChecker::FlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	if flipCheck( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
		flipDoIt( flipChecker, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	end
end

getFlipCheckerName( flipCheckerType::Type{<:FlipChecker} ) = throwFlipCheckerUndefined();

function getFlipCheckerName( flipChecker::FlipChecker )
	return getFlipCheckerName( typeof(flipChecker) );
end

getFlipCheckerAttrLst( flipCheckerType::Type{<:FlipChecker} ) = throwFlipCheckerUndefined();
getFlipCheckerAttrLst( flipChecker::FlipChecker ) = getFlipCheckerAttrLst( typeof(flipChecker) );

getFlipCheckerValLst( flipChecker::FlipChecker; rndDigs = rndDigsLpsMC ) = throwFlipCheckerUndefined();

# function getFlipCheckerAttrLst( flipChecker::FlipChecker )
	# return getFlipCheckerAttrLst( typeof(flipChecker) );
# end

function getFlipCheckerAttrValLst( flipChecker::FlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = getFlipCheckerValLst( flipChecker; rndDigs = rndDigs );
	
	return attrLst, valLst;
end



abstract type AbstractSwitchingFlipChecker <: FlipChecker end

struct SwitchingFlipChecker{T_tuple} <: AbstractSwitchingFlipChecker
	flipCheckerTup::Tuple{Vararg{FlipChecker}};
	
	function SwitchingFlipChecker( flipCheckerTup::Tuple{Vararg{FlipChecker}} )
		new{typeof(flipCheckerTup)}(flipCheckerTup);
	end
end

SwitchingFlipChecker( flipCheckers::FlipChecker... ) = SwitchingFlipChecker( flipCheckers );

function SwitchingFlipChecker{T_tuple}( cArea, cPerim, cFerro = 0 ) where {T_tuple}
	flipCheckerTup = ntuple( ii -> T_tuple.parameters[ii](cArea, cPerim, cFerro), length(T_tuple.parameters) );
	
	return SwitchingFlipChecker( flipCheckerTup );
end

abstract type AbstractNeighborFlipChecker <: FlipChecker end

struct NeighborFlipChecker <: AbstractNeighborFlipChecker
	pFlipLst::Array{Float64};
	cParamLst::Array{Real};
	
	function NeighborFlipChecker( cArea, cPerim, cFerroSigned = 0 )
		pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
		cParamLst = Real[cArea, cPerim, cFerroSigned];
		
		new(pFlipLst,cParamLst);
	end
end

struct CubeFlipChecker <: FlipChecker
	pFlipLst::Vector{Float64};
	cArea::Real;
	
	function CubeFlipChecker( cArea )
		pFlipLst = genPFlipLstCube( cArea = cArea );
		
		new(pFlipLst,cArea);
	end
end

CubeFlipChecker( cArea, cPerim, cFerro = 0 ) = CubeFlipChecker( cArea );

struct IsingFlipChecker <: FlipChecker
	pFlipLst::Array{Float64};
	pCubeFlipLst::Vector{Float64};
	cParamLst::Vector{Real};
	
	function IsingFlipChecker( cArea, cPerim, cFerroSigned = 0 )
		pFlipLst = genPFlipLst( cArea = cArea, cPerim = cPerim, cFerroSigned = cFerroSigned );
		cParamLst = Real[cArea, cPerim, cFerroSigned];
		pCubeFlipLst = genPFlipLstCube( cArea = cArea );
		
		new(pFlipLst,pCubeFlipLst,cParamLst);
	end
end





throwAuxDataUndefined() = error( "Loops_MC: AuxData not defined" );

genAuxData( auxDataType::Type{<:AuxData}, params::ParamsLoops, flipChecker::FlipChecker, itController::ItController ) = genAuxData( auxDataType, params, flipChecker, getItNumLst( itController )... );

genAuxData( auxDataType::Type{<:AuxData}, params::ParamsLoops, itController::ItController ) = genAuxData( auxDataType, params, getItNumLst( itController )... );

genAuxData( auxDataType::Type{<:AuxData}, flipChecker::FlipChecker, itController::ItController ) = genAuxData( auxDataType, flipChecker, getItNumLst( itController )... );

genAuxData( auxDataType::Type{<:AuxData}, params::ParamsLoops, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) = auxDataType( params, itNum, itNumSample, itNumStartSample );
genAuxData( auxDataType::Type{<:AuxData}, params::ParamsLoops, flipChecker::FlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) = auxDataType( params, flipChecker, itNum, itNumSample, itNumStartSample );
genAuxData( auxDataType::Type{<:AuxData}, flipChecker::FlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) = auxDataType( flipChecker, itNum, itNumSample, itNumStartSample );

# genAuxData( auxDataType::Type{<:AuxData}, params::ParamsLoops, itController::ItController ) = genAuxData( auxDataType, params, getItNumLst( itController )... );

getAuxDataSummaryName( auxDataType::Type{<:AuxData} ) = throwAuxDataUndefined();
getAuxDataSummaryName( auxData::AuxData ) = getAuxDataSummaryName( typeof(auxData) );

function getAuxDataSummaryItSampleLstName( auxDataType::Type{<:AuxData} )
	return "itSampleLst";
end
getAuxDataSummaryItSampleLstName( auxData::AuxData ) = getAuxDataSummaryItSampleLstName( typeof(auxData) );

function getAuxDataSummarySampleName( auxDataType::Type{<:AuxData} )
	return getAuxDataSummaryName( auxDataType ) * "Sample";
end
getAuxDataSummarySampleName( auxData::AuxData ) = getAuxDataSummarySampleName( typeof(auxData) );

function getAuxDataSummaryStartSampleName( auxDataType::Type{<:AuxData} )
	return getAuxDataSummaryName( auxDataType ) * "StartSample";
end
getAuxDataSummaryStartSampleName( auxData::AuxData ) = getAuxDataSummaryStartSampleName( typeof(auxData) );

function getAuxDataSummaryNumName( auxDataType::Type{<:AuxData} )
	return getAuxDataSummaryName( auxDataType ) * "Num";
end
getAuxDataSummaryNumName( auxData::AuxData ) = getAuxDataSummaryNumName( typeof(auxData) );

getAuxDataNameLst( auxDataType::Type{<:AuxData} ) = throwAuxDataUndefined();
getAuxDataNameLst( auxData::AuxData ) = getAuxDataNameLst( typeof(auxData) );

function getAuxDataSampleNameLst( auxDataType::Type{<:AuxData} )
	return getAuxDataNameLst(auxDataType) .*= "Sample";
end
getAuxDataSampleNameLst( auxData::AuxData ) = getAuxDataSampleNameLst( typeof(auxData) );

function getAuxDataStartSampleNameLst( auxDataType::Type{<:AuxData} )
	return getAuxDataNameLst(auxDataType) .*= "StartSample";
end
getAuxDataStartSampleNameLst( auxData::AuxData ) = getAuxDataStartSampleNameLst( typeof(auxData) );

function getAuxDataNumNameLst( auxDataType::Type{<:AuxData} )
	return getAuxDataNameLst(auxDataType) .*= "Num";
end
getAuxDataNumNameLst( auxData::AuxData ) = getAuxDataNumNameLst( typeof(auxData) );

renewAuxDataOutLst!( auxData::AuxData ) = throwAuxDataUndefined();

function renewAuxJldVarLst!( auxData::AuxData )
	renewAuxDataOutLst!(auxData);
	lstSplicing!( auxData.jldVarSampleLst, getAuxDataSampleNameLst( auxData ), auxData.dataSampleOutLst );
	lstSplicing!( auxData.jldVarStartSampleLst, getAuxDataStartSampleNameLst( auxData ), auxData.dataStartSampleOutLst );
	lstSplicing!( auxData.jldVarNumLst, getAuxDataNumNameLst( auxData ), auxData.dataNumOutLst );
	lstSplicing!( auxData.jldVarItSampleLst, ["itSampleLst"], [auxData.itSampleLst] )
end

function calcAuxData!( auxData::AuxData, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	throwAuxDataUndefined();
end

function flipAuxData!( auxData::AuxData, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	throwAuxDataUndefined();
end

function storeAuxDataSample( auxData::AuxData, itSample::Int64, it::Int64 )
	checkDataSampleAndItLstExtend!( auxData, itSample, it );
	
	storeAuxDataSampleNoBndCheck( auxData, itSample );
end

function checkDataSampleAndItLstExtend!( auxData::AuxData, itSample::Int64, it::Int64 )
	if itSample > length(auxData.itSampleLst)
		resize!( auxData.itSampleLst, itSample );
	end
	auxData.itSampleLst[itSample] = it;
	
	checkDataSampleExtend!( auxData::AuxData, itSample::Int64 );
end

function checkDataExtendHelper!( dataSampleLst::Vector, dataLst::Vector, itSample::Int64 )
	if isempty(dataSampleLst)
		return;
	elseif itSample > length( dataSampleLst[1] )
		for iD = 1 : length(dataSampleLst)
			# @infiltrate length(dataLst) == 2
			push!(dataSampleLst[iD], deepcopy(dataLst[iD]));
		end
	end
end

function checkDataSampleExtend!( auxData::AuxData, itSample::Int64 )
	checkDataExtendHelper!( auxData.dataSampleLst, auxData.dataLst, itSample )
end

storeAuxDataSampleNoBndCheck( auxData::AuxData, itSample::Int64 ) = throwAuxDataUndefined();

function storeAuxDataStartSample( auxData::AuxData, itStartSample::Int64 )
	checkDataStartSampleExtend!( auxData, itStartSample );
	storeAuxDataStartSampleNoBndCheck( auxData, itStartSample );
end
storeAuxDataStartSampleNoBndCheck( auxData::AuxData, it::Int64 ) = throwAuxDataUndefined();

function checkDataStartSampleExtend!( auxData::AuxData, itStartSample::Int64 )
	checkDataExtendHelper!( auxData.dataStartSampleLst, auxData.dataLst, itStartSample );
end

function storeAuxDataNum( auxData::AuxData, it::Int64 )
	if ~isempty(auxData.dataNumLst)
		checkDataNumExtend!( auxData, it );
		storeAuxDataNumNoBndCheck( auxData, it );
	end
end
storeAuxDataNumNoBndCheck( auxData::AuxData, it::Int64 ) = throwAuxDataUndefined();

function checkDataNumExtend!( auxData::AuxData, it::Int64 )
	checkDataExtendHelper!( auxData.dataNumLst, auxData.dataNumSnapLst, it );
end

function getJldVarSampleLst( auxData::AuxData )
	return auxData.jldVarSampleLst;
end

function getJldVarStartSampleLst( auxData::AuxData )
	return auxData.jldVarStartSampleLst;
end

function getJldVarNumLst( auxData::AuxData )
	return auxData.jldVarNumLst;
end

function getJldVarItSampleLst( auxData::AuxData )
	return auxData.jldVarItSampleLst;
end

function saveAuxDataAll( auxData::AuxData, attrLst::Vector{String}, valLst::Vector; fMod::String = "" )
	renewAuxJldVarLst!( auxData );
	funFNameAux = ( f -> fNameFunc( f( auxData ), attrLst, valLst, jld2Type; fMod = fMod ) );
	fNameLst = funFNameAux.( [getAuxDataSummarySampleName, getAuxDataSummaryStartSampleName, getAuxDataSummaryNumName, getAuxDataSummaryItSampleLstName] );
	
	funSaveAux = (fName, varFunc) -> save( fName, varFunc( auxData )... );
	funSaveAux.( fNameLst, [getJldVarSampleLst, getJldVarStartSampleLst, getJldVarNumLst, getJldVarItSampleLst] );
end

struct NoAuxData <: AuxData

end

NoAuxData( arg... ) = NoAuxData();

struct ZakArrAuxData <: AuxData
	itSampleLst::Vector{Int64};

	zakArr::Array{Bool,3};
	zakArrSampleLst::Vector{Array{Bool,3}};
	zakArrStartSampleLst::Vector{Array{Bool,3}};
	zakMeanLst::Vector{Vector{Float64}};
	
	dataLst::Vector{Array{Bool,3}};
	dataNumSnapLst::Vector{Vector{Float64}}
	dataSampleLst::Vector{Vector{Array{Bool,3}}};
	dataStartSampleLst::Vector{Vector{Array{Bool,3}}};
	dataNumLst::Vector{Vector{Vector{Float64}}};
	
	# dataOutLst::Vector{Array{Bool,3}};
	dataSampleOutLst::Vector{Array{Bool,4}};
	dataStartSampleOutLst::Vector{Array{Bool,4}};
	dataNumOutLst::Vector{Array{Float64,2}};
	
	jldVarSampleLst::Vector{Any};
	jldVarStartSampleLst::Vector{Any};
	jldVarNumLst::Vector{Any};
	jldVarItSampleLst::Vector{Any};
	
	nDim::Int64;
	
	function ZakArrAuxData( params::ParamsLoops, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 )
		zakArr = zeros( Bool, ntuple( i -> params.divNum, 2 )..., params.nDim );
		# zakArrSampleLst = zeros( Bool, size(zakArr)..., itNumSample );
		# zakArrStartSampleLst = zeros( Bool, size(zakArr)..., itNumStartSample );;
		# zakMeanLst = zeros(Float64, params.nDim, itNum);
		itSampleLst = zeros(Int64, itNumSample);
		
		zakArrSampleLst = [ similar(zakArr) for it = 1 : itNumSample ];
		zakArrStartSampleLst = [ similar(zakArr) for it = 1 : itNumStartSample ];
		zakMeanLst = [ zeros(params.nDim) for it = 1 : itNum ];
		
		dataLst = [zakArr];
		dataNumSnapLst = [zeros(params.nDim)];
		dataSampleLst = [zakArrSampleLst];
		dataStartSampleLst = [zakArrStartSampleLst];
		dataNumLst = [zakMeanLst];
		
		dataSampleOutLst = Vector{Array{Bool,4}}(undef,length(dataSampleLst));
		dataStartSampleOutLst = Vector{Array{Bool,4}}(undef,length(dataStartSampleLst));
		dataNumOutLst = Vector{Array{Float64,2}}(undef,length(dataNumLst));
		
		jldVarSampleLst = Vector{Any}(undef,0);
		jldVarStartSampleLst = similar(jldVarSampleLst);
		jldVarNumLst = similar(jldVarSampleLst);
		jldVarItSampleLst = similar(jldVarSampleLst);
		
		nDim = params.nDim;
		
		auxData = new( itSampleLst, zakArr, zakArrSampleLst, zakArrStartSampleLst, zakMeanLst, dataLst, dataNumSnapLst, dataSampleLst, dataStartSampleLst, dataNumLst, dataSampleOutLst, dataStartSampleOutLst, dataNumOutLst, jldVarSampleLst, jldVarStartSampleLst, jldVarNumLst, jldVarItSampleLst, nDim );
		
		# renewAuxJldVarLst!( auxData );
		
		return auxData;
	end
end

ZakArrAuxData( params::ParamsLoops ) = ZakArrAuxData( params, 0, 0, 0 );



struct BLinkAuxData{D} <: AuxData
	itSampleLst::Vector{Int64};
	
	BfieldLst::Vector{Array{Bool,D}};
	linkLst::Vector{Array{Bool,D}};
	linkFerroLst::Array{Array{Bool,D}};
	
	dataLst::Vector{Array{Array{Bool,D}}};
	dataNumSnapLst::Vector{Array{Int64}};
	dataSampleLst::Vector{Vector{Array{Array{Bool,D}}}};
	dataStartSampleLst::Vector{Vector{Array{Array{Bool,D}}}};
	dataNumLst::Vector{Vector{Array{Int64}}};
	
	dataSampleOutLst::Vector{Array{Bool}};
	dataStartSampleOutLst::Vector{Array{Bool}};
	dataNumOutLst::Vector{Array{Int64}};
	
	jldVarSampleLst::Vector{Any};
	jldVarStartSampleLst::Vector{Any};
	jldVarNumLst::Vector{Any};
	jldVarItSampleLst::Vector{Any};
	
	nDim::Int64;
	nDimB::Int64;
	nDimLayer::Int64;
	
	params::ParamsLoops{D};
	
	function BLinkAuxData( params::ParamsLoops, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 )
		itSampleLst = zeros(Int64, itNumSample);
		
		divTup = Tuple( params.divLst );
		BfieldLst, linkLst, linkFerroLst = genBfieldLinkArr( params );
		dataLst = Array{Array{Bool,params.nDim}}[BfieldLst, linkLst, linkFerroLst];
		dataSampleLst = [ [ deepcopy( dataLst[iD] ) for it = 1 : itNumSample ] for iD = 1 : length(dataLst) ];
		dataStartSampleLst = [ [ deepcopy( dataLst[iD] ) for it = 1 : itNumStartSample ] for iD = 1 : length(dataLst) ];
		# BfieldSampleLst, linkSampleLst, linkFerroSampleLst = genBfieldLinkFerroArrSample( params, itNumSample );
		# dataSampleLst = Vector{Array{Array{Bool,params.nDim}}}[BfieldSampleLst, linkSampleLst, linkFerroSampleLst];
		# BfieldStartSampleLst, linkStartSampleLst, linkFerroStartSampleLst = genBfieldLinkFerroArrSample( params, itNumStartSample );
		# dataStartSampleLst = Vector{Array{Array{Bool,params.nDim}}}[ BfieldStartSampleLst, linkStartSampleLst, linkFerroStartSampleLst ];
		
		# numBfieldLst = [zeros(Int64, params.nDimB) for it = 1 : itNum];
		# numLinkLst = [zeros(Int64, params.nDim) for it = 1 : itNum];
		# numLinkFerroLst = [zeros(Int64, params.nDimLayer, params.nDimB) for it = 1 : itNum];
		numBfieldSnap = zeros(Int64, params.nDimB);
		numLinkSnap = zeros(Int64, params.nDim);
		numLinkFerroSnap = zeros(Int64, params.nDimLayer, params.nDimB);
		dataNumSnapLst = Array{Int64}[numBfieldSnap, numLinkSnap, numLinkFerroSnap];
		dataNumLst = [ Array{Int64}[ deepcopy( dataNumSnapLst[iD] ) for it = 1 : itNum ] for iD = 1 : length( dataNumSnapLst ) ];
		# dataNumLst = Vector{Array{Int64}}[numBfieldLst, numLinkLst, numLinkFerroLst];
		
		# numBfieldOutLst = zeros(itNum, params.nDimB);
		# numLinkOutLst = zeros( itNum, params.nDim );
		# numLinkFerroOutLst = zeros( itNum, params.nDimB, params.nDimLayer );
		# dataNumOutLst = Array{Float64}[numBfieldOutLst, numLinkOutLst, numLinkFerroOutLst];
		
		dataSampleOutLst = Vector{Array}(undef,length(dataLst));
		dataStartSampleOutLst = Vector{Array}(undef,length(dataLst));
		dataNumOutLst = Vector{Array}(undef,length(dataLst));
		
		jldVarSampleLst = Vector{Any}(undef,0);
		jldVarStartSampleLst = similar(jldVarSampleLst);
		jldVarNumLst = similar(jldVarSampleLst);
		jldVarItSampleLst = similar(jldVarSampleLst);
		
		new{params.nDim}( itSampleLst, BfieldLst, linkLst, linkFerroLst, dataLst, dataNumSnapLst, dataSampleLst, dataStartSampleLst, dataNumLst, dataSampleOutLst, dataStartSampleOutLst, dataNumOutLst, jldVarSampleLst, jldVarStartSampleLst, jldVarNumLst, jldVarItSampleLst, params.nDim, params.nDimB, params.nDimLayer, params );
	end
end




abstract type LoopsUpdater end

function updateLoops( updater::LoopsUpdater, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	error( "Loops_MC: updater not defined yet" );
end

function updateLoops( updater::LoopsUpdater, flipChecker::FlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	error( "Loops_MC: updater or flipChecker not defined yet" );
end

function updateLoops( updater::LoopsUpdater, flipChecker::FlipChecker, flipProposer::FlipProposer, auxData::AuxData, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	updateLoops( updater, flipChecker, flipProposer, BfieldLst, linkLst, linkFerroLst, params );
	calcAuxData!( auxData, params, BfieldLst, linkLst, linkFerroLst );
end

function getUpdaterFMod( updaterType::(Type{T} where T <: LoopsUpdater) )
	error( "Loops_MC: updater not defined yet" );
end

abstract type AbstractSwitchingUpdater <: LoopsUpdater end

struct SwitchingUpdater{T_tuple} <: AbstractSwitchingUpdater
	updaterTup::Tuple{Vararg{LoopsUpdater}};
	counterRef::Ref{Int64};
	
	function SwitchingUpdater( updaterTup::Tuple{Vararg{LoopsUpdater}} )
		counterVal = 1;
		
		new{typeof(updaterTup)}( updaterTup, Ref(counterVal) );
	end
end

struct SingleUpdater <: LoopsUpdater
	;
end

SingleUpdater( params::ParamsLoops ) = SingleUpdater();

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

abstract type AbstractCubeUpdater <:LoopsUpdater end

struct CubeUpdater <: AbstractCubeUpdater
	
end

CubeUpdater( params::ParamsLoops ) = CubeUpdater();

struct CubeStaggeredCubeUpdater{N,Nplus1} <: AbstractCubeUpdater
	posLstSh0::CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}};
	posLstAdvOrNot::Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}} }};
	posShOrNotLst::Vector{Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}}}};
	posStagCubeLst::Vector{Array{CartesianIndex{N},Nplus1}};
	idStagLst::CartesianIndices{Nplus1,NTuple{Nplus1,Base.OneTo{Int64}}};
	
	iDimLst::UnitRange{Int64};
	iIsShLst::UnitRange{Int64};
	randIDimLst::Array{Int64,Nplus1};
	randIShLst::Array{Int64,Nplus1};
	
	function CubeStaggeredCubeUpdater( params::ParamsLoops )
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





function getFModLoopsMC( fMod::String, updaterType::Type{<:LoopsUpdater}, isInit0::Bool=false; isFModMethod = true, flipChecker::Union{FlipChecker,Nothing} = nothing, flipProposer::Union{FlipProposer,Nothing} = nothing )
	fModOut = fMod;
	if isFModMethod
		fModOut = Utils.strAppendWith_( fModOut, getUpdaterFMod(updaterType) );
	end
	if !isnothing( flipChecker )
		fModOut = Utils.strAppendWith_( fModOut, getFlipCheckerName(flipChecker) );
	end
	if !isnothing( flipProposer )
		fModOut = Utils.strAppendWith_( fModOut, getFlipProposerName(flipProposer) );
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

function genAttrLstLttcFullUpdater( divNum::Int64, nDim::Int64, flipChecker::FlipChecker, initializer::BLinkInitializer; itController::Union{ItController,Nothing} = nothing, rndDigs = rndDigsLpsMC )
	attrLstLttc, valLstLttc = genAttrValLstLttc( divNum, nDim );
	attrLstItCtrl, valLstItCtrl = getAttrValLstItController( itController );
	attrLstFlip, valLstFlip = getFlipCheckerAttrValLst( flipChecker );
	attrLstInit, valLstInit = getAttrValInitializer( initializer );
	
	attrLstLst = [attrLstLttc, attrLstItCtrl, attrLstFlip, attrLstInit];
	valLstLst = [valLstLttc, valLstItCtrl, valLstFlip, valLstInit];
	
	attrLst = append!( attrLstLst... );
	valLst = append!( valLstLst... );
	
	return attrLst, valLst;
end

function loops_MC_methods_cALFCube( divNum = 64, itNum = 10000; updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, fMod = "", cArea = 1, cPerim = 1, cFerro = 0, itNumSample = 100, itStartSample = 50, nDim = 3, probInit = nothing, isFileNameOnly::Bool = false, fMainOutside::String = "" )
	flipChecker = CubeFlipChecker( cArea );
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

function loops_MC_methods_cALF( divNum = 64, itNum = 10000; updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, fMod = "", cArea = 1, cPerim = 1, cFerro = 0, itNumSample = 100, itStartSample = 50, nDim = 3, probInit = nothing, isFileNameOnly::Bool = false, fMainOutside::String = "", flipCheckerType::Type{<:FlipChecker} = NeighborFlipChecker, flipProposerType::Union{Nothing,Type{<:FlipProposer}} = nothing, itControllerType::Union{Type{<:ItController},Nothing} = nothing )
	# if flipCheckerType == NeighborFlipChecker
		# flipChecker = NeighborFlipChecker( cArea, cPerim, cFerro );
	# elseif flipCheckerType == CubeFlipChecker
		# flipChecker = CubeFlipChecker( cArea );
	# end
	flipChecker = flipCheckerType( cArea, cPerim, cFerro );
	if initializerType == BinomialInitializer
		if isnothing(probInit)
			initializer = genMeanFieldInitializer( cArea );
		else
			initializer = genProbInitializer( probInit );
		end
	elseif initializerType == ConstantInitializer
		initializer = ConstantInitializer();
	end
	
	auxDataType = ZakArrAuxData;
	
	if isnothing(flipProposerType)
		flipProposer = nothing;
	else
		flipProposer = flipProposerType();
	end
	
	if isnothing(itControllerType)
		fNameOutside = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, initializer = initializer, flipProposer = flipProposer, auxDataType = auxDataType, fMod = fMod, itNumSample = itNumSample, itStartSample = itStartSample, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside, nDim = nDim );
	else
		itController = itControllerType( itNum, itNumSample, itStartSample );
		fNameOutside = loops_MC_methods_Base( divNum; updaterType = updaterType, flipChecker = flipChecker, initializer = initializer, flipProposer = flipProposer, auxDataType = auxDataType, itController = itController, fMod = fMod, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside, nDim = nDim );
	end
	
	return fNameOutside;
end

# function loops_MC_methods_cALF( divNum = 64, itNum = 10000; updaterType::Type{<:LoopsUpdater}, initializerType::Type{<:BLinkInitializer}, fMod = "", cArea = 1, cPerim = 1, cFerro = 0, itNumSample = 100, itStartSample = 50, nDim = 3, probInit = nothing, isFileNameOnly::Bool = false, fMainOutside::String = "", flipCheckerType::Type{<:FlipChecker} = NeighborFlipChecker, flipProposerType::Union{Nothing,Type{<:FlipProposer}} = nothing )
	# # if flipCheckerType == NeighborFlipChecker
		# # flipChecker = NeighborFlipChecker( cArea, cPerim, cFerro );
	# # elseif flipCheckerType == CubeFlipChecker
		# # flipChecker = CubeFlipChecker( cArea );
	# # end
	# flipChecker = flipCheckerType( cArea, cPerim, cFerro );
	# if initializerType == BinomialInitializer
		# if isnothing(probInit)
			# initializer = genMeanFieldInitializer( cArea );
		# else
			# initializer = genProbInitializer( probInit );
		# end
	# elseif initializerType == ConstantInitializer
		# initializer = ConstantInitializer();
	# end
	
	# auxDataType = ZakArrAuxData;
	
	# if isnothing(flipProposerType)
		# flipProposer = nothing;
	# else
		# flipProposer = flipProposerType();
	# end
	
	# fNameOutside = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, initializer = initializer, flipProposer = flipProposer, auxDataType = auxDataType, fMod = fMod, itNumSample = itNumSample, itStartSample = itStartSample, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside, nDim = nDim );
	
	# return fNameOutside;
# end

function loops_MC_methods_Base_detailed_output( divNum = 64; updaterType::Type{<:LoopsUpdater}, flipChecker::FlipChecker, flipProposer::Union{FlipProposer,Nothing} = nothing, auxDataType::Type{<:AuxData} = NoAuxData, initializer::BLinkInitializer, itController::ItController, fMod = "", nDim = 3, isFileNameOnly::Bool = false, fMainOutside::Union{String, Nothing}= "", attrLstMod::Vector{String} = Vector{String}(undef,0), valLstMod::Vector = [] )
	fModOut = getFModLoopsMC( fMod, updaterType; flipChecker = flipChecker, flipProposer = flipProposer );
	
	fMain = fMainLoopsMC;
	attrLst, valLst = genAttrLstLttcFullUpdater( divNum, nDim, flipChecker, initializer; itController = itController );
	append!( attrLst, attrLstMod );
	append!( valLst, valLstMod );
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
	
	if isFileNameOnly
		fNameOutside = fNameFunc( fMainOutside, attrLst, valLst, jld2Type; fMod = fModOut );
		return fNameOutside, nothing;
	end
	
	params = ParamsLoops( divNum, nDim );
	
	updater = updaterType( params );
	
	auxData = genAuxData( auxDataType, params, flipChecker, itController );
	
	# itSample = 1;
	# while( testItNotDone( itController ) )
		# it = itController.itRef[];
		# print( "it = ", it, "         \r" )
		
		# if testItDoSample( itController )
			# storeAuxDataSample( bLinkData, itSample, it );
			# storeAuxDataSample( auxData, itSample, it );
			# itSample += 1;
		# end
		
		# if testItDoStartSample( itController )
			# storeAuxDataStartSample( bLinkData, it );
			# storeAuxDataStartSample( auxData, it );
		# end
		
		# storeAuxDataNum( bLinkData, it );
		# storeAuxDataNum( auxData, it );
		
		# if isnothing(flipProposer)
			# updateLoops( updater, flipChecker, BfieldLst, linkLst, linkFerroLst, params );
		# else
			# updateLoops( updater, flipChecker, flipProposer, auxData, BfieldLst, linkLst, linkFerroLst, params );
		# end
		
		# advanceItControl( itController );
	# end
	
	# if testItDoSample( itController )
		# it = itController.itRef[];
		# storeAuxDataSample( bLinkData, itSample, it );
		# storeAuxDataSample( auxData, itSample, it );
		# itSample += 1;
	# end
	
	# save( fName, "divNum", divNum, "itNum", itNum );
	
	# saveAuxDataAll( bLinkData, attrLst, valLst; fMod = fModOut );
	# saveAuxDataAll( auxData, attrLst, valLst; fMod = fModOut );
	
	# if isFileNameOnly
		# return fNameOutside;
	# end
	
	# return fName, auxData;
	fName = loops_MC_NoPrefabHelper_Base( params = params, updater = updater, flipChecker = flipChecker, flipProposer = flipProposer, auxData = auxData, initializer = initializer, itController = itController, fMod = fMod, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	
	return fName, auxData;
end

function loops_MC_NoPrefabHelper_Base( ; params::ParamsLoops, updater::LoopsUpdater, flipChecker::FlipChecker, flipProposer::FlipProposer = OneFlipProposer, auxData::AuxData, initializer::BLinkInitializer, itController::ItController, fMod = "", isFileNameOnly::Bool = false, fMainOutside::Union{String, Nothing}= "" )
	fModOut = getFModLoopsMC( fMod, typeof( updater ); flipChecker = flipChecker, flipProposer = flipProposer );
	
	fMain = fMainLoopsMC;
	attrLst, valLst = genAttrLstLttcFullUpdater( params.divNum, params.nDim, flipChecker, initializer; itController = itController );
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
	
	if isFileNameOnly
		fNameOutside = fNameFunc( fMainOutside, attrLst, valLst, jld2Type; fMod = fModOut );
		return fNameOutside;
	end
	
	bLinkData = genAuxData( BLinkAuxData, params, itController );
	BfieldLst, linkLst, linkFerroLst = bLinkData.dataLst;
	
	initializeBL( initializer, BfieldLst, params );
	updateLinkFrom0ByBAllDims( BfieldLst, linkLst, linkFerroLst, params );
	
	calcAuxData!( auxData, params, BfieldLst, linkLst, linkFerroLst );
	resetItControl( itController );
	
	itNum, itNumSample, itNumStartSample = getItNumLst( itController );
	
	itSample = 1;
	while( testItNotDone( itController ) )
		it = itController.itRef[];
		# if Threads.threadid() == 1
			# print( "it = ", it, "         \r" )
		# end
		
		if testItDoSample( itController )
			storeAuxDataSample( bLinkData, itSample, it );
			storeAuxDataSample( auxData, itSample, it );
			itSample += 1;
		end
		
		if testItDoStartSample( itController )
			storeAuxDataStartSample( bLinkData, it );
			storeAuxDataStartSample( auxData, it );
		end
		
		storeAuxDataNum( bLinkData, it );
		storeAuxDataNum( auxData, it );
		
		if isnothing(flipProposer)
			updateLoops( updater, flipChecker, BfieldLst, linkLst, linkFerroLst, params );
		else
			updateLoops( updater, flipChecker, flipProposer, auxData, BfieldLst, linkLst, linkFerroLst, params );
		end
		
		advanceItControl( itController );
	end
	
	if testItDoSample( itController )
		it = itController.itRef[];
		storeAuxDataSample( bLinkData, itSample, it );
		storeAuxDataSample( auxData, itSample, it );
		itSample += 1;
	end
	
	# save( fName, "divNum", divNum, "itNum", itNum );
	
	saveAuxDataAll( bLinkData, attrLst, valLst; fMod = fModOut );
	saveAuxDataAll( auxData, attrLst, valLst; fMod = fModOut );
	
	if isFileNameOnly
		return fNameOutside;
	end
	
	return fName;
end

# function loops_MC_methods_Base( divNum = 64; updaterType::Type{<:LoopsUpdater}, flipChecker::FlipChecker, flipProposer::Union{FlipProposer,Nothing} = nothing, auxDataType::Type{<:AuxData} = NoAuxData, initializer::BLinkInitializer, itController::ItController, fMod = "", nDim = 3, isFileNameOnly::Bool = false, fMainOutside::Union{String, Nothing}= "" )
	# fName, auxData = loops_MC_methods_Base_detailed_output( divNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, auxDataType = auxDataType, initializer = initializer, itController = itController, fMod = fMod, nDim = nDim, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	
	# return fName;
# end

function loops_MC_methods_Base( args...; kwargs... )
	fName, auxData = loops_MC_methods_Base_detailed_output( args...; kwargs... );
	
	return fName;
end

function loops_MC_methods_Base_detailed_output( divNum, itNum; updaterType::Type{<:LoopsUpdater}, flipChecker::FlipChecker, flipProposer::Union{FlipProposer,Nothing} = nothing, auxDataType::Type{<:AuxData} = NoAuxData, initializer::BLinkInitializer, fMod = "", itNumSample = 100, itStartSample = 50, nDim = 3, isFileNameOnly::Bool = false, fMainOutside::Union{String, Nothing}= "" )
	fModOut = getFModLoopsMC( fMod, updaterType; flipChecker = flipChecker, flipProposer = flipProposer );
	
	fMain = fMainLoopsMC;
	attrLst, valLst = genAttrLstLttcFlipInit( divNum, itNum, nDim, flipChecker, initializer );
	fName = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
	
	if isFileNameOnly
		fNameOutside = fNameFunc( fMainOutside, attrLst, valLst, jld2Type; fMod = fModOut );
		return fNameOutside;
	end
	
	params = ParamsLoops( divNum, nDim );
	
	itStep = max( Int64( floor(itNum / itNumSample) ), 1 );
	lnSample = Int64( floor( itNum / itStep ) );
	itStartSample = min( itStartSample, itNum );
	
	BfieldLst, linkLst, linkFerroLst = genBfieldLinkArr( params );
	
	numBfieldLst, numLinkLst = genBfieldLinkNumLst( params, itNum );
	BfieldSampleLst, linkSampleLst = genBfieldLinkArrSample( params, lnSample );
	BfieldStartSampleLst, linkStartSampleLst = genBfieldLinkArrSample( params, itStartSample );

	initializeBL( initializer, BfieldLst, params );
	updateLinkFrom0ByBAllDims( BfieldLst, linkLst, linkFerroLst, params );
	updater = updaterType( params );
	
	auxData = auxDataType( params, flipChecker, itNum, lnSample, itStartSample );
	calcAuxData!( auxData, params, BfieldLst, linkLst, linkFerroLst );
	
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
			storeAuxDataSample( auxData, itSample, it );
			itSample += 1;
		end
		
		if it <= itStartSample
			Threads.@threads for dim = 1 : params.nDim
				linkStartSampleLst[it][dim] .= linkLst[dim];
			end
			Threads.@threads for dim = 1 : params.nDimB
				BfieldStartSampleLst[it][dim] .= BfieldLst[dim];
			end
			storeAuxDataStartSample( auxData, it );
		end
		
		if isnothing(flipProposer)
			updateLoops( updater, flipChecker, BfieldLst, linkLst, linkFerroLst, params );
		else
			updateLoops( updater, flipChecker, flipProposer, auxData, BfieldLst, linkLst, linkFerroLst, params );
		end
		
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
	
	saveAuxDataAll( auxData, attrLst, valLst; fMod = fModOut );
	
	if isFileNameOnly
		return fNameOutside;
	end
	
	return fName, auxData;
end

function loops_MC_methods_Base( divNum, itNum; updaterType::Type{<:LoopsUpdater}, flipChecker::FlipChecker, flipProposer::Union{FlipProposer,Nothing} = nothing, auxDataType::Type{<:AuxData} = NoAuxData, initializer::BLinkInitializer, fMod = "", itNumSample = 100, itStartSample = 50, nDim = 3, isFileNameOnly::Bool = false, fMainOutside::Union{String, Nothing}= "" )
	fName, auxData = loops_MC_methods_Base_detailed_output( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, auxDataType = auxDataType, initializer = initializer, fMod =fMod, itNumSample = itNumSample, itStartSample = itStartSampe, nDim = nDim, isFilenameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	
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
	ELst = - cArea .* Float64[2*nDimB : -2 : -2*nDimB;] ;
	dELst = -2 .* ELst;
	pFlipLst = 1 ./ (1 .+ exp.(dELst) );
	
	return pFlipLst;
end

function boolToOnePN( varBool::Bool )
	# return -(-1).^varBool;
	return varBool ? -1 : 1;
end

function boolToOnePN( varBoolNum::Real; cnt::Int64 = 1 )
	# return -(-1).^varBool;
	return cnt - 2*varBoolNum;
end

function onePNToBool( varOnePN::Real; cnt::Int64 = 1 )
	# return -(-1).^varBool;
	return Int64( ( cnt - varOnePN ) / 2 );
end

function boolToFlipChange( sBool::Bool )
	return boolToOnePN( sBool )
end

function boolToFlipOnePNChange( sBool::Bool )
	return -2 * boolToOnePN( sBool );
end

function ratioToBinId( ratio::Number, divNum::Int64 )
	return Int64( min( floor( ratio * divNum ) + 1, divNum  ) );
end

include("loops_MC_initializer.jl")

include("loops_MC_flipProposer.jl")

include("loops_MC_flipChecker.jl")

include("loops_MC_itController.jl")

include("loops_MC_auxArr.jl")

include("loops_MC_updater.jl")

include("loops_MC_legacy.jl")

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
