#Loops_MC module

function getSwitchingFlipTypeLst( flipCheckerType::Type{<:SwitchingFlipChecker} )
	return (flipCheckerType.parameters[1]).parameters;
end

function getFlipCheckerName( flipCheckerType::Type{<:SwitchingFlipChecker} )
	return strLstJoinWith_( getFlipCheckerName.( getSwitchingFlipTypeLst( flipCheckerType ) ) );
end

function getFlipCheckerAttrLst( flipChecker::SwitchingFlipChecker )
	return append!( getFlipCheckerAttrLst.( flipChecker.flipCheckerTup )... );
end

function getFlipCheckerValLst( flipChecker::SwitchingFlipChecker; rndDigs = rndDigsLpsMC )
	return append!( getFlipCheckerValLst.( flipChecker.flipCheckerTup; rndDigs = rndDigs )... );
end


function flipCheck( flipChecker::SwitchingFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	error( "SwitchingFlipChecker not to be flipped directly" );
end

function getSwitchLstLen( flipChecker::SwitchingFlipChecker )
	return length(flipChecker.flipCheckerTup);
end



function getFlipCheckerAttrLst( flipChecker::AbstractNeighborFlipChecker )
	return ["cArea", "cPerim", "cFerro"];
end

function flipDoIt( flipChecker::AbstractNeighborFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = dim );
end

function getFlipCheckerName( flipCheckerType::Type{NeighborFlipChecker} )
	return "cAcLcFFlip";
end

# function getFlipCheckerAttrLst( flipChecker::NeighborFlipChecker )
	# return ["cArea", "cPerim", "cFerro"];
# end

function getFlipCheckerAttrLst( flipChecker::NeighborFlipChecker )
	attrLst = ["cArea", "cPerim", "cFerro"];
	if flipChecker.cParamLst[3] == 0
		attrLst = attrLst[1:2];
	end
	
	return attrLst;
end

function getFlipCheckerValLst( flipChecker::AbstractNeighborFlipChecker; rndDigs = rndDigsLpsMC )
	valLst = roundKeepInt.( deepcopy( flipChecker.cParamLst ); digits = rndDigs );
	if flipChecker.cParamLst[3] == 0
		valLst = valLst[1:2];
	end
	
	return valLst;
end

# function getFlipCheckerAttrValLst( flipChecker::AbstractNeighborFlipChecker; rndDigs = rndDigsLpsMC )
	# attrLst = getFlipCheckerAttrLst( flipChecker );
	# valLst = roundKeepInt.( deepcopy( flipChecker.cParamLst ); digits = rndDigs );
	# if flipChecker.cParamLst[3] == 0
		# attrLst = attrLst[1:2];
		# valLst = valLst[1:2];
	# end
	
	# return attrLst, valLst;
# end

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




function flipDoIt( flipChecker::CubeFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	flipBfieldCubeAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos );
end

function getFlipCheckerName( flipCheckerType::Type{CubeFlipChecker} )
	return "cAcLcFCubeFlip";
end

function getFlipCheckerAttrLst( flipChecker::CubeFlipChecker )
	return ["cArea"];
end

function getFlipCheckerValLst( flipChecker::CubeFlipChecker; rndDigs = rndDigsLpsMC )
	valLst = roundKeepInt.( [flipChecker.cArea]; digits = rndDigs );
	
	return valLst;
end

# function getFlipCheckerAttrValLst( flipChecker::CubeFlipChecker; rndDigs = rndDigsLpsMC )
	# attrLst = getFlipCheckerAttrLst( flipChecker );
	# valLst = roundKeepInt.( [flipChecker.cArea]; digits = rndDigs );
	
	# return attrLst, valLst;
# end

function flipCheck( flipChecker::CubeFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	iArea = 1;
	for dimB = 1 : params.nDimB
		iArea += BfieldLst[dimB][pos];
		iArea += BfieldLst[dimB][params.posLstShLst[dimB,1][pos]];
	end
	
	pSwitchRand = rand();
	return pSwitchRand < flipChecker.pFlipLst[iArea];
end




function flipBfieldCubeAtPos( params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}; pos::CartesianIndex{D} ) where {D}
	for dimB = 1 : params.nDimB
		BfieldLst[dimB][pos] = !BfieldLst[dimB][pos];
		BfieldLst[dimB][params.posLstShLst[dimB,1][pos]] = !BfieldLst[dimB][params.posLstShLst[dimB,1][pos]];
	end
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
