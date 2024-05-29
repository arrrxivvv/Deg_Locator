#Loops_MC module

function getFlipProposerName( flipProposerType::Type{<:SwitchingFlipProposer} )
	return join( getFlipProposerName.((flipProposerType.parameters[1]).parameters), "_" )
end



function getFlipProposerName( flipProposerType::Type{OneFlipProposer} )
	return "oneFlip";
end

function flipDoIt( flipProposer::OneFlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = dim );
end




function getFlipProposerName( flipProposerType::Type{CubeFlipProposer} )
	return "cubeFlip"
end

function flipDoIt( flipProposer::CubeFlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	flipBfieldCubeAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos );
end
