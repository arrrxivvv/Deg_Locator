#Loops_MC module

function flipDoIt( flipProposer::OneFlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = dim );
end
