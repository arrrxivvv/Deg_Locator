struct DegParams
	divLst::Vector{Int64};
	nDim::Int64;
	N::Int64;
	nonPeriodic::Bool;
	posLst::AbstractArray{CartesianIndex{N}, N} where N;
	
	minLst::Vector{Float64};
	maxLst::Vector{Float64};
	stepLst::Vector{Float64};
	gridLst::Vector{ Vector{Float64} };
	mesh::Array{ Vector{Float64} };
end

function degParamsNonPeriodic( N, divLst, minLst, maxLst, nDim )
	posLst = CartesianIndices(Tuple(divLst));
	if !isa(minLst, Array)
		minLst = fill( minLst, nDim );
	end
	if !isa(maxLst, Array)
		maxLst = fill( maxLst, nDim );
	end
	nonPeriodic = false;
	stepLst = ( maxLst .- minLst ) ./ divLst;
	gridLst = [ collect( range( minLst[iDim], maxLst[iDim] - stepLst[iDim], length = divLst[iDim] ) ) for iDim = 1:nDim ];
	mesh = [ [ gridLst[j][ind[j]] for j = 1:nDim ] for ind in posLst ];
	
	return DegParams( divLst, nDim, N, nonPeriodic, posLst, minLst, maxLst, stepLst, gridLst, mesh );
end

function makeArrOverGrid( type::DataType, params::DegParams )
	return Array{type,params.nDim}(fillVal,params.divLst...);
end

function linIdFromIdVec( idVec::Vector{Int64}, params::DegParams )
	id = idVec[end]-1;
	for iDim = (params.nDim-1):-1:1
		id = id * params.divLst[iDim] + idVec[iDim]-1;
	end
	id = id+1;
end

function wrapIdVec!( idVec::Vector{Int64}, params::DegParams )
	idVec .-= 1;
	idVec .= mod.( idVec, params.divLst );
	idVec .+= 1;
end
