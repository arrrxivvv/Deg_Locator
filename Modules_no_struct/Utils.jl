module Utils

using LinearAlgebra
using EllipsisNotation
using Logging; using LoggingExtras;

const fileType = ".pdf";
const jldType = ".jld";
const npyType = ".npy";

export fileType, jldType, npyType

export printlnLogFile, printLogFile, SortABbyA, EigenSort!, dotEachCol, wrap, arrLstToArr, arrToArrLst, cartIndLstToArr, gridLst, wrapCartInd!

function printlnLogFile( file, msg... )
	open(f -> (println(f, msg...); println(stdout, msg...)), file, "a");
end

function printLogFile( file, msg... )
	open(f -> (print(f, msg...); print(stdout, msg...)), file, "a");
end

function SortABbyA( key, object )
	indSort = sortperm( key );
	key .= key[indSort];
	object .= object[indSort,..];
	return key, object;
end

# function EigenSort!( valLSt, vecLst )
	# indSort = sortperm( valLSt );
	# valLSt .= valLSt[indSort];
	# vecLst .= vecLst[..,indSort];
	# # return valLSt, vecLst
# end

function wrapCartInd!( ind, arr )
	if !checkbounds( Bool, arr, ind )
		@debug "not in bound"
		bndLst = size(arr);
		indArr = [ ind[ii] for ii = 1 : length(ind) ];
		ind = CartesianIndex( Tuple( mod.( indArr .- 1, bndLst ) .+ 1 ) );
		@debug "ind: " ind
	end
	return ind;
end

function EigenSort!( valLSt, vecLst )
	indSort = sortperm( valLSt );
	valLSt = valLSt[indSort];
	vecLst = vecLst[..,indSort];
	# return valLSt, vecLst
end

function dotEachCol( vec1, vec2 )
	# @debug( "dotEachCol", size(vec1), size(vec2) );
	return dot.( eachcol(vec1), eachcol(vec2) );
end

function wrap( point, bnd )
	isWrappedArr = point .> (bnd/2);
	point = -bnd .* isWrappedArr + point;
end

function arrLstToArr( arrLst )
	return [ arrLst[idOut][idIn] for idOut in CartesianIndices(arrLst), idIn in CartesianIndices(arrLst[1]) ];
end

function arrLstToArrCat( arrLst )
	numDim = ndims(arrLst[1]);
	arr = cat( arrLst..., dims = numDim+1 );
	arr = permutedims( arr, circshift( [1:numDim], 1 ) );
end

function arrToArrLst( arr )
	return [ arr[it,:] for it = 1 : size(arr)[1] ];
end

function cartIndLstToArr( indLst, dim )
	return [ indLst[idOut][idIn] for idOut in CartesianIndices(indLst), idIn = 1:dim ];
end

function gridLst( space_dim, space_divide )
	xLst = Vector{ Array{ Complex{Float64} } }(undef,space_dim);
	for dim = 1 : space_dim
		ln_lst = ones(Int64, space_dim);
		ln_lst[dim] = space_divide;	
		xLstTmp = collect( 0 : space_divide-1 );
		xLst[dim] = reshape( xLstTmp, ln_lst... );
	end
	return xLst;
end

end
