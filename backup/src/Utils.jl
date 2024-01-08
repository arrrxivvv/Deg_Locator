module Utils

# using Infiltrator
using LinearAlgebra
using EllipsisNotation
using Logging; using LoggingExtras;
using ShiftedArrays

const fileType = ".pdf";
const jldType = ".jld";
const jld2Type = ".jld2";
const npyType = ".npy";

@enum MatchMethod mX = 1 mExtended mSameN mPropagate;

export error_only_logger

export fileType, jldType, jld2Type, npyType, MatchMethod, mX, mExtended, mSameN, mPropagate

export printlnLogFile, printLogFile, SortABbyA, EigenSort!, dotEachCol, dotEachCol!, wrap, arrLstToArr, arrLstToArrDepth, arrToArrLst, cartIndLstToArr, gridLst, wrapCartInd!, selectDimLst, wrapIntInd, getLinInd, wrapIndArr!, wrapCoorArr!, wrapDiffArr!, funArrDims, shIdVec!

export thrId

export myreverse!, myreverse

export colabIOreset

export structAssign!

error_only_logger = MinLevelLogger( current_logger(), Logging.Error );

include("Utils_DimDrops.jl")
export maxDropDims, minDropDims, sumDropDims

include("Utils_time_print.jl")
export timeMemStr, timeInfo

include("Utils_sort_perms.jl")
export permuteArr!, permuteCol2d!, permute1d!

include("Utils_arrManip.jl")
export arrSlcLst, arrShAllLst, assignArrOfArrs!

function colabIOreset()
	if isdefined(Main, :IJulia)
		Main.IJulia.stdio_bytes[] = 0
	end	
	flush(stdout);
	flush(stderr);
end

function selectDimLst( arr, dimLst, indLst )
	sliceLst = repeat(Any[:],ndims(arr));
	sliceLst[dimLst] .= indLst;
	# @infiltrate	
	return @view( arr[sliceLst...] );
end

function printlnLogFile( file, msg... )
	open(f -> (println(f, msg...); println(stdout, msg...)), file, "a");
end

function printLogFile( file, msg... )
	open(f -> (print(f, msg...); print(stdout, msg...)), file, "a");
end

function thrId()
	return Threads.threadid();
end


function getLinInd( indArr, szLst, nDim::Int64 )
	ind = indArr[end]-1;
	for id = nDim-1 : -1 : 1
		ind = ind * szLst[id] + indArr[id]-1;
	end
	ind = ind+1;
	return ind;
end

function wrapIntInd( ind, bnd )
	return mod( ind - 1, bnd ) + 1;
end

function wrapIndArr!( indArr, bndArr )
	indArr .= mod.( indArr .- 1, bndArr ) .+ 1;
end

function wrapCoorArr!( coorArr, bndArr )
	coorArr .= mod.( coorArr, bndArr );
end

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

function wrapDiffArr!( indArr, bndArr )
	indArr .= ( indArr .> bndArr./2 ? indArr : indArr .- bndArr );
end

function dotEachCol( vec1, vec2 )
	# @debug( "dotEachCol", size(vec1), size(vec2) );
	return dot.( eachcol(vec1), eachcol(vec2) );
end

function dotEachCol!( dotted, vec1, vec2 )
	# @debug( "dotEachCol", size(vec1), size(vec2) );
	N = size(vec1, 1);
	for n = 1 : N
		dotted[n] = dot( view(vec1, :,n), view(vec2,:,n) );
	end
end

function wrap( point, bnd )
	isWrappedArr = point .> (bnd/2);
	point = -bnd .* isWrappedArr + point;
end

function arrLstToArr( arrLst )
	return length(arrLst) != 0 ? [ arrLst[idOut][idIn] for idOut in CartesianIndices(arrLst), idIn in CartesianIndices(arrLst[1]) ] : [];
end

function arrLstToArrDepth( arrLst, depth )
	arrSub = arrLst;
	szLst = Vector{Any}(undef,depth);
	indLst = Vector{Any}(undef,depth);
	for dNow = 1 : depth
		szLst[dNow] = size(arrSub);
		indLst[dNow] = CartesianIndices(arrSub);
		# if dNow < depth
		arrSub = arrSub[1];
		# end
	end
	arrArr = zeros( typeof(arrSub), (szLst...)... );
	# for itFull in CartesianIndices(arrArr)
	for itFull in Iterators.product( indLst... )
		arrSub = arrLst;
		for dNow = 1 : depth
			arrSub = arrSub[itFull[dNow]];
		end
		# arrArr[itFull...] = arrSub;
		arrArr[itFull...] = 0;
		# @infiltrate
	end
	return arrArr;
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

myreverse(A::AbstractArray; dims=:) = myreverse!(copy(A); dims)
myreverse!(A::AbstractArray; dims=:) = _myreverse!(A, dims)
_myreverse!(A::AbstractArray{<:Any,N}, ::Colon) where {N} = _myreverse!(A, ntuple(identity, Val{N}()))
_myreverse!(A, dim::Integer) = _myreverse!(A, (Int(dim),))
_myreverse!(A, dims::NTuple{M,Integer}) where {M} = _myreverse!(A, Int.(dims))

function _myreverse!(A::AbstractArray{<:Any,N}, dims::NTuple{M,Int}) where {N,M}
    if N < M || !allunique(dims) || !all(d -> 1 ≤ d ≤ N, dims)
        throw(ArgumentError("invalid dimensions $dims to reverse"))
    end
    M == 0 && return A # nothing to reverse
    
    idx = Vector{Int}(undef, N) # temporary array for index calculations
    
    # compute sizes for half of reversed array
    ntuple(k -> @inbounds(idx[k] = size(A,k)), Val{N}())
    @inbounds idx[dims[1]] = (idx[dims[1]] + 1) >> 1
    
    last1 = ntuple(k -> lastindex(A,k)+1, Val{N}())
    for i in CartesianIndices(ntuple(k -> firstindex(A,k):firstindex(A,k)-1+@inbounds(idx[k]), Val{N}()))
        ntuple(k -> @inbounds(idx[k] = i[k]), Val{N}())
        ntuple(Val{M}()) do k
            @inbounds j = dims[k]
            @inbounds idx[j] = last1[j] - idx[j]
        end
        iᵣ = CartesianIndex(ntuple(k -> @inbounds(idx[k]), Val{N}()))
        @inbounds A[iᵣ], A[i] = A[i], A[iᵣ]
    end
    return A
end


function funArrDims( fun, arr, dims )
	return dropdims( fun( arr; dims = dims ); dims = dims );
end

function structAssign!( dest, src )
	fldLst = fieldnames( typeof(src) );
	for fld in fldLst
		if isa( getfield( src, fld ), Array )
			getfield( dest, fld ) .= getfield( src, fld );
		else
			setfield!( dest, fld, getfield(src, fld) );
		end
	end
end

function shIdVec!( idVec, nDim, iSh )
	idVec .= 0;
	iSh -= 1;
	for iDim = 1 : nDim
		idVec[iDim] += (iSh & 1);
		iSh = iSh >> 1;
	end
end

end
