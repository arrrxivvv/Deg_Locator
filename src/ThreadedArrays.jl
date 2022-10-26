module ThreadedArrays

export ThrArray, ThrStruct, thrStructFill, thrStructCpyTheRest, threaded_zeros, threaded_ones, threaded_fill, thrArr_empty, getThrInst, broadcastAssig;

struct ThrArray{T,N} # <: AbstractArray{T,N}
	data::Vector{Array{T,N}};
end

struct ThrStruct{Tstruct}
	data::Vector{Tstruct};
end

function thrStructFill( T::DataType, args... )
	ThrStruct{T}( [ T(args...) for ii = 1 : Threads.nthreads() ] );
end

function thrStructCpyTheRest( obj )
	T = typeof(obj);
	thrObj = ThrStruct{T}( Vector{T}(undef,Threads.nthreads()) );
	thrObj.data[1] = obj;
	for ii = 2 : Threads.nthreads()
		thrObj.data[ii] = deepcopy(obj);
	end
	
	return thrObj;
end

function getThrInst( thrStruct::ThrStruct{T} ) where{T}
	if length(thrStruct.data) == 0
		return nothing;
	else
		return thrStruct.data[Threads.threadid()];
	end
end

ThrArray{T}( data::Vector{Array{T,N}} ) where{T,N} = ThrArray{T,ndims(data[1])}( data );

function threaded_zeros( dims... )
	threaded_zeros( Float64, dims... );
end

function threaded_zeros( T::DataType, dims... )
	ThrArray{T}( [ zeros( T, dims... ) for ii = 1 : Threads.nthreads() ] );
end

function threaded_ones( T::DataType=Float64, dims... )
	ThrArray{T}( [ ones(T, dims...) for ii = 1 : Threads.nthreads() ] );
end

function threaded_fill( T, elements, dims )
	ThrArray{T}( [ [ T(elements...) for id in CartesianIndices(dims) ] for iTh = 1 : Threads.nthreads() ] );
end

function thrArr_empty(T::DataType=Float64)
	ThrArray{T,1}(Vector{Vector{T}}(undef,0));
end

function getThrInst( arrTh::ThrArray{T} ) where{T}
	if length(arrTh.data) == 0
		return nothing;
	else
		return arrTh.data[Threads.threadid()];
	end
end

function broadcastAssign!( arrTh::ThrArray{T}, src ) where{T}
	getThrInst( arrTh ) .= src;
end

# function Base.getindex( arrTh::ThrArray{T}, ind...  ) where{T}
	# arrTh.data[Threads.threadid()][ind];
# end
Base.getindex( arrTh::ThrArray{T}, ind... ) where{T} = Base.getindex( getThrInst(arrTh), ind... );
Base.setindex!( arrTh::ThrArray{T}, X, ind... ) where{T} = Base.setindex!( getThrInst(arrTh), X, ind... );
Base.view( arrTh::ThrArray{T}, ind... ) where{T} = Base.view( getThrInst(arrTh), ind... );
Base.similar( arrTh::ThrArray{T} ) where{T} = ThrArray{T}( [ Base.similar( arrTh.data[ii] ) for ii = 1 : Threads.nthreads() ] );
Base.ndims( arrTh::ThrArray{T} ) where{T} = Base.ndims( arrTh.data[Threads.threadid()] );
Base.size( arrTh::ThrArray{T,N} ) where{T,N} = Base.size( arrTh.data[Threads.threadid()] );
# Base.BroadcastStyle(::Type{<:ThrArray{T}}) wherer {T} = 
# Base.broadcast!( f, arrTh::ThrArray{T}, As... ) = Base.broadcast!( f, arrTh.data[Threads.threadid()], As... );

end
