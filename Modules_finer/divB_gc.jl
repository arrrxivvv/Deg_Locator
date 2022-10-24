using LinearAlgebra
# using ArgParse
using EllipsisNotation
using Distributions
using Random
using JLD
using Logging; using LoggingExtras;
# using JLD2, FileIO
include("utils.jl")
include("Umat.jl")
include("Chern.jl")

function locator_div(N,param_divide_num;Hmat_seed = 1234, polarity = +1)
	param_dim = 3;
	dimLst = collect(1:param_dim);
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	# param_divide_num = 50;
	param_ln_num = param_divide_num+1;
	param_divide = param_divide_num * ones(param_dim);
	param_ln = Int64.( param_divide .+ 1 );
	param_step = ( param_max - param_min ) ./ param_divide;
	param_grids = hcat( collect( range( param_min, length=param_ln_num, stop=param_max ) ) );
	
	# N=5;
	
	# Hmat_seed = 1234;
	H_GUElst = fill( Array{Complex{Float64}}(undef,(N,N)), (2,param_dim) );
	Random.seed!(Hmat_seed);
	for it = 1:length(H_GUElst)
		H_GUElst[it] = H_GUE(N);
	end
	HmatFun = (x,y,z)->Hmat_3comb([x,y,z],H_GUElst);
	
	arrForInd = zeros( Bool, param_ln... );
	
	HmatLst = [ zeros( Complex{Float64}, N,N ) for it in CartesianIndices(arrForInd) ];
	
	param_arr = zeros(param_dim);
	for ind in CartesianIndices(HmatLst)
		for j = 1 : param_dim
			param_arr[j] = param_grids[ind[j]][j];
		end
		HmatLst[ind] = HmatFun( param_arr... );
	end
	
	@info( "eigen" );
	
	eigenLst = eigen.( HmatLst );
	@debug(size(eigenLst));
	@debug(size(HmatLst));
	Elst = (p-> Real.( p.values )).( eigenLst );
	@debug(typeof(Elst), size(Elst));
	vecLst = (p->p.vectors).( eigenLst );
	EigenSort!.( Elst, vecLst );
	
	HmatLst = nothing;
	GC.gc();
	
	@info( "linkLst" ); 
	
	linkLst = Vector{ Array{Vector{Complex{Float64}},param_dim} }(undef,param_dim);
	for dim = 1:param_dim
		sh_lst = zeros(param_dim);
		sh_lst[dim] = -1;
		vecLst_sh = circshift(vecLst,Tuple(sh_lst));
		@debug(size(vecLst_sh), size(vecLst));
		linkLst[dim] = dotEachCol.(vecLst_sh,vecLst);
		linkLst[dim] = ( x -> x ./ abs.(x) ).(linkLst[dim]);
		vecLst_sh = nothing;
	end
	
	eigenLst = nothing;
	Elst = nothing;
	vecLst = nothing;
	GC.gc();
	
	@info( "Bfield" ); 
	Bfield_ln = Int64( round( param_dim * (param_dim-1) / 2 ) );
	# typeOfBvec = Vector{Complex(Float64)}(undef,Bfield_ln);
	BfieldLst = Vector{ Array{Vector{Complex{Float64}},param_dim} }(undef,Bfield_ln);
	id_B = 0;
	for dim1 = 1 : param_dim
		for dim2 = dim1+1 : param_dim
			dimLst = [dim1,dim2];
			dimOtherLst = [dim2,dim1];
			id_B += 1;
			linkRatio = Vector{ Array{Vector{Complex{Float64}},param_dim} }(undef,2);
			for it = 1 : 2
				dim = dimLst[it];
				dimSh = dimOtherLst[it];
				sh_lst = zeros(param_dim);
				# @debug( "zeros, ", typeof(sh_lst) );
				sh_lst[dimSh] = -1;
				# @debug( "circshift, ", typeof(sh_lst) );
				sh_lst = Int64.(sh_lst);
				linkLstDim_sh = circshift(linkLst[dim],sh_lst);
				linkRatio[it] = broadcast.( /, linkLst[dim], linkLstDim_sh );
				linkLstDim_sh = nothing;
				GC.gc();
			end
			parity = (-1)^(dim1+dim2-1);
			BfieldLst[id_B] = parity /1im * broadcast.( log, broadcast.( /, linkRatio[1], linkRatio[2] ) );
			linkRatio = nothing;
			GC.gc()
		end
	end
	
	linkLst = nothing;
	linkRatio = nothing;
	GC.gc();
	
	@info( "divB" );
	
	# divBlst = Array{ Vector{Complex{Float64}}, 2 }(undef, param_ln);
	divBlst = [ zeros(Float64,N) for it in CartesianIndices(arrForInd) ];
	dimLstLeviCivita = [3,2,1];
	for Bdim = 1 : Bfield_ln
		vecDim = dimLstLeviCivita[Bdim];
		sh_lst = zeros(Bfield_ln);
		sh_lst[vecDim] = -1;
		@debug(typeof(divBlst));
		@debug(typeof(divBlst[1]));
		divBlst = divBlst .+ circshift( BfieldLst[Bdim],sh_lst ) .- BfieldLst[Bdim];
	end
	
	BfieldLst = nothing;
	GC.gc();
	
	threshold = 0.00001;
	divBlstReal = (p-> real.(p)).(divBlst);
	divBposInd = broadcast.( >, divBlstReal, threshold );
	divBnegInd = broadcast.( <, divBlstReal, -threshold );
	# numDegLst = count.( broadcast.( >, divBlstReal, threshold ) );
	posCount = sum( divBposInd );
	negCount = sum( divBnegInd );
	GC.gc();
	for n = 1:N
		println( "level ", n, ": ", "N_+ = ", posCount[n], ", N_- = ", negCount[n] );
	end
	println( "total N: N_+ = ", sum(posCount), ", N_- = ", sum(negCount) );
	# @debug(divBlstReal);
	@debug(size(divBlstReal[1]));
	# println( "nonzero div B count: ", sum(numDegLst) );
end
