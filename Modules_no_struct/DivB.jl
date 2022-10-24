module DivB

using LinearAlgebra
using EllipsisNotation
using Distributions
using Random
using JLD
using Logging; using LoggingExtras;

using Utils
using UMat
# include("Chern.jl")

export divB, linkLstFun, BfieldFromLink, Bfield_dA

function linkLstFun( param_dim, vecLst )
	# linkLst = Vector{ Array{Vector{Complex{Float64}},param_dim} }(undef,param_dim);	
	xyzLst = CartesianIndices(vecLst);
	N = size( vecLst[1], 1 );
	@debug "N: " N
	param_divide = size( vecLst );
	linkLst = [ [ zeros( Complex{Float64}, N ) for pos in xyzLst ] for dim = 1:param_dim ];
	for dim = 1:param_dim
		sh_lst = zeros(Int64,param_dim);
		# sh_lst[dim] = -1;
		sh_lst[dim] = +1;
		shCartInd = CartesianIndex( Tuple(sh_lst) );
		# vecLst_sh = circshift(vecLst,Tuple(sh_lst));
		Threads.@threads for pos in CartesianIndices(vecLst)
			nextPos = pos + shCartInd;
			nextPos = wrapCartInd!( nextPos, vecLst );
			linkLst[dim][pos] = dotEachCol( vecLst[nextPos], vecLst[pos] );
		end
		# linkLst[dim] = dotEachCol.(vecLst_sh,vecLst);
		linkLst[dim] = ( x -> x ./ abs.(x) ).(linkLst[dim]);
	end
	return linkLst;
end

function BfieldFromLink( param_dim, linkLst )
	Bfield_ln = Int64( round( param_dim * (param_dim-1) / 2 ) );
	# BfieldLst = Vector{ Array{Vector{Complex{Float64}},param_dim} }(undef,Bfield_ln);
	posLst = CartesianIndices( linkLst[1] );
	N = size( linkLst[1][1], 1 );
	BfieldLst = [ [ zeros( Complex{Float64}, N ) for pos in posLst ] for l = 1 : Bfield_ln ];
	id_B = 0;
	linkRatio = [ [ zeros( Complex{Float64}, N ) for pos in posLst ] for id = 1:2 ];
	linkLstSh = [ zeros(Complex{Float64}, N) for pos in posLst ];
	for dim1 = 1 : param_dim
		for dim2 = dim1+1 : param_dim
			dimLst = [dim1,dim2];
			dimOtherLst = [dim2,dim1];
			id_B += 1;
			# linkRatio = Vector{ Array{Vector{Complex{Float64}},param_dim} }(undef,2);
			for it = 1 : 2
				dim = dimLst[it];
				dimSh = dimOtherLst[it];
				sh_lst = zeros(Int64, param_dim);
				sh_lst[dimSh] = -1;
				# sh_lst[dimSh] = 1;
				# shCartInd = CartesianIndex( Tuple(sh_lst) );
				# Threads.@threads for pos in posLst
					# nextPos = pos + shCartInd;
					# nextPos = wrapCartInd!( nextPos, linkLst[1] );
					# linkRatio[it] = broadcast.( /, linkLst[dim], circshift(linkLst[dim],sh_lst) );
					linkRatio[it] = broadcast.( /, linkLst[dim], circshift!(linkLstSh,linkLst[dim],sh_lst) );
					# linkRatio[it][pos] = linkLst[dim][pos] ./ linkLst[dim][nextPos];
				# end
			end
			parity = (-1)^(dim1+dim2-1);
			BfieldLst[id_B] = parity /1im * broadcast.( log, broadcast.( /, linkRatio[1], linkRatio[2] ) );
		end
	end
	return Bfield_ln, BfieldLst;
end

function Bfield_dA( param_dim, param_divide, linkLst )
	Bln, Blst = BfieldFromLink( param_dim, linkLst );
	id_B = 0;
	for dim1 = 1 : param_dim
		for dim2 = dim1+1 : param_dim
			id_B += 1;
			dA = param_divide(dim1) * param_divide(dim2);
			Blst[id_B] = broadcast.( /, Blst[id_B], dA );
		end
	end
	return Bln, Blst;
end

function divB( BfieldLst, Bfield_ln, N )
	dimLstLeviCivita = [3,2,1];
	
	divBlst = [ zeros(Float64,N) for it in CartesianIndices(BfieldLst[1]) ];
	BfieldSh = [ zeros(Complex{Float64},N) for it in CartesianIndices(BfieldLst[1]) ];
	for Bdim = 1 : Bfield_ln
		vecDim = dimLstLeviCivita[Bdim];
		sh_lst = zeros(Int64, Bfield_ln);
		sh_lst[vecDim] = -1;
		@debug(typeof(divBlst));
		@debug(typeof(divBlst[1]));
		divBlst = divBlst .+ circshift!( BfieldSh, BfieldLst[Bdim],sh_lst ) .- BfieldLst[Bdim];
	end
	return divBlst;
end

end
