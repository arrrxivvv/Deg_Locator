module DivB

using LinearAlgebra
using EllipsisNotation
using Distributions
using Random
using JLD
using Logging; using LoggingExtras;
using ShiftedArrays

using Utils
using UMat
# include("Chern.jl")

export divB!, linkLstFun!, BfieldFromLink!, Bfield_dA

function linkLstFun!( degObj )
	for dim = 1:degObj.param_dim
		sh_lst = zeros(Int64,degObj.param_dim);
		sh_lst[dim] = +1;
		shCartInd = CartesianIndex( Tuple(sh_lst) );
		Threads.@threads for pos in degObj.posLst
			nextPos = pos + shCartInd;
			nextPos = wrapCartInd!( nextPos, degObj.vecLst );
			dotEachCol!( degObj.linkLst[dim][pos], degObj.vecLst[nextPos], degObj.vecLst[pos] );
		end
		Threads.@threads for pos in degObj.posLst
			degObj.linkLst[dim][pos] .= degObj.linkLst[dim][pos] ./ abs.(degObj.linkLst[dim][pos]);
		end
	end
	return degObj.linkLst;
end

function BfieldFromLink!( degObj )
	id_B = 0;
	linkRatio = degObj.linkRatio;
	linkLstSh = degObj.linkLstSh;
	for dim1 = 1 : degObj.param_dim
		for dim2 = dim1+1 : degObj.param_dim
			dimLst = [dim1,dim2];
			dimOtherLst = [dim2,dim1];
			id_B += 1;
			for it = 1 : 2
				dim = dimLst[it];
				dimSh = dimOtherLst[it];
				sh_lst = zeros(Int64, degObj.param_dim);
				linkLstShTmp = ShiftedArrays.circshift( degObj.linkLst[dim], ntuple( i -> i==dimSh ? -1 : 0, degObj.param_dim ) );
				Threads.@threads for pos in degObj.posLst
					linkRatio[it][pos] .= degObj.linkLst[dim][pos] ./ linkLstShTmp[pos];
				end
			end
			parity = (-1)^(dim1+dim2-1);
			Threads.@threads for ii in degObj.posLst
				degObj.BfieldLst[id_B][ii] .= parity ./ 1im .* log.( linkRatio[1][ii] ./ linkRatio[2][ii] );
			end
		end
	end
	return degObj.Bfield_ln, degObj.BfieldLst;
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

function divB!( degObj )
	dimLstLeviCivita = [3,2,1];
	
	for ii in eachindex(degObj.divBLst)
		degObj.divBLst[ii] .= 0;
	end
	# BfieldSh = [ zeros(Complex{Float64},degObj.N) for it in degObj.posLst ];
	BfieldSh = degObj.BfieldSh;
	for Bdim = 1 : degObj.Bfield_ln
		vecDim = dimLstLeviCivita[Bdim];
		sh_lst = zeros(Int64, degObj.Bfield_ln);
		sh_lst[vecDim] = -1;
		circshift!( BfieldSh, degObj.BfieldLst[Bdim],sh_lst );
		Threads.@threads for ii in eachindex(degObj.divBLst)
			degObj.divBLst[ii] .+= BfieldSh[ii] .- degObj.BfieldLst[Bdim][ii];
		end
		# degObj.divBLst .+= circshift!( BfieldSh, degObj.BfieldLst[Bdim],sh_lst ) - degObj.BfieldLst[Bdim];
	end
	return degObj.divBLst;
end

end
