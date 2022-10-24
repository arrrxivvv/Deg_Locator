using LinearAlgebra
# using ArgParse
using EllipsisNotation
using Distributions
using Random
using JLD
# using JLD2, FileIO
include("Umat.jl")
include("Chern.jl")

function dotEachCol( vec1, vec2 )
	return dot.( eachcol(vec1), eachcol(vec2) );
end

function wrap( point, bnd )
	isWrappedArr = point .> (bnd/2);
	point = -bnd .* isWrappedArr + point;
	
end

function Bfield( point_origin, HmatFun, dr_incre, N, dim; n=0 )
	incrDimLst = Int64.( 2 * ones(dim) );
	# vecLst_DimLst = vcat( [N], incrDimLst, dim );
	# DimLst = vcat( [N], incrDimLst );
	typeOfNvec = Array{Complex{Float64},2};
	# typeOfNvec = Vector{typeOfVec}(undef,N);
	vecLst = Array{ typeOfNvec,dim }( undef, Tuple(incrDimLst) );
	typeOfNval = Vector{Complex{Float64}};
	valLst = Array{ typeOfNval,dim }(undef, Tuple(incrDimLst));
	# block = zeros( Tuple(incrDimLst) );
	
	for it = 1:length(valLst)
		it_cart_arr = collect( Tuple( CartesianIndices(valLst)[it] ) );
		it_cart_arr_incr = it_cart_arr .- 1;
		point = point_origin + dr_incre * it_cart_arr_incr;
		
		Hmat = HmatFun( Tuple(point)... );
		valTmp, vecTmp = eigen( Hmat );
		# vecTmp = transpose(vecTmp);
		valTmp = real(valTmp);
		for n = 1:N
			vecTmp[:,n] = vecTmp[:,n] / norm(vecTmp[:,n]);
		end
		# print( valTmp );
		EigenSort( valTmp, vecTmp );
		valLst[it] = valTmp;
		vecLst[it] = vecTmp;
	end
	
	typeOfNlinks = Vector{Complex{Float64}};
	typeOfLinkArr = Array{typeOfNlinks,dim};
	linkLst = Vector{typeOfLinkArr}(undef,dim);
	
	for ln_dim = 1:dim
		sh_lst = zeros(dim);
		sh_lst[ln_dim] = 1;
		vecLst_sh = circshift(vecLst,Tuple(sh_lst));
		# linkLst[ln_dim] = dot.(eachcol(vecLst_sh),eachcol(vecLst));
		linkLst[ln_dim] = dotEachCol.(vecLst_sh,vecLst);
		linkLst[ln_dim] = ( x -> x ./ abs.(x) ).(linkLst[ln_dim]);
	end
	
	Bfield_ln = Int64( round( dim * (dim-1) / 2 ) );
	# typeOfBvec = Vector{Complex(Float64)}(undef,Bfield_ln);
	BfieldLst = Array{Complex{Float64}}(undef,(N,Bfield_ln));
	id_B = 0;
	for dim1 = 1:dim
		for dim2 = dim1+1:dim
			id_B += 1;
			id_org = Int64.( ones(dim) );
			id_dim1 = copy(id_org);
			id_dim2 = copy(id_org);
			id_dim1[dim1] = 2;
			id_dim2[dim2] = 2;
			parity = (-1)^(dim1+dim2-1);
			BfieldLst[:,id_B] = parity * 
				log.( ( linkLst[dim1][Tuple(id_org)...] ./ linkLst[dim1][Tuple(id_dim2)...] ) ./
				(linkLst[dim2][Tuple(id_org)...] ./ linkLst[dim2][Tuple(id_dim1)...]) ) / 1im;
		end
	end
	return BfieldLst;
end

function main()
	num_ini = 20;
	dim = 3;
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	num_it = 4000;
	dr_incre = 0.001;
	Hmat_seed = 1234;
	
	N=2;
	# HmatFun = Hmat_3sin;
	H_GUElst = fill( Array{Complex{Float64}}(undef,(N,N)), (2,dim) );
	Random.seed!(Hmat_seed);
	for it = 1:length(H_GUElst)
		H_GUElst[it] = H_GUE(N);
	end
	HmatFun = (x,y,z)->Hmat_3comb([x,y,z],H_GUElst);
	
	ini_arr = zeros( num_ini, dim );
	degen_arr = [ Vector{Vector{Float64}}(undef,0) for n = 1:N ];
	Random.seed!();
	for d = 1:dim
		ini_arr[:,d] = rand( Uniform(param_min[d],param_max[d]), num_ini, 1 );
	end
	println( 1, round.( ini_arr; digits=2 ) );
	
	for it_try = 1:num_ini
		for n = 1:N
			point = ini_arr[it_try,:];
			for it = 1:num_it
				BvecLst = Bfield( point, HmatFun, dr_incre, N, dim );
				Bvec = BvecLst[n,..];
				Bvec = Bvec[[3,2,1]];
				dr = real.( Bvec / norm(Bvec) * dr_incre );
				point = mod.( point+dr, 2*pi );
				# if mod(it,100) == 0
					# println( it, round.( real.( Bvec ); digits=2 ),  round.( real.( point ); digits=2 ) );
				# end
			end
			if isempty(degen_arr[n])
				push!( degen_arr[n], wrap( point, param_max ) );
				print( "pushed", round.( wrap( point, param_max ); digits=2 ), "\n" );
			else
				isDuplicate = false;
				for it = 1:length(degen_arr[n])
					dr_pt = abs.( point - degen_arr[n][it] );
					isWrappedArr = dr_pt .> pi;
					dr_pt = 2*pi * isWrappedArr + (-1) .^ isWrappedArr .* dr_pt;
					print( "\n", it_try, ", n ", n, " it ", it, ", " );
					print( round.( dr_pt; digits=5 ) );
					print( round.( point; digits=2 ) );
					print( round.( degen_arr[n][it]; digits=2 ), "\n" );
					if norm( dr_pt ) <= 10*dr_incre
						isDuplicate = true;
						println("caught ", round.( wrap( point, param_max );digits=2 ) );
						break;
					end
				end
				if !isDuplicate
					print( "point raw ", round.( point; digits=2 ), "\n" );
					# isWrappedArr = point .> pi;
					# point = -2*pi * isWrappedArr + point;
					point = wrap(point, param_max);
					push!( degen_arr[n], point );
					print( "pushed", round.( point; digits=2 ), "\n" );
				end
			end
		end
	end
	print("\n");
	println("matrix dimension: ", N);
	println("trace iteration count: ", num_it);
	println("number of trial points: ", num_ini);
	for n = 1 : N
		println( "state ", n, " :" );
		for it = 1 : length(degen_arr[n])
			println( "\t point ", it, ' ',
			round.( real.(degen_arr[n][it]); digits=2 ));
		end
	end
	
	# oFileName = "/output/deg_" + "N_" + string(N) + "_seed_" + string(Hmat_seed) + "_iniPoints_" + string(num_ini) + "_it_" + string(num_ini);
	# save( oFileName, "N", N, "seed", Hmat_seed, "iniPoints", num_ini, "iterations", num_ini, "degLst", degen_arr );
	
	# println( "state ", 2, ' ', round.( real.(degen_arr[:,2,:]); digits=2 ));
	# println(round.( imag.(degen_arr[:,1,:]); digits=2 ));
end