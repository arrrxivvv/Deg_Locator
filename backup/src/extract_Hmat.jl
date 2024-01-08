using DegLocatorDiv
using JLD
using UMat
using Utils
using LinearAlgebra

function readLocLst( N, param_divide, Mmat_seed, avgNum )
	param_dim = 3;
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	
	# H_GUElst = fill( Array{Complex{Float64}}(undef,(N,N)), (2,param_dim) );
	# if Hmat_seed > 0
		# Random.seed!(Hmat_seed);
	# end
	# for it = 1:length(H_GUElst)
		# H_GUElst[it] = H_GUE(N);
	# end
	# HmatFun = (x,y,z)->Hmat_3comb([x,y,z],H_GUElst);
	
	dimLst = collect(1:param_dim);
	# param_ln_num = param_divide_num;	
	param_divide = Int64.( param_divide );
	if !( param_divide isa Array{Int64} )
		param_divide_num = param_divide;
		param_divide = param_divide_num * ones(param_dim);
		param_divide = Int64.( param_divide );
	end
	param_step = ( param_max - param_min ) ./ param_divide;
	param_grids = [ collect( range( param_min[iDim], param_max[iDim] - param_step[iDim], length = param_divide[iDim] ) ) for iDim = 1:param_dim ];
	dummyMat = zeros( Tuple(param_divide) );
	
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum,  Mmat_seed );
	degFileName = string( "deg", "_", fileNameAttr );
	varFileName = string( degFileName, jldType );
	data = load(varFileName);
	
	posLocLst = data["posLocLst"];
	negLocLst = data["negLocLst"];
	H_GUE_lst = data["H_GUE_lst"];
	
	HmatFun = [ (x,y,z)->Hmat_3comb([x,y,z],H_GUE_lst[it]) for it = 1:avgNum ];
	
	return posLocLst, negLocLst, HmatFun, H_GUE_lst, param_grids, param_step;
end

function getH1Lst( posLoc1Lst, negLoc1Lst, HmatFun1, param_dim, param_grids, param_step )
	locLstPol = [posLoc1Lst, negLoc1Lst];
	N = size(posLoc1Lst)[1];
	HmatLstPol = [];
	paramLstPol = [];
	for idPol = 1:2
		locLst = locLstPol[idPol];
		paramLst = broadcast( x->Float64.(x), locLst );
		for n = 1 : N
			for iDeg = 1 : size( locLst[n] )[1]
				# println( size(paramLst[n]) );
				paramLst[n][iDeg,:] = nToParam( locLst[n][iDeg,:], param_dim, param_grids, param_step );
			end
		end
		HmatLst = [ [ HmatFun1( paramLst[n][iDeg,:]... ) for iDeg = 1 : size(locLst[n])[1] ] for n = 1 : N ];
		append!(HmatLstPol, [HmatLst]);
		append!(paramLstPol,[paramLst]);
	end
	
	return HmatLstPol, paramLstPol;
end

function nToParam( nLst, param_dim, param_grids, param_step )
	param_arr = zeros(param_dim);
	for j = 1 : param_dim
		param_arr[j] = param_grids[j][nLst[j]] + param_step[j]/2 ;
	end
	
	# param_arr = param_arr .+ 
	return param_arr;
end

function printHmatLst( HmatLstPol, N, n )
	# for idPol = 1 : 2
		# for it = 1 : length( HmatLstPol[1][n] )
			# display( round.( HmatLstPol[1][n][it]; digits=2 ) );
		# end
	# end
	for it = 1 : length( HmatLstPol[1][n] )
		display( hcat( round.( HmatLstPol[1][n][it]; digits=2 ), zeros(N), round.( HmatLstPol[2][n][it]; digits=2 ) ) );
	end
end

function printHmatLstEigen( HmatLstPol, N, n )
	# for idPol = 1 : 2
		# for it = 1 : length( HmatLstPol[1][n] )
			# display( round.( HmatLstPol[1][n][it]; digits=2 ) );
		# end
	# end
	for it = 1 : length( HmatLstPol[1][n] )
		ElstPos = eigen( HmatLstPol[1][n][it] ).values;
		ElstNeg = eigen( HmatLstPol[2][n][it] ).values;
		display( hcat( round.( ElstPos; digits=2 ), zeros(N), round.( ElstNeg; digits=2 ) ) );
	end
end
