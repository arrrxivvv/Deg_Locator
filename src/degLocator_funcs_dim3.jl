using LinearAlgebra
using EllipsisNotation
using Random
using Logging
using LoggingExtras

function locateNon0Div( divBlst, threshold, N, param_dim )
	divBlstReal = (p-> real.(p)).(divBlst);
	divBposInd = broadcast.( >, divBlstReal, threshold );
	divBnegInd = broadcast.( <, divBlstReal, -threshold );
	divBposArr = arrLstToArr( divBposInd );
	divBnegArr = arrLstToArr( divBnegInd );
	posLocLst = Vector{ Array{ Int64,2 } }(undef,N);
	negLocLst = Vector{ Array{ Int64,2 } }(undef,N);
	for n = 1:N
		posLocLst[n] = cartIndLstToArr( findall( divBposArr[..,n] ), param_dim );
		negLocLst[n] = cartIndLstToArr( findall( divBnegArr[..,n] ), param_dim );
	end	
	posCount = sum( divBposInd );
	negCount = sum( divBnegInd );
	
	return posCount, negCount, posLocLst, negLocLst;
end

function locLstDistill( N, param_dim, posLocLst, negLocLst )
	posLocLstPure = Vector{ Array{ Array{Int64} } }(undef,N);
	negLocLstPure = Vector{ Array{ Array{Int64} } }(undef,N);
	locLstPols = broadcast( (x->arrToArrLst.(x)), [posLocLst, negLocLst] );
	locLstPolsPure = [posLocLstPure, negLocLstPure];
	idPolsOppLst = [2,1];
	for ip = 1 : 2, n in [1,N]
		locLstPolsPure[ip][n] = locLstPols[ip][n];
	end
	for n = 2 : N-1
		for idPols = 1 : 2
			idPolsOpp = idPolsOppLst[idPols];
			locLstPolsPure[idPols][n] = filter( !(x->x in locLstPolsPure[idPolsOpp][n-1]), locLstPols[idPols][n] );
			if length( locLstPolsPure[idPols][n] ) != length( locLstPols[idPols][n] ) - length(locLstPolsPure[idPolsOpp][n-1])
				@warn( "imprecise exclusion: ", string( length( locLstPols[idPols][n] ) - length(locLstPolsPure[idPolsOpp][n-1]) - length( locLstPolsPure[idPols][n] ) ) );
			end
		end
	end
	locLstPolsPure = broadcast( (x->arrLstToArr.(x)), locLstPolsPure );
	
	return locLstPolsPure[1], locLstPolsPure[2];
end

function locator_div(N, HmatFun, param_dim, param_divide, param_max, param_min)
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
		
	@time HmatLst = [ zeros( Complex{Float64}, N, N ) for i1 in CartesianIndices(dummyMat) ];
	
	@time begin
		param_arr = zeros(param_dim);
		for ind in CartesianIndices(dummyMat)
			for j = 1 : param_dim
				param_arr[j] = param_grids[j][ind[j]];
			end
			HmatLst[ind] = HmatFun( param_arr... );
		end
	end
	
	@time eigenLst = eigen.( HmatLst );
	Elst = (p-> Real.( p.values )).( eigenLst );
	vecLst = (p->p.vectors).( eigenLst );
	@time EigenSort!.( Elst, vecLst );
	
	HmatLst = nothing;
	GC.gc();

	@time linkLst = linkLstFun( param_dim, vecLst );
	@time Bfield_ln, BfieldLst = BfieldFromLink( param_dim, linkLst );
	@time divBlst = divB( BfieldLst, Bfield_ln, N );
	
	threshold = 0.00001;
	divBlstReal = (p-> real.(p)).(divBlst);
	posCount, negCount, posLocLst, negLocLst = locateNon0Div( divBlst, threshold, N, param_dim );
	
	return posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal;
end

function locator_div_GUE( N, param_divide_num; Hmat_seed = 1234, polarity = +1 )
	param_dim = 3;
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	
	H_GUElst = fill( Array{Complex{Float64}}(undef,(N,N)), (2,param_dim) );
	if Hmat_seed > 0
		Random.seed!(Hmat_seed);
	end
	for it = 1:length(H_GUElst)
		H_GUElst[it] = H_GUE(N);
	end
	HmatFun = (x,y,z)->Hmat_3comb([x,y,z],H_GUElst);
	
	posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min ); 
	return posCount, negCount, posLocLst, negLocLst, H_GUElst, BfieldLst, divBlstReal;
end

function locator_div_GUE_scale( N, param_divide_num, scale; Hmat_seed = 1234, polarity = +1 )
	param_dim = 3;
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	
	H_GUElst = fill( Array{Complex{Float64}}(undef,(N,N)), (2,param_dim) );
	if Hmat_seed > 0
		Random.seed!(Hmat_seed);
	end
	for it = 1:length(H_GUElst)
		H_GUElst[it] = H_GUE(N);
	end
	HmatFun = (x,y,z)-> scale .* Hmat_3comb([x,y,z],H_GUElst);
	
	return locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min );
end

function locator_div_GUE_ratio( N, param_divide_num; ratio=[1,1,1], Hmat_seed = 1234, polarity = +1 )
	param_dim = 3;
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	
	H_GUElst = fill( Array{Complex{Float64}}(undef,(N,N)), (2,param_dim) );
	if Hmat_seed > 0
		Random.seed!(Hmat_seed);
	end
	for it = 1:length(H_GUElst)
		H_GUElst[it] = H_GUE(N);
	end
	HmatFun = (x,y,z)-> Hmat_3comb_ratio([x,y,z],H_GUElst, ratio);
	
	posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min ); 
	return posCount, negCount, posLocLst, negLocLst, H_GUElst, BfieldLst, divBlstReal;
end

function locator_div_GOE_ratio( N, param_divide_num; ratio=[1,1], Hmat_seed = 1234, polarity = +1 )
	param_dim = 2;
	param_max = [2*pi,2*pi];
	param_min = [0,0];
	
	H_GUElst = fill( Array{Complex{Float64}}(undef,(N,N)), (2,param_dim) );
	if Hmat_seed > 0
		Random.seed!(Hmat_seed);
	end
	for it = 1:length(H_GOElst)
		H_GOElst[it] = H_GOE(N);
	end
	HmatFun = (x,y)-> Hmat_3comb_ratio([x,y],H_GOElst, ratio);
	
	posCount, negCount, posLocLst, negLocLst, BfieldLst, divBlstReal = locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min ); 
	return posCount, negCount, posLocLst, negLocLst, H_GUElst, BfieldLst, divBlstReal;
end

function locator_div_sin3( param_divide_num )
	N=2;	
	HmatFun = Hmat_3sin;
	
	param_dim = 3;
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	
	return locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min );
end