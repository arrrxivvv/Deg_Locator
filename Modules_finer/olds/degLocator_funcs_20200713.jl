using LinearAlgebra
using EllipsisNotation
using Random
using Logging
using LoggingExtras

function locateNon0Div( divBlst, threshold, N, param_dim; clipBnd = 0 )
	divBlstReal = (p-> real.(p)).(divBlst);
	if clipBnd > 0
		divBlstReal = divBlstReal[1:end-clipBnd,1:end-clipBnd,1:end-clipBnd];
	end
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

function locator_div(N, HmatFun, param_dim, param_divide, param_max, param_min; maxDepth = 0, depth = 0, nSpecific = 0)
	dimLst = collect(1:param_dim);
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
	if depth > 0
		Nlimited = 1;
		vecLst = (x->x[:,nSpecific]).(vecLst);
		Elst = (x->[x[nSpecific]]).(Elst);
		clipBnd = 1;
		# Nlimited = N;
	else
		Nlimited = N;
		clipBnd = 0;
	end
	@time EigenSort!.( Elst, vecLst );
	
	HmatLst = nothing;
	GC.gc();

	@time linkLst = linkLstFun( param_dim, vecLst );
	@time Bfield_ln, BfieldLst = BfieldFromLink( param_dim, linkLst );
	@time divBlst = divB( BfieldLst, Bfield_ln, Nlimited );
	
	threshold = 0.00001;
	divBlstReal = (p-> real.(p)).(divBlst);
	posCount, negCount, posLocLst, negLocLst = locateNon0Div( divBlst, threshold, Nlimited, param_dim; clipBnd = clipBnd );
	
	locToParamHere = x -> locToParamArr( x, param_grids, param_dim );
	
	posParamLst = locToParamHere.( posLocLst );
	negParamLst = locToParamHere.( negLocLst );
	
	paramLstPol = [posParamLst, negParamLst];
	locLstPol = [posLocLst, negLocLst];
	countPol = [posCount, negCount];
	
	# println("depth: ", depth);
	
	if maxDepth > 0 && depth < maxDepth
		for iPol = 1 : 2
			paramLst = paramLstPol[iPol];
			locLst = locLstPol[iPol];
			degCount = countPol[iPol];
			for n = 1 : Nlimited
				itStop = degCount[n];
				if depth > 0
					itStop = min(1,itStop);
				end
				for it = 1 : itStop
					param = paramLst[n][it,:];
					param_min_next = param;
					param_max_next = param .+ 1.5 * param_step;
					param_divide_next = 3;
					depthNext = depth + 1;
					if depth > 0
						nSpecificNext = nSpecific;
					else
						nSpecificNext = n;
					end
					# if depth >0
						# display( divBlst[ locLst[n][it,:]... ][n] );
						# display(param_min);
						# display(param_step / (2*pi));
						# println("it: ", it, ", n", n);
						# display(locLst[n]);
						# display(degCount[n]);
					# end
					trNum = 3;
					for iTr = 1 : trNum
						posCountFine,negCountFine,posParamLstFine,negParamLstFine, BfieldLstFine, divBlstFine = locator_div( N, HmatFun, param_dim, param_divide_next, param_max_next, param_min_next; maxDepth=maxDepth, depth=depthNext,nSpecific=nSpecificNext );
						paramLstFinePol = [posParamLstFine, negParamLstFine];
						countFinePol = [posCountFine,negCountFine];
						countFine = countFinePol[iPol];
						paramLstFine = paramLstFinePol[iPol];
						nLimInd = 1;
						# println(posCountFine, ' ', negCountFine);
						# println( paramLstFine[nLimInd] );
						# if depth>0
							# display(divBlstFine);
						# end
						if countFine[nLimInd] > 0
							paramLst[n][it,:] = paramLstFine[nLimInd][1,:];
							break;
						else 
							param_divide_next += 1;
						end
					end
				end
			end
		end
	end
	
	return posCount, negCount, posParamLst, negParamLst, BfieldLst, divBlstReal, posLocLst, negLocLst;
end

function locator_div_GUE( N, param_divide_num; Hmat_seed = 1234, maxDepth = 0, polarity = +1)
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
	
	posCount, negCount, posParamLst, negParamLst, BfieldLst, divBlstReal, posLocLst, negLocLst = locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min; maxDepth = maxDepth ); 
	return posCount, negCount, posParamLst, negParamLst, H_GUElst, BfieldLst, divBlstReal, posLocLst, negLocLst;
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
	HmatFun = (x,y,z)-> scale .* Hmat_3comb_ratio([x,y,z],H_GUElst, ratio);
	
	return locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min );
end

function locator_div_sin3( param_divide_num )
	N=2;	
	HmatFun = Hmat_3sin;
	
	param_dim = 3;
	param_max = [2*pi,2*pi,2*pi];
	param_min = [0,0,0];
	
	return locator_div( N, HmatFun, param_dim, param_divide_num, param_max, param_min );
end

function locToParam( loc, param_grids, param_dim )
	paramArr = zeros( param_dim );
	for iDim = 1 : param_dim
		paramArr[iDim] = param_grids[iDim][loc[iDim]]
	end
	
	return paramArr;
end

function locToParamArr( locArr, param_grids, param_dim )
	paramArr = zeros( size(locArr) );
	for it = 1 : size(locArr)[1]
		paramArr[it,:] = locToParam( locArr[it,:], param_grids, param_dim );
	end
	return paramArr;
end
