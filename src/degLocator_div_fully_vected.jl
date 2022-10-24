include("divB.jl")

function locator_div(N, HmatFun, param_dim, param_divide_num, param_max, param_min)
	dimLst = collect(1:param_dim);
	param_ln_num = param_divide_num;
	param_divide = param_divide_num * ones(param_dim);
	param_ln = Int64.( param_divide );
	param_step = ( param_max - param_min ) ./ param_divide;
	param_grids = hcat( collect( range( param_min, length=param_ln_num, stop=param_max-param_step ) ) );
	
	HmatLst = Array{ Array{Complex{Float64},2}, param_dim }(undef,Tuple(param_ln));
	
	param_arr = zeros(param_dim);
	for ind in CartesianIndices(HmatLst)
		for j = 1 : param_dim
			param_arr[j] = param_grids[ind[j]][j];
		end
		HmatLst[ind] = HmatFun( param_arr... );
	end
	
	@time eigenLst = eigen.( HmatLst );
	@debug(size(eigenLst));
	@debug(size(HmatLst));
	Elst = (p-> Real.( p.values )).( eigenLst );
	@debug(typeof(Elst), size(Elst));
	vecLst = (p->p.vectors).( eigenLst );
	@time EigenSort!.( Elst, vecLst );
	
	HmatLst = nothing;
	GC.gc();

	@time linkLst = linkLstFun( param_dim, vecLst );
	
	
	@time Bfield_ln, BfieldLst = BfieldFromLink( param_dim, linkLst );
	
	divBlst = [ zeros(Float64,N) for it in CartesianIndices(Elst) ];
	dimLstLeviCivita = [3,2,1];
	for Bdim = 1 : Bfield_ln
		vecDim = dimLstLeviCivita[Bdim];
		sh_lst = zeros(Bfield_ln);
		sh_lst[vecDim] = -1;
		@debug(typeof(divBlst));
		@debug(typeof(divBlst[1]));
		divBlst = divBlst .+ circshift( BfieldLst[Bdim],sh_lst ) .- BfieldLst[Bdim];
	end
	
	threshold = 0.00001;
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
	@debug(size(divBlstReal[1]));
	
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