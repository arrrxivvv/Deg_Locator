include("utils.jl")
include("divB.jl")
using FFTW
using Logging
using EllipsisNotation

function structure_factor_FFT( N, param_dim, param_divide, param_step, BfieldLst )
	tempB = BfieldLst[1];
	BfieldLst[1] = BfieldLst[3];
	BfieldLst[3] = tempB;
	BfieldArr = arrLstToArr.( BfieldLst );
	BfieldFFT = ( field -> fft( field, collect(1:param_dim) ) ).(BfieldArr);
	
	@debug( "start kLst construction" );
	kLst = gridLst( param_dim, param_divide );
	
	dvecLst = ( x-> ( exp.(+1im * x * param_step) .- 1 ) / param_step ).( kLst );
	
	@debug( "end kLst construction" );
	structFactFFT = sum( broadcast.( *, dvecLst, BfieldFFT ) );
	
	return structFactFFT;
end

function structure_factor_true( N, param_dim, param_divide, param_step, posLocLst, negLocLst )
	locLst = [ posLocLst, negLocLst ] .* param_step;
	kLst = gridLst( param_dim, param_divide );
	
	ln_dim = ones( Int64, param_dim ) * param_divide;
	full_ln = ln_dim;
	push!(full_ln, N);
	expArr = zeros( ln_dim... );
	structFact = zeros( Complex{Float64}, full_ln... );
	for n = 1 : N
		for idPol = 1:2
			polarity = (-1)^(idPol-1);
			for idLoc = 1 : size( locLst[idPol][n] )[1]
				expArr = 0;
				loc = locLst[idPol][n][idLoc,:];
				for dim = 1 : param_dim
					expArr = ( expArr .+ kLst[dim] * loc[dim] );
				end
				structFact[..,n] += polarity * exp.(-1im * expArr);
			end
		end
	end
	return structFact;
end

function structure_factor_true( N, param_dim, param_divide, param_step, posLocLst, negLocLst )
	locLst = [ posLocLst, negLocLst ] .* param_step;
	kLst = gridLst( param_dim, param_divide );
	
	ln_dim = ones( Int64, param_dim ) * param_divide;
	full_ln = ln_dim;
	push!(full_ln, N);
	expArr = zeros( ln_dim... );
	structFact = zeros( Complex{Float64}, full_ln... );
	for n = 1 : N
		for idPol = 1:2
			polarity = (-1)^(idPol-1);
			for idLoc = 1 : size( locLst[idPol][n] )[1]
				expArr = 0;
				loc = locLst[idPol][n][idLoc,:];
				for dim = 1 : param_dim
					expArr = ( expArr .+ kLst[dim] * loc[dim] );
				end
				structFact[..,n] += polarity * exp.(-1im * expArr);
			end
		end
	end
	return structFact;
end

function structure_factor_pol( N, param_dim, param_divide, param_step, locLstRaw, polarity )
	locLst = locLstRaw .* param_step;
	kLst = gridLst( param_dim, param_divide );
	
	ln_dim = ones( Int64, param_dim ) * param_divide;
	full_ln = ln_dim;
	push!(full_ln, N);
	expArr = zeros( ln_dim... );
	structFact = zeros( Complex{Float64}, full_ln... );
	for n = 1 : N
		for idLoc = 1 : size( locLst[n] )[1]
			expArr = 0;
			loc = locLst[n][idLoc,:];
			for dim = 1 : param_dim
				expArr = ( expArr .+ kLst[dim] * loc[dim] );
			end
			structFact[..,n] += polarity * exp.(-1im * expArr);
		end
	end
	return structFact;
end