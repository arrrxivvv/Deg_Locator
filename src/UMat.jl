module UMat
using LinearAlgebra
using CUDA

export Hmat_3sin, H_GUE, H_GOE, Hmat_3GUE, Hmat_3comb, U_mat_QWZ_raw, Hmat_3comb_ratio, Hmat_3comb_ratio!, Hmat_3comb_ratio_Hoffset!, H_GUElst_orth!

function U_element_fun( N, tau, V1, V2, ang1, ang2, k1, k2 )
	U_ans = 0;
	for n = 1 : N
		U_ans += ( 1/N * exp(1im * 2 * pi * (k1-k2) * n / N ) 
			 * exp(-1im * tau * V2 * cos( (2*pi*n + ang2)/N ) ) 
			 * exp(-1im * tau * V1 * cos( (2*pi*k2 + ang1) / N )) );
	end
	return U_ans
end

function U_mat_flq_raw(  N, tau, V1, V2, ang1, ang2 )
	k_num = N;
	U_mat = zeros( Complex{Float64}, (k_num, k_num) );
	for k1 = 1:k_num
		for k2 = 1 : k_num
			U_mat[k1,k2] = U_element_fun( N, tau, V1, V2, ang1, ang2, k1, k2 );
		end
	end
	return U_mat
end

pauliX = [0 1;1 0];
pauliY = [0 -1im;1im 0];
pauliZ = [1 0;0 -1];

function U_mat_QWZ_raw( ang1, ang2, u )
	return sin(ang1) * pauliX + sin(ang2) * pauliY + (u + cos(ang1) + cos(ang2)) * pauliZ
end

function Hmat_3sin(x1,x2,x3)
	sin(x1) * pauliX + sin(x2) * pauliY + sin(x3) * pauliZ;
end

function H_GUE( N )
	matTmp = randn( N,N ) + im * rand( N,N );
	return 1/sqrt(2) * ( matTmp + matTmp' );
end

function H_GOE( N )
	matTmp = randn( N,N );
	return 1/sqrt(2) * ( matTmp + matTmp' );
end

function Hmat_3GUE(N, x1,x2,x3)
	xLst = [x1,x2,x3];
	Hmat = zeros(N,N);
	for it = 1:length(xLst)
		Hmat += cos(xLst[it]) * H_GUE(N) + sin(xLst[it]) * H_GUE(N);
	end
	return Hmat
end

function Hmat_3comb( xLst, Hlst )
	Hmat = zeros( size(Hlst[1]) );
	for it = 1:length(xLst)
		Hmat += cos(xLst[it]) * Hlst[1,it] + sin(xLst[it]) * Hlst[2,it];
	end
	return Hmat;
end

function paramGridSetup!( paramGrid, param_divide, param_dim, param_min )
	ln = prod(param_divide);
	indx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stridex = blockDim().x * gridDim().x;
	
end

function Hmat_3comb_CUDA!( HmatLst, param_divide, HbaseLst )
	numDim = length(param_divide);
	indx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stridex = blockDim().x * gridDim().x;
	indy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    stridey = blockDim().y * gridDim().y;
	indz = (blockIdx().z - 1) * blockDim().z + threadIdx().z
    stridez = blockDim().z * gridDim().z;
	
	for kk = indz : stridez : size( HmatLst, 3 )
		paramLst = kk;
		for jj = indy : stridey : size( HmatLst, 2 )
			for ii = indx : stridex : size( HmatLst, 1 )
				HmatLst[ii,jj,kk] = HmatFunCu([ii,jj,paramLst]);
				# for dim = 1 : numDim
					# theta = thetaLst[kk][dim];
					# HmatLst[ii,jj,kk] += CUDA.cos(theta) * HbaseLst[ii,jj,2*dim-1] + CUDA.sin(theta) * HbaseLst[ii,jj,2*dim];
				# end
			end
		end
	end
end

function HmatGUE_CUDA( ii, jj, paramLst, Hlst )
	numH = size()
end

function H_GUElst_orth!( Hlst, N )
	for ii in eachindex(Hlst)
		Hlst[ii] = H_GUE(N);
		for jj = 1 : ii-1
			@debug "jj" jj
			@debug "ii" ii
			@debug "Hlst ii", Hlst[ii]
			@debug "Hlst jj", Hlst[jj]
			Hlst[ii] -= Hlst[jj] * tr( Hlst[ii] * Hlst[jj] ) / tr( Hlst[jj] * Hlst[jj] );
		end
	end
end

function Hmat_3comb_ratio( xLst, Hlst, ratio )
	Hmat = zeros( size(Hlst[1]) );
	for it = 1:length(xLst)
		Hmat += ratio[it] .* ( cos(xLst[it]) * Hlst[1,it] + sin(xLst[it]) * Hlst[2,it] );
	end
	return Hmat;
end

function Hmat_3comb_ratio!( Hmat, xLst, Hlst, ratio=nothing )
	if isnothing(ratio)
		ratio = fill(1,length(xLst));
	end
	Hmat .= 0;
	for it = 1:length(xLst)
		for ii in eachindex(Hmat)
			Hmat[ii] += ratio[it] * ( cos(xLst[it]) * Hlst[1,it][ii] + sin(xLst[it]) * Hlst[2,it][ii] );
		end
		# Hmat .+= ratio[it] * ( cos(xLst[it]) * Hlst[1,it] + sin(xLst[it]) * Hlst[2,it] );
	end
end

function Hmat_3comb_ratio_Hoffset!( Hmat, xLst, Hlst, ratio, Hoffset; c1 = 1, cOff = 1 )
	Hmat .= 0;
	for it = 1:length(xLst)
		for ii in eachindex(Hmat)
			Hmat[ii] += c1 * ratio[it] * ( cos(xLst[it]) * Hlst[1,it][ii] + sin(xLst[it]) * Hlst[2,it][ii] );
		end
		# Hmat .+= ratio[it] * ( cos(xLst[it]) * Hlst[1,it] + sin(xLst[it]) * Hlst[2,it] );
	end
	
	Hmat .+= cOff * Hoffset;
end
	
end
