module UMat

export Hmat_3sin, H_GUE, Hmat_3GUE, Hmat_3comb, U_mat_QWZ_raw

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

end