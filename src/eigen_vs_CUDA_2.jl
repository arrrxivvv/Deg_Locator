using LinearAlgebra
using CUDA

N = 20;
matNum = 1000;

matReLst = [ Symmetric( rand(Float32, N,N) ) for it = 1 : matNum ];
matReArr = zeros( N, N, matNum );
for n = 1 : matNum
	matReArr[:,:,n] = matReLst[n];
end
matReArrCu = cu(matReArr);

function testEigenRe()
	Threads.@threads for ii in eachindex(matReLst)
		eigen(matReLst[ii]);
	end
end

function testCuSolRe()
	sols = CUDA.CUSOLVER.syevjBatched!('V','U',matReArrCu);
end
