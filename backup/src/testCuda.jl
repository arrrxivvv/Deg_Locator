using CUDA
using LinearAlgebra

M = 32;
matNum = 50000;

matCxLst = [rand( Complex{Float64}, M, M ) for ii = 1:matNum ];
matCxLst .+= conj.( transpose.(matCxLst) );

# matCxArr = cat( matCxLst..., dims=3 );
# matCxArrCu = cu( matCxArr );

matCxArr = rand( Complex{Float64}, M, M, matNum );
matCxArr += conj.( permutedims( matCxArr, [2,1,3] ) );

matReLst = [rand( M, M ) for ii = 1:matNum ];
matReLst .+= transpose.( matReLst );

matReArr = cat( matReLst..., dims=3 );
matReArrCu = cu( matReArr );

matRe2Lst = [ rand( 2*M, 2*M ) for ii = 1:matNum ];
matRe2Lst .+= transpose.( matRe2Lst );

# matRe2Arr = cat( matRe2Lst..., dims=3 );
matRe2Arr = rand( 2*M, 2*M, matNum );

matRe2ArrCu = cu( matRe2Arr );