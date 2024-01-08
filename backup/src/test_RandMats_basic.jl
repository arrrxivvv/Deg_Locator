using UMat
using EigCustom
using JLD2
using DegLocatorDiv

itNum = 1000;
nDim = 3;
N = 50;

HLstLst = [ Hmat_3comb_HLstGen(N, nDim, GOE) for it = 1 : itNum ];

H2NLst = [ zeros() for dim = 1 : nDim, iCos = 1:2 ];

eigWork = eigWorkStructFromNum!(N);

valLstLst = zeros(N, itNum);
valLst = zeros(N);
# vecLst = zeros( ComplexF64, N, N );
vecLst = zeros( Float64, N, N );

for it = 1 : itNum
	Hmat = H_GOE(N);
	
	eigenZheevrStruct!( Hmat, valLst, vecLst, eigWork );
	
	@view(valLstLst[:,it]) .= valLst;
end

save( "GOEValLst.jld2", "valLstLst", valLstLst );