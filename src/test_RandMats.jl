using UMat
using EigCustom

itNum = 100;
N = 10;

eigWork = eigPreworkStructFromNum(N);

valLstLst = zeros(N, itNum);
valLst = zeros(N);
vecLst = zeros( ComplexF64, N );

for it = 1 : itNum
	Hmat = H_GUE(N);
	
	eigenZheevrStruct!( Hmat, valLst, vecLst, eigWork );
	
	@view(valLstLst[:,it]) .= valLst;
end


