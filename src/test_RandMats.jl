using UMat
using EigCustom
using JLD2
using DegLocatorDiv
# using Infiltrator

itNum = 1000;
nDim = 3;
N = 50;
elT = Float64;

divX = 6;
xMax = 2*pi;
dx = xMax / divX;
xMaxPrev = xMax - dx;

HLstLst = [ H_GOE(N) for iCos = 1:2, iDim = 1:nDim, it = 1 : itNum ];

H2NLst = [ zeros(elT, 2*N, 2*N) for iCos = 1:2, iDim = 1 : nDim ];

eigWork = eigWorkStructFromNum!(2*N);

valLstLst = zeros(2*N, itNum, divX, divX, divX);
valLst = zeros(2*N);
# vecLst = zeros( ComplexF64, N, N );
vecLst = zeros( elT, 2*N, 2*N );

xLst = zeros(nDim);

Hmat = zeros(elT, 2*N, 2*N);

for it = 1 : itNum
	print("iteration: $it / $itNum                       \r");

	for iCos = 1:2, iDim = 1:nDim
		H2NLst[iCos, iDim] .= 0;
	end
	
	iCos = 1;
	for iDim = 1:nDim
		@view( H2NLst[iCos,iDim][1:N,N+1:2*N] ) .= HLstLst[iCos, iDim, it];
		@view( H2NLst[iCos,iDim][N+1:2*N,1:N] ) .= HLstLst[iCos, iDim, it];
	end
	
	iCos = 1;
	for iDim = 1:nDim
		@view( H2NLst[iCos,iDim][1:N,N+1:2*N] ) .= HLstLst[iCos, iDim, it];
		@view( H2NLst[iCos,iDim][N+1:2*N,1:N] ) .= HLstLst[iCos, iDim, it];
	end
	
	iCos = 2;
	for iDim = 1:nDim
		@view( H2NLst[iCos,iDim][1:N,1:N] ) .= HLstLst[iCos, iDim, it];
		@view( H2NLst[iCos,iDim][N+1:2*N,N+1:2*N] ) .= - HLstLst[iCos, iDim, it];
	end
	
	# @infiltrate
	
	for i1 = 1:divX, i2 = 1:divX, i3 = 1:divX
		xLst[1] = i1; 
		xLst[2] = i2;
		xLst[3] = i3;
		
		xLst .-= 1;
		xLst .*= dx;
		
		Hmat_3comb!( Hmat, xLst, H2NLst );
		
		eigenZheevrStruct!( Hmat, valLst, vecLst, eigWork );
		@view(valLstLst[:,it,i1,i2,i3]) .= valLst;
	end
end

save( "GOE_RSym_Lst.jld2", "valLstLst", valLstLst );