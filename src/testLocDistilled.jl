locNCum = zeros( avgNum, N-1 );

for iA = 1 : avgNum
	for n = 1 : N-1
		locNCum[iA,n] = size( locCumPol[1][iA][n], 1 ) - size( locCumPol[2][iA][n], 1 );
	end
end

