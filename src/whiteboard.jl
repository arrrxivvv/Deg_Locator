locLstCumPol =[ 
	if n == 1 
		zeros(2);
	elseif n == 3
		zeros(N-1);
	else 
		ones(3);
	end for n = 1 : 3 ];