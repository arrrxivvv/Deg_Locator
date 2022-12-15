using DegLocatorDiv
using Utils

nlst = [10:5:60;];

nDim = 2;

divLstBase = [64,64];
divLst = copy( divLstBase );

resBase = 16;
nBase = 10;

itNum = 10;

seed = 1000;

for n in nlst
	divLst .= Int64.( floor.( sqrt( n / nBase ) .* divLstBase ) );
	
	println( "n = $n, res = $(divLst)" );
	with_logger( errLogger )do 
		@time divB_profile_flux( n, divLst, itNum, seed; enumSaveMem = memEig, nDim = nDim );
	end
end
