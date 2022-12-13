using DegLocatorDiv
using Utils

nlst = [10:5:60;];

divLstBase = [80, 32, 32];
divLst = copy( divLstBase );

resBase = 16;
nBase = 10;

itNum = 100;

seed = 1000;

for n in nlst
	divLst .= Int64.( floor.( sqrt( n / nBase ) .* divLstBase ) );
	
	println( "n = $n, res = $(divLst[end])" );
	with_logger( errLogger )do 
		@time divB_profile_flux( n, divLst, itNum, seed; enumSaveMem = memEig );
	end
end
