using DegLocatorDiv
using Utils

nlst = [10:10:60;];

divLstBase = [80, 16, 16];
divLst = copy( divLstBase );

resBase = 16;
nBase = 10;

itNum = 100;

seed = 1000;

for n in nlst
	divLst .= Int64.( floor.( n ./ nBase .* divLstBase ) );
	
	println( "n = $n, res = $(divLst[end])" );
	with_logger( errLogger )do 
		@time divB_profile_flux( n, divLst, itNum, seed; enumSaveMem = memEig );
	end
end
