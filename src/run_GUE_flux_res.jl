using DegLocatorDiv
using Utils

mSz = 10;

resRat = 2;
numRes = 4;

divLstBase = [16,16,16];
divLst = copy( divLstBase );

resBase = 16;
nBase = 10;

itNum = 1;

seed = 1000;

for iRes = 1:numRes
	divLst .= resRat^(iRes-1) .* divLstBase;
	
	println( "n = $mSz, res = $(divLst)" );
	with_logger( errLogger )do 
		@time divB_profile_flux( mSz, divLst, itNum, seed; enumSaveMem = memEig );
	end
end
