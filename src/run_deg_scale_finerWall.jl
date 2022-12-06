using DegLocatorDiv
using Utils

nStep = 10;
nLst = collect(10:nStep:60);

ratFine = 4;

nMin = 10;
nMax = 50;
dNres = nMax - nMin;

resMin = 16;
resMax = 32;
dRes = resMax - resMin;

divLst = [80,16,16];

itNum = 100;
seed = 1000;

dim = 3;

for n in nLst
	res = 16 + Int64(floor( ( n - nMin ) / dNres * dRes ));
	println("n = $n, res = $res, ratFine = $ratFine");
	
	divLst .= [80,res,res];
	with_logger(errLogger) do 
		@time divB_profile_flux_cell( n, divLst, itNum, seed; ratFine = ratFine );
	end
	GC.gc();
end
