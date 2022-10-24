using DegLocatorDiv
using Infiltrator

function testDegGC()
	degObjTest = DegObj( 10, 3, [64,64,64], fill(0.0,3), fill(2*pi,3); );
	@infiltrate
	degObjTest=nothing;
	@infiltrate
	GC.gc();
	@infiltrate;
end

