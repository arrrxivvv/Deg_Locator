using DegLocatorDiv
using UMat

N=30;
num_it = 5;

param_dim = 3;
param_min = fill(0.0,param_dim);
param_max = fill(2*pi,param_dim);
param_divide = [80,16,16];
ratio = [1,1,1];

degObj = DegObj( N, param_dim, param_divide, param_min, param_max );

H_GUElst = DegLocatorDiv.HlstFunc( H_GUE, degObj.param_dim, degObj.N, false );

HmatFun = (H,xlst)-> Hmat_3comb_ratio!(H,xlst,H_GUElst, ratio);

@time begin
	Threads.@threads for ind in degObj.posLst
		HmatFun( degObj.HmatLst[ind], degObj.param_mesh[ind] );
	end
end
