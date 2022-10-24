N = 10;
param_dim = 3;
divNum = 20;
param_divide_tmp_raw = [divNum,1,1];
param_divide_tmp = param_divide_tmp_raw + ones(Int64, param_dim);
param_step_tmp = degObj.param_step ./ param_divide_tmp_raw;
ratio = ones(Int64, param_dim);

locStart = [6,10,2];
# sh = [0.5,0,0.5] .* degObj.param_step;
sh = zeros(Int64, 3);
# sh = [1,0,0] .* degObj.param_step;
shEndScale = 1;
shEnd = shEndScale * fill(1,3) .* degObj.param_step;
param_min_tmp = degObj.param_mesh[ locStart... ] + sh;
param_max_tmp = param_min_tmp + shEnd + param_step_tmp;

degObjTmp = DegObj( N, param_dim, param_divide_tmp, param_min_tmp, param_max_tmp, nonPeriodic = true );
# HmatFun = (H,xlst)-> Hmat_3comb_ratio!(H,xlst,H_GUElstCurrent, ratio);

posNlstTmp, negNlstTmp, posLocLstTmp, negLocLstTmp, H_GUE_lstTmp = locator_div( degObjTmp, HmatFun ); 