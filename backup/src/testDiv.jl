using JLD
using Utils
using DegLocatorDiv
using UMat
using Random

avgNum = 1000; N = 10; param_divide = [80, 16, 16]; seedFedStart = 1000; fileNameMod = "_degObj"; param_dim = 3;

param_min = fill(0.0,param_dim);
param_max = fill(2*pi,param_dim);

fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum, seedFedStart; dim = param_dim );
fileNameAttr = string( fileNameAttr, fileNameMod );

seedFed = seedFedStart;
if seedFed > 0
	Random.seed!(seedFed);
end

fileNameMain = "deg";
fileNameMainDistill = "locDistilled";

degFileName = string( fileNameMain, "_", fileNameAttr );
varFileName = string( degFileName, jldType );

distillFileName = string( fileNameMainDistill, "_", fileNameAttr );
distillFileNameFull = string( distillFileName, jldType );
locLstPolsDistilled = load( distillFileNameFull, "locDistilledLst" );

posLocLst = load(varFileName, "posLocLst");
negLocLst = load(varFileName, "negLocLst");
posNlst = load(varFileName, "posNlst");
negNlst = load(varFileName, "negNlst");
H_GUE_lst = load(varFileName, "H_GUE_lst");
locLstPols = [ posLocLst, negLocLst ];

locLstPols = [ [ [ sortslices( locLstPols[iPol][it][n],dims=1 ) for n = 1 : N ] for it = 1 : avgNum ] for iPol = 1 : 2 ] ;

it = 1;
ratio = ones(Int64, 3);

H_GUElstCurrent = H_GUE_lst[it];
HmatFun = (H,xlst)-> Hmat_3comb_ratio!(H,xlst,H_GUElstCurrent, ratio);

degObj = DegObj( N, param_dim, param_divide, param_min, param_max );

posCount, negCount, posLocs, negLocs, BfieldLst, divBlstReal = locator_div( degObj, HmatFun );

locsPols = [posLocs, negLocs];
locsPols = [ [ sortslices(locsPols[iPol][n], dims=1) for n = 1 : N ] for iPol = 1 : 2 ];

locLstPolsPureTest, pureNlstPolsTest, locLstPolsPure, locLstPolsDirty = locLstPurify_detailedOutput( posLocLst, negLocLst, posNlst, param_divide );
locLstPolsPureRevTest, pureNlstPolsRevTest, locLstPolsPureRev, locLstPolsDirtyRev = locLstPurify_detailedOutput( posLocLst, negLocLst, posNlst, param_divide, isReversed = true );