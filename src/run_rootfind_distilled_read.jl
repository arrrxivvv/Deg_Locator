ratLst = [10:10:200;];
ratVal = 10;
ratLn = length(ratLst);
itNum = 100;
fMain = "NLstRootDistilled";
nDim = 3;
seed = 1000;
mSz = 10;
res = 15;
divLst = fill(res, nDim);
thresVal = 1e-9;
thresSz = 1e-9;
fExt = jld2Type;

NLstLst = zeros(Int64, itNum, ratLn);

for iRat = 1:ratLn
	ratSz = ratLst[iRat];
	attrMoreLst = ["thresVal", "thresSz", "thresValRatio", "thresSzRatio"];
	valMoreLst = Any[thresVal, thresSz, ratVal, ratSz];
	
	attrLst, valLst = fAttrOptLstFunc( mSz, divLst, itNum, seed; attrMoreLst = attrMoreLst, valMoreLst = valMoreLst, dim = nDim );
	fName = fNameFunc( fMain, attrLst, valLst, fExt );
	
	NLst = load( fName, "NLstDistilledLst" );
	@view( NLstLst[:,iRat] ) .= NLst;
end

NmeanLst = mean( NLstLst; dims=1 );
save( "NLstRootfind_N_$(mSz)_res_$(res).jld2", "NMeanLst", NmeanLst, "thresSz", thresSz, "szRatLst", ratLst );