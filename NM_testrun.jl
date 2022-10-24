using NelderMeadTest

dim = 2;
numPt = dim+1;

valLst = zeros( numPt );
ptLst = zeros( dim, numPt );

xSt = 1;
ySt = 1;
dx = 0.1;
dy = 0.1;

@view(ptLst[1,:]) .= xSt;
@view(ptLst[2,:]) .= ySt;
ptLst[1,2] += dx;
ptLst[2,3] += dy;

valLstTmp = similar(valLst);
ptLstTmp = similar(ptLst);
ptLstOrg = deepcopy(ptLst);

ixLst = similar( valLst, Int64 );

cnt = nmOpt!( cosXY, dim, ptLstOrg, valLst, ixLst, ptLst, ptLstTmp, valLstTmp );