
# zakCorrLstCmplx = similar( zakCorrLst1Pass, ComplexF64 );

# for iN = 1 : lnNCirc, iR = 1 : lnRCirc, it = 1 : itNum
	# RandomCircle.calcZakCorr!( @view( zakCorrLst1Pass[:,:,it,iR,iN] ), @view( zakCorrLstCmplx[:,:,it,iR,iN] ), @view( zakArrLst1Pass[:,:,it,iR,iN] ) );
# end

# divNum1PassHalf = div( divNum1Pass, 2 );
# xLstZakHalf = xLstZak[1:divNum1PassHalf];

# zakCorrLstAvg1Pass = mean( zakCorrLst1Pass; dims = 3 );
# zakCorrLstAvg1d1Pass = zakCorrLstAvg1Pass[:,1,1,:,:];
# zakCorrLstAvg1d1PassArrLst = [ @view( zakCorrLstAvg1d1Pass[1:divNum1PassHalf,iR,iN] ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

# fittedModelLst = [ curve_fit( modelExp, xLstZakHalf, zakCorrLstAvg1d1PassArrLst[iR,iN], p0Fit ) for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];

# expScaleLst = (x->x.param[1]).( fittedModelLst );
# expShLst = (x->x.param[3]).( fittedModelLst );
# corrLenLst = (x->x.param[2]).( fittedModelLst );

# jldsave( "zakCorrTest"; zakCorrAvg1d = zakCorrLstAvg1d1Pass, expScaleLst, expShLst, corrLenLst );

dataTest = randCircDataLst[3];

RandomCircle.setZakArr!( dataTest, @view zakArrLst1Pass[:,:,1,1,3] );

RandomCircle.calcZakCorr!( dataTest );

zakCorrTest = RandomCircle.getZakCorr( dataTest );

zakCorr1dTest = zakCorrTest[:,1];

fittedModelTest = curve_fit( modelExp, xLstZakHalf, @view( zakCorr1dTest[1:divNum1PassHalf] ), p0Fit );

epSc = fittedModelTest.param[1];
corrLen = fittedModelTest.param[2];
epSh = fittedModelTest.param[3];

jldsave( "zakCorrTest.jld2"; zakArr = RandomCircle.getZakArr( dataTest ), zakCorr = RandomCircle.getZakCorr( dataTest ), epSc, corrLen, epSh );
