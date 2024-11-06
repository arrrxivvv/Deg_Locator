
function genRandPt3dLst( numPts::Int64 )
	[ rand(nDim3)  for ii = 1 : numPts ];
end

function genRandPt2dFrom3dLst( pt3dLst::Vector{<:AbstractVector{Float64}} )
	[ pt3dLst[iPt][1:2] for iPt = 1 : length(pt3dLst) ];
end

function genRandQuartLst( numQuarts::Int64 )
	quartLst = [ OffsetArray( ( @MVector randn(4) ), idQuatRng ) for ii = 1 : numQuarts ];
	
	for ii = 1 : numQuarts
		quartLst[ii] .= quartLst[ii] ./ norm( quartLst[ii] );
	end
	
	return quartLst;
end

function genRotMatFromQuart( quart::OffsetVector{Float64,MVector{nDimQuat,Float64}} )
	rotMat = Utils.genStaticIdentityMat(nDim3, Float64);
	normQuart = norm(quart)^2;
	normSqV = 0;
	for ii = 1 : nDim3
		normSqV += quart[ii]^2;
	end
	
	for ii = 1 : nDim3
		rotMat[ii,ii] = rotMat[ii,ii] - 2*normSqV + 2 * quart[ii]^2;
	end
	
	for ii = 1 : nDim3, jj = 1 : nDim3
		if ii == jj 
			continue;
		end
		rotMat[ii,jj] = rotMat[ii,jj] + 2 * quart[ii] * quart[jj] - 2 * quart[0] * quart[leviCivita3rdIdMat[ii,jj]] * leviCivita3rdSgnMat[ii,jj];
	end
	
	return rotMat
end

function scaleRotMatLst!( rotMatLst::Vector{<:AbstractMatrix{Float64}}, scale::Float64 )
	(x -> x .= x .* scale).(rotMatLst);
end

function genRotMatLstFromQuartLst( quartLst::Vector{<:OffsetVector{Float64,<:MVector{nDimQuat,Float64}}} )
	return genRotMatFromQuart.( quartLst );
end

function genRotMatLst( numRots::Int64 )
	quartLst = genRandQuartLst( numRots );
	
	return genRotMatLstFromQuartLst( quartLst );
end

function genRotMat2dLst( rotMatLst::Vector{<:AbstractMatrix{Float64}} )
	rotMat2dLst = [ zeros(2, 2) for ii = 1 : length(rotMatLst) ];
	
	for iMat = 1 : length(rotMatLst), ii = 1 : 2, jj = 1 : 2
		rotMat2dLst[iMat][ii,jj] = rotMatLst[iMat][ii,jj];
	end
	
	return rotMat2dLst;
end

function genRotMat2dInvLst( rotMat2dLst::Vector{<:AbstractMatrix{Float64}} )
	rotMat2dInvLst = inv.( rotMat2dLst );
	
	return rotMat2dInvLst;
end

function genEllipseSampleLst( rotMat2dLst::Vector{<:AbstractMatrix{Float64}}, pt2dLst::Vector{<:AbstractVector{Float64}}, numSamples::Int64 )
	thetaLst = [1:numSamples;]';
	xyLst = vcat( cos.(thetaLst), sin.(thetaLst) );
	
	xySampleLst = [ rotMat2dLst[ii] * xyLst for ii = 1 : length(pt2dLst) ];
	((x,y) -> x .= x .+ y).(xySampleLst, pt2dLst);
	
	return xySampleLst;
end

function genEllipseWidthLst( rotMat2dLst::Vector{<:AbstractMatrix{Float64}}, rotMat2dInvLst::Vector{<:AbstractMatrix{Float64}} )
	areaLst = abs.( det.(rotMat2dLst) );
	xySqEqCoeffLst = dropdims.( sum.( x->x^2, rotMat2dInvLst; dims = 1 ); dims = 1 );
	
	interceptXYLst = reverse.( xySqEqCoeffLst );
	( x -> x .= 1 ./ sqrt.(x) ).(interceptXYLst);
	xyWidthLst = ( (x,y) -> x ./ y ).(areaLst, interceptXYLst);
	
	# @infiltrate
	
	return xyWidthLst;
end

function genEllipseBndFromWidthLst( xyWidthLst::Vector{<:AbstractVector{Float64}}, pt2dLst::Vector{<:AbstractVector{Float64}} )
	# bndLst = [ [ pt2dLst[iCirc][iXY] - xyWidthLst[iCirc][iXY] for iXY = 1 : 2 ] for iCirc = 1 : length(pt2dLst) ];
	
	bndLst = ( (pLst, wLst)->( (p,w)->[p-w,p+w] ).(pLst, wLst) ).(pt2dLst, xyWidthLst);
	
	return bndLst;
end

function calcShOutBndLeft( coord::Float64 )
	return (coord < 0 ? 1 : 0);
end

function calcShOutBndRight( coord::Float64 )
	return ( coord > 1 ? -1 : 0 );
end

const calcShOutBndFuncLst = [calcShOutBndLeft, calcShOutBndRight];

function genEllipseBndModExtendedLst( bndLst::Vector{<:AbstractVector{<:AbstractVector{Float64}}}, pt2dLst::AbstractVector{<:AbstractVector{Float64}} )
	# bndModXYLst = [ [ [ bndLst[iBnd][iXY] .+ calcShOutBndLst[iStartEnd](bndLst[iBnd][iXY][iStartEnd]) for iBnd = 1 : length(bndLst) ] for iStartEnd = 1 : 2 ] for iXY = 1 : 2 ];
	
	bndModXYLst = [ [ [ copy( bndLst[iBnd][iXY] ) for iBnd = 1 : length(bndLst) ] for iStartEnd = 1 : 2 ] for iXY = 1 : 2 ];
	
	bndExtendedXYLst = [ [ copy( bndLst[iBnd][iXY] ) for iBnd = 1 : length(bndLst) ] for iXY = 1 : 2 ];
	
	pt2dModXYLst = [ [ deepcopy(pt2dLst) for iStartEnd = 1 : 2 ] for iXY = 1 : 2 ];
	pt2dExtendedXYLst = [ deepcopy(pt2dLst) for iXY = 1 : 2 ];
	
	idExtendedXYLst = [ [1:length(bndLst);] for iXY = 1 : 2 ];
	
	sh = 0;
	for iXY = 1 : 2, iStartEnd = 1 : 2, iBnd = 1 : length(bndLst)
		sh = calcShOutBndFuncLst[iStartEnd](bndLst[iBnd][iXY][iStartEnd]);
		bndModXYLst[iXY][iStartEnd][iBnd] .+= sh;
		pt2dModXYLst[iXY][iStartEnd][iBnd][iXY] += sh;
		# if iStartEnd == 1 && 
		if sh != 0
			push!( bndExtendedXYLst[iXY], bndModXYLst[iXY][iStartEnd][iBnd] );
			push!( pt2dExtendedXYLst[iXY], pt2dModXYLst[iXY][iStartEnd][iBnd] );
			push!( idExtendedXYLst[iXY], iBnd );
		end
	end
	
	return bndModXYLst, bndExtendedXYLst, idExtendedXYLst, pt2dModXYLst, pt2dExtendedXYLst;
end

function sortBndLst( bndModLst::Vector{<:AbstractVector{<:AbstractVector{<:AbstractVector{Float64}}}}, bndExtendedLst::Vector{<:AbstractVector{<:AbstractVector{Float64}}} )
	idSortedBndModLst = [ [ sortperm( bndModLst[iXY][iStartEnd] ) for iStartEnd = 1 :2 ] for iXY = 1 : 2 ];
	idSortedBndExtendedLst = sortperm.( bndExtendedLst );
	
	return idSortedBndModLst, idSortedBndExtendedLst;
end

function genEllipseEqCoeffLst( rotMatInvLst::Vector{<:AbstractMatrix{Float64}} )
	eqCoeffLst = [ zeros(3) for ii = 1 : length(rotMatInvLst) ];
	
	for ii = 1 : length(rotMatInvLst)
		eqCoeffLst[ii][1] = rotMatInvLst[ii][1,1]^2 + rotMatInvLst[ii][2,1]^2;
		eqCoeffLst[ii][3] = rotMatInvLst[ii][1,2]^2 + rotMatInvLst[ii][2,2]^2;
		eqCoeffLst[ii][2] = 2 * ( rotMatInvLst[ii][1,1] * rotMatInvLst[ii][1,2] + rotMatInvLst[ii][2,1] * rotMatInvLst[ii][2,2] );
	end
	
	return eqCoeffLst;
end

function solveOtherXY( eqCoeffs::AbstractVector{Float64}, pt::AbstractVector{Float64}, xyVal::Real; xySolved = 'X' )
	if xySolved == 'X'
		eqCoeffs1Var = copy( eqCoeffs );
		xyVal -= pt[2];
		solSh = pt[1];
	elseif xySolved == 'Y'
		eqCoeffs1Var = reverse( eqCoeffs );
		xyVal -= pt[1];
		solSh = pt[2];
	end
	
	for iPow = 1 : 3
		eqCoeffs1Var[iPow] *= xyVal^(iPow-1);
	end
	eqCoeffs1Var[3] -= 1;
	# @infiltrate
	
	discr = sqrt( eqCoeffs1Var[2]^2 - 4 * eqCoeffs1Var[1] * eqCoeffs1Var[3] ) / (2 * eqCoeffs1Var[1]);
	b2a = - eqCoeffs1Var[2] / (2 * eqCoeffs1Var[1]);
	
	xySolLst = [ b2a - discr, b2a + discr ];
	xySolLst .= xySolLst .+ solSh;
	
	return xySolLst;
end	

function genZakArr( divNum::Int64, eqCoeffLst::Vector{<:AbstractVector{Float64}}, pt2dLst::AbstractVector{<:AbstractVector{Float64}}, pt2dModXYLst::AbstractVector{<:AbstractVector{<:AbstractVector{<:AbstractVector{Float64}}}}, pt2dExtendedXYLst::AbstractVector{<:AbstractVector{<:AbstractVector{Float64}}}, bndModLst::Vector{<:AbstractVector{<:AbstractVector{<:AbstractVector{Float64}}}}, bndExtendedLst::Vector{<:AbstractVector{<:AbstractVector{Float64}}}, idExtendedXYLst::Vector{<:AbstractVector{Int64}}, idSortedBndModLst::Vector{<:AbstractVector{<:AbstractVector{Int64}}}, idSortedBndExtendedLst::Vector{<:AbstractVector{Int64}} )
	zakArr = zeros(Bool, divNum, divNum);
	zakXLst = zeros(Bool, divNum);
	
	yVal = 0;
	iYPast0 = Utils.searchSortedFirstByVal( idSortedBndModLst[2][2], yVal ; by = x -> bndModLst[2][2][x][1] );
	iXStart = 1;
	xBndLst = zeros(0);
	# @infiltrate
	for iYBnd = 1 : iYPast0 - 1
		id = idSortedBndModLst[2][2][iYBnd];
		xSolLst = solveOtherXY( eqCoeffLst[id], pt2dModXYLst[2][2][id], yVal; xySolved = 'X' );
		for iEnd = 1 : 2
			push!(xBndLst, xSolLst[iEnd]);
		end
	end
	xBndLst .= mod.( xBndLst, 1 );
	sort!(xBndLst);
	
	iXStart = 1;
	zakVal = false;
	for iBnd = 1 : length(xBndLst)
		iXEnd = Int64(floor( xBndLst[iBnd] * divNum )) + 1;
		for iX = iXStart : iXEnd
			zakXLst[iX] = zakVal;
		end
		iXStart = iXEnd + 1;
		zakVal = !zakVal;
	end
	for iX = iXStart : divNum
		zakXLst[iX] = zakVal;
	end
	
	bndExtendedXLst = bndExtendedLst[1];
	idSortedBndExtendedXLst = idSortedBndExtendedLst[1];
	
	iBndX = 1;
	
	xLst = [0:1/divNum:(1-1/divNum);];
	
	iXBndQueue = BinaryHeap{Int64}(Base.By( ii -> bndExtendedXLst[ii][2] ));
	
	lnBndExt = length(bndExtendedLst[1]);
	
	yBndLst = zeros(0);
	
	iBnd = 1;
	for iX  = 1 : divNum
		while iBnd <= lnBndExt && bndExtendedXLst[idSortedBndExtendedXLst[iBnd]][1] < xLst[iX]
			push!( iXBndQueue, idSortedBndExtendedXLst[iBnd] );
			iBnd += 1;
		end
		while !isempty(iXBndQueue) && bndExtendedXLst[first(iXBndQueue)][2] < xLst[iX]
			pop!(iXBndQueue);
		end
		
		empty!(yBndLst);
		for iXBnd in iXBndQueue
			idNoExtend = idExtendedXYLst[1][iXBnd];
			
			yBndSolved = solveOtherXY( eqCoeffLst[idNoExtend], pt2dExtendedXYLst[1][iXBnd], xLst[iX]; xySolved = 'Y' );
			for iEnd = 1 : 2
				push!( yBndLst, yBndSolved[iEnd] );
			end
			yBndLst .= mod.( yBndLst, 1 );
			sort!(yBndLst);
		end
		
		iYStart = 1;
		zakVal = false;
		for iYBnd = 1 : length(yBndLst)
			iYBndEnd = Int64( floor( yBndLst[iYBnd] * divNum ) ) + 1;
			for iY = iYStart : iYBndEnd
				zakArr[iX,iY] = xor( zakVal, zakXLst[iX] );
			end
			iYStart = iYBndEnd + 1;
			zakVal = !zakVal;
		end
		for iY = iYStart : divNum
			zakArr[iX,iY] = xor( zakVal, zakXLst[iX] );
		end
	end
	
	return zakArr;
end
