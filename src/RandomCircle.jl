module RandomCircle

using StaticArrays
using OffsetArrays
using Random
using Utils
using LinearAlgebra
using DataStructures
using FFTW

using Infiltrator

const nDim3 = 3;
const nDimQuat = nDim3 + 1;

const uVecLst3 = Utils.genStaticIdentityMat( nDim3 );
const idQuatRng = 0:nDim3;

function boolToIntPosNeg( valBool::Bool )
	return valBool ? 1 : -1;
end

struct QuartSolHelper
	solLst::MVector{2,Float64};
	coeffs1Var::MVector{3,Float64};
	
	function QuartSolHelper()
		solLst = @MVector zeros(2);
		coeffs1Var = @MVector zeros(3);
		
		new(solLst, coeffs1Var);
	end
end

struct RandCircData{N_circ}
	rCirc::Base.RefValue{Float64};
	nCirc::Int64;

	ptLst::MVector{N_circ,<:AbstractVector{Float64}};
	pt2dLst::MVector{N_circ,<:AbstractVector{Float64}};
	quatLst::MVector{N_circ,<:AbstractVector{Float64}};
	quatNormLst::MVector{N_circ,Float64};
	quatNormSqLst::MVector{N_circ,Float64};
	
	rotMatLst::MVector{N_circ,<:AbstractMatrix{Float64}};
	rotMat2dLst::MVector{N_circ,<:AbstractMatrix{Float64}};
	rotMat2dInvLst::MVector{N_circ,<:AbstractMatrix{Float64}};
	
	sampleThetaLst::Vector{Float64};
	sampleCosSinLst::Matrix{Float64};
	
	sampleCircLst::Vector{Matrix{Float64}};
	
	eqCoeffLst::MVector{N_circ,MVector{3,Float64}};
	solHelper::QuartSolHelper;
	
	areaLst::MVector{N_circ, Float64};
	widthXYLst::MVector{N_circ, MVector{2,Float64}};
	
	bndLst::MVector{N_circ,<:MVector{2,<:AbstractVector{Float64}}};
	bndModLst::MMatrix{2,2,<:MVector{N_circ,<:AbstractVector{Float64}},4};
	bndExtendedLst::MVector{2,<:AbstractVector{<:AbstractVector{Float64}}};	
	pt2dModLst::MMatrix{2,2,MVector{N_circ,MVector{2,Float64}},4};
	pt2dExtendedLst::MVector{2,Vector{MVector{2,Float64}}};
	bndExtendedXLst::AbstractVector{MVector{2,Float64}};
	
	idExtendedLst::MVector{2,Vector{Int64}};
	idSortedBndModLst::MMatrix{2,2,MVector{N_circ,Int64},4};
	idSortedBndExtendedLst::MVector{2,Vector{Int64}};
	idSortedBndExtendedXLst::Vector{Int64};
	
	circHeap::BinaryHeap{Int64};
	
	xLst::Vector{Float64};
	xMidLst::Vector{Float64};
	
	zakArrRef::Base.RefValue{Matrix{Bool}};
	zakXLst::Vector{Bool};
	
	bndXOnY0Lst::Vector{Float64};
	bndYOnXValLst::Vector{Float64};
	
	zakArrFloatRef::Base.RefValue{Matrix{Float64}};
	zakCorrArrCmplxRef::Base.RefValue{Matrix{ComplexF64}};
	zakCorrArrRef::Base.RefValue{Matrix{Float64}};
	
	function RandCircData( nCirc::Int64, rCirc::Real; nSample = 1000, divNum = 128 )
		ptLst = Utils.@MVectorCompr [ @MVector( zeros(3) ) for ii = 1 : nCirc ];
		pt2dLst = Utils.@MVectorCompr [ @MVector zeros(2) for ii = 1 : nCirc ];
		quatLst = Utils.@MVectorCompr [ OffsetArray( @MVector( zeros(4) ), idQuatRng ) for ii = 1 : nCirc ];
		quatNormLst = @MVector zeros(nCirc);
		quatNormSqLst = @MVector zeros(nCirc);
		
		rotMatLst = Utils.@MVectorCompr [ @MMatrix zeros(3,3) for ii = 1 : nCirc ];
		rotMat2dLst = Utils.@MVectorCompr [ @MMatrix zeros(2,2) for ii = 1 : nCirc ];
		rotMat2dInvLst = Utils.@MVectorCompr [ @MMatrix zeros(2,2) for ii = 1 : nCirc ];
		
		sampleThetaLst = zeros(nSample);
		sampleCosSinLst = zeros(2,nSample);
		sampleCircLst = [ zeros( 2, nSample ) for iCirc = 1 : nCirc ];
		
		eqCoeffLst = Utils.@MVectorCompr [ @MVector zeros(3) for ii = 1 : nCirc ];
		quartHelper = QuartSolHelper();
		
		areaLst = @MVector zeros(nCirc);
		widthXYLst = Utils.@MVectorCompr [ @MVector zeros(2) for ii = 1 : nCirc ];
		
		bndLst = Utils.@MVectorCompr [ Utils.@MVectorCompr [ @MVector zeros(2) for iXY = 1 : 2 ] for ii = 1 : nCirc ];
		bndModLst = Utils.@MMatrixCompr [ Utils.@MVectorCompr [ @MVector zeros(2) for ii = 1 : nCirc ] for iXY = 1 : 2, iStartEnd = 1 : 2 ];
		bndExtendedLst = Utils.@MVectorCompr [ [ @MVector zeros(2) for ii = 1 : nCirc ] for iXY = 1 : 2 ];
 		pt2dModLst = Utils.@MMatrixCompr [ Utils.@MVectorCompr [ @MVector zeros(2) for iCirc = 1 : nCirc ] for iXY = 1:2, iStartEnd = 1:2];
		pt2dExtendedLst = Utils.@MVectorCompr [ [ @MVector zeros(2) for iCirc = 1 : nCirc ] for iXY = 1 : 2 ];
		bndExtendedXLst = bndExtendedLst[1];
		
		idExtendedLst = [ [1:nCirc;] for iXY = 1 : 2 ];
		idSortedBndModLst = Utils.@MMatrixCompr [ @MVector zeros(Int64, nCirc) for iXY = 1 : 2, iStartEnd = 1 : 2 ];
		idSortedBndExtendedLst = Utils.@MVectorCompr [ zeros(Int64, nCirc) for iXY = 1 : 2 ];
		idSortedBndExtendedXLst = idSortedBndExtendedLst[1];
		
		solHelper = QuartSolHelper();
		
		circHeap = BinaryHeap{Int64}( Base.By( ii -> bndExtendedLst[1][ii][2] ) );
		
		xLst = [0.0:divNum-1;];
		xLst .= xLst ./ divNum;
		xMidLst = copy(xLst);
		xMidLst .+= 0.5/divNum;
		
		zakArrRef = Ref( zeros(Bool, divNum, divNum) );
		zakXLst = zeros(Bool, divNum);
		
		bndXOnY0Lst = zeros(0);
		bndYOnXValLst = zeros(0);
		
		zakArrFloatRef = Ref( similar( zakArrRef[], Float64 ) );
		zakCorrArrRef = Ref( zeros( divNum, divNum ) );
		zakCorrArrCmplxRef = Ref( zeros( ComplexF64, divNum, divNum ) );
		
		data = new{nCirc}( Ref(Float64(rCirc)), nCirc, ptLst, pt2dLst, quatLst, quatNormLst, quatNormSqLst, rotMatLst, rotMat2dLst, rotMat2dInvLst, sampleThetaLst, sampleCosSinLst, sampleCircLst, eqCoeffLst, solHelper, areaLst, widthXYLst, bndLst, bndModLst, bndExtendedLst, pt2dModLst, pt2dExtendedLst, bndExtendedXLst, idExtendedLst, idSortedBndModLst, idSortedBndExtendedLst, idSortedBndExtendedXLst, circHeap, xLst, xMidLst, zakArrRef, zakXLst, bndXOnY0Lst, bndYOnXValLst, zakArrFloatRef, zakCorrArrCmplxRef, zakCorrArrRef );
		
		refreshSampleBaseLst!( data );
		
		
		return data;
	end
end

function getRCirc( data::RandCircData )
	return data.rCirc[];
end

function setRCirc!( data::RandCircData, rCirc::Float64 )
	if rCirc != data.rCirc[]
		
		scaleRotMat!( data, rCirc / data.rCirc[] );
		setRCircNoRotUpdate( data, rCirc );
		refreshRotMat2dInv!( data );
	end
end

function setRCircNoRotUpdate!( data::RandCircData, rCirc::Float64 )
	data.rCirc[] = rCirc;
end

function setDivNum!( data::RandCircData, divNum::Int64 )
	resize!( data.zakXLst, divNum );
	data.zakArrRef[] = zeros( Bool, divNum, divNum );
	data.zakCorrArrRef[] = similar( data.zakArrRef[], Float64 );
	data.zakCorrArrCmplxRef[] = similar( data.zakArrRef[], ComplexF64 );
	
	
	
	GC.gc()
end

function getZakArr( data::RandCircData )
	return data.zakArrRef[];
end

function setZakArr!( data::RandCircData, zakArr::AbstractMatrix )
	data.zakArrRef[] .= zakArr;
end

function getZakCorr( data::RandCircData )
	return data.zakCorrArrRef[];
end

function getZakCorrCmplx( data::RandCircData )
	return data.zakCorrArrCmplxRef[];
end

function refreshXLst!( data )
	;
end

function refreshPtLst!( randCircData::RandCircData )
	rand!.( randCircData.ptLst );
	for ii = 1 : randCircData.nCirc, iXY = 1 : 2
		randCircData.pt2dLst[ii][iXY] = randCircData.ptLst[ii][iXY]
	end
end

function refreshQuatLst!( randCircData::RandCircData )
	randn!.( randCircData.quatLst );
	randCircData.quatNormLst .= norm.(randCircData.quatLst);
	( (x,y) -> x .= x ./ y ).( randCircData.quatLst, randCircData.quatNormLst );
	randCircData.quatNormSqLst .= 0;
	for iMat = 1 : randCircData.nCirc
		for ii = 1 : 3
			randCircData.quatNormSqLst[iMat] += randCircData.quatLst[iMat][ii]^2;
		end
	end
	randCircData.quatNormLst .= sqrt.( randCircData.quatNormSqLst );
end

function refreshRotMatFull!( data::RandCircData )
	refreshRotMat!( data );
	scaleRotMat!( data );
	refreshRotMat2dInv!( data );
end

function refreshRotMat!( randCircData::RandCircData )
	for iMat = 1 : randCircData.nCirc
		for ii = 1 : 3
			randCircData.rotMatLst[iMat][ii,ii] = 1 - 2 * randCircData.quatNormSqLst[iMat] + 2 * randCircData.quatLst[iMat][ii]^2;
		end
		for ii = 1 : 3, jj = 1 : 3
			if ii == jj 
				continue;
			end
			randCircData.rotMatLst[iMat][ii,jj] = 2 * randCircData.quatLst[iMat][ii] * randCircData.quatLst[iMat][jj] - 2 * randCircData.quatLst[iMat][0] * randCircData.quatLst[iMat][leviCivita3rdIdMat[ii,jj]] * leviCivita3rdSgnMat[ii,jj];
		end
	end
end

function scaleRotMat!( data::RandCircData )
	rCirc = getRCirc( data );
	# ( x -> x .*= rCirc ).( data.rotMatLst );
	scaleRotMat!( data, rCirc );
end

function scaleRotMat!( data::RandCircData, scale::Float64 )
	( x -> x .*= scale ).( data.rotMatLst );
end

function refreshRotMat2dInv!( randCircData::RandCircData )
	for iMat = 1 : randCircData.nCirc
		for ii = 1 : 2, jj = 1 : 2
			randCircData.rotMat2dLst[iMat][ii,jj] = randCircData.rotMatLst[iMat][ii,jj];
		end
		Utils.inv!( randCircData.rotMat2dInvLst[iMat], randCircData.rotMat2dLst[iMat] );
	end
end

function refreshEqCoeff!( data::RandCircData )
	for iCirc = 1 : data.nCirc
		matInv = data.rotMat2dInvLst[iCirc];
		data.eqCoeffLst[iCirc][1] = matInv[1,1]^2 + matInv[2,1]^2;
		data.eqCoeffLst[iCirc][3] = matInv[1,2]^2 + matInv[2,2]^2;
		data.eqCoeffLst[iCirc][2] = 2 * ( matInv[1,1] * matInv[1,2] + matInv[2,1] * matInv[2,2] );
	end
end

function refreshBndLst!( data::RandCircData )
	data.areaLst .= abs.( det.( data.rotMat2dLst ) );
	
	for iCirc = 1 : data.nCirc
		data.widthXYLst[iCirc][1] = data.eqCoeffLst[iCirc][3];
		data.widthXYLst[iCirc][2] = data.eqCoeffLst[iCirc][1];
		data.widthXYLst[iCirc] .= data.areaLst[iCirc] .* sqrt.( data.widthXYLst[iCirc] );
		for iXY = 1 : 2, iStartEnd = 1 : 2
			data.bndLst[iCirc][iXY][iStartEnd] = data.pt2dLst[iCirc][iXY] + ( iStartEnd == 1 ? -data.widthXYLst[iCirc][iXY] : data.widthXYLst[iCirc][iXY] );
		end
	end
end

function refreshBndPtsModExtendLst!( data::RandCircData )
	resize!.( data.bndExtendedLst, data.nCirc );
	resize!.( data.pt2dExtendedLst, data.nCirc );
	resize!.( data.idExtendedLst, data.nCirc );
	
	for iCirc = 1 : data.nCirc
		for iXY = 1 : 2
			data.bndExtendedLst[iXY][iCirc] .= data.bndLst[iCirc][iXY];
			data.pt2dExtendedLst[iXY][iCirc] .= data.pt2dLst[iCirc];
		end
		for iXY = 1 : 2, iStartEnd = 1 : 2
			data.bndModLst[iXY,iStartEnd][iCirc] .= data.bndLst[iCirc][iXY];
			data.pt2dModLst[iXY,iStartEnd][iCirc] .= data.pt2dLst[iCirc];
			sh = 0; 
			if iStartEnd == 1 && data.bndLst[iCirc][iXY][iStartEnd] < 0
				sh = 1;
			elseif iStartEnd == 2 && data.bndLst[iCirc][iXY][iStartEnd] > 1
				sh = -1;
			end
			data.bndModLst[iXY,iStartEnd][iCirc] .+= sh;
			data.pt2dModLst[iXY,iStartEnd][iCirc][iXY] += sh;
			if sh != 0
				push!( data.bndExtendedLst[iXY], data.bndModLst[iXY,iStartEnd][iCirc] );
				push!( data.pt2dExtendedLst[iXY], data.pt2dModLst[iXY,iStartEnd][iCirc] );
				push!( data.idExtendedLst[iXY], iCirc );
			end
		end
	end
end

function backupRotMat!( rotMatLst, rotMat2dLst, rotMat2dInvLst, data::RandCircData )
	for iCirc = 1 : data.nCirc
		rotMatLst[iCirc] .= data.rotMatLst[iCirc];
		rotMat2dLst[iCirc] .= data.rotMat2dLst[iCirc];
		rotMat2dInvLst[iCirc] .= data.rotMat2dInvLst[iCirc];
	end
end

function restoreRotMatOnly!( data, rotMatLst, rotMat2dLst, rotMat2dInvLst )
	for iCirc = 1 : data.nCirc
		data.rotMatLst[iCirc] .= rotMatLst[iCirc];
		data.rotMat2dLst[iCirc] .= rotMat2dLst[iCirc];
		data.rotMat2dInvLst[iCirc] .= rotMat2dInvLst[iCirc];
	end
end

function restoreRotMat!( data, rotMatLst, rotMat2dLst, rotMat2dInvLst )
	restoreRotMatOnly!( data, rotMatLst, rotMat2dLst, rotMat2dInvLst );
	
	refreshEqCoeff!( data );
	refreshBndLstFull!( data );
end

function refreshSortBndLst!( data::RandCircData )
	((x,y) -> resize!(x, length(y))).( data.idSortedBndExtendedLst, data.bndExtendedLst );
	
	sortperm!.( data.idSortedBndModLst, data.bndModLst );
	sortperm!.( data.idSortedBndExtendedLst, data.bndExtendedLst );
end

function refreshBndLstFull!( data::RandCircData )
	refreshBndLst!( data );
	refreshBndPtsModExtendLst!( data );
	refreshSortBndLst!( data );
end

function resizeSampleLn!( data::RandCircData, nSample::Int64 )
	if nSample != length(data.sampleThetaLst)
		resize!( data.sampleThetaLst, nSample );
		resize!( data.sampleCosSinLst, 2, nSample );
		refreshSampleBaseLst!( data );
	end
end

function refreshSampleBaseLst!( data::RandCircData )
	nSample = length( data.sampleThetaLst );
	for ii = 1 : nSample
		data.sampleThetaLst[ii] = (ii-1)/nSample * 2*pi;
	end
	data.sampleCosSinLst[1,:] .= cos.( data.sampleThetaLst );
	data.sampleCosSinLst[2,:] .= sin.( data.sampleThetaLst );
end

function calcSampleLst!( data::RandCircData, nSample::Int64 )
	resizeSampleLn!( data, nSample );
	
	calcSampleLst!( data::RandCircData );
end

function calcSampleLst!( data::RandCircData )
	for iCirc = 1 : data.nCirc
		mul!( data.sampleCircLst[iCirc], data.rotMat2dLst[iCirc], data.sampleCosSinLst );
		data.sampleCircLst[iCirc] .+= data.pt2dLst[iCirc];
	end
end

function backupSampleLst!( sampleCircLst, data::RandCircData )
	( (x, y) -> x .= y ).( sampleCircLst, data.sampleCircLst );
end

function solveQuartEqXY!( solHelper::QuartSolHelper, coeffs::AbstractVector{Float64}, xyVal::Real; xySolved = 'X' )
	coeffs1Var = solHelper.coeffs1Var;
	coeffs1Var .= coeffs;
	if xySolved == 'Y'
		reverse!( coeffs1Var );
	end
	
	for ii = 1 : 3
		coeffs1Var[ii] *= xyVal^(ii-1);
	end
	
	coeffs1Var[3] -= 1;
	
	discr = sqrt( coeffs1Var[2]^2 - 4 * coeffs1Var[1] * coeffs1Var[3] ) / (2 * coeffs1Var[1]);
	b2a = - coeffs1Var[2] / (2 * coeffs1Var[1]);
	
	solHelper.solLst .= b2a;
	solHelper.solLst[1] -= discr;
	solHelper.solLst[2] += discr;
end

function solveQuartEqXYSh!( solHelper::QuartSolHelper, coeffs::AbstractVector{Float64}, pt2d::AbstractVector{Float64}, xyVal::Real; xySolved = 'X' )
	if xySolved == 'X'
		shInput = pt2d[2];
		shOutput = pt2d[1];
	else
		shInput = pt2d[1];
		shOutput = pt2d[2];
	end
	
	xyVal = xyVal - shInput;
	solveQuartEqXY!( solHelper, coeffs, xyVal; xySolved = xySolved );
	solHelper.solLst .+= shOutput;
end

function xValToId( data::RandCircData, xVal::Float64 )
	divNum = length( data.xLst );
	id = Int64( floor( xVal * divNum ) ) + 1;
	
	return id;
end

function calcZakArr!( data::RandCircData )
	zakArr = getZakArr( data );

	divNum = length(data.zakXLst);
	yVal = 0;
	iYUpTo0 = Utils.searchSortedFirstByVal( data.idSortedBndModLst[2,2], yVal; by = id -> data.bndModLst[2,2][id][1] ) - 1;
	
	if iYUpTo0 > 0
		resize!( data.bndXOnY0Lst, 2*iYUpTo0 );
	end
	
	i1 = 1;
	i2 = 1;
	for iSorted = 1 : iYUpTo0
		i2 = 2*iSorted;
		i1 = i2 - 1;
		iCirc = data.idSortedBndModLst[2,2][iSorted];
		solveQuartEqXYSh!( data.solHelper, data.eqCoeffLst[iCirc], data.pt2dModLst[2,2][iCirc], yVal; xySolved = 'X' );
		data.solHelper.solLst .= mod.( data.solHelper.solLst, 1 );
		data.bndXOnY0Lst[i1] = data.solHelper.solLst[1];
		data.bndXOnY0Lst[i2] = data.solHelper.solLst[2];
	end
	sort!(data.bndXOnY0Lst);
	
	data.zakXLst .= false;
	iXStart = 1;
	zakVal = false;
	for iXBnd = 1 : length(data.bndXOnY0Lst)
		iXEnd = xValToId( data, data.bndXOnY0Lst[iXBnd] );
		for iX = iXStart : iXEnd
			data.zakXLst[iX] = zakVal;
		end
		iXStart = iXEnd + 1;
		zakVal = !zakVal;
	end
	for iX = iXStart : divNum
		data.zakXLst[iX] = zakVal;
	end
	
	bndExtendedXLst = data.bndExtendedLst[1];
	empty!( data.circHeap );
	
	iSorted = 1;
	for iX = 1 : divNum
		xVal = data.xLst[iX];
		while iSorted <= length(data.bndExtendedXLst) && data.bndExtendedXLst[data.idSortedBndExtendedXLst[iSorted]][1] < xVal
			iCirc = data.idSortedBndExtendedXLst[iSorted];
			push!( data.circHeap, iCirc );
			iSorted += 1;
		end
		
		while !isempty(data.circHeap) && data.bndExtendedXLst[first(data.circHeap)][2] < xVal
			pop!(data.circHeap);
		end
		
		resize!( data.bndYOnXValLst, 2 * length( data.circHeap ) );
		i1 = 1;
		i2 = 2;
		for ii = 1 : length(data.circHeap)
			i2 = 2*ii;
			i1 = i2 - 1;
			iCirc = data.circHeap[ii];
			iCircNoExtend = data.idExtendedLst[1][iCirc];
			solveQuartEqXYSh!( data.solHelper, data.eqCoeffLst[iCircNoExtend], data.pt2dExtendedLst[1][iCirc], xVal; xySolved = 'Y' );
			data.solHelper.solLst .= mod.( data.solHelper.solLst, 1 );
			
			data.bndYOnXValLst[i1] = data.solHelper.solLst[1];
			data.bndYOnXValLst[i2] = data.solHelper.solLst[2];
		end
		sort!( data.bndYOnXValLst );
		
		iYStart = 1;
		zakVal = false;
		for iBnd = 1 : length(data.bndYOnXValLst)
			iYEnd = xValToId( data, data.bndYOnXValLst[iBnd] );
			for iY = iYStart : iYEnd
				zakArr[iX,iY] = xor( zakVal, data.zakXLst[iX] );
			end
			iYStart = iYEnd + 1;
			zakVal = !zakVal;
		end
		for iY = iYStart : divNum
			zakArr[iX,iY] = xor( zakVal, data.zakXLst[iX] );
		end
	end
end

function calcZakCorr!( data::RandCircData )
	zakCorrArrCmplx = getZakCorrCmplx( data );
	zakCorrArr = getZakCorr( data );
	zakArr = getZakArr( data );
	
	# dArea = 1 / length(data.zakXLst)^2;
	# zakCorrArrCmplx .= (x -> x ? 1 : 0).( zakArr );
	# zakCorrArrCmplx .= boolToIntPosNeg.( zakArr );
	# fft!( zakCorrArrCmplx );
	# zakCorrArrCmplx .= abs.( zakCorrArrCmplx ).^2 .* dArea;
	# ifft!( zakCorrArrCmplx );
	# calcCorrCmplx!( zakCorrArrCmplx );
	# zakCorrArr .= real.( zakCorrArrCmplx );
	
	calcZakCorr!( zakCorrArr, zakCorrArrCmplx, zakArr );
end

function calcZakCorr!( corrArr::AbstractArray, corrArrCmplx::AbstractArray, zakArr::AbstractArray )
	corrArrCmplx .= boolToIntPosNeg.(zakArr);
	calcCorrCmplx!( corrArrCmplx );
	corrArr .= real.( corrArrCmplx );
end

function calcCorrCmplx!( arrCmplx )
	# zakCorrArrCmplx .= boolToIntPosNeg.( zakArr );
	dArea = 1 / size(arrCmplx,1)^2;
	fft!( arrCmplx );
	arrCmplx .= abs.( arrCmplx ).^2 .* dArea;
	ifft!( arrCmplx );
end

function backupZakArrCorr!( zakArr, zakCorrArr, data::RandCircData )
	zakArr .= getZakArr( data );
	zakCorrArr .= getZakCorr( data );
end

include("randomCircleFunc_nonStruct.jl");

end # endmodule
