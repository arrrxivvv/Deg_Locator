
const modelExp( x, p ) = p[1] .* exp.( - p[2] .* x ) .+ p[3];
p0FitModelExp = [1, 0.5, 0];

const nDim2 = 2;

struct ZakCorrArrTmpData
	tmpArrRef::Base.RefValue{<:Array{Float64}};
	meanArrRef::Base.RefValue{<:Array{Float64}};
	mean1dLstRef::Base.RefValue{<:AbstractVector{Float64}};
	mean1dHalfLstRef::Base.RefValue{<:AbstractVector{Float64}};
end

getTmpArr( data::ZakCorrArrTmpData ) = data.tmpArrRef[];
getMeanArr( data::ZakCorrArrTmpData ) = data.meanArrRef[];
getMean1dLst( data::ZakCorrArrTmpData ) = data.mean1dLstRef[];
getMean1dHalfLst( data::ZakCorrArrTmpData ) = data.mean1dHalfLstRef[];

function ZakCorrArrTmpData()
	tmpArrRef = Ref( Array{Float64}(undef,0,0,0) );
	meanArrRef = Ref( Array{Float64}(undef,0,0) );
	mean1dLstRef = Base.RefValue{SubArray{Float64,1}}();
	mean1dHalfLstRef = Base.RefValue{SubArray{Float64,1}}();
	
	ZakCorrArrTmpData( tmpArrRef, meanArrRef, mean1dLstRef, mean1dHalfLstRef );
end

function setDivNum!( data::ZakCorrArrTmpData, divNum::Int64, divNumHalf::Int64, itNum::Int64 )
	data.tmpArrRef[] = zeros( divNum, divNum, itNum );
	data.meanArrRef[] = zeros( divNum, divNum );
	data.mean1dLstRef[] = @view data.meanArrRef[][:,1];
	data.mean1dHalfLstRef[] = @view data.mean1dLstRef[][1:divNumHalf];
end

function meanZakCorr!( data::ZakCorrArrTmpData );
	mean!( getMeanArr( data ), getTmpArr( data ) );
end

struct RunRandCircData
	nCircLst::Vector{Int64};
	rCircLst::Vector{Float64};
	itNum1PassRef::Base.RefValue{Int64};
	itNumFineRef::Base.RefValue{Int64};
	itNumNowRef::Base.RefValue{Int64};
	divNum1Pass::Int64;
	
	isStoreCorrFull::Bool
	
	randCircDataLst::Vector{RandCircData};
	
	rotMatBackupLst::Array{Vector{Vector{Matrix{Float64}}},2};
	rotMat2dBackupLst::Array{Vector{Vector{Matrix{Float64}}},2};
	rotMat2dInvBackupLst::Array{Vector{Vector{Matrix{Float64}}},2};
	
	sampleLstLst::Array{Vector{Vector{Matrix{Float64}}},2};
	
	# zakArrLst1Pass::Array{Bool};
	zakCorr1dLst1Pass::Array{Float64};
	
	zakCorr1dMeanLst1PassSingleton::Array{Float64};
	zakCorr1dMeanLst1Pass::Array{Float64};
	zakCorrMean1dHalf1Pass::Array{<:AbstractVector{Float64}};
	
	fitModel1PassLst::Array{LsqFit.LsqFitResult};
	xLst1Pass::Vector{Float64};
	xLstHalf1Pass::AbstractVector{Float64};
	corrLenLst1Pass::Array{Float64};
	corrLenInvLst1Pass::Array{Float64};
	
	# zakCorrTmpLstRef::Base.RefValue{<:Array{Float64}};
	# zakCorrMeanArrRef::Base.RefValue{<:Array{Float64}};
	# zakCorrMean1dRef::Base.RefValue{<:AbstractVector{Float64}};
	# zakCorrMean1dHalfRef::Base.RefValue{<:AbstractVector{Float64}};
	zakTmpData::ZakCorrArrTmpData;
	zakCorrMeanFullStoreLst::Array{Array{Float64}};
	zakCorr1dStoreLst::Array{Array{Float64}};
	zakCorrMean1dStoreLst::Array{Array{Float64}};
	zakArrAvgLst::Array{Float64};
	zakArrAvgMeanLstSingleton::Array{Float64};
	zakArrAvgMeanLst::Array{Float64};
	xLstFineLst::Array{Vector{Float64}};
	xHalfLstFineLst::Array{AbstractVector{Float64}};
	corrLenFineLst::Array{Float64};
	expScaleFineLst::Array{Float64};
	expShFineLst::Array{Float64};
	
	divNumNxtLst::Matrix{Int64};
	divNumNxtHalfLst::Matrix{Int64};
	
	numCornerLst::Array{Int64};
	idCornerLstLst::Array{Vector{MVector{nDim2,Float64}}};
	numCornerMeanSingletonLst::Array{Float64};
	numCornerMeanLst::AbstractArray{Float64};
end

function RunRandCircData( nCircLst::Vector{Int64}, rCircLst::Vector{Float64}, itNum1Pass::Int64, itNumFine::Int64, nSample::Int64, divNum1Pass::Int64; isStoreCorrFull::Bool = false )
	divNumHalf1Pass = div( divNum1Pass, 2 );
	randCircDataLst = RandCircData.(nCircLst);
	
	lnNCirc, lnRCirc = length.( (nCircLst, rCircLst) );
	
	rotMatBackupLst = [ [ [ zeros(3,3) for iCirc = 1 : nCircLst[iN] ] for it = 1 : itNum1Pass ] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
	rotMat2dBackupLst = [ [ [ zeros(2,2) for iCirc = 1 : nCircLst[iN] ] for it = 1 : itNum1Pass ] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
	rotMat2dInvBackupLst = deepcopy( rotMat2dBackupLst );
	
	sampleLstLst = [ [ [ zeros(2, nSample) for iCirc = 1 : nCircLst[iNCirc] ] for it = 1 : itNum1Pass ] for iRCirc = 1 : lnRCirc, iNCirc = 1 : lnNCirc ];
	
	# zakArrLst1Pass = zeros(Bool, divNum1Pass, divNum1Pass, itNum1Pass, lnRCirc, lnNCirc);
	zakCorr1dLst1Pass = zeros( Float64, divNum1Pass, itNum1Pass, lnRCirc, lnNCirc );
	
	zakCorr1dMeanLst1PassSingleton = zeros( Float64, divNum1Pass, 1, lnRCirc, lnNCirc );
	zakCorr1dMeanLst1Pass = dropdims( zakCorr1dMeanLst1PassSingleton; dims = 2 );
	zakCorrMean1dHalf1Pass = [ @view zakCorr1dMeanLst1Pass[1:divNumHalf1Pass,iR,iN] for iR = 1 : lnRCirc, iN = 1 : lnNCirc ];
	
	fitModel1PassLst = Array{LsqFit.LsqFitResult}(undef, lnRCirc, lnNCirc);
	xLst1Pass = [0.0:divNum1Pass-1;];
	xLst1Pass .= xLst1Pass ./ divNum1Pass;
	divNum1PassHalf = div( divNum1Pass, 2 );
	xLstHalf1Pass = @view( xLst1Pass[1:divNum1PassHalf] );
	corrLenLst1Pass = zeros( lnRCirc, lnNCirc );
	corrLenInvLst1Pass = similar( corrLenLst1Pass );
	
	# zakCorrTmpLstRef = Ref( Array{Float64}(undef, 0,0,0) );
	# zakCorrMeanArrRef = Ref( Array{Float64}(undef, 0,0) );
	# zakCorrMean1dRef = Base.RefValue{SubArray{Float64,1}}();
	# zakCorrMean1dHalfRef = Base.RefValue{SubArray{Float64,1}}();
	zakTmpData = ZakCorrArrTmpData();
	zakCorrMeanFullStoreLst = Array{Array{Float64}}(undef, lnRCirc, lnNCirc);
	zakCorr1dStoreLst = [ Array{Float64}(undef,0) for iR = 1 : lnRCirc, iR = 1 : lnNCirc ]
	zakCorrMean1dStoreLst = similar( zakCorr1dStoreLst );
	zakArrAvgLst = zeros( itNumFine, lnRCirc, lnNCirc );
	zakArrAvgMeanLstSingleton = zeros( 1, lnRCirc, lnNCirc );
	zakArrAvgMeanLst = dropdims( zakArrAvgMeanLstSingleton; dims = 1 );
	xLstFineLst = Array{Vector{Float64}}(undef, lnRCirc, lnNCirc);
	xHalfLstFineLst = Array{AbstractVector{Float64}}(undef, lnRCirc, lnNCirc);
	corrLenFineLst = zeros( lnRCirc, lnNCirc );
	expScaleFineLst = similar( corrLenFineLst );
	expShFineLst = similar( corrLenFineLst );
	
	divNumNxtLst = ones( Int64, lnRCirc, lnNCirc );
	divNumNxtHalfLst = similar(divNumNxtLst);
	
	numCornerLst = zeros( Int64, itNumFine, lnRCirc, lnNCirc );
	idCornerLstLst = Array{Array{Vector{MVector{nDim2,Float64}}}}( undef, itNumFine, lnRCirc, lnNCirc );
	numCornerMeanSingletonLst = zeros( Float64, 1, lnRCirc, lnNCirc );
	numCornerMeanLst = dropdims( numCornerMeanSingletonLst; dims = 1 );
	
	return RunRandCircData( nCircLst, rCircLst, Ref(itNum1Pass), Ref(itNumFine), Ref(itNum1Pass), divNum1Pass,  isStoreCorrFull, randCircDataLst, rotMatBackupLst, rotMat2dBackupLst, rotMat2dInvBackupLst, sampleLstLst, zakCorr1dLst1Pass, zakCorr1dMeanLst1PassSingleton, zakCorr1dMeanLst1Pass, zakCorrMean1dHalf1Pass, fitModel1PassLst, xLst1Pass, xLstHalf1Pass, corrLenLst1Pass, corrLenInvLst1Pass, zakTmpData, zakCorrMeanFullStoreLst, zakCorr1dStoreLst, zakCorrMean1dStoreLst, zakArrAvgLst, zakArrAvgMeanLstSingleton, zakArrAvgMeanLst, xLstFineLst, xHalfLstFineLst, corrLenFineLst, expScaleFineLst, expShFineLst, divNumNxtLst, divNumNxtHalfLst, numCornerLst, idCornerLstLst, numCornerMeanSingletonLst, numCornerMeanLst );
end

getItNum1Pass( runData::RunRandCircData ) = runData.itNum1PassRef[];
getItNumFine( runData::RunRandCircData ) = runData.itNumFineRef[];
getItNum( runData::RunRandCircData ) = runData.itNumNowRef[];
getLnNCirc( runData::RunRandCircData ) = length( runData.nCircLst );
getLnRCirc( runData::RunRandCircData ) = length( runData.rCircLst );
getDivNum1Pass( runData::RunRandCircData ) = runData.divNum1Pass;
getZakCorrTmpLst( runData::RunRandCircData ) = runData.zakCorrTmpLstRef[];
getZakCorrMeanArr( runData::RunRandCircData ) = runData.zakCorrMeanArrRef[];
getZakCorrMean1d( runData::RunRandCircData ) = runData.zakCorrMean1dRef[];
getZakCorrMean1dHalf( runData::RunRandCircData ) = runData.zakCorrMean1dHalfRef[];

function setDivNum!( runData::RunRandCircData, divNum::Int64, divNumHalf::Int64 )
	runData.zakCorrTmpLstRef[] = zeros( divNum, divNum, getItNumFine(runData) );
	runData.zakCorrMeanArrRef[] = zeros( divNum, divNum );
	runData.zakCorrMean1dRef[] = @view getZakCorrMeanArr( runData )[:,1];
	runData.zakCorrMean1dHalfRef[] = @view runData.zakCorrMean1dRef[][1:divNumHalf];
	GC.gc();
end

function setItNum!( runData::RunRandCircData, itNum::Int64 );
	itNumOld = getItNum( runData );
	runData.itNumNowRef[] = itNum;
	if itNum != itNumOld
		nSample = size( runData.sampleLstLst[1,1][1][1], 2 );
		for iN = 1 : getLnNCirc( runData ), iR = 1 : getLnRCirc( runData )
			runData.rotMatBackupLst[iR,iN] = [ [ zeros(3,3) for iCirc = 1 : runData.nCircLst[iN] ] for it = 1 : itNum ];
			runData.rotMat2dBackupLst[iR,iN] = [ [ zeros(2,2) for iCirc = 1 : runData.nCircLst[iN] ] for it = 1 : itNum ];
			runData.rotMat2dInvBackupLst[iR,iN] = [ [ zeros(2,2) for iCirc = 1 : runData.nCircLst[iN] ] for it = 1 : itNum ];
			runData.sampleLstLst[iR,iN] = [ [ zeros(2, nSample) for iCirc = 1 : runData.nCircLst[iN] ] for it = 1 : itNum ];
		end
		GC.gc();
		runBaseInfo!( runData );
	end
end

function storeZakCorrInRun!( runData::RunRandCircData, data::RandCircData, it::Int64 )
	getZakCorrTmpLst( runData )[:,:,it] .= getZakCorr( data );
end

function storeZakCorrInRun!( tmpData::ZakCorrArrTmpData, data::RandCircData, it::Int64 )
	getTmpArr( tmpData )[:,:,it] .= getZakCorr( data );
end

function meanZakCorr!( runData::RunRandCircData )
	mean!( getZakCorrMeanArr( runData ), getZakCorrTmpLst( runData ) );
end

function meanZakCorr1Pass!( runData::RunRandCircData )
	mean!( runData.zakCorr1dMeanLst1PassSingleton, runData.zakCorr1dLst1Pass );
end

function meanNumCornerLst!( runData::RunRandCircData )
	mean!( runData.numCornerMeanSingletonLst, runData.numCornerLst );
end

function meanZakAvgLst!( runData::RunRandCircData )
	mean!( runData.zakArrAvgMeanLstSingleton, runData.zakArrAvgLst );
end

function calcFitExp1Pass!( runData::RunRandCircData )
	runData.fitModel1PassLst .= ( corr -> LsqFit.curve_fit( modelExp, runData.xLstHalf1Pass, corr, p0FitModelExp ) ).( runData.zakCorrMean1dHalf1Pass );
	runData.corrLenInvLst1Pass .= ( x -> x.param[2] ).( runData.fitModel1PassLst );
	runData.corrLenLst1Pass .= 1 ./ runData.corrLenInvLst1Pass;
end

function calcDivNumNxt!( runData::RunRandCircData )
	runData.divNumNxtLst .= ( x -> Int64( floor( ( x < runData.corrLenLst1Pass[1] ? runData.divNum1Pass * runData.corrLenLst1Pass[1] / x : runData.divNum1Pass ) ) ) ).( runData.corrLenLst1Pass );
	runData.divNumNxtHalfLst .= div.( runData.divNumNxtLst, 2 );
	for iN = 1 : getLnNCirc( runData ), iR = 1 : getLnRCirc( runData )
		runData.xLstFineLst[iR,iN] = [0:runData.divNumNxtLst[iR,iN]-1;];
		runData.xLstFineLst[iR,iN] ./= runData.divNumNxtLst[iR,iN];
		runData.xHalfLstFineLst[iR,iN] = runData.xLstFineLst[iR,iN][1:runData.divNumNxtHalfLst[iR,iN]];
		runData.zakCorr1dStoreLst[iR,iN] = zeros( runData.divNumNxtLst[iR,iN], getItNumFine( runData ) );
		runData.zakCorrMean1dStoreLst[iR,iN] = zeros( runData.divNumNxtLst[iR,iN] );
	end
end

# function calcFitExpFine!( runData::RunRandCircData, iR::Int64, iN::Int64 )
	# # @infiltrate
	# fitModelFine = curve_fit( modelExp, runData.xHalfLstFineLst[iR,iN], getZakCorrMean1dHalf( runData ), p0FitModelExp );
	# runData.corrLenFineLst[iR,iN] = 1 / fitModelFine.param[2];
	# runData.expScaleFineLst[iR,iN] = fitModelFine.param[1];
	# runData.expShFineLst[iR,iN] = fitModelFine.param[3];
# end

function calcFitExpFine!( runData::RunRandCircData, iR::Int64, iN::Int64 )
	calcFitExpFine!( runData, getZakCorrMean1dHalf( runData ), iR, iN );
end

function calcFitExpFine!( runData::RunRandCircData, zakTmpData::ZakCorrArrTmpData, iR::Int64, iN::Int64 )
	calcFitExpFine!( runData, getMean1dHalfLst( zakTmpData ), iR, iN );
end

function calcFitExpFine!( runData::RunRandCircData, corrMean1dHalf::AbstractVector, iR::Int64, iN::Int64 )
	fitModelFine = curve_fit( modelExp, runData.xHalfLstFineLst[iR,iN], corrMean1dHalf, p0FitModelExp );
	runData.corrLenFineLst[iR,iN] = 1 / fitModelFine.param[2];
	runData.expScaleFineLst[iR,iN] = fitModelFine.param[1];
	runData.expShFineLst[iR,iN] = fitModelFine.param[3];
end

function runBaseInfo!( runData::RunRandCircData )
	runBaseInfo!( runData, runData.rCircLst );
end

function runBaseInfo!( runData::RunRandCircData, rCircLst::Vector{Float64} )
	rotMatBackupLst = runData.rotMatBackupLst;
	rotMat2dBackupLst = runData.rotMat2dBackupLst;
	rotMat2dInvBackupLst = runData.rotMat2dInvBackupLst;
	sampleLstLst = runData.sampleLstLst;
	
	for iN = 1 : getLnNCirc( runData )
		data = runData.randCircDataLst[iN];
		for iR = 1 : getLnRCirc( runData )
			setRCircNoRotUpdate!( data, rCircLst[iR] );
			for it = 1 : getItNum( runData )
				refreshPtLst!( data );
				refreshQuatLst!( data );
				refreshRotMatFull!( data );
				refreshEqCoeff!( data );
				refreshBndLstFull!( data );
				calcZakArr!( data );
				calcZakCorr!( data );
				
				backupRotMat!( rotMatBackupLst[iR, iN][it], rotMat2dBackupLst[iR, iN][it], rotMat2dInvBackupLst[iR, iN][it], data );
				# backupZakArrCorr!( @view( runData.zakArrLst1Pass[:, :, it, iR, iN] ), @view( runData.zakCorrLst1Pass[:, :, it, iR, iN] ), data );
				
				calcSampleLst!( data );
				backupSampleLst!( sampleLstLst[iR, iN][it], data );
			end
		end
	end
end

function run1Pass!( runData::RunRandCircData )
	rotMatBackupLst = runData.rotMatBackupLst;
	rotMat2dBackupLst = runData.rotMat2dBackupLst;
	rotMat2dInvBackupLst = runData.rotMat2dInvBackupLst;
	
	for iN = 1 : getLnNCirc( runData )
		data = runData.randCircDataLst[iN];
		for iR = 1 : getLnRCirc( runData ), it = 1 : getItNum1Pass( runData )
			restoreRotMat!( data, rotMatBackupLst[iR,iN][it], rotMat2dBackupLst[iR,iN][it], rotMat2dInvBackupLst[iR,iN][it] );
			calcZakArr!( data );
			calcZakCorr!( data );
			# backupZakArrCorr!( @view( runData.zakArrLst1Pass[:, :, it, iR, iN] ), @view( runData.zakCorrLst1Pass[:, :, it, iR, iN] ), data );
			runData.zakCorr1dLst1Pass[:,it,iR,iN] .= @view getZakCorr( data )[:,1];
		end
	end
	
	meanZakCorr1Pass!( runData );
	calcFitExp1Pass!( runData );
	calcDivNumNxt!( runData );
end

function runFine!( runData::RunRandCircData; isCornerDetect = true )
	setItNum!( runData, getItNumFine(runData) );
	for iN = 1 : getLnNCirc( runData )
		data = runData.randCircDataLst[iN];
		for iR = 1 : getLnRCirc( runData )
			GC.gc();
			divNum = runData.divNumNxtLst[iR,iN];
			setDivNum!( data, divNum );
			# setDivNum!( runData, divNum, runData.divNumNxtHalfLst[iR,iN] );
			setDivNum!( runData.zakTmpData, divNum, runData.divNumNxtHalfLst[iR,iN], getItNumFine( runData ) );
			cornerHelper = CornerDetector.JointFiltHelper( divNum );
			for it = 1 : getItNumFine( runData )
				restoreRotMat!( data, runData.rotMatBackupLst[iR,iN][it], runData.rotMat2dBackupLst[iR,iN][it], runData.rotMat2dInvBackupLst[iR,iN][it] );
				calcZakArr!( data );
				calcZakCorr!( data );
				# storeZakCorrInRun!( runData, data, it );
				storeZakCorrInRun!( runData.zakTmpData, data, it );
				# runData.zakCorr1dStoreLst[iR,iN][:,it] .= @view getZakCorrTmpLst( runData )[:,1,it];
				runData.zakCorr1dStoreLst[iR,iN][:,it] .= @view getTmpArr( runData.zakTmpData )[:,1,it];
				if isCornerDetect
					idCornerLst, numCorner = CornerDetector.genCornerIdNumFromJoint!( cornerHelper, getZakArr( data ) );
					runData.numCornerLst[it, iR, iN] = numCorner;
					runData.idCornerLstLst[it, iR, iN] = idCornerLst;
				end
				runData.zakArrAvgLst[it, iR, iN] = calcZakAvgAbs( data );
			end
			# meanZakCorr!( runData );
			meanZakCorr!( runData.zakTmpData );
			if runData.isStoreCorrFull
				# runData.zakCorrMeanFullStoreLst[iR,iN] = getZakCorrMeanArr( runData );
				runData.zakCorrMeanFullStoreLst[iR,iN] = getMeanArr( runData.zakTmpData );
			end
			# runData.zakCorrMean1dStoreLst[iR,iN] .= getZakCorrMean1d( runData );
			runData.zakCorrMean1dStoreLst[iR,iN] .= getMean1dLst( runData.zakTmpData );
			# calcFitExpFine!( runData, iR, iN );
			calcFitExpFine!( runData, runData.zakTmpData, iR, iN );
		end
	end
	meanNumCornerLst!( runData );
	meanZakAvgLst!( runData );
end

function exportData( runData::RunRandCircData )
	return runData.corrLenFineLst, runData.numCornerMeanLst;
end

function exportDataDetailed( runData::RunRandCircData )
	dataOutLst = [ runData.corrLenFineLst, runData.expScaleFineLst, runData.expShFineLst, runData.zakCorrMean1dStoreLst, runData.zakCorr1dStoreLst, runData.zakArrAvgLst, runData.zakArrAvgMeanLst ];
	return dataOutLst;
end

function exportDataCorner( runData::RunRandCircData )
	return runData.numCornerLst, runData.numCornerMeanLst, runData.idCornerLstLst;
end

function exportParams( runData::RunRandCircData )
	return runData.nCircLst, runData.rCircLst, getItNumFine( runData ), runData.divNum1Pass, runData.divNumNxtLst;
end

function exportFullCorrMean( runData::RunRandCircData )
	return runData.zakCorrMeanFullStoreLst;
end
