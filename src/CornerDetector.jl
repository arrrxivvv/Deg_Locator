module CornerDetector

using ImageFiltering
using OffsetArrays
using ImgProcessing
using DataStructures
using StaticArrays
using Utils

using Infiltrator

function genSteerFilts( sz::Int64 )
	xLst = [-sz:sz;];
	yLst = deepcopy( xLst )';
	
	steerFiltLstRaw = [ xLst .* yLst, ( xLst.^2 .- (yLst.^2) ) ./ 2 ];
	
	filtFact = ImgProcessing.gaussFiltLstNoOffset[sz+1];
	@infiltrate
	println(size(steerFiltLstRaw[1]));
	println(size(filtFact));
	
	steerFiltLst1 = similar(steerFiltLstRaw[1]);
	steerFiltLst1 .= steerFiltLstRaw[1] .* filtFact;
	
	steerFiltLst = ( (x) -> (x .* filtFact) ).( steerFiltLstRaw );
	
	return steerFiltLst;
end

struct PrefabSteerFilts
	steerFiltLstLst::Vector{Vector{OffsetMatrix{Float64,Matrix{Float64}}}};
	
	function PrefabSteerFilts( lnMax::Int64 = 3 )
		steerFiltLstLstRaw = genSteerFilts.( [1:lnMax;] );
		steerFiltLstLst = (x->centered.(x)).( steerFiltLstLstRaw );
		@infiltrate
		println( typeof(steerFiltLstLstRaw) )
		println( typeof(steerFiltLstLst) )
		
		new( steerFiltLstLst );
	end
end

const prefabSteerFilts = PrefabSteerFilts();
const steerFiltLstLst = prefabSteerFilts.steerFiltLstLst;

abstract type AbstractFiltHelperData end

function throwUndefined( filtDataType::Type{<:AbstractFiltHelperData} )
	error( "FiltData not defined" );
end
throwUndefined( filtData::AbstractFiltHelperData ) = throwUndefined( typeof(filtData) );

getMaxedFiltArr( filtData::AbstractFiltHelperData ) = throwUndefined(typeof(filtData));
getIsMaxedArr( filtData::AbstractFiltHelperData ) = throwUndefined(typeof(filtData));
getIsVisitedArr( filtData::AbstractFiltHelperData ) = throwUndefined(typeof(filtData));
getFiltedArr( filtData::AbstractFiltHelperData ) = throwUndefined(typeof(filtData));


struct SteerFiltHelperData <: AbstractFiltHelperData
	steerFiltedArr::Matrix{Float64};
	steerFiltedXYLst::Vector{Matrix{Float64}};
	maxedFiltArr::Matrix{Float64};
	isMaxedArr::Matrix{Bool};
	isVisitedArr::Matrix{Bool};
	
	function SteerFiltHelperData( sz::Int64 )
		steerFiltedArr = zeros( sz, sz );
		steerFiltedXYLst = [ similar(steerFiltedArr) for iXY = 1 : 2 ];
		maxedFiltArr = similar( steerFiltedArr );
		isMaxedArr = similar( maxedFiltArr, Bool );
		isVisitedArr = similar(isMaxedArr);
		
		new( steerFiltedArr, steerFiltedXYLst, maxedFiltArr, isMaxedArr, isVisitedArr );
	end
end

function getMaxedFiltArr( steerFiltData::SteerFiltHelperData )
	return steerFiltData.maxedFiltArr
end

function getIsMaxedArr( steerFiltData::SteerFiltHelperData )
	return steerFiltData.isMaxedArr;
end

function getIsVisitedArr( steerFiltData::SteerFiltHelperData )
	return steerFiltData.isVisitedArr;
end

function getFiltedArr( steerFiltData::SteerFiltHelperData )
	return steerFiltData.steerFiltedArr;
end


struct HarrisFiltHelperData <: AbstractFiltHelperData
	filtedArr::Matrix{Float64};
	maxedFiltArr::Matrix{Float64};
	isMaxedArr::Matrix{Bool};
	isVisitedArr::Matrix{Bool};
	
	diffArrLst::Vector{Matrix{Float64}};
	covArrMat::Matrix{Matrix{Float64}};
	covFiltedArrMat::Matrix{Matrix{Float64}};
	
	function HarrisFiltHelperData( sz::Int64 )
		filtedArr = zeros( sz, sz );
		maxedFiltArr = similar( filtedArr );
		isMaxedArr = similar( maxedFiltArr, Bool );
		isVisitedArr = similar(isMaxedArr);
		
		diffArrLst = [ similar(filtedArr) for iXY = 1 : 2 ];
		covArrMat = [ similar(filtedArr) for iX = 1 : 2, iY = 1:2 ];
		covFiltedArrMat = deepcopy( covArrMat );
		
		new( filtedArr, maxedFiltArr, isMaxedArr, isVisitedArr, diffArrLst, covArrMat, covFiltedArrMat );
	end
end

function getFiltedArr( harrisFiltData::HarrisFiltHelperData )
	return harrisFiltData.filtedArr;
end

function getMaxedFiltArr( harrisFiltData::HarrisFiltHelperData )
	return harrisFiltData.maxedFiltArr;
end

function getIsMaxedArr( harrisFiltData::HarrisFiltHelperData )
	return harrisFiltData.isMaxedArr
end

function getIsVisitedArr( harrisFiltData::HarrisFiltHelperData )
	return harrisFiltData.isVisitedArr;
end




function genCovMat!( harrisFiltData::HarrisFiltHelperData, imgArr::AbstractArray{<:Number} )
	return genCovMat!( harrisFiltData.covArrMat, harrisFiltData.diffArrLst, imgArr );
end

function genCovMat!( covMat::Matrix{<:AbstractArray{<:Number}}, diffLst::Vector{<:AbstractArray{<:Number}}, imgArr::AbstractArray{<:Number} )
	for iDiff = 1 : 2
		ImgProcessing.convolute!( diffLst[iDiff], imgArr, ImgProcessing.sobelFiltLst[iDiff] );
	end
	
	diffX, diffY = diffLst;
	
	covMat[1,1] .= diffX.^2;
	covMat[2,2] .= diffY.^2;
	covMat[1,2] .= diffX .* diffY;
	covMat[2,1] .= diffY .* diffX;	
end

function genHarrisCornerFiltFromCovMat!( harrisFiltData::HarrisFiltHelperData, filt::AbstractMatrix{<:Number} )
	return genHarrisCornerFiltFromCovMat!( getFiltedArr( harrisFiltData ), harrisFiltData.covFiltedArrMat, harrisFiltData.covArrMat, filt );
end

function genHarrisCornerFiltFromCovMatBox!( harrisFiltData::HarrisFiltHelperData; filtLen::Int64 = 1 )
	return genHarrisCornerFiltFromCovMatBox!( getFiltedArr( harrisFiltData ), harrisFiltData.covFiltedArrMat, harrisFiltData.covArrMat; filtLen = filtLen );
end

function genHarrisCornerFiltFromCovMatGauss!( harrisFiltData::HarrisFiltHelperData; filtLen::Int64 = 1 )
	return genHarrisCornerFiltFromCovMatGauss!( getFiltedArr( harrisFiltData ), harrisFiltData.covFiltedArrMat, harrisFiltData.covArrMat; filtLen = filtLen );
end

function genHarrisCornerFiltFromCovMat!( harrisArr::Matrix{<:Number}, covMatFilted::Matrix{<:AbstractMatrix{<:Number}}, covMat::Matrix{<:AbstractMatrix{<:Number}}, filt::AbstractMatrix{<:Number}; kappa::Real = 0.05 )
	for ii in eachindex(covMat)
		ImgProcessing.convolute!( covMatFilted[ii], covMat[ii], filt );
	end
	
	harrisArr .= ( covMatFilted[1,1] .* covMatFilted[2,2] .- covMatFilted[1,2] .* covMatFilted[2,1] ) .- kappa .* (covMatFilted[1,1] .+ covMatFilted[2,2]).^2;
end

function genHarrisCornerFiltFromCovMatBox!( harrisArr::Matrix{<:Number},covMatFilted::Matrix{<:AbstractMatrix{<:Number}}, covMat::Matrix{<:AbstractMatrix{<:Number}}; filtLen::Int64 = 1 )
	boxFilt = ImgProcessing.boxFiltLst[filtLen];
	
	genHarrisCornerFiltFromCovMat!( harrisArr, covMatFilted, covMat, boxFilt );
end

function genHarrisCornerFiltFromCovMatGauss!( harrisArr::Matrix{<:Number},covMatFilted::Matrix{<:AbstractMatrix{<:Number}}, covMat::Matrix{<:AbstractMatrix{<:Number}}; filtLen::Int64 = 1 )
	gaussFilt = ImgProcessing.gaussFiltLst[filtLen];
	
	genHarrisCornerFiltFromCovMat!( harrisArr, covMatFilted, covMat, gaussFilt );
end

function genSteerCornerFilt!( steerArr::Matrix{<:Number}, steerFiltedLst::AbstractVector{<:AbstractMatrix{<:Number}}, imgArr::AbstractMatrix{<:Number}, lnFilt::Int64 )
	steerFiltLst = CornerDetector.steerFiltLstLst[lnFilt];
	
	for ii = 1 : 2
		ImgProcessing.convolute!( steerFiltedLst[ii], imgArr, steerFiltLst[ii] );
		steerFiltedLst[ii] .= steerFiltedLst[ii].^2;
	end
	
	# steerArr .= sum.( steerFiltedLst );
	steerArr .= 0;
	for ii = 1 : 2
		steerArr .= steerArr .+ steerFiltedLst[ii];
	end
end


function genSteerCornerFilt!( steerData::SteerFiltHelperData, imgArr::AbstractMatrix{<:Number}, lnFilt::Int64 )
	return genSteerCornerFilt!( steerData.steerFiltedArr, steerData.steerFiltedXYLst, imgArr, lnFilt );
end



function nonMaxSuppress!( filtData::AbstractFiltHelperData; lnMax::Int64 = 1 )
	return ImgProcessing.nonMaxSuppress!( getMaxedFiltArr( filtData ), getFiltedArr( filtData ); ln = lnMax );
end

function extractMaxIdWithConnectedComp!( filtData::AbstractFiltHelperData, thres::Float64, lnFilt::Int64 )
	return extractMaxIdWithConnectedComp!( getIsMaxedArr( filtData ), getIsVisitedArr( filtData ), getMaxedFiltArr( filtData ), thres, lnFilt );
end

function extractMaxIdWithConnectedComp!( isMaxArr, isVisitedArr, maxedArr::AbstractArray{<:Number}, thres::Float64, lnFilt::Int64 )
	nDim = 2;
	sz = size( maxedArr, 1 );
	idLst = findall( x -> x > thres, maxedArr );
	# isMaxArr = zeros(Bool, size(maxedArr) );
	isVisitedArr .= false;
	isMaxArr .= false;
	for ii in idLst
		isMaxArr[ii] = true;
	end
	iSh0 = ntuple(x->0,nDim);
	
	idSearchQueue = Queue{MVector{nDim,Int64}}();
	idSearchLst = Vector{MVector{nDim,Int64}}(undef,0);
	idMergedLst = Vector{MVector{nDim,Float64}}(undef,0);;
	
	rngWind1d = -lnFilt : lnFilt;
	rngWind = Iterators.product( rngWind1d, rngWind1d );
	idTmp = ones(MVector{nDim,Int64});
	idThis = copy( idTmp );
	idMerged = ones(MVector{nDim,Float64});
	
	for iId = 1 : length(idLst)
		idCart = idLst[iId];
		if !isVisitedArr[idCart]
			empty!(idSearchLst);
			for ii = 1 : length(idThis)
				idThis[ii] = idCart[ii];
			end
			# @infiltrate
			enqueue!( idSearchQueue, copy(idThis) );
			push!( idSearchLst, last(idSearchQueue) );
			isVisitedArr[idThis...] = true;
			# @time begin
			while !isempty( idSearchQueue )
				id = dequeue!( idSearchQueue );
				for iSh in rngWind
					if iSh == iSh0
						continue;
					end
					idTmp .= id;
					idTmp .= idTmp .+ iSh;
					idTmp .= wrapIntInd.( idTmp, sz );
					if isMaxArr[idTmp...] && !isVisitedArr[idTmp...]
						enqueue!( idSearchQueue,  copy( idTmp ) );
						push!( idSearchLst, last( idSearchQueue ) );
						isVisitedArr[idTmp...] = true;
					end
				end
				# @infiltrate
			end
			# end
			# @infiltrate
			
			if length( idSearchQueue ) == 1
				idMerged .= idSearchLst[1];
			else
				idMerged .= 0;
				for iSearch = 1 : length(idSearchLst)
					idMerged .= idMerged .+ idSearchLst[iSearch];
				end
				idMerged .= idMerged ./ length(idSearchLst);
			end
			push!( idMergedLst, copy(idMerged) );
		end
	end
	
	return idMergedLst;
end

# function extractMaxIdWithConnectedComp!( isMaxArr, isVisitedArr, maxedArr::AbstractArray{<:Number}, thres::Float64, lnFilt::Int64 )
	# nDim = 2;
	# sz = size( maxedArr, 1 );
	# idLst = findall( x -> x > thres, maxedArr );
	# # isMaxArr = zeros(Bool, size(maxedArr) );
	# isVisitedArr .= false;
	# isMaxArr .= false;
	# for ii in idLst
		# isMaxArr[ii] = true;
	# end
	# iSh0 = ntuple(x->0,nDim);
	
	# idSearchQueue = Queue{MVector{nDim,Int64}}();
	# idSearchLst = Vector{MVector{nDim,Int64}}(undef,0);
	# idMergedLst = Vector{MVector{nDim,Float64}}(undef,0);;
	
	# rngWind1d = -lnFilt : lnFilt;
	# rngWind = Iterators.product( rngWind1d, rngWind1d );
	# idTmp = ones(MVector{nDim,Int64});
	# idThis = copy( idTmp );
	# idMerged = ones(MVector{nDim,Float64});
	
	# for iId = 1 : length(idLst)
		# idCart = idLst[iId];
		# if !isVisitedArr[idCart]
			# empty!(idSearchLst);
			# for ii = 1 : length(idThis)
				# idThis[ii] = idCart[ii];
			# end
			# # @infiltrate
			# enqueue!( idSearchQueue, copy(idThis) );
			# push!( idSearchLst, last(idSearchQueue) );
			# isVisitedArr[idThis...] = true;
			# # @time begin
			# while !isempty( idSearchQueue )
				# id = dequeue!( idSearchQueue );
				# for iSh in rngWind
					# if iSh == iSh0
						# continue;
					# end
					# idTmp .= id;
					# idTmp .= idTmp .+ iSh;
					# idTmp .= wrapIntInd.( idTmp, sz );
					# if isMaxArr[idTmp...] && !isVisitedArr[idTmp...]
						# enqueue!( idSearchQueue,  copy( idTmp ) );
						# push!( idSearchLst, last( idSearchQueue ) );
						# isVisitedArr[idTmp...] = true;
					# end
				# end
				# # @infiltrate
			# end
			# # end
			# # @infiltrate
			
			# if length( idSearchQueue ) == 1
				# idMerged .= idSearchLst[1];
			# else
				# idMerged .= 0;
				# for iSearch = 1 : length(idSearchLst)
					# idMerged .= idMerged .+ idSearchLst[iSearch];
				# end
				# idMerged .= idMerged ./ length(idSearchLst);
			# end
			# push!( idMergedLst, copy(idMerged) );
		# end
	# end
	
	# return idMergedLst, isMaxArr;
# end

function genCornerIdIdArrFromFilt( filtedArr::AbstractMatrix{<:Number}, thres::Float64 )
	;
end

function genCornerIdSteerHarrisMerged( idSteerMergedLst, isHarrisMaxArr::Array{Bool}, lnFilt::Int64 )
	nDim = 2;
	idSteerHarrisMergedLst = Vector{MVector{nDim,Float64}}(undef,0);
	sz = size( isHarrisMaxArr, 1 );
	
	idFloored = zeros( MVector{nDim,Int64} );
	idTmp = copy(idFloored);
	rng1d = -lnFilt : lnFilt;
	rng1dFloat = -lnFilt : lnFilt+1;
	rng2d = Iterators.product( rng1d, rng1d );
	rng2dFloat = Iterators.product( rng1dFloat, rng1dFloat );
	
	isIntegerLst = zeros(Bool, nDim);
	isAllInteger = false;
	
	for id in idSteerMergedLst
		idFloored .= floor.( id );
		
		isIntegerLst .= isinteger.(id);
		isAllInteger = reduce(&,isIntegerLst);
		if isAllInteger
			rng = rng2d;
		else
			rng = rng2dFloat;
		end
		
		isNearHarris = false;
		for iSh in rng 
			idTmp .= idFloored .+ iSh;
			idTmp .= wrapIntInd.( idTmp, sz );
			if isHarrisMaxArr[idTmp...]
				isNearHarris = true;
				push!(idSteerHarrisMergedLst,id);
				break;
			end
		end
		# @infiltrate
		# if isNearHarris		
	end
	
	return idSteerHarrisMergedLst;
end

end # endmodule
