# Loops_MC module

fNameFileLstWL = "fNameFileLstWL.txt";

function calcDL( params::ParamsLoops, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D}
	dL = 0;
	for iLnkDim in 1 : params.nDimLayer
		dimLink = params.linkDimLst[dim][iLnkDim];
		dimLinkSh = params.linkDimShLst[dim][iLnkDim];
		dL += boolToFlipChange( linkLst[dimLink][pos] );
		dL += boolToFlipChange( linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] );
	end
	
	return dL;
end

abstract type AbstractWLHistDos end

throwWLHistDosUndefined() = error( "Lops_MC: WLHistDos undefined" );

getWLHistDosName( wlHistDosType::Type{<:AbstractWLHistDos} ) = throwWLHistDosUndefined();

getAttrLst( wlHistDosType::Type{<:AbstractWLHistDos} ) = throwWLHistDosUndefined();
getAttrLst( wlHistDos::AbstractWLHistDos ) = getAttrLst(typeof(wlHistDos));

getValLst( wlHistDos::AbstractWLHistDos ) = throwWLHistDosUndefined();

function getAttrValLst( wlHistDos::AbstractWLHistDos ) 
	attrLst = getAttrLst(wlHistDos);
	valLst = getValLst(wlHistDos);
	
	return attrLst, valLst;
end

function getHistArr( histDos::AbstractWLHistDos )
	if getIsHistOutBackup( histDos )
		return histDos.histArrFlatBackup;
	else
		return histDos.histArr;
	end
end

function getDosArr( histDos::AbstractWLHistDos )
	return histDos.dosArr;
end

getIsHistOutBackup( histDos::AbstractWLHistDos ) = histDos.isHistOutBackupRef[];

function setHistOutBackup!( histDos::AbstractWLHistDos )
	histDos.isHistOutBackupRef[] = true;
end

function unsetHistOutBackup!( histDos::AbstractWLHistDos )
	histDos.isHistOutBackupRef[] = false;
end

getHistUpdateId( wlHistDos::AbstractWLHistDos, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D} = throwWLHistDosUndefined();

getHistUpdateId( wlHistDos::AbstractWLHistDos, params::ParamsLoops, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D} = throwWLHistDosUndefined();

getHistThisId( wlHistDos::AbstractWLHistDos, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}} ) where {D} = throwWLHistDosUndefined();

getAndTestHistThisId( wlHistDos::AbstractWLHistDos, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}} ) where {D} = throwWLHistDosUndefined();

getHistThisId( wlHistDos::AbstractWLHistDos, params::ParamsLoops, linkLst::Vector{Array{Bool,D}} ) where {D} = throwWLHistDosUndefined();

getLinkHistId( wlHistDos::AbstractWLHistDos, lnkVal::Int64 ) = throwWLHistDosUndefined();

function histUpdateExchange( histDos1::AbstractWLHistDos, histDos2::AbstractWLHistDos, params::ParamsLoops, bLinkData1::BLinkAuxData, bLinkData2::BLinkAuxData, dosIncr1::Float64, dosIncr2::Float64 )
	isFlip = histUpdateExchangeNoBackupSet( histDos1, histDos2, params, bLinkData1, bLinkData2, dosIncr1, dosIncr2 )
	unsetHistOutBackup!( histDos1 );
	unsetHistOutBackup!( histDos2 );
	
	return isFlip;
end

function histUpdateExchangeNoBackupSet( histDos1::AbstractWLHistDos, histDos2::AbstractWLHistDos, params::ParamsLoops, bLinkData1::BLinkAuxData, bLinkData2::BLinkAuxData, dosIncr1::Float64, dosIncr2::Float64 )
	histDosLst = [histDos1, histDos2];
	bLinkDataLst = [bLinkData1, bLinkData2];
	idLst = zeros(eltype(histDos1.posLst), length(histDosLst), length(bLinkDataLst));
	isOutBndLst = similar( idLst, Bool );
	for iHist = 1 : length(histDosLst), iBLink = 1 : length(bLinkDataLst)
		histDos = histDosLst[iHist];
		bLinkData = bLinkDataLst[iBLink];
		idLst[iHist, iBLink], isOutBndLst[iHist, iBLink] = getAndTestHistThisId( histDos, params, getBfieldLst( bLinkData ), getLinkLst( bLinkData ) )
	end
	if any(isOutBndLst)
		isFlip = false;
		return isFlip;
	end
	# idLst = [ getHistThisId( histDos, params, bLinkData.BfieldLst, bLinkData.linkLst ) for histDos in histDosLst, bLinkData in bLinkDataLst ];
	# id1 = getHistThisId( histDos1, params, bLinkData1.BfieldLst, bLinkData1.linkLst );
	# id2 = getHistThisId( histDos1, params, bLinkData2.BfieldLst, bLinkData2.linkLst );
	
	p = exp( histDos1.dosArr[idLst[1,1]] - histDos1.dosArr[idLst[1,2]] - histDos2.dosArr[idLst[2,1]] + histDos2.dosArr[idLst[2,2]] );
	isFlip = rand() < p;
	
	if isFlip
		idEnd1 = idLst[1,2];
		idEnd2 = idLst[2,1];
	else
		idEnd1 = idLst[1,1] ;
		idEnd2 = idLst[2,2];
	end
	
	histUpdateWithId( histDos1, dosIncr1, idEnd1 );
	histUpdateWithId( histDos2, dosIncr2, idEnd2 );
		
	return isFlip;
end

function histUpdate( histDos::AbstractWLHistDos, params::ParamsLoops, flipProposer::FlipProposer, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D}, dosIncr::Float64 ) where {D}
	# id, idNxt = getHistUpdateId( histDos, params, linkLst, dim, pos );
	# isFlip = histUpdateWithId( histDos, dosIncr, id, idNxt );
	isFlip = histUpdateNoBackupSet( histDos, params, flipProposer, BfieldLst, linkLst, dim, pos, dosIncr )
	unsetHistOutBackup!( histDos );
	
	return isFlip;
end

function histUpdateNoBackupSet( histDos::AbstractWLHistDos, params::ParamsLoops, flipProposer::FlipProposer, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D}, dosIncr::Float64 ) where {D}
	id, idNxt = getHistUpdateId( histDos, params, BfieldLst, linkLst, dim, pos );
	isFlip = histUpdateWithId( histDos, dosIncr, id, idNxt );
	
	return isFlip;
end

function histUpdateWithId( histDos::AbstractWLHistDos, dosIncr::Float64, lId::CartesianIndex{histD}, lIdNxt::CartesianIndex{histD} ) where {histD}
	p = exp( histDos.dosArr[lId] - histDos.dosArr[lIdNxt] );
	isFlip = rand() < p;
	
	if isFlip
		histUpdateWithId( histDos, dosIncr, lIdNxt );
	else
		histUpdateWithId( histDos, dosIncr, lId );
	end
	
	
	return isFlip;
end

function histUpdateWithId( histDos::AbstractWLHistDos, dosIncr::Float64, lId::CartesianIndex{D_hist} ) where {D_hist}
	histDos.histArr[lId] += 1;
	histDos.dosArr[lId] += dosIncr;
end

function histResetZero!( histDos::AbstractWLHistDos )
	histDos.histArrFlatBackup .= histDos.histArr;
	histDos.histArr .= 0;
	
	setHistOutBackup!( histDos );
end

function getLinkHistIdFull( wlHistDos::AbstractWLHistDos, lnkVal::Int64 )
	lnkId = Int64( lnkVal / 2 );
	if lnkId == wlHistDos.numLnkVal
		lnkId -= 2;
	elseif lnkId > 0
		lnkId -= 1;
	end
	
	lnkId += 1;
	
	return lnkId;
end



struct WLHistDosFull{nDim,D_hist} <: AbstractWLHistDos
	histArr::Array{Int64,D_hist};
	histArrFlatBackup::Array{Int64,D_hist};
	dosArr::Array{Float64,D_hist};
	
	isHistOutBackupRef::Ref{Bool};
	
	numLnkVal::Int64;
	
	posLst::CartesianIndices{D_hist,Tuple{Vararg{Base.OneTo{Int64},D_hist}}};
	
	function WLHistDosFull{nDim}( histArr::Array{Int64,D_hist}, numLnkVal::Int64 ) where {nDim, D_hist}
		histArrFlatBackup = similar(histArr);
		# D_hist = ndims(histArr);
		dosArr = similar(histArr, Float64);
		dosArr .= 0;
		
		isHistOutBackupVal = false;
		
		posLst = CartesianIndices(histArr);
		
		new{nDim,D_hist}( histArr, histArrFlatBackup, dosArr, isHistOutBackupVal, numLnkVal, posLst );
	end
end

function getLinkHistId( wlHistDos::WLHistDosFull, lnkVal::Int64 )
	# lnkId = Int64( lnkVal / 2 );
	# if lnkId == wlHistDos.grdNum
		# lnkId -= 2;
	# elseif lnkId > 0
		# lnkId -= 1;
	# end
	
	# lnkId += 1;
	
	return getLinkHistIdFull( wlHistDos, lnkVal );
end

const WLHistDosFull1dDos = WLHistDosFull{2,1};
const WLHistDosFull2dDos = WLHistDosFull{2,2};

function WLHistDosFull1dDos( divNum )
	D_hist = 1;
	
	nDim = 2;
	numLnkVal = divNum^nDim;
	histArr = zeros( Int64, numLnkVal+1-2 );
	
	return WLHistDosFull{nDim}( histArr, numLnkVal );
end

function getWLHistDosName( wlHistDosType::Type{<:WLHistDosFull} )
	return "WLHistDosFull" * "_" * string(wlHistDosType.parameters[2]) * "dDos";
end

function getAttrLst( wlHistDosType::Type{<:WLHistDosFull} )
	return [ getWLHistDosName( wlHistDosType ) * "_grdNum" ]
end

function getValLst( wlHistDos::WLHistDosFull )
	return [ wlHistDos.numLnkVal ];
end

function getHistUpdateId( wlHistDos::WLHistDosFull1dDos, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D}
	dL = calcDL( params, linkLst, dim, pos );
	
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
	lId = getLinkHistId( wlHistDos, lTotal );
	lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	return wlHistDos.posLst[lId], wlHistDos.posLst[lIdNxt];
end

function getHistUpdateId( wlHistDos::WLHistDosFull{nDim,1}, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D,nDim}
	dL = calcDL( params, linkLst, dim, pos );
	
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
	lId = getLinkHistId( wlHistDos, lTotal );
	lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	return wlHistDos.posLst[lId], wlHistDos.posLst[lIdNxt];
end

# function getLinkHistId( wlHistDos::WLHistDosFull1dDos, lnkVal::Int64 )
	# # lnkId = Int64( lnkVal / 2 );
	# # if lnkId == wlHistDos.grdNum
		# # lnkId -= 2;
	# # elseif lnkId > 0
		# # lnkId -= 1;
	# # end
	
	# # lnkId += 1;
	
	# return getLinkHistIdFull( wlHistDos, lnkVal );
# end




function WLHistDosFull2dDos( divNum )
	D_hist = 2;
	nDim = 2;
	numLnkVal = divNum^nDim;
	histArr = zeros( Int64, numLnkVal+1, numLnkVal+1-2 );
	
	WLHistDosFull{nDim}( histArr, numLnkVal );
end

function getHistUpdateId( wlHistDos::WLHistDosFull{nDim,2}, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D, nDim}
	dS = boolToFlipChange( BfieldLst[dim][pos] );
	dL = calcDL( params, linkLst, dim, pos );
	
	sTotal = sum( sum.( BfieldLst ) );
	lTotal = sum( sum.( linkLst ) );
	
	sNxt = sTotal + dS;
	lNxt = lTotal + dL;
	
	sId = sTotal + 1;
	sIdNxt = sNxt + 1;
	
	lId = getLinkHistId( wlHistDos, lTotal );
	lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	return wlHistDos.posLst[sId,lId], wlHistDos.posLst[sIdNxt,lIdNxt];
end

# function getHistUpdateId( wlHistDos::WLHistDosFull2dDos, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D}
	# dS = boolToFlipChange( BfieldLst[dim][pos] );
	# dL = calcDL( params, linkLst, dim, pos );
	
	# sTotal = sum( sum.( BfieldLst ) );
	# lTotal = sum( sum.( linkLst ) );
	
	# sNxt = sTotal + dS;
	# lNxt = lTotal + dL;
	
	# sId = sTotal + 1;
	# sIdNxt = sNxt + 1;
	
	# lId = getLinkHistId( wlHistDos, lTotal );
	# lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	# return wlHistDos.posLst[sId,lId], wlHistDos.posLst[sIdNxt,lIdNxt];
# end



const WL3dHistDosFull1dDos = WLHistDosFull{3,1};
const WL3dHistDosFull2dDos = WLHistDosFull{3,2};


function getWLHistDosName( wlHistDosType::Type{<:WLHistDosFull{3}} )
	 return "WLHistDosFull" * "_" * string(wlHistDosType.parameters[2]) * "dDos" * "_" * string(wlHistDosType.parameters[1]) * "d" ;
end



function WL3dHistDosFull1dDos( divNum )
	D_hist = 1;
	
	nDim = 3;
	grdNum = nDim * divNum^nDim;
	numLnkVal = Int64( floor( grdNum/2 ) );
	histArr = zeros( Int64, numLnkVal + 1 - 2 );
	
	return WLHistDosFull{nDim}( histArr, numLnkVal );
end

# function getHistUpdateId( wlHistDos::WL3dHistDosFull1dDos, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D}
	# dL = calcDL( params, linkLst, dim, pos );
	
	# lTotal = sum( sum.(linkLst) );
	# lNxt = lTotal + dL;
	
	# lId = getLinkHistId( wlHistDos, lTotal );
	# lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	# return wlHistDos.posLst[lId], wlHistDos.posLst[lIdNxt];
# end




function WL3dHistDosFull2dDos( divNum )
	D_hist = 2;
	nDim = 3;
	grdNum = nDim * divNum^nDim;
	numLnkVal = Int64( floor( grdNum/2 ) );
	histArr = zeros( Int64, grdNum + 1, numLnkVal + 1 - 2 );
	
	WLHistDosFull{nDim}( histArr, numLnkVal );
end

# function getHistUpdateId( wlHistDos::WLHistDosFull2dDos, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D}
	# dS = boolToFlipChange( BfieldLst[dim][pos] );
	# dL = calcDL( params, linkLst, dim, pos );
	
	# sTotal = sum( sum.( BfieldLst ) );
	# lTotal = sum( sum.( linkLst ) );
	
	# sNxt = sTotal + dS;
	# lNxt = lTotal + dL;
	
	# sId = sTotal + 1;
	# sIdNxt = sNxt + 1;
	
	# lId = getLinkHistId( wlHistDos, lTotal );
	# lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	# return wlHistDos.posLst[sId,lId], wlHistDos.posLst[sIdNxt,lIdNxt];
# end





struct WLHistDos2DFull <: AbstractWLHistDos
	histArr::Vector{Int64};
	histArrFlatBackup::Vector{Int64};
	dosArr::Vector{Float64};
	
	isHistOutBackupRef::Ref{Bool};
	
	grdNum::Int64;
	
	posLst::CartesianIndices{1,Tuple{Vararg{Base.OneTo{Int64},1}}};
	
	function WLHistDos2DFull( divNum )
		nDim = 2;
		grdNum = divNum^nDim;
		histArr = zeros( Int64, grdNum+1-2 );
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr, Float64);
		dosArr .= 0;
		
		isHistOutBackupVal = false;
		
		posLst = CartesianIndices(histArr);
		
		new( histArr, histArrFlatBackup, dosArr, isHistOutBackupVal, grdNum, posLst );
	end
end

function getWLHistDosName( wlHistDosType::Type{<:WLHistDos2DFull} )
	return "WLHist2DFull";
end

function getAttrLst( wlHistDosType::Type{<:WLHistDos2DFull} )
	return [ getWLHistDosName( wlHistDosType ) * "_grdNum" ]
end

function getValLst( wlHistDos::WLHistDos2DFull )
	return [ wlHistDos.grdNum ];
end

function getHistUpdateId( wlHistDos::WLHistDos2DFull, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D}
	dL = calcDL( params, linkLst, dim, pos );
	
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
	lId = getLinkHistId( wlHistDos, lTotal );
	lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	return wlHistDos.posLst[lId], wlHistDos.posLst[lIdNxt];
end

function getLinkHistId( wlHistDos::WLHistDos2DFull, lnkVal::Int64 )
	lnkId = Int64( lnkVal / 2 );
	if lnkId == wlHistDos.grdNum
		lnkId -= 2;
	elseif lnkId > 0
		lnkId -= 1;
	end
	
	lnkId += 1;
	
	return lnkId;
end




struct WLHistDosJoint2DFull <: AbstractWLHistDos
	histArr::Matrix{Int64};
	histArrFlatBackup::Matrix{Int64};
	dosArr::Matrix{Float64};
	
	isHistOutBackupRef::Ref{Bool};
	
	grdNum::Int64;
	
	posLst::CartesianIndices{2,Tuple{Vararg{Base.OneTo{Int64},2}}};
	
	function WLHistDosJoint2DFull( divNum )
		nDim = 2;
		grdNum = divNum^nDim;
		histArr = zeros( Int64, grdNum+1, grdNum+1-2 );
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr, Float64);
		dosArr .= 0;
		
		isHistOutBackupVal = false;
		
		posLst = CartesianIndices(histArr);
		
		new( histArr, histArrFlatBackup, dosArr, isHistOutBackupVal, grdNum, posLst );
	end
end

function getWLHistDosName( wlHistDosType::Type{<:WLHistDosJoint2DFull} )
	return "WLHistJoint2DFull";
end

function getAttrLst( wlHistDosType::Type{<:WLHistDosJoint2DFull} )
	return [ getWLHistDosName( wlHistDosType ) * "_grdNum" ]
end

function getValLst( wlHistDos::WLHistDosJoint2DFull )
	return [ wlHistDos.grdNum ];
end

function getHistUpdateId( wlHistDos::WLHistDosJoint2DFull, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D}
	dS = boolToFlipChange( BfieldLst[dim][pos] );
	dL = calcDL( params, linkLst, dim, pos );
	
	sTotal = sum( sum.( BfieldLst ) );
	lTotal = sum( sum.( linkLst ) );
	
	sNxt = sTotal + dS;
	lNxt = lTotal + dL;
	
	sId = sTotal + 1;
	sIdNxt = sNxt + 1;
	
	lId = getLinkHistId( wlHistDos, lTotal );
	lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	return wlHistDos.posLst[sId,lId], wlHistDos.posLst[sIdNxt,lIdNxt];
end

function getLinkHistId( wlHistDos::WLHistDosJoint2DFull, lnkVal::Int64 )
	lnkId = Int64( lnkVal / 2 );
	if lnkId == wlHistDos.grdNum
		lnkId -= 2;
	elseif lnkId > 0
		lnkId -= 1;
	end
	
	lnkId += 1;
	
	return lnkId;
end




struct WLHistDos2DHalf <: AbstractWLHistDos
	histArr::Vector{Int64};
	histArrFlatBackup::Vector{Int64};
	dosArr::Vector{Float64};
	
	isHistOutBackupRef::Ref{Bool};
	
	grdNum::Int64;
	grdNumHalf::Int64;
	
	posLst::CartesianIndices{1,Tuple{Vararg{Base.OneTo{Int64},1}}};
	
	function WLHistDos2DHalf( divNum )
		nDim = 2;
		grdNum = divNum^nDim;
		grdNumHalf = Int64(grdNum / 2);
		histArr = zeros( Int64, Int64( grdNum / 2) );
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr, Float64);
		dosArr .= 0;
		
		isHistOutBackupVal = false;
		
		posLst = CartesianIndices(histArr);
		
		new( histArr, histArrFlatBackup, dosArr, isHistOutBackupVal, grdNum, grdNumHalf, posLst );
	end
end

function getWLHistDosName( wlHistDosType::Type{<:WLHistDos2DHalf} )
	return "WLHist2DHalf";
end

function getAttrLst( wlHistDosType::Type{<:WLHistDos2DHalf} )
	return [ getWLHistDosName( wlHistDosType ) * "_grdNum" ]
end

function getValLst( wlHistDos::WLHistDos2DHalf )
	return [ wlHistDos.grdNum ];
end

function getHistUpdateId( wlHistDos::WLHistDos2DHalf, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D}
	dL = calcDL( params, linkLst, dim, pos );
	
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
	lId = getLinkHistId( wlHistDos, lTotal );
	lIdNxt = getLinkHistId( wlHistDos, lNxt );
	
	return wlHistDos.posLst[lId], wlHistDos.posLst[lIdNxt];
end

function getLinkHistId( wlHistDos::WLHistDos2DHalf, lnkVal::Int64 )
	lnkId = Int64( lnkVal / 2 );
	if lnkId > wlHistDos.grdNumHalf
		lnkId = wlHistDos.grdNum - lnkId;
	end
	if lnkId > 0
		lnkId -= 1;
	end
	
	lnkId += 1;
	
	return lnkId;
end




abstract type AbstractWLHistDosZoned <: AbstractWLHistDos end

function getLinkHistIdFull( grdNum::Int64, lnkVal::Int64 )
	lnkId = Int64( lnkVal / 2 );
	if lnkId == grdNum
		lnkId -= 2;
	elseif lnkId > 0
		lnkId -= 1;
	end
	
	lnkId += 1;
	
	return lnkId;
end

function eRatioToLnk01ValFlt( eRatio::Real, numLnkVal::Int64 )
	return ( eRatio/2 * numLnkVal + numLnkVal );
end

function roundLnkVal( lnkValFlt::Float64, numLnkVal::Int64, funRnd::Union{typeof(floor),typeof(ceil)} )
	return Int64( funRnd( lnkValFlt / 2 ) * 2 );
end

function eRatioToLnk01ValRnd( eRatio::Real, numLnkVal::Int64, funRnd::Union{typeof(floor),typeof(ceil)} )
	return roundLnkVal( eRatioToLnk01ValFlt( eRatio, numLnkVal ), numLnkVal, funRnd );
end

function histUpdateNoBackupSet( histDos::AbstractWLHistDosZoned, params::ParamsLoops, flipProposer::FlipProposer, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D}, dosIncr::Float64 ) where {D}
	# @infiltrate
	id, idNxt, isInBnd = getHistUpdateId( histDos, params, BfieldLst, linkLst, dim, pos );
	if isInBnd
		isFlip = histUpdateWithId( histDos, dosIncr, id, idNxt );
	else
		isFlip = false;
		histUpdateWithId( histDos, dosIncr, id );
	end
	# isFlip = isFlip && isInBnd;
	
	return isFlip;
end




abstract type AbstractWLHistDosZonedInE <: AbstractWLHistDosZoned end




struct WLHistDosZonedInE{nDim,D_hist} <: AbstractWLHistDosZonedInE
	histArr::Array{Int64,D_hist};
	histArrFlatBackup::Array{Int64,D_hist};
	dosArr::Array{Float64,D_hist};
	
	isHistOutBackupRef::Ref{Bool};
	
	numLnkVal::Int64;
	idMin::Int64;
	idMax::Int64;
	
	posLst::CartesianIndices{D_hist,Tuple{Vararg{Base.OneTo{Int64},D_hist}}};
	
	sSumDimLst::Vector{Int64};
	lSumDimLst::Vector{Int64};
	
	sTotalArr::Array{Int64};
	lTotalArr::Array{Int64};
	
	
	tmpArr::Matrix{Int64};
	
	function WLHistDosZonedInE{nDim,D_hist}( idMin::Int64, idMax::Int64, histArr::Array{Int64,D_hist}, numLnkVal::Int64 ) where {nDim,D_hist}
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr, Float64);
		dosArr .= 0;
		
		posLst = CartesianIndices( histArr );
		
		isHistFlatVal = false;
		
		nDimB = Int64( nDim * (nDim-1) / 2 );
		
		lSumDimLst = zeros(Int64, nDim);
		sSumDimLst = zeros(Int64, nDimB);
		
		sTotalArr = zeros(Int64,1);
		lTotalArr = zeros(Int64,1);
		
		tmpArr = zeros(Int64, 1, 1);
		
		new{nDim,D_hist}( histArr, histArrFlatBackup, dosArr, Ref(isHistFlatVal), numLnkVal, idMin, idMax, posLst, sSumDimLst, lSumDimLst, sTotalArr, lTotalArr, tmpArr );
	end
end

function getWLHistDosName( wlHistDosType::Type{<:WLHistDosZonedInE} )
	return "WLHistDosZonedInE" * "_" * string(wlHistDosType.parameters[1]) * "dDos";
end

function getAttrLst( wlHistDosType::Type{<:WLHistDosZonedInE} )
	return [ getWLHistDosName( wlHistDosType ) * "_grdNum", "idMin", "idMax" ];
end

function getValLst( wlHistDos::WLHistDosZonedInE )
	return [ wlHistDos.numLnkVal, wlHistDos.idMin, wlHistDos.idMax ];
end

function getGrdNum_NumLnkVal( typeHist::Type{<:WLHistDosZonedInE}, divNum::Int64 )
	throwWLHistDosUndefined();
end

function getIdMinMaxWLHist( histType::Type{<:WLHistDosZonedInE}, divNum::Int64, EMinRatio::Real, EMaxRatio::Real )
	_, numLnkVal = getGrdNum_NumLnkVal( histType, divNum );
	
	lnkValMin = eRatioToLnk01ValRnd( EMinRatio, numLnkVal, floor );
	lnkValMax = eRatioToLnk01ValRnd( EMaxRatio, numLnkVal, ceil );
	idMin = getLinkHistIdFull( numLnkVal, lnkValMin );
	idMax = getLinkHistIdFull( numLnkVal, lnkValMax );
	
	return idMin, idMax;
end

function WLHistDosZonedInE{nDim,D_hist}( divNum::Int64, EMinRatio::Real, EMaxRatio::Real ) where {nDim, D_hist}
	histType = WLHistDosZonedInE{nDim,D_hist};
	
	idMin, idMax = getIdMinMaxWLHist( histType, divNum, EMinRatio, EMaxRatio );
	
	histType( divNum, idMin, idMax );
end

function WLHistDosZonedInE{nDim,D_hist}( divNum::Int64, idMin::Int64, idMax::Int64 ) where {nDim, D_hist}
	throwWLHistDosUndefined();
end

function WLHistDosZonedInE{nDim,1}( divNum::Int64, idMin::Int64, idMax::Int64 ) where {nDim}
	histType = WLHistDosZonedInE{nDim,1};
	grdNum, numLnkVal = getGrdNum_NumLnkVal( histType, divNum );
	
	histArr = zeros(Int64, idMax - idMin + 1);
	
	histType( idMin, idMax, histArr, numLnkVal );
end

function WLHistDosZonedInE{nDim,2}( divNum::Int64, idMin::Int64, idMax::Int64 ) where {nDim}
	histType = WLHistDosZonedInE{nDim,2};
	grdNum, numLnkVal = getGrdNum_NumLnkVal( histType, divNum );
	
	histArr = zeros(Int64, grdNum+1, idMax - idMin + 1);
	
	histType( idMin, idMax, histArr, numLnkVal );
end





const WLHistDosZonedInE1dDos = WLHistDosZonedInE{2,1};
const WLHistDosZonedInE2dDos = WLHistDosZonedInE{2,2};

function getGrdNum_NumLnkVal( typeHist::Type{<:WLHistDosZonedInE{2}}, divNum::Int64 )
	nDim = 2;
	grdNum = divNum ^ nDim;
	numLnkVal = grdNum;
	
	return grdNum, numLnkVal;
end

function getHistUpdateId( wlHistDos::WLHistDosZonedInE{nDim,1}, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {nDim,D}
	dL = calcDL( params, linkLst, dim, pos );
	
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
	# lId = getLinkHistIdFull( wlHistDos.grdNum, lTotal );
	# lIdNxt = getLinkHistIdFull( wlHistDos.grdNum, lNxt );
	lId = getLinkHistIdFull( wlHistDos, lTotal );
	lIdNxt = getLinkHistIdFull( wlHistDos, lNxt );
	
	if testHistIdOutBnd( wlHistDos, lId ) 
		error( "Loops_MC.WLHistDos2DZoned: id out of bound" )
	end
	
	isInBnd = true;
	if testHistIdOutBnd( wlHistDos, lIdNxt )
		lIdNxt = lId;
		isInBnd = false;
	end
	
	lIdZoned = getZonedIdFromFull( wlHistDos, lId );
	lIdNxtZoned = getZonedIdFromFull( wlHistDos, lIdNxt );
	
	return wlHistDos.posLst[lIdZoned], wlHistDos.posLst[lIdNxtZoned], isInBnd;
end

function getAndTestHistThisId( wlHistDos::WLHistDosZonedInE{nDim,1}, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}} ) where {nDim,D}
	lTotal = sum( sum.(linkLst) );
	
	lId = getLinkHistIdFull( wlHistDos.grdNum, lTotal );
	
	isOutBnd = testHistIdOutBnd( wlHistDos, lId );
	if isOutBnd
		# error( "Loops_MC.WLHistDos2DZoned: id out of bound" )
		lIdZoned = 1;
	else
		lIdZoned = getZonedIdFromFull( wlHistDos, lId );
	end
	
	return wlHistDos.posLst[lIdZoned], isOutBnd;
end

function testHistIdOutBnd( wlHistDos::WLHistDosZonedInE, id::Int64 )
	if id < wlHistDos.idMin || id > wlHistDos.idMax
		return true;
	else
		return false;
	end
end

function getZonedIdFromFull( wlHistDos::WLHistDosZonedInE, idFull::Int64 )
	idZoned = idFull - wlHistDos.idMin + 1;
	
	return idZoned;
end

# function genWLHistDos2DZonedFull( divNum::Int64 )
	# return WLHistDos2DZoned( divNum, -2, 2 );
# end




function getHistUpdateId( wlHistDos::WLHistDosZonedInE{nDim,2}, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D,nDim}
	dL = calcDL( params, linkLst, dim, pos );
	
	# for iDim = 1 : nDim
		# sum!( wlHistDos.tmpArr, linkLst[iDim] );
		# wlHistDos.lSumDimLst[iDim] = wlHistDos.tmpArr[1];
	# end
	wlHistDos.lSumDimLst .= sum.(linkLst);
	# lTotal = sum( wlHistDos.lSumDimLst );
	sum!( wlHistDos.lTotalArr, wlHistDos.lSumDimLst );
	# lTotal = sum( sum.(linkLst) );
	# lTotal = sum( sum.(linkLst) );
	# wlHistDos.lTotalArr[1] = lTotal;
	# lNxt = lTotal + dL;
	lNxt = wlHistDos.lTotalArr[1] + dL;
	
	# lId = getLinkHistIdFull( wlHistDos.grdNum, lTotal );
	# lIdNxt = getLinkHistIdFull( wlHistDos.grdNum, lNxt );
	lId = getLinkHistIdFull( wlHistDos, wlHistDos.lTotalArr[1] );
	lIdNxt = getLinkHistIdFull( wlHistDos, lNxt );
	
	if testHistIdOutBnd( wlHistDos, lId ) 
		error( "Loops_MC.WLHistDos2DZoned: id out of bound" )
	end
	
	# for iDim = 1 : params.nDimB
		# sum!( wlHistDos.tmpArr, BfieldLst[iDim] );
		# wlHistDos.sSumDimLst[iDim] = wlHistDos.tmpArr[1];
	# end
	wlHistDos.sSumDimLst .= sum.(BfieldLst);
	sum!( wlHistDos.sTotalArr, wlHistDos.sSumDimLst );
	# sTotal = sum( wlHistDos.sSumDimLst );
	# sTotal = sum( sum.(BfieldLst) );
	dS = boolToFlipChange( BfieldLst[dim][pos] );
	sTotalNxt = wlHistDos.sTotalArr[1] + dS;
	
	sId = wlHistDos.sTotalArr[1] + 1;
	sIdNxt = sTotalNxt + 1;
	
	isInBnd = true;
	if testHistIdOutBnd( wlHistDos, lIdNxt )
		lIdNxt = lId;
		sIdNxt = sId;
		isInBnd = false;
	end
	
	lIdZoned = getZonedIdFromFull( wlHistDos, lId );
	lIdNxtZoned = getZonedIdFromFull( wlHistDos, lIdNxt );
	
	return wlHistDos.posLst[sId, lIdZoned], wlHistDos.posLst[sIdNxt, lIdNxtZoned], isInBnd;
end

function getAndTestHistThisId( wlHistDos::WLHistDosZonedInE{nDim,2}, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}} ) where {nDim,D}
	lTotal = sum( sum.(linkLst) );
	sTotal = sum( sum.(BfieldLst) );
	
	lId = getLinkHistIdFull( wlHistDos.grdNum, lTotal );
	sId = sTotal + 1;
	
	isOutBnd = testHistIdOutBnd( wlHistDos, lId );
	if isOutBnd
		lIdZoned = 1;
	else
		lIdZoned = getZonedIdFromFull( wlHistDos, lId );
	end
	
	return wlHistDos.posLst[sId,lIdZoned], isOutBnd;
end




const WL3dHistDosZonedInE1dDos = WLHistDosZonedInE{3,1};
const WL3dHistDosZonedInE2dDos = WLHistDosZonedInE{3,2};




function getGrdNum_NumLnkVal( histType::Type{<:WLHistDosZonedInE{3}}, divNum::Int64 )
	nDim = 3;
	grdNum = nDim * divNum ^ nDim;
	numLnkVal = Int64( floor( grdNum/2 ) );
	
	return grdNum, numLnkVal;
end





struct WLHistDos2DZoned <: AbstractWLHistDosZonedInE
	histArr::Vector{Int64};
	histArrFlatBackup::Vector{Int64};
	dosArr::Vector{Float64};
	
	isHistOutBackupRef::Ref{Bool};
	
	grdNum::Int64;
	idMin::Int64;
	idMax::Int64;
	
	posLst::CartesianIndices{1,Tuple{Vararg{Base.OneTo{Int64},1}}};
	
	function WLHistDos2DZoned( divNum::Int64, EMinRatio::Real, EMaxRatio::Real )
		nDim = 2;
		grdNum = divNum ^ nDim;
		lnkValMin = Int64( floor( eRatioToLnk01ValFlt( EMinRatio, grdNum )/2 ) )*2;
		lnkValMax = Int64( ceil( eRatioToLnk01ValFlt( EMaxRatio, grdNum )/2 ) )*2;
		idMin = getLinkHistIdFull( grdNum, lnkValMin );
		idMax = getLinkHistIdFull( grdNum, lnkValMax );
		
		histArr = zeros(Int64, idMax - idMin + 1);
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr, Float64);
		dosArr .= 0;
		
		posLst = CartesianIndices( histArr );
		
		isHistFlatVal = false;
		
		new( histArr, histArrFlatBackup, dosArr, Ref(isHistFlatVal), grdNum, idMin, idMax, posLst );
	end
end

function genWLHistDos2DZonedFull( divNum::Int64 )
	return WLHistDos2DZoned( divNum, -2, 2 );
end

function getWLHistDosName( wlHistDosType::Type{<:WLHistDos2DZoned} )
	return "WLHist2Zoned";
end

function getAttrLst( wlHistDosType::Type{<:WLHistDos2DZoned} )
	return [ getWLHistDosName( wlHistDosType ) * "_grdNum", "idMin", "idMax" ];
end

function getValLst( wlHistDos::WLHistDos2DZoned )
	return [ wlHistDos.grdNum, wlHistDos.idMin, wlHistDos.idMax ];
end

function getHistUpdateId( wlHistDos::WLHistDos2DZoned, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, dim::Int64, pos::CartesianIndex{D} ) where {D}
	dL = calcDL( params, linkLst, dim, pos );
	
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
	lId = getLinkHistIdFull( wlHistDos.grdNum, lTotal );
	lIdNxt = getLinkHistIdFull( wlHistDos.grdNum, lNxt );
	
	if testHistIdOutBnd( wlHistDos, lId ) 
		error( "Loops_MC.WLHistDos2DZoned: id out of bound" )
	end
	
	isInBnd = true;
	if testHistIdOutBnd( wlHistDos, lIdNxt )
		lIdNxt = lId;
		isInBnd = false;
	end
	
	lIdZoned = getZonedIdFromFull( wlHistDos, lId );
	lIdNxtZoned = getZonedIdFromFull( wlHistDos, lIdNxt );
	
	return wlHistDos.posLst[lIdZoned], wlHistDos.posLst[lIdNxtZoned], isInBnd;
end

function getAndTestHistThisId( wlHistDos::WLHistDos2DZoned, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}} ) where {D}
	lTotal = sum( sum.(linkLst) );
	
	lId = getLinkHistIdFull( wlHistDos.grdNum, lTotal );
	
	isOutBnd = testHistIdOutBnd( wlHistDos, lId );
	if isOutBnd
		# error( "Loops_MC.WLHistDos2DZoned: id out of bound" )
		lIdZoned = 1;
	else
		lIdZoned = getZonedIdFromFull( wlHistDos, lId );
	end
	
	return wlHistDos.posLst[lIdZoned], testHistIdOutBnd( wlHistDos, lId ) ;
end

function testHistIdOutBnd( wlHistDos::WLHistDos2DZoned, id::Int64 )
	if id < wlHistDos.idMin || id > wlHistDos.idMax
		return true;
	else
		return false;
	end
end

function getZonedIdFromFull( wlHistDos::WLHistDos2DZoned, idFull::Int64 )
	idZoned = idFull - wlHistDos.idMin + 1;
	
	return idZoned;
end





abstract type AbstractFlipCheckerWithProposer <: FlipChecker end

function flipCheckDoIt( flipChecker::AbstractFlipCheckerWithProposer, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	if flipCheck( flipChecker, flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst )
		flipDoIt( flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	end
end

function flipCheck( flipChecker::AbstractFlipCheckerWithProposer, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	error("Loops_MC: flipChecker or flipProposer not defined yet");
end




abstract type AbstractWangLandauFlipChecker <: AbstractFlipCheckerWithProposer end


struct WangLandauItController <: ItController
	dosIncrMin::Float64;
	flipChecker::AbstractWangLandauFlipChecker;
	
	itRef::Ref{Int64};
	
	function WangLandauItController( dosIncrMin::Float64, flipChecker::AbstractWangLandauFlipChecker )
		itVal = 1;
		
		new( dosIncrMin, flipChecker, Ref(itVal) );
	end
end

function getAttrLstItController( itControllerType::Type{WangLandauItController} )
	return ["dosIncrMin"];
end

function getValLstItController( itController::WangLandauItController )
	return [itController.dosIncrMin];
end

function getItNumLst( itController::WangLandauItController )
	return zeros(Int64, 3);
end

function testItNotDone( wlItController::WangLandauItController )
	return getDosIncr( wlItController.flipChecker ) > wlItController.dosIncrMin;
end

function testItDoSample( wlItController::WangLandauItController )
	return getIsHistFlat( wlItController.flipChecker );
end

testItDoStartSample( wlItController::WangLandauItController ) = false;



abstract type AbstractWLReplicaItController <: ItController end

function getItExchange( itController::AbstractWLReplicaItController )
	return itController.itExchange;
end

function getValLstItController( itController::AbstractWLReplicaItController )
	return [itController.dosIncrMin];
end

function getAttrLstItController( itControllerType::Type{<:AbstractWLReplicaItController} )
	return [getItControllerName(itControllerType) * "dosIncrMin"];
end

function testItNotDone( wlItController::AbstractWLReplicaItController )
	wlItController.dosIncrLst .= getDosIncr.( wlItController.flipCheckerLst );
	return maximum( wlItController.dosIncrLst ) > wlItController.dosIncrMin;
end

testItDoStartSample( wlItController::AbstractWLReplicaItController ) = false;

testItDoExchange( wlItController::AbstractWLReplicaItController ) = mod( wlItController.itRef[], wlItController.itExchange ) == 0;

testItDoExchange( jointItController::JointItController{<:Tuple{ItController,AbstractWLReplicaItController}} ) = testItDoExchange( jointItController.itControllerTup[2] );

getItExchange( jointItController::JointItController{<:Tuple{ItController,AbstractWLReplicaItController}} ) = getItExchange( jointItController.itControllerTup[2] );

struct WangLandauReplicaItController <: AbstractWLReplicaItController
	dosIncrMin::Float64;
	itExchange::Int64;
	flipCheckerLst::Array{<:AbstractWangLandauFlipChecker};
	
	dosIncrLst::Array{Float64};
	dosIncrMaxRef::Ref{Float64};
	
	itRef::Ref{Int64};
	
	function WangLandauReplicaItController( flipCheckerLst::Array{<:AbstractWangLandauFlipChecker}, dosIncrMin::Float64, itExchange::Int64 )
		itVal = 1;
		dosIncrLst = similar(flipCheckerLst, Float64)
		dosIncrMax = 1.0;
		
		new( dosIncrMin, itExchange, flipCheckerLst, dosIncrLst, Ref(dosIncrMax), Ref(itVal) );
	end
end

getItControllerName( itControllerType::Type{WangLandauReplicaItController} ) = "WLReplicaItCtrl";

function getItNumLst( itController::WangLandauReplicaItController )
	return zeros(Int64, 3);
end

function testItDoSample( wlItController::WangLandauReplicaItController )
	# return reduce( |, getIsHistFlat.( wlItController.flipCheckerLst ) );
	wlItController.dosIncrLst .= getDosIncr.( wlItController.flipCheckerLst );
	dosIncrMaxNow = maximum( wlItController.dosIncrLst );
	# isDoSample = wlItController.dosIncrMax > dosIncrMaxNow;
	if wlItController.dosIncrMaxRef[] > dosIncrMaxNow;
		wlItController.dosIncrMaxRef[] = dosIncrMaxNow;
		return true;
	else
		return false;
	end
	
	return false;
end




struct WangLandauReplicaFineSampleItController <: AbstractWLReplicaItController
	dosIncrMin::Float64;
	itExchange::Int64;
	flipCheckerLst::Array{<:AbstractWangLandauFlipChecker};
	
	dosIncrLst::Array{Float64};
	dosIncrPrevLst::Array{Float64};
	isDosIncrUpdatedLst::Array{Bool};
	dosIncrMaxRef::Ref{Float64};
	
	itRef::Ref{Int64};
	
	function WangLandauReplicaFineSampleItController( flipCheckerLst::Array{<:AbstractWangLandauFlipChecker}, dosIncrMin::Float64, itExchange::Int64 )
		itVal = 1;
		dosIncrLst = similar(flipCheckerLst, Float64)
		dosIncrPrevLst = similar(dosIncrLst);
		isDosIncrUpdatedLst = similar(dosIncrLst, Bool);
		dosIncrMax = 1.0;
		dosIncrPrevLst .= dosIncrMax;
		
		new( dosIncrMin, itExchange, flipCheckerLst, dosIncrLst, dosIncrPrevLst, isDosIncrUpdatedLst, Ref(dosIncrMax), Ref(itVal) );
	end
end

getItControllerName( itControllerType::Type{WangLandauReplicaFineSampleItController} ) = "WLReplicaItCtrlFineSample";

function testItDoSample( wlItController::WangLandauReplicaFineSampleItController )
	# return reduce( |, getIsHistFlat.( wlItController.flipCheckerLst ) );
	wlItController.dosIncrLst .= getDosIncr.( wlItController.flipCheckerLst );
	wlItController.isDosIncrUpdatedLst .= wlItController.dosIncrLst .< wlItController.dosIncrPrevLst;
	
	isAnyUpdated = reduce(|, wlItController.isDosIncrUpdatedLst);
	
	if isAnyUpdated
		wlItController.dosIncrPrevLst .= wlItController.dosIncrLst;
	end
	
	return isAnyUpdated;
end




abstract type ErrCheckItController <: ItController end

struct NanCheckWLItController <: ErrCheckItController
	flipCheckerLst::Array{<:AbstractWangLandauFlipChecker};
	# wlHistDosArr::Array{<:AbstractWLHistDos};
	# histArrArr::Array;
	# dosArrArr::Array;
	histDosJointArrArr::Vector{<:Array};
	itRef::Ref{Int64};
	arrIdLst::( CartesianIndices{N,NTuple{N,Base.OneTo{Int64}}} where {N} );
	
	isNanHistDosArr::Vector{<:Array{<:Array{Bool}}};
	isNanReducedHistDosArr::Vector{Array{Bool}};
	
	# function NanCheckWLItController( flipCheckerLst::Array{<:AbstractWangLandauFlipChecker}; )
		# itVal = 1;
		# histArrArr = getHistArr.( flipCheckerLst );
		# dosArrArr = getDosArr.( flipCheckerLst );
		
		# new( flipCheckerLst, histArrArr, dosArrArr, Ref(itVal) );
	# end
end

getAttrLstItController( itControllerType::Type{NanCheckWLItController} ) = String[];
getValLstItController( itController::NanCheckWLItController ) = [];

function NanCheckWLItController( flipCheckerLst::Array{<:AbstractWangLandauFlipChecker}; )
	itVal = 1;
	histArrArr = getHistArr.( flipCheckerLst );
	dosArrArr = getDosArr.( flipCheckerLst );
	histDosJointArrArr = [histArrArr, dosArrArr];
	arrIdLst = CartesianIndices(histArrArr);
	isNanHistDosArr = [ [ zeros(Bool, size(histArrArr[ii])) for ii in CartesianIndices(arrIdLst) ] for iHistDos = 1 : 2 ];
	isNanReducedHistDosArr = [ zeros(Bool, size(histArrArr)) for iHistDos = 1 : 2 ];
	
	NanCheckWLItController( flipCheckerLst, histDosJointArrArr, Ref(itVal), arrIdLst, isNanHistDosArr, isNanReducedHistDosArr );
end

# function NanCheckWLItController( flipCheckerArr::Array{<:AbstractWangLandauFlipChecker} )
	# wlHistDosArr = getHistDos.( flipCheckerArr );
	# isNanHistArr = zeros( Bool, size(histDosArr[1])..., size(histDosArr).. );
	# isNanDosArr = zeros( Bool, size(histDosArr) );
	
	# NanCheckWLIdController( histDosArr );
# end

function testItNotDone( itController::NanCheckWLItController )
	for iHistDos = 1 : 2
		for ii in itController.arrIdLst
			itController.isNanHistDosArr[iHistDos][ii] .= isnan.( itController.histDosJointArrArr[iHistDos][ii] );
			itController.isNanReducedHistDosArr[iHistDos][ii] = reduce( |, itController.isNanHistDosArr[iHistDos][ii] );
		end
	end
	
	isNanFound = reduce( |, reduce.(|, itController.isNanReducedHistDosArr) );
	return !isNanFound;
end

testItDoSample( itController::NanCheckWLItController ) = false;
testItDoStartSample( itController::NanCheckWLItController ) = false;
getItNumLst( itController::NanCheckWLItController ) = zeros(Int64, 3);
testItDoExchange( itController::NanCheckWLItController ) = false;





function flipCheck( flipChecker::AbstractWangLandauFlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	isFlip = wangLandauUpdateHistDos( flipChecker, flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	
	if wangLandauHistResetCheck( flipChecker )
		setIsHistFlat!( flipChecker );
		wangLandauUpdateDosIncr( flipChecker );
	else
		unsetIsHistFlat!(flipChecker);
	end
	
	return isFlip;
end

function getDosIncr( flipChecker::AbstractWangLandauFlipChecker )
	return flipChecker.dosIncrRef[];
end

function getIsHistFlat( flipChecker::AbstractWangLandauFlipChecker )
	return flipChecker.isHistFlatRef[];
end

function setIsHistFlat!( flipChecker::AbstractWangLandauFlipChecker )
	flipChecker.isHistFlatRef[] = true;
	flipChecker.histArrFlatBackup .= flipChecker.histArr;
end

function unsetIsHistFlat!( flipChecker::AbstractWangLandauFlipChecker )
	flipChecker.isHistFlatRef[] = false;
end

function getHistArr( flipChecker::AbstractWangLandauFlipChecker )
	if getIsHistFlat( flipChecker )
		return flipChecker.histArrFlatBackup;
	else
		return flipChecker.histArr;
	end
end

function getHistDos( flipChecker::AbstractWangLandauFlipChecker )
	return flipChecker.histDos;
end

function getDosArr( flipChecker::AbstractWangLandauFlipChecker )
	return flipChecker.dosArr;
end

function getHistDosArr( flipChecker::AbstractWangLandauFlipChecker )
	return [getHistArr( flipChecker ), getDosArr( flipChecker )];
end

function wangLandauUpdateHistDos( flipChecker::AbstractWangLandauFlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} )where {D}
	error("WangLandau Flipper or FlipProser not defined")
end

function wangLandauHistResetCheck( flipChecker::AbstractWangLandauFlipChecker )
	error("WangLandau Flipper not defined")
end

function wangLandauUpdateDosIncr( flipChecker::AbstractWangLandauFlipChecker )
	flipChecker.histArr .= 0
	flipChecker.dosIncrRef[] /= 2;
end


abstract type AbstractWLHistStructFlipChecker <: AbstractWangLandauFlipChecker end

function flipCheckExchange( flipChecker1::AbstractWLHistStructFlipChecker, flipChecker2::AbstractWLHistStructFlipChecker, params::ParamsLoops, bLinkData1::BLinkAuxData, bLinkData2::BLinkAuxData )
	isFlip = histUpdateExchange(  getHistDos( flipChecker1 ), getHistDos( flipChecker2 ), params, bLinkData1, bLinkData2, getDosIncr( flipChecker1 ), getDosIncr( flipChecker2 ) );
	@infiltrate isnothing(isFlip)
	if isFlip
		swapBLinkData!( bLinkData1, bLinkData2 );
	end
	
	return isFlip;
end

function swapBLinkData!( bLinkData1::BLinkAuxData, bLinkData2::BLinkAuxData )
	dataLstTmp = copy(bLinkData2.dataLst);
	bLinkData2.dataLst .= bLinkData1.dataLst;
	bLinkData1.dataLst .= dataLstTmp;
end

function flipCheck( flipChecker::AbstractWLHistStructFlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	isFlip = wangLandauUpdateHistDos( flipChecker, flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
	
	wangLandauHistReset( flipChecker );
	
	return isFlip;
end

function wangLandauUpdateHistDos( flipChecker::AbstractWLHistStructFlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} )where {D}
	return histUpdate( getHistDos( flipChecker ), params, flipProposer, BfieldLst, linkLst, dim, pos, flipChecker.dosIncrRef[] );
end

function wangLandauHistReset( flipChecker::AbstractWLHistStructFlipChecker )
	if wangLandauHistResetCheck( flipChecker )
		wangLandauHistResetNoCheck( flipChecker );
	end
end

function wangLandauHistResetNoCheck( flipChecker::AbstractWLHistStructFlipChecker )
	wangLandauUpdateDosIncr( flipChecker );
	histResetZero!( flipChecker.histDos );
end

function histResetZero!( flipChecker::AbstractWLHistStructFlipChecker );
	histResetZero!( flipChecker.histDos );
end

function getDosIncr( flipChecker::AbstractWLHistStructFlipChecker )
	return flipChecker.dosIncrRef[];
end

function getIsHistFlat( flipChecker::AbstractWLHistStructFlipChecker )
	return getIsHistOutBackup( flipChecker.histDos );
end

function setIsHistFlat!( flipChecker::AbstractWLHistStructFlipChecker )
	nothing;
end

function unsetIsHistFlat!( flipChecker::AbstractWLHistStructFlipChecker )
	nothing;
end

function getHistArr( flipChecker::AbstractWLHistStructFlipChecker )
	return getHistArr( flipChecker.histDos );
end

function getDosArr( flipChecker::AbstractWLHistStructFlipChecker )
	return getDosArr( flipChecker.histDos );
end

function getHistDosArr( flipChecker::AbstractWLHistStructFlipChecker )
	return [getHistArr( flipChecker ), getDosArr( flipChecker )];
end

# function wangLandauUpdateHistDos( flipChecker::AbstractWLHistStructFlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} )where {D}
	# error("WangLandau Flipper or FlipProser not defined")
# end

function wangLandauHistResetCheck( flipChecker::AbstractWLHistStructFlipChecker )
	error("WangLandau Flipper not defined")
end

function wangLandauUpdateDosIncr( flipChecker::AbstractWLHistStructFlipChecker )
	flipChecker.dosIncrRef[] /= 2;
end





struct WLFriendsFlipChecker{T_flip<:AbstractWangLandauFlipChecker} <: AbstractWLHistStructFlipChecker
	selfChecker::T_flip;
	friendCheckerLst::Vector{WLFriendsFlipChecker{T_flip}};
	
	function WLFriendsFlipChecker{T_flip}( args...; kwargs... ) where {T_flip <: AbstractWangLandauFlipChecker}
		selfChecker = T_flip(args...; kwargs...);
		friendCheckerLst = Vector{T_flip}(undef,0);
		
		new{T_flip}( selfChecker, friendCheckerLst );
	end
end

function setFriendCheckerLst!( friendChecker::WLFriendsFlipChecker, checkers... )
	resize!( friendChecker.friendCheckerLst, length(checkers) )
	for iCh = 1 : length(checkers)
		friendChecker.friendCheckerLst[iCh] = checkers[iCh];
	end
end

function getFlipCheckerName( flipCheckerType::Type{<:WLFriendsFlipChecker} )
	return getFlipCheckerName( flipCheckerType.parameters[1] );
end

function getFlipCheckerAttrLst( flipCheckerType::Type{<:WLFriendsFlipChecker} )
	return getFlipCheckerAttrLst( flipCheckerType.parameters[1] );
end

function getFlipCheckerAttrLst( flipChecker::WLFriendsFlipChecker )
	return getFlipCheckerAttrLst( typeof(flipChecker) );
end

function getFlipCheckerValLst( flipChecker::WLFriendsFlipChecker; rndDigs::Int64 )
	return getFlipCheckerValLst( flipChecker.selfChecker; rndDigs );
end

function getIsHistFlat( checker::WLFriendsFlipChecker )
	return reduce( &, ( c -> getIsHistFlat( c.selfChecker ) ).(checker.friendCheckerLst) );
end

function wangLandauUpdateHistDos( flipChecker::WLFriendsFlipChecker, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} )where {D}
	return wangLandauUpdateHistDos( flipChecker.selfChecker, flipProposer, params, dim, pos, BfieldLst, linkLst, linkFerroLst );
end

function wangLandauHistReset( flipChecker::WLFriendsFlipChecker )
	nothing;
end

function wangLandauHistResetNoCheck( flipChecker::WLFriendsFlipChecker )
	wangLandauHistResetNoCheck( flipChecker.selfChecker );
end

function wangLandauHistResetSynchronized( flipCheckerLst::AbstractVector{<:WLFriendsFlipChecker} )
	if wangLandauHistResetCheck( flipCheckerLst )
		# @time begin
		wangLandauHistResetNoCheck.( flipCheckerLst );
		avgDosAllFriends!( flipCheckerLst );
		# end
		
		# @infiltrate
	end
end

function wangLandauHistResetCheck( flipChecker::WLFriendsFlipChecker )
	return wangLandauHistResetCheck( flipChecker.friendCheckerLst );
end

function wangLandauHistResetCheck( flipCheckerLst::AbstractVector{<:WLFriendsFlipChecker} )
	return reduce( &, ( c -> wangLandauHistResetCheck( c.selfChecker ) ).(flipCheckerLst) );
end

function wangLandauUpdateDosIncr( flipChecker::WLFriendsFlipChecker )
	# avgDosFriends!( flipChecker );
	
	return wangLandauUpdateDosIncr( flipChecker.selfChecker );
end

function avgDosFriends!( flipChecker::WLFriendsFlipChecker, flipCheckerLst::AbstractVector{<:WLFriendsFlipChecker} )
	getDosArr( flipChecker ) .= getDosArr( flipCheckerLst[1] );
	
	for iF = 2 : length(flipCheckerLst)
		getDosArr( flipChecker ) .+= getDosArr( flipCheckerLst[iF] );
	end
	
	getDosArr( flipChecker ) ./= length(flipCheckerLst);
end

function avgDosFriends!( flipChecker::WLFriendsFlipChecker )
	avgDosFriends!( flipChecker, flipChecker.friendCheckerLst );
end

function avgDosAllFriends!( flipCheckerLst::AbstractVector{<:WLFriendsFlipChecker} )
	avgDosFriends!( flipCheckerLst[1], flipCheckerLst );
	
	for iF = 2 : length(flipCheckerLst)
		getDosArr( flipCheckerLst[iF] ) .= getDosArr( flipCheckerLst[1] );
	end
end

function avgDosAllFriends!( flipChecker::WLFriendsFlipChecker )
	avgDosFriends!( flipChecker );
	
	for iF = 1 : length(flipChecker.friendCheckerLst)
		getDosArr( flipChecker.friendCheckerLst[iF] ) .= getDosArr( flipChecker.selfChecker );
	end
end

function getHistDos( flipChecker::WLFriendsFlipChecker )
	return getHistDos( flipChecker.selfChecker );
end

function getDosIncr( flipChecker::WLFriendsFlipChecker )
	return getDosIncr( flipChecker.selfChecker );
end

function setIsHistFlat!( flipChecker::WLFriendsFlipChecker )
	return setIsHistFlat!( flipChecker.selfChecker );
end

function unsetIsHistFlat!( flipChecker::WLFriendsFlipChecker )
	return unsetIsHistFlat!( flipChecker.selfChecker );
end

function getHistArr( flipChecker::WLFriendsFlipChecker )
	return getHistArr( flipChecker.selfChecker );
end

function getDosArr( flipChecker::WLFriendsFlipChecker )
	return getDosArr( flipChecker.selfChecker );
end

function getHistDosArr( flipChecker::WLFriendsFlipChecker )
	return getHistDosArr( flipChecker.selfChecker );
end







struct WL2dHistStructFlipChecker{T_histDos} <: AbstractWLHistStructFlipChecker
	divNum::Int64;
	grdNum::Int64;
	
	histDos::T_histDos;
	dosIncrRef::Ref{Float64};
	
	histCutoffThres::Float64;
	histStdRatioThres::Float64;	
	
	histIsOverArr::Array{Bool};
	histIsOverArrPrev::Array{Bool};
	histMaskedArr::Array{Int64};
	histStdMaskedArr::Array{Float64};
	
	wlCounter::Ref{Int64};
	wlResetInterval::Int64;
	
	function WL2dHistStructFlipChecker{T_histDos}( divNum::Int64, nDim::Int64; dosIncrVal::Float64 = 1.0, histCutoffThres = 0.5, histStdRatioThres::Float64 = 0.15, wlResetInterval::Int64 = 100, wlHistDosArgs = () ) where {T_histDos <: AbstractWLHistDos}
		grdNum = divNum^nDim;
		histDos = T_histDos( divNum, wlHistDosArgs... );
		D_hist = ndims( getHistArr(histDos) );
		histIsOverArr = similar( getHistArr(histDos), Bool );
		histIsOverArrPrev = similar( histIsOverArr );
		histIsOverArr .= false;
		histIsOverArrPrev .= false;
		histMaskedArr = similar( getHistArr(histDos) );
		histStdMaskedArr = similar( histMaskedArr, Float64 );
		wlCounterVal = 1;
		
		new{T_histDos}( divNum, grdNum, histDos, Ref(dosIncrVal), histCutoffThres, histStdRatioThres, histIsOverArr, histIsOverArrPrev, histMaskedArr, histStdMaskedArr, Ref(wlCounterVal), wlResetInterval );
	end
end

function WL2dHistStructFlipChecker( divNum::Int64, nDim::Int64; dosIncrVal::Float64 = 1.0, histCutoffThres = 0.5, histStdRatioThres::Float64 = 0.15, wlResetInterval::Int64 = 100, wlHistDosType::Type{<:AbstractWLHistDos} = WLHistDos2DFull, wlHistDosArgs = () )
	return WL2dHistStructFlipChecker{wlHistDosType}( divNum, nDim; dosIncrVal = dosIncrVal, histCutoffThres = histCutoffThres, histStdRatioThres = histStdRatioThres, wlResetInterval = wlResetInterval, wlHistDosArgs = wlHistDosArgs );
end

function getFlipCheckerName( flipCheckerType::Type{<:WL2dHistStructFlipChecker} ) 
	return "WL2dExploredFlatFlip";
end

function getFlipCheckerAttrLst( flipCheckerType::Type{<:WL2dHistStructFlipChecker} )
	return append!( ["histStdRatioThres", "histCutoffThres", "wlResetIntvl"], getAttrLst( flipCheckerType.parameters[1] ) );
end

function getFlipCheckerValLst( flipChecker::WL2dHistStructFlipChecker; rndDigs = rndDigsLpsMC )
	return append!( roundKeepInt.( [flipChecker.histStdRatioThres, flipChecker.histCutoffThres, flipChecker.wlResetInterval]; digits = rndDigs ), getValLst(flipChecker.histDos) );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WL2dHistStructFlipChecker )
	isHistReset = false;
	flipChecker.wlCounter[] += 1;
	# @infiltrate typeof( getHistDos( flipChecker ) ) <: AbstractWLHistDosZoned
	if mod( flipChecker.wlCounter[], flipChecker.wlResetInterval ) == 0
		# @time begin
		maxCnt, cntId = findmax( getHistArr( flipChecker.histDos ) );
		cntCutoff = maxCnt * flipChecker.histCutoffThres;
		@strided flipChecker.histIsOverArr .= getHistArr( flipChecker.histDos ) .>= cntCutoff .|| flipChecker.histIsOverArrPrev;
		@strided flipChecker.histMaskedArr .= getHistArr( flipChecker.histDos ) .* flipChecker.histIsOverArr;
		
		cntOver = ThreadsX.sum(flipChecker.histIsOverArr);
		meanHistOver = ThreadsX.sum( flipChecker.histMaskedArr ) / cntOver;
		# stdHistOver = sqrt( sum( ( flipChecker.histMaskedArr .- meanHistOver ).^2 .* flipChecker.histIsOverArr ) / (cntOver - 1) );
		@strided flipChecker.histStdMaskedArr .= ( flipChecker.histMaskedArr .- meanHistOver ).^2 .* flipChecker.histIsOverArr;
		stdHistOver = sqrt( ThreadsX.sum( flipChecker.histStdMaskedArr ) / (cntOver - 1) );
		# stdHistOver = sqrt( sum( ( flipChecker.histMaskedArr .- meanHistOver ).^2 .* flipChecker.histIsOverArr ) / (cntOver - 1) );
		
		if stdHistOver < meanHistOver * flipChecker.histStdRatioThres
			isHistReset = true;
			flipChecker.histIsOverArrPrev .= flipChecker.histIsOverArr;
		end
		# end
		# @infiltrate typeof( getHistDos( flipChecker ) ).parameters[1] == 3 && typeof( getHistDos( flipChecker ) ).parameters[2] == 2
	end
	
	return isHistReset
end





abstract type AbstractWangLandau2dLinkOnlyFlipChecker <: AbstractWangLandauFlipChecker end

function getLinkHistId( flipChecker::AbstractWangLandau2dLinkOnlyFlipChecker, lnkVal::Int64 )
	# lnkId = Int64( ( lnkVal + 2*flipChecker.grdNum ) / 4 );
	lnkId = Int64( lnkVal / 2 );
	if lnkId == flipChecker.grdNum
		lnkId -= 2;
	elseif lnkId > 0
		lnkId -= 1;
	end
	
	lnkId += 1;
	
	return lnkId;
end

function wangLandauUpdateHistDos( flipChecker::AbstractWangLandau2dLinkOnlyFlipChecker, flipProposer::OneFlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} )where {D}
	dL = calcDL( params, linkLst, dim, pos );
	
	lTotal = sum( sum.(linkLst) );
	lNxt = lTotal + dL;
	
	lId = getLinkHistId( flipChecker, lTotal );
	lIdNxt = getLinkHistId( flipChecker, lNxt );
	
	p = exp( flipChecker.dosArr[lId] - flipChecker.dosArr[lIdNxt] );
	isFlip = rand() < p;
	
	if isFlip
		flipChecker.histArr[lIdNxt] += 1;
		flipChecker.dosArr[lIdNxt] += flipChecker.dosIncrRef[];
	else
		flipChecker.histArr[lId] += 1;
		flipChecker.dosArr[lId] += flipChecker.dosIncrRef[];
	end
	
	return isFlip;
end

struct WL2dLinkFlatFlipChecker <: AbstractWangLandau2dLinkOnlyFlipChecker
	divNum::Int64;
	grdNum::Int64;
	
	histArr::Vector{Int64};
	histArrFlatBackup::Vector{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	isHistFlatRef::Ref{Bool};
	
	histMinRatioThres::Float64;
	
	wlCounter::Ref{Int64};
	wlResetInterval::Int64;
	
	function WL2dLinkFlatFlipChecker( divNum::Int64, nDim::Int64; dosIncrVal::Float64 = 1.0, histMinRatioThres = 0.8, wlResetInterval::Int64 = 100 )
		grdNum = divNum^nDim;
		histArr = zeros( Int64, grdNum+1-2 );
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr);
		dosArr .= 0;
		wlCounterVal = 1;
		isHistFlatVal = false;
		
		new( divNum, grdNum, histArr, histArrFlatBackup, dosArr, Ref(dosIncrVal), Ref(isHistFlatVal), histMinRatioThres, Ref(wlCounterVal), wlResetInterval );
	end
end

function getFlipCheckerName( flipCheckerType::Type{WL2dLinkFlatFlipChecker} )
	return "WangLandau2dFlatHist";
end

function getFlipCheckerAttrLst( flipChecker::WL2dLinkFlatFlipChecker )
	return ["histMinRatio"];
end

function getFlipCheckerAttrValLst( flipChecker::WL2dLinkFlatFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.([flipChecker.histMinRatioThres]; digits = rndDigs);
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WL2dLinkFlatFlipChecker )
	isHistRest = false;
	flipChecker.wlCounter[] += 1;
	if mod( flipChecker.wlCounter[], flipChecker.wlResetInterval ) == 0
		meanHist = mean(flipChecker.histArr);
		minHist = minimum(flipChecker.histArr);
		
		if minHist > meanHist * flipChecker.histMinRatioThres
			isHistRest = true;
		end
	end
	
	return isHistRest;
end

struct WL2dLinkPartFlatFlipChecker <: AbstractWangLandau2dLinkOnlyFlipChecker
	divNum::Int64;
	grdNum::Int64;
	
	histArr::Vector{Int64};
	histArrFlatBackup::Vector{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	isHistFlatRef::Ref{Bool};
	
	histCutoffThres::Float64;
	histStdRatioThres::Float64;	
	
	histIsOverArr::Vector{Bool};
	histIsOverArrPrev::Vector{Bool};
	histMaskedArr::Vector{Int64};
	
	wlCounter::Ref{Int64};
	wlResetInterval::Int64;
	
	function WL2dLinkPartFlatFlipChecker( divNum::Int64, nDim::Int64; dosIncrVal::Float64 = 1.0, histCutoffThres = 0.5, histStdRatioThres::Float64 = 0.15, wlResetInterval::Int64 = 100 )
		grdNum = divNum^nDim;
		histArr = zeros( Int64, grdNum+1-2 );
		histArrFlatBackup = similar(histArr);
		dosArr = similar(histArr);
		dosArr .= 0;
		histIsOverArr = similar( histArr, Bool );
		histIsOverArrPrev = similar( histIsOverArr );
		histIsOverArr .= false;
		histIsOverArrPrev .= false;
		histMaskedArr = similar( histArr );
		wlCounterVal = 1;
		isHistFlatVal = false;
		
		new( divNum, grdNum, histArr, histArrFlatBackup, dosArr, Ref(dosIncrVal), Ref(isHistFlatVal), histCutoffThres, histStdRatioThres, histIsOverArr, histIsOverArrPrev, histMaskedArr, Ref(wlCounterVal), wlResetInterval );
	end
end

function getFlipCheckerName( flipCheckerType::Type{WL2dLinkPartFlatFlipChecker} ) 
	return "WL2dExploredFlatFlip"
end

function getFlipCheckerAttrLst( flipChecker::WL2dLinkPartFlatFlipChecker )
	return ["histStdRatioThres", "histCutoffThres", "wlResetIntvl"];
end

function getFlipCheckerAttrValLst( flipChecker::WL2dLinkPartFlatFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.( [flipChecker.histStdRatioThres, flipChecker.histCutoffThres, flipChecker.wlResetInterval]; digits = rndDigs );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WL2dLinkPartFlatFlipChecker )
	isHistReset = false;
	flipChecker.wlCounter[] += 1;
	if mod( flipChecker.wlCounter[], flipChecker.wlResetInterval ) == 0
		maxCnt, cntId = findmax( flipChecker.histArr );
		# maxDos = flipChecker.dosArr[ cntId ];
		cntCutoff = maxCnt * flipChecker.histCutoffThres;
		# dosCutoff = maxDos * flipChecker.histCutoffThres;
		# flipChecker.histIsOverArr .= flipChecker.histArr .>= cntCutoff .|| flipChecker.dosArr .>= dosCutoff;
		flipChecker.histIsOverArr .= flipChecker.histArr .>= cntCutoff .|| flipChecker.histIsOverArrPrev;
		flipChecker.histMaskedArr .= flipChecker.histArr .* flipChecker.histIsOverArr;
		
		cntOver = sum(flipChecker.histIsOverArr);
		meanHistOver = sum( flipChecker.histMaskedArr ) / cntOver;
		stdHistOver = sqrt.( sum( ( flipChecker.histMaskedArr .- meanHistOver ).^2 .* flipChecker.histIsOverArr ) / (cntOver - 1) );
		
		if stdHistOver < meanHistOver .* flipChecker.histStdRatioThres
			isHistReset = true;
			flipChecker.histIsOverArrPrev .= flipChecker.histIsOverArr;
		end
	end
	
	return isHistReset
end

abstract type AbstractWangLandau3dFlipChecker <: AbstractWangLandauFlipChecker end

function wangLandauUpdateHistDos( flipChecker::AbstractWangLandau3dFlipChecker, flipProposer::OneFlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	dS = boolToFlipChange( BfieldLst[dim][pos] );
	dL = calcDL( params, linkLst, dim, pos );
	
	sTotal = sum( sum.( BfieldLst ) );
	lTotal = sum( sum.( linkLst ) );
	
	sTotalNxt = sTotal + dS;
	lTotalNxt = lTotal + dL;
	
	sRatio = sTotal / params.nDim / params.grdNum;
	lRatio = lTotal / params.nDim / params.grdNum;
	
	sRatioNxt = sTotalNxt / params.nDim / params.grdNum;
	lRatioNxt = lTotalNxt / params.nDim / params.grdNum;
	
	sId = ratioToBinId( sRatio, flipChecker.histDivNum );
	lId = ratioToBinId( lRatio, flipChecker.histDivNum );
	sIdNxt = ratioToBinId( sRatioNxt, flipChecker.histDivNum );
	lIdNxt = ratioToBinId( lRatioNxt, flipChecker.histDivNum );
	
	p = exp( flipChecker.dosArr[sId,lId] - flipChecker.dosArr[sIdNxt,lIdNxt] );
	isFlip = rand() < p;
	if isFlip
		flipChecker.histArr[sIdNxt,lIdNxt] += 1;
		flipChecker.dosArr[sIdNxt,lIdNxt] += flipChecker.dosIncrRef[];
	else
		flipChecker.histArr[sId,lId] += 1;
		flipChecker.dosArr[sId,lId] += flipChecker.dosIncrRef[];
	end
	
	return isFlip;
end

struct WangLandauNoResetFlipChecker <: AbstractWangLandau3dFlipChecker
	histDivNum::Int64;
	histArr::Array{Int64};
	histArrFlatBackup::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	isHistFlatRef::Ref{Bool};
	
	function WangLandauNoResetFlipChecker( histDivNum::Int64; dosIncrVal::Float64 = 1.0 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		histArrFlatBackup = similar( histArr );
		dosArr = zeros(Int64, histDivNum, histDivNum);
		isHistFlatVal = false;
		
		new( histDivNum, histArr, histArrFlatBackup, dosArr, Ref(dosIncrVal), Ref(isHistFlatVal) );
	end
end

function getFlipCheckerName( flipType::Type{WangLandauNoResetFlipChecker} )
	return "WangLandauNoResetFlip";
end

function getFlipCheckerAttrLst( flipChecker::WangLandauNoResetFlipChecker )
	return ["histDiv"];
end

function getFlipCheckerAttrValLst( flipChecker::WangLandauNoResetFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst(flipChecker);
	valLst = roundKeepInt.( [flipChecker.histDivNum]; digits = rndDigs );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WangLandauNoResetFlipChecker )
	return false;
end

struct WangLandauFlipChecker <: AbstractWangLandau3dFlipChecker
	histDivNum::Int64;
	histArr::Array{Int64};
	histArrFlatBackup::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	histCoverThreshold::Float64;
	histLowRatioThreshold::Float64;
	isHistFlatRef::Ref{Bool};
	
	function WangLandauFlipChecker( histDivNum::Int64; histLowRatioThreshold::Float64 = 0.8, histCoverThreshold::Float64 = 0.8, dosIncrVal::Float64 = 1.0 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		histArrFlatBackup = similar(histArr);
		dosArr = zeros(Int64, histDivNum, histDivNum);
		
		isHistFlatVal = false;
		
		new( histDivNum, histArr, histArrFlatBackup, dosArr, Ref(dosIncrVal), Ref(isHistFlatRef), histCoverThreshold, histLowRatioThreshold );
	end
end

function getFlipCheckerName( flipType::Type{WangLandauFlipChecker} )
	return "WangLandauFlip";
end

function getFlipCheckerAttrLst( flipChecker::WangLandauFlipChecker )
	attrLst = ["histCoverThres", "histLowRatioThres"];
end

function getFlipCheckerAttrValLst( flipChecker::WangLandauFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.( [ flipChecker.histCoverThreshold, flipChecker.histLowRatioThreshold ]; digits = rndDigs );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WangLandauFlipChecker )
	isHistRest = false;
	coverNum = sum( x-> x>0, flipChecker.histArr );
	coverRatio = coverNum / flipChecker.histDivNum^2;
	if coverRatio > flipChecker.histCoverThreshold
		histAvg = sum( flipChecker.histArr ) / coverNum;
		histMin = minimum( x -> x == 0 ? Inf : x, flipChecker.histArr );
		isHistRest = histMin > flipChecker.histLowRatioThreshold * histAvg
	end
	
	return isHistRest;
end

function wangLandau_dSdL( flipChecker::WangLandauFlipChecker, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	dS = boolToFlipChange( BfieldLst[dim][pos] );
	dL = 0;
	for iLnkDim in 1 : params.nDimLayer
		dimLink = params.linkDimLst[dim][iLnkDim];
		dimLinkSh = params.linkDimShLst[dim][iLnkDim];
		dL += boolToFlipChange( linkLst[dimLink][pos] );
		dL += boolToFlipChange( linkLst[dimLink][params.posLstShLst[dimLinkSh,1][pos]] );
	end
	
	return dS, dL;
end

function wangLandau_flipCheck_dSdL( dS::Int64, dL::Int64, flipChecker::WangLandauFlipChecker, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	sTotal = sum( sum.( BfieldLst ) );
	lTotal = sum( sum.( linkLst ) );
	
	sTotalNxt = sTotal + dS;
	lTotalNxt = lTotal + dL;
	
	sRatio = sTotal / params.nDim / params.grdNum;
	lRatio = lTotal / params.nDim / params.grdNum;
	
	sRatioNxt = sTotalNxt / params.nDim / params.grdNum;
	lRatioNxt = lTotalNxt / params.nDim / params.grdNum;
	
	sId = ratioToBinId( sRatio, flipChecker.histDivNum );
	lId = ratioToBinId( lRatio, flipChecker.histDivNum );
	sIdNxt = ratioToBinId( sRatioNxt, flipChecker.histDivNum );
	lIdNxt = ratioToBinId( lRatioNxt, flipChecker.histDivNum );
	
	p = exp( flipChecker.dosArr[sId,lId] - flipChecker.dosArr[sIdNxt,lIdNxt] );
	isFlip = rand() < p;
	if isFlip
		flipChecker.histArr[sIdNxt,lIdNxt] += 1;
		flipChecker.dosArr[sIdNxt,lIdNxt] += flipChecker.dosIncrRef[];
	else
		flipChecker.histArr[sId,lId] += 1;
		flipChecker.dosArr[sId,lId] += flipChecker.dosIncrRef[];
	end
	
	return isFlip;
end

struct WangLandauStdFlipChecker <: AbstractWangLandau3dFlipChecker
	histDivNum::Int64;
	histArr::Array{Int64};
	histArrFlatBackup::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	isHistFlatRef::Ref{Bool};
	
	histMaxCntThres::Int64;
	histCntRatioCutoff::Float64;
	histStdRatioThres::Float64;
	histCoverCntThres::Int64;
	
	histIsOverThresArr::Array{Bool};
	histMaskedOverArr::Array{Int64};
	
	function WangLandauStdFlipChecker( histDivNum::Int64; dosIncrVal::Float64 = 1.0, histMaxCntThres::Int64 = 20, histCntRatioCutoff::Float64 = 0.5, histStdRatioThres::Float64 = 0.1, histCoverCntThres::Int64 = 5 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		histArrFlatBackup = similar(histArr);
		dosArr = zeros(Float64, histDivNum, histDivNum);
		
		histIsOverThresArr = similar( histArr, Bool );
		histMaskedOverArr = similar( histArr );
		
		isHistFlatVal = false;
		
		new( histDivNum, histArr, dosArr, Ref(dosIncrVal), Ref(isHistFlatVal), histMaxCntThres, histCntRatioCutoff, histStdRatioThres, histCoverCntThres, histIsOverThresArr, histMaskedOverArr );
	end
end

function getFlipCheckerName( flipCheckerType::Type{WangLandauStdFlipChecker} )
	return "WangLandauStdFlip";
end

function getFlipCheckerAttrLst( flipChecker::WangLandauStdFlipChecker )
	return ["histDiv", "histMaxCntThres","histCntRatioCut","histStdRatio","histCoverCntThres"];
end

function getFlipCheckerAttrValLst( flipChecker::WangLandauStdFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	
	valLst = roundKeepInt.( [ flipChecker.histDivNum, flipChecker.histMaxCntThres, flipChecker.histCntRatioCutoff, flipChecker.histStdRatioThres, flipChecker.histCoverCntThres ]; digits = rndDigs );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WangLandauStdFlipChecker )
	isReset = false;
	
	maxCnt, idMax = findmax( flipChecker.histArr );
	# maxCnt = flipChecker.histArr[idMax];
	maxDos = flipChecker.dosArr[idMax];
	if maxCnt > flipChecker.histMaxCntThres
		cntThres = maxCnt * flipChecker.histCntRatioCutoff;
		dosThres = maxDos * flipChecker.histCntRatioCutoff;
		flipChecker.histIsOverThresArr .= flipChecker.histArr .> cntThres .|| flipChecker.dosArr .> dosThres;
		
		cntOverThres = sum( flipChecker.histIsOverThresArr );
		if cntOverThres > flipChecker.histCoverCntThres
			flipChecker.histMaskedOverArr .= flipChecker.histArr .* flipChecker.histIsOverThresArr;
			histMean = sum( flipChecker.histMaskedOverArr ) / cntOverThres;
			histStd = sqrt.( sum( flipChecker.histIsOverThresArr.* ( flipChecker.histMaskedOverArr .- histMean ).^2 ) / (cntOverThres - 1) );
			flipChecker.dosIncrRef[] < 1
			if histStd < histMean * flipChecker.histStdRatioThres
				isReset = true;
				# @infiltrate
			end
		end
	end
	
	
	return isReset;
end

struct WL3dPartFlatStdFlipChecker <: AbstractWangLandau3dFlipChecker
	histDivNum::Int64
	histArr::Array{Int64};
	histArrFlatBackup::Array{Int64};
	dosArr::Array{Float64};
	dosIncrRef::Ref{Float64};
	isHistFlatRef::Ref{Bool};
	
	wlCounterRef::Ref{Int64};
	wlResetInterval::Int64;
	
	histStdRatioThres::Float64;
	histCutoffThres::Float64;
	
	histIsOverArr::Array{Bool};
	histIsOverArrPrev::Array{Bool};
	histMaskedArr::Array{Int64};
	
	function WL3dPartFlatStdFlipChecker( histDivNum::Int64; dosIncrVal::Float64 = 1.0, wlResetInterval = 1000, histStdRatioThres = 0.15, histCutoffThres = 0.5 )
		histArr = zeros(Int64, histDivNum, histDivNum);
		histArrFlatBackup = similar( histArr );
		dosArr = similar(histArr, Float64);
		dosArr .= 0;
		
		wlCounterVal = 1;
		
		histIsOverArr = similar(histArr, Bool);
		histIsOverArrPrev = similar(histArr, Bool);
		histIsOverArr .= false;
		histIsOverArrPrev .= false;
		histMaskedArr = similar(histArr);
		
		new( histDivNum, histArr, dosArr, Ref(dosIncrVal), Ref(isHistFlatVal), Ref(wlCounterVal), wlResetInterval, histStdRatioThres, histCutoffThres, histIsOverArr, histIsOverArrPrev, histMaskedArr );
	end
end

function getFlipCheckerName( flipCheckerType::Type{WL3dPartFlatStdFlipChecker} )
	return "WL3dPartFlatStdFlip";
end

function getFlipCheckerAttrLst( flipChecker::WL3dPartFlatStdFlipChecker )
	return ["histStdRatioThres", "histCutoff", "wlResetIntrvl"];
end

function getFlipCheckerAttrValLst( flipChecker::WL3dPartFlatStdFlipChecker; rndDigs = rndDigsLpsMC )
	attrLst = getFlipCheckerAttrLst( flipChecker );
	valLst = roundKeepInt.( [flipChecker.histStdRatioThres, flipChecker.histCutoffThres, flipChecker.wlResetInterval]; digits = rndDigsLpsMC );
	
	return attrLst, valLst;
end

function wangLandauHistResetCheck( flipChecker::WL3dPartFlatStdFlipChecker )
	isHistReset = false;
	flipChecker.wlCounterRef[] += 1;
	if mod( flipChecker.wlCounterRef[], flipChecker.wlResetInterval ) == 0
		maxCnt, idMaxCnt = findmax( flipChecker.histArr );
		histCutoff = maxCnt * flipChecker.histCutoffThres;
		flipChecker.histIsOverArr .= flipChecker.histArr .> histCutoff .|| flipChecker.histIsOverArrPrev;
		flipChecker.histMaskedArr .= flipChecker.histArr .* flipChecker.histIsOverArr;
		cntOver = sum(flipChecker.histIsOverArr);
		
		meanHistOver = sum( flipChecker.histMaskedArr ) / cntOver;
		stdHistOver = sqrt( sum( ( flipChecker.histMaskedArr .- meanHistOver ).^2 .* flipChecker.histIsOverArr ) ) / (cntOver - 1);
		
		if stdHistOver < meanHistOver * flipChecker.histStdRatioThres
			isHistReset = true;
			flipChecker.histIsOverArrPrev .= flipChecker.histIsOverArr;
		end
	end
	
	return isHistReset;
end





struct WangLandauAuxData <: AuxData
	itSampleLst::Vector{Int64};
	
	flipChecker::FlipChecker;
	
	dataLst::Vector{Array};
	dataSampleLst::Vector{<:Vector{<:Array}};
	dataStartSampleLst::Vector{<:Vector{<:Array}};
	dataNumLst::Vector{Vector};
	
	dataSampleOutLst::Vector{Array};
	dataStartSampleOutLst::Vector{Array};
	dataNumOutLst::Vector{Array};
	
	jldVarSampleLst::Vector{Any};
	jldVarStartSampleLst::Vector{Any};
	jldVarNumLst::Vector{Any};
	jldVarItSampleLst::Vector{Any};
	
	function WangLandauAuxData( flipChecker::AbstractWangLandauFlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 )
		itSampleLst = zeros(Int64, itNumSample);
		
		dataLst = Array[getHistArr(flipChecker), getDosArr(flipChecker)];
		nDimHist = ndims(getHistArr(flipChecker));
		dataSampleLst = [ Array{<:Real,nDimHist}[ similar( dataLst[ii] ) for itSample = 1 : itNumSample ] for ii = 1 : length(dataLst) ];
		# dataSampleLst = Vector{Array}[histArrSampleLst, dosArrSampleLst];
		# dataStartSampleLst = Vector{Array}[histArrStartSampleLst, dosArrStartSampleLst];
		dataStartSampleLst = [ Array{<:Real,nDimHist}[ similar( dataLst[ii] ) for itSample = 1 : itNumStartSample ] for ii = 1 : length(dataLst) ];
		dataNumLst = Vector{Vector}[];
		# dataNumLst = Vector{Vector}[ [] for ii = 1 : length(dataLst) ];
		
		dataSampleOutLst = Vector{Array}(undef, length(dataSampleLst));
		dataStartSampleOutLst = Vector{Array}(undef, length(dataStartSampleLst));
		dataNumOutLst = deepcopy(dataNumLst);
		
		jldVarSampleLst = Vector{Any}(undef,0);
		jldVarStartSampleLst = similar(jldVarSampleLst);
		jldVarNumLst = similar(jldVarSampleLst);
		jldVarItSampleLst = similar(jldVarSampleLst);
		
		new( itSampleLst, flipChecker, dataLst, dataSampleLst, dataStartSampleLst, dataNumLst, dataSampleOutLst, dataStartSampleOutLst, dataNumOutLst, jldVarSampleLst, jldVarStartSampleLst, jldVarNumLst, jldVarItSampleLst );
	end
end

WangLandauAuxData( flipChecker::AbstractWangLandauFlipChecker ) = WangLandauAuxData( flipChecker, 0, 0, 0 );
WangLandauAuxData( params::ParamsLoops, flipChecker::AbstractWangLandauFlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) = WangLandauAuxData( flipChecker, itNum, itNumSample, itNumStartSample );
WangLandauAuxData( params::ParamsLoops, flipChecker::AbstractWangLandauFlipChecker ) = WangLandauAuxData( flipChecker );
genAuxData( auxDataType::Type{WangLandauAuxData}, params::ParamsLoops, flipChecker::AbstractWangLandauFlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) = auxDataType( params, flipChecker, itNum, itNumSample, itNumStartSample );

function getAuxDataSummaryName( wlAuxDataType::Type{WangLandauAuxData} )
	return "WLHistDos";
end

function getAuxDataNameLst( wlAuxDataType::Type{WangLandauAuxData} )
	return ["histArr", "dosArr"];
end

getAuxDataNumNameLst( wlAuxDataType::Type{WangLandauAuxData} ) = String[];

function storeAuxDataSampleNoBndCheck( wlAuxData::WangLandauAuxData, itSample::Int64 )
	for ii = 1 : length(wlAuxData.dataLst)
		wlAuxData.dataSampleLst[ii][itSample] .= wlAuxData.dataLst[ii];
	end
end

function storeAuxDataStartSampleNoBndCheck( wlAuxData::WangLandauAuxData, itStartSample::Int64 )
	for ii = 1 : length(wlAuxData.dataLst)
		wlAuxData.dataSampleLst[ii][itStartSample] .= wlAuxData.dataLst[ii];
	end
end

storeAuxDataNum( auxData::WangLandauAuxData, it::Int64 ) = nothing;
storeAuxDataNumNoBndCheck( wlAuxData::WangLandauAuxData, it::Int64 ) = nothing;

function flipAuxData!( wlAuxData::WangLandauAuxData, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	nothing;
end

function calcAuxData!( wlAuxData::WangLandauAuxData, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	# nothing;
	wlAuxData.dataLst[1] = getHistArr( wlAuxData.flipChecker );
	wlAuxData.dataLst[2] = getDosArr( wlAuxData.flipChecker );
end

function renewAuxDataOutLst!( wlAuxData::WangLandauAuxData )
	funStackArrs = ( arr -> cat( arr...; dims = ndims(wlAuxData.dataLst[1])+1 ) );
	wlAuxData.dataSampleOutLst .= funStackArrs.(wlAuxData.dataSampleLst);
	wlAuxData.dataStartSampleOutLst .= funStackArrs.(wlAuxData.dataStartSampleLst);
	
	GC.gc();
end




struct BConfigOfLnkValAuxData{D} <: AuxData
	itSampleLst::Vector{Int64};
	
	configArr::Array{Vector{Array{Bool,D}}};
	
	flipChecker::FlipChecker;
	
	dataLst::Vector{Array};
	dataSampleLst::Vector{<:Vector{<:Array}};
	dataStartSampleLst::Vector{<:Vector{<:Array}};
	dataNumLst::Vector{Vector};
	
	dataSampleOutLst::Vector{Array};
	dataStartSampleOutLst::Vector{Array};
	dataNumOutLst::Vector{Array};
	
	jldVarSampleLst::Vector{Any};
	jldVarStartSampleLst::Vector{Any};
	jldVarNumLst::Vector{Any};
	jldVarItSampleLst::Vector{Any};
	
	foundCntRef::Ref{Int64};
	
	function BConfigOfLnkValAuxData{D}( flipChecker::AbstractWangLandauFlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) where {D}
		itSampleLst = zeros(Int64, itNumSample);
		
		configArr = similar( getHistArr(flipChecker), Vector{Array{Bool,D}} )
		dataLst = Array[configArr];
		dataSampleLst = [ Vector[ similar( dataLst[ii] ) for itSample = 1 : itNumSample ] for ii = 1 : length(dataLst) ];
		dataStartSampleLst = [ Vector[ similar( dataLst[ii] ) for itSample = 1 : itNumStartSample ] for ii = 1 : length(dataLst) ];
		# dataNumLst = Vector{Vector}[ [] for ii = 1 : length(dataLst) ];
		dataNumLst = Vector{Vector}[];
		
		dataSampleOutLst = dataSampleLst;
		dataStartSampleOutLst = dataStartSampleLst;
		dataNumOutLst = dataNumLst;
		
		jldVarSampleLst = Vector{Any}(undef,0);
		jldVarStartSampleLst = similar(jldVarSampleLst);
		jldVarNumLst = similar(jldVarSampleLst);
		jldVarItSampleLst = similar(jldVarSampleLst);
		
		foundCntVal = 0;
		
		new{D}( itSampleLst, configArr, flipChecker, dataLst, dataSampleLst, dataStartSampleLst, dataNumLst, dataSampleOutLst, dataStartSampleOutLst, dataNumOutLst, jldVarSampleLst, jldVarStartSampleLst, jldVarNumLst, jldVarItSampleLst, Ref(foundCntVal) );
	end
end

BConfigOfLnkValAuxData{D}( flipChecker::AbstractWangLandauFlipChecker ) where {D} = BConfigOfLnkValAuxData{D}( flipChecker, 0,0,0 );

function getAuxDataSummaryName( wlAuxDataType::Type{<:BConfigOfLnkValAuxData} )
	return "BConfigOfLnkVal";
end

function getAuxDataNameLst( wlAuxDataType::Type{<:BConfigOfLnkValAuxData} )
	return ["BfieldOfLnkVal"];
end

getAuxDataNumNameLst( wlAuxDataType::Type{<:BConfigOfLnkValAuxData} ) = String[];

function flipAuxData!( wlAuxData::BConfigOfLnkValAuxData{D}, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	nothing;
end

function calcAuxData!( wlAuxData::BConfigOfLnkValAuxData{D}, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	lnkVal = sum( sum.(linkLst) );
	lId = getLinkHistIdFull( params.grdNum, lnkVal );
	if !isassigned( wlAuxData.configArr, lId )
		wlAuxData.configArr[lId] = deepcopy( BfieldLst );
		wlAuxData.foundCntRef[] += 1;
	end
end

function getIsConfigAllFnd( wlAuxData::BConfigOfLnkValAuxData )
	# isFnd = true;
	for it = 1 : length(wlAuxData.configArr)
		if !isassigned( wlAuxData.configArr[it] )
			return false;
		end
	end
	# return wlAuxData.foundCntRef[] >= length(wlAuxData.configArr);
	return true;
end

function renewAuxDataOutLst!( wlAuxData::BConfigOfLnkValAuxData )
	nothing;
end

storeAuxDataNum( auxData::BConfigOfLnkValAuxData, it::Int64 ) = nothing;



abstract type AbstractConfigAuxData <: AuxData end

throwConfigAuxDataUndefined() = error( "Loops_MC: config AuxData undefined" );

getIsConfigAllFnd( wlAuxData::AbstractConfigAuxData ) = throwConfigAuxDataUndefined();





# struct BConfigOfLnkValAuxData{D} <: AbstractConfigAuxData
	# itSampleLst::Vector{Int64};
	
	# configArr::Array{Vector{Array{Bool,D}}};
	
	# flipChecker::FlipChecker;
	
	# dataLst::Vector{Array};
	# dataSampleLst::Vector{<:Vector{<:Array}};
	# dataStartSampleLst::Vector{<:Vector{<:Array}};
	# dataNumLst::Vector{Vector};
	
	# dataSampleOutLst::Vector{Array};
	# dataStartSampleOutLst::Vector{Array};
	# dataNumOutLst::Vector{Array};
	
	# jldVarSampleLst::Vector{Any};
	# jldVarStartSampleLst::Vector{Any};
	# jldVarNumLst::Vector{Any};
	# jldVarItSampleLst::Vector{Any};
	
	# foundCntRef::Ref{Int64};
	
	# function BConfigOfLnkValAuxData{D}( flipChecker::AbstractWangLandauFlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) where {D}
		# itSampleLst = zeros(Int64, itNumSample);
		
		# configArr = similar( getHistArr(flipChecker), Vector{Array{Bool,D}} )
		# dataLst = Array[configArr];
		# dataSampleLst = [ Vector[ similar( dataLst[ii] ) for itSample = 1 : itNumSample ] for ii = 1 : length(dataLst) ];
		# dataStartSampleLst = [ Vector[ similar( dataLst[ii] ) for itSample = 1 : itNumStartSample ] for ii = 1 : length(dataLst) ];
		# dataNumLst = Vector{Vector}[];
		
		# dataSampleOutLst = dataSampleLst;
		# dataStartSampleOutLst = dataStartSampleLst;
		# dataNumOutLst = dataNumLst;
		
		# jldVarSampleLst = Vector{Any}(undef,0);
		# jldVarStartSampleLst = similar(jldVarSampleLst);
		# jldVarNumLst = similar(jldVarSampleLst);
		# jldVarItSampleLst = similar(jldVarSampleLst);
		
		# foundCntVal = 0;
		
		# new{D}( itSampleLst, configArr, flipChecker, dataLst, dataSampleLst, dataStartSampleLst, dataNumLst, dataSampleOutLst, dataStartSampleOutLst, dataNumOutLst, jldVarSampleLst, jldVarStartSampleLst, jldVarNumLst, jldVarItSampleLst, Ref(foundCntVal) );
	# end
# end

# BConfigOfLnkValAuxData{D}( flipChecker::AbstractWangLandauFlipChecker ) where {D} = BConfigOfLnkValAuxData{D}( flipChecker, 0,0 ,0 );

# function getAuxDataSummaryName( wlAuxDataType::Type{<:BConfigOfLnkValAuxData} )
	# return "BConfigOfLnkVal";
# end

# function getAuxDataNameLst( wlAuxDataType::Type{<:BConfigOfLnkValAuxData} )
	# return ["BfieldOfLnkVal"];
# end

# getAuxDataNumNameLst( wlAuxDataType::Type{<:BConfigOfLnkValAuxData} ) = String[];

# function flipAuxData!( wlAuxData::BConfigOfLnkValAuxData{D}, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	# nothing;
# end

# function calcAuxData!( wlAuxData::BConfigOfLnkValAuxData{D}, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	# lnkVal = sum( sum.(linkLst) );
	# lId = getLinkHistIdFull( params.grdNum, lnkVal );
	# if !isassigned( wlAuxData.configArr, lId )
		# wlAuxData.configArr[lId] = deepcopy( BfieldLst );
		# wlAuxData.foundCntRef[] += 1;
	# end
# end

# function getIsConfigAllFnd( wlAuxData::BConfigOfLnkValAuxData )
	# return wlAuxData.foundCntRef[] >= length(wlAuxData.configArr);
# end

# function renewAuxDataOutLst!( wlAuxData::BConfigOfLnkValAuxData )
	# nothing;
# end




struct BConfigOfLnkValZonedAuxData{D} <: AbstractConfigAuxData
	itSampleLst::Vector{Int64};
	
	configArr::Array{Vector{Array{Bool,D}}};
	
	# flipChecker::FlipChecker;
	
	dataLst::Vector{Array};
	dataSampleLst::Vector{<:Vector{<:Array}};
	dataStartSampleLst::Vector{<:Vector{<:Array}};
	dataNumLst::Vector{Vector};
	
	dataSampleOutLst::Vector{Array};
	dataStartSampleOutLst::Vector{Array};
	dataNumOutLst::Vector{Array};
	
	jldVarSampleLst::Vector{Any};
	jldVarStartSampleLst::Vector{Any};
	jldVarNumLst::Vector{Any};
	jldVarItSampleLst::Vector{Any};
	
	foundCntRef::Ref{Int64};
	
	idMinLst::Vector{Int64};
	idMaxLst::Vector{Int64};
	
	function BConfigOfLnkValZonedAuxData{D}( idMinLst::Vector{Int64}, idMaxLst::Vector{Int64}, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) where {D}
		itSampleLst = zeros(Int64, itNumSample);
		
		configArr = similar( idMinLst, Vector{Array{Bool,D}} )
		dataLst = Array[configArr];
		dataSampleLst = [ Vector[ similar( dataLst[ii] ) for itSample = 1 : itNumSample ] for ii = 1 : length(dataLst) ];
		dataStartSampleLst = [ Vector[ similar( dataLst[ii] ) for itSample = 1 : itNumStartSample ] for ii = 1 : length(dataLst) ];
		dataNumLst = Vector[];
		# dataNumLst = Vector{Vector}[ [] for ii = 1 : length(dataLst) ];
		
		dataSampleOutLst = dataSampleLst;
		dataStartSampleOutLst = dataStartSampleLst;
		dataNumOutLst = deepcopy( dataNumLst );
		
		jldVarSampleLst = Vector{Any}(undef,0);
		jldVarStartSampleLst = similar(jldVarSampleLst);
		jldVarNumLst = similar(jldVarSampleLst);
		jldVarItSampleLst = similar(jldVarSampleLst);
		
		foundCntVal = 0;		
		
		new{D}( itSampleLst, configArr, dataLst, dataSampleLst, dataStartSampleLst, dataNumLst, dataSampleOutLst, dataStartSampleOutLst, dataNumOutLst, jldVarSampleLst, jldVarStartSampleLst, jldVarNumLst, jldVarItSampleLst, Ref(foundCntVal), idMinLst, idMaxLst );
	end
end

function BConfigOfLnkValZonedAuxData{D}( histDosLst::AbstractWLHistDosZoned, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) where {D}
	idMinLst = (h -> h.idMin).( histDosLst );
	idMaxLst = (h -> h.idMax).( histDosLst );
	
	BConfigOfLnkValZonedAuxData{D}( idMinLst, idMaxLst, itNum, itNumSample, itNumStartSample )
end

BConfigOfLnkValZonedAuxData{D}( histDosLst::AbstractWLHistDosZoned ) where {D} = BConfigOfLnkValAuxData{D}( histDosLst, 0,0 ,0 );

function getAuxDataSummaryName( wlAuxDataType::Type{<:BConfigOfLnkValZonedAuxData} )
	return "BConfigOfLnkValZoned";
end

function getAuxDataNameLst( wlAuxDataType::Type{<:BConfigOfLnkValZonedAuxData} )
	return ["BfieldOfLnkVal"];
end

getAuxDataNumNameLst( wlAuxDataType::Type{<:BConfigOfLnkValZonedAuxData} ) = String[];

function flipAuxData!( wlAuxData::BConfigOfLnkValZonedAuxData{D}, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	nothing;
end

function calcAuxData!( wlAuxData::BConfigOfLnkValZonedAuxData{D}, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	lnkVal = sum( sum.(linkLst) );
	_, numLnkVal = getGrdNum_NumLnkVal( WLHistDosZonedInE{params.nDim,1}, params.divNum );
	lId = getLinkHistIdFull( numLnkVal, lnkVal );
	idHistHigh = searchsortedlast( wlAuxData.idMinLst, lId );
	idHistLow = searchsortedfirst( wlAuxData.idMaxLst, lId );
	for idHist = idHistLow : idHistHigh
		if !isassigned( wlAuxData.configArr, idHist )
			wlAuxData.configArr[idHist] = deepcopy( BfieldLst );
			wlAuxData.foundCntRef[] += 1;
		end
	end
end

function getIsConfigAllFnd( wlAuxData::BConfigOfLnkValZonedAuxData )
	for it = 1 : length(wlAuxData.configArr)
		if !isassigned( wlAuxData.configArr, it )
			return false;
		end
	end
	return true;
	# return wlAuxData.foundCntRef[] >= length(wlAuxData.configArr);
end

function renewAuxDataOutLst!( wlAuxData::BConfigOfLnkValZonedAuxData )
	nothing;
end

storeAuxDataSampleNoBndCheck( auxData::BConfigOfLnkValZonedAuxData, itSample::Int64 ) = nothing;
storeAuxDataStartSampleNoBndCheck( auxData::BConfigOfLnkValZonedAuxData, itStartSample::Int64 ) = nothing;
storeAuxDataNum( auxData::BConfigOfLnkValZonedAuxData, it::Int64 ) = nothing;
storeAuxDataNumNoBndCheck( auxData::BConfigOfLnkValZonedAuxData, it::Int64 ) = nothing;





struct FindConfigItController <: ItController
	configAuxData::AbstractConfigAuxData;
	
	itRef::Ref{Int64}
	
	function FindConfigItController( configAuxData::AbstractConfigAuxData )
		itVal = 1;
		
		new(configAuxData, Ref(itVal));
	end
end

getItControllerName( itControllerType::Type{FindConfigItController} ) = "FindConfigController";

function getAttrLstItController( itControllerType::Type{FindConfigItController} )
	return [getItControllerName(itControllerType)];
end

function getValLstItController( itController::FindConfigItController )
	return [length( itController.configAuxData.configArr )];
end

function testItNotDone( wlItController::FindConfigItController )
	return !getIsConfigAllFnd( wlItController.configAuxData );
end

function getItNumLst( wlItController::FindConfigItController )
	return zeros(Int64, 3);
end

testItDoSample( wlItController::FindConfigItController ) = false;

testItDoStartSample( wlItController::FindConfigItController ) = false;




struct FindConfigReplicaItController <: ItController
	findConfigItController::FindConfigItController;
	wlReplicaItController::WangLandauReplicaItController;
	
	itRef::Ref{Int64};
	
	function FindConfigReplicaItController( configAuxData::AbstractConfigAuxData, flipCheckerLst::Array{<:AbstractWangLandauFlipChecker}, itExchange::Int64 )
		findConfigItCtrl = FindConfigItController( configAuxData );
		dosIncrMin = 0.1;
		wlReplicaItCtrl = WangLandauReplicaItController( flipCheckerLst, dosIncrMin, itExchange );
		
		itVal = 1;
		
		new( findConfigItCtrl, wlReplicaItCtrl, Ref(itVal) );
	end
end

getItControllerName( itControllerType::Type{FindConfigReplicaItController} ) = "FindConfigREplicaController";

function getAttrLstItController( itControllerType::Type{FindConfigReplicaItController} )
	return append!( getAttrLstItController( FindConfigItController ), getAttrLstItController( WangLandauReplicaItController ) );
end

function getValLstItController( itController::FindConfigReplicaItController )
	return append!( Vector{Any}( getValLstItController( itController.findConfigItController ) ), getValLstItController( itController.wlReplicaItController ) );
end

function testItNotDone( wlItController::FindConfigReplicaItController )
	return testItNotDone( wlItController.findConfigItController );
end

function getItNumLst( wlItController::FindConfigReplicaItController )
	return zeros(Int64, 3);
end

function getItExchange( wlItController::FindConfigReplicaItController )
	return getItExchange( wlItController.wlReplicaItController );
end

testItDoSample( wlItController::FindConfigReplicaItController ) = false;

testItDoStartSample( wlItController::FindConfigReplicaItController ) = false;

testItDoExchange( wlItController::FindConfigReplicaItController ) = testItDoExchange( wlItController.wlReplicaItController );




struct StaggeredCubeUpdaterWangLandau{N,Nplus1} <: LoopsUpdater
	posLstSh0::CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}};
	posLstAdvOrNot::Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}} }};
	posShOrNotLst::Vector{Vector{CircShiftedArray{CartesianIndex{N}, N, CartesianIndices{N,Tuple{Vararg{Base.OneTo{Int64},N}}}}}};
	posStagCubeLst::Vector{Array{CartesianIndex{N},Nplus1}};
	idStagLst::CartesianIndices{Nplus1,NTuple{Nplus1,Base.OneTo{Int64}}};
	randIsFlipStagLst::Array{Bool};
	
	iDimLst::UnitRange{Int64};
	iIsShLst::UnitRange{Int64};
	randIDimLst::Array{Int64,Nplus1};
	randIShLst::Array{Int64,Nplus1};
	
	function StaggeredCubeUpdaterWangLandau( params::ParamsLoops; )
		posLstSh0 = ShiftedArrays.circshift( params.posLst, ntuple(x->0,params.nDim) );
		posLstAdvanced = [ ShiftedArrays.circshift( params.posLst, ntuple( ( i -> ( i == dim ? -1 : 0 ) ), params.nDim ) ) for dim = 1 : params.nDim ];
		posLstAdvOrNot = pushfirst!( posLstAdvanced, ShiftedArrays.circshift( params.posLst, ntuple(x->0,params.nDim) ) );
		posShOrNotLst = [ [ posLstSh0, posLstAdvanced[dim] ] for dim = 1 : params.nDim ];
		posShOrNotArr = [ (iSh == 1 ? posLstSh0 : posLstAdvanced[dim]) for dim = 1 : params.nDim, iSh = 1:2 ];
		
		coordsStagA = ntuple( iDim -> 1:2:params.divLst[iDim], params.nDim );
		coordsStagB = ntuple( iDim -> 2:2:params.divLst[iDim], params.nDim );
		posStagCubeLstA = @view( params.posLst[coordsStagA...] );
		posStagCubeLstB = @view( params.posLst[coordsStagB...] );
		posStagCubeLst = [ cat( @view( posLstAdvOrNot[iAdv][coordsStagA...]), @view(posLstAdvOrNot[iAdv][coordsStagB...]); dims = params.nDim+1 ) for iAdv = 1 : params.nDim+1 ];
		idStagLst = CartesianIndices( posStagCubeLst[1] );
		randIsFlipStagLst = similar( idStagLst, Bool );
		
		iDimLst = 1:params.nDim;
		iIsShLst = 1:2;
		randIDimLst = similar( posStagCubeLst[1], Int64 );
		randIShLst = similar(randIDimLst);
		
		new{params.nDim,params.nDim+1}( posLstSh0, posLstAdvOrNot, posShOrNotLst, posStagCubeLst, idStagLst, randIsFlipStagLst, iDimLst, iIsShLst, randIDimLst, randIShLst );
	end
end

function getUpdaterFMod( updaterType::Type{StaggeredCubeUpdaterWangLandau} )
	return "upStagCube";
end

function updateLoops( updater::StaggeredCubeUpdaterWangLandau, flipChecker::WangLandauFlipChecker, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}}, params::ParamsLoops ) where {D}
	for iAdv = 1 : params.nDim+1
		rand!( updater.randIDimLst, updater.iDimLst );
		rand!( updater.randIShLst, updater.iIsShLst );
		rand!( updater.randIsFlipStagLst );
		dSSumRef = Threads.Atomic{Int}(0);
		dLSumRef = Threads.Atomic{Int}(0);
		Threads.@threads for idStag in updater.idStagLst
			posCube = updater.posStagCubeLst[iAdv][idStag];
			pos = updater.posShOrNotLst[updater.randIDimLst[idStag]][updater.randIShLst[idStag]][posCube];
			
			if updater.randIsFlipStagLst[idStag]
				dS, dL = wangLandau_dSdL( flipChecker, params, updater.randIDimLst[idStag], pos, BfieldLst, linkLst, linkFerroLst );
				Threads.atomic_add!( dSSumRef, dS );
				Threads.atomic_add!( dLSumRef, dL );
			end
		end
		if wangLandau_flipCheck_dSdL( dSSumRef[], dLSumRef[], flipChecker, params, BfieldLst, linkLst, linkFerroLst )
			Threads.@threads for idStag in updater.idStagLst
				if updater.randIsFlipStagLst[idStag]
					posCube = updater.posStagCubeLst[iAdv][idStag];
					pos = updater.posShOrNotLst[updater.randIDimLst[idStag]][updater.randIShLst[idStag]][posCube];
					
					flipBLinkAtPos( params, BfieldLst, linkLst, linkFerroLst; pos = pos, dim = updater.randIDimLst[idStag] );
				end
			end
		end
	end
end


function loops_MC_methods_WangLandauStaggered( divNum = 64, itNum = 10000; histDivNum = 64, itNumSample = 100, itStartSample = 50, isInit0 = false, dosIncrInit = 1, cAreaInit = 0 )
	flipChecker = WangLandauFlipChecker( histDivNum; dosIncrVal = dosIncrInit );
	updaterType = StaggeredCubeUpdaterWangLandau;
	
	numBfieldLst, numLinkLst, BfieldSampleLst, linkSampleLst, BfieldStartSampleLst, linkStartSampleLst, zakMeanLst, zakLstSampleLst = loops_MC_methods_inJulia( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, itNumSample = itNumSample, itStartSample = itStartSample, isInit0 = isInit0, cAreaInit = cAreaInit );
	
	fMain = "loops_WL_staggered";
	attrLst = ["divNum", "itNum", "histDivNum", "cAreaInit", "dosIncrInit"];
	valLst = [divNum, itNum, histDivNum, cAreaInit, dosIncrInit];
	fNameWL = fNameFunc( fMain, attrLst, valLst, jld2Type );
	
	save( fNameWL, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	
	return fNameWL;
end

function loops_MC_methods_WangLandau( divNum = 64, itNum = 10000; histDivNum = 64, itNumSample = 100, itStartSample = 50, isInit0 = false, cAreaInit = 0, dosIncrInit = 1, nDim = 3 )
	# flipChecker = WangLandauFlipChecker( histDivNum; dosIncrVal = dosIncrInit );
	# flipChecker = WangLandauNoResetFlipChecker( histDivNum );
	# flipChecker = WangLandauStdFlipChecker( histDivNum; histMaxCntThres = 10, histStdRatioThres = 0.15 );
	flipChecker = WL3dPartFlatStdFlipChecker( histDivNum; wlResetInterval = 1000, histStdRatioThres = 0.15 );
	flipProposer = OneFlipProposer();
	initializer = genMeanFieldInitializer( cAreaInit );
	updaterType = SingleUpdater;
	
	fName = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, itNumSample = itNumSample, itStartSample = itStartSample );
	
	fMain = "loops_WL_single";
	attrLst, valLst = genAttrLstLttcFlipInit( divNum, itNum, nDim, flipChecker, initializer );
	fNameWL = fNameFunc( fMain, attrLst, valLst, jld2Type );
	
	println( "f = ", flipChecker.dosIncrRef[] );
	
	save( fNameWL, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	
	open( dirLog * fNameFileLstWL, "w" ) do io
		println( io, fNameWL );
	end
	
	return fNameWL;
end

function loops_MC_methods_WL2d( divNum, itNum; itNumSample = 100, itStartSample = 50, isInit0 = false, cAreaInit = 0, dosIncrInit = 1, dosIncrMin = 0.01, nDim = 2, isFileNameOnly = false, fMainOutside = "", wlHistDosType = WLHistDos2DFull, wlHistDosArgs = () )
	# flipChecker = WangLandauFlipChecker( histDivNum; dosIncrVal = dosIncrInit );
	# flipChecker = WangLandauNoResetFlipChecker( histDivNum );
	# flipChecker = WL2dLinkFlatFlipChecker( divNum, nDim; histMinRatioThres = 0.8 );
	flipChecker = WL2dLinkPartFlatFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = 1000 );
	flipProposer = OneFlipProposer();
	initializer = genMeanFieldInitializer( cAreaInit );
	updaterType = SingleUpdater;
	auxDataType = WangLandauAuxData;
	
	itController = ItNumItController( itNum, itNumSample, itStartSample );
	
	# fName = loops_MC_methods_Base( divNum, itNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, auxDataType = auxDataType, itNumSample = itNumSample, itStartSample = itStartSample, nDim = nDim, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	
	fName = loops_MC_methods_Base( divNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, auxDataType = auxDataType, itController = itController, nDim = nDim, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	
	fMain = "loops_WL2d";
	attrLst, valLst = genAttrLstLttcFlipInit( divNum, itNum, nDim, flipChecker, initializer );
	fNameWL = fNameFunc( fMain, attrLst, valLst, jld2Type );
	
	println( "f = ", flipChecker.dosIncrRef[] );
	
	save( fNameWL, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	
	open( dirLog * fNameFileLstWL, "w" ) do io
		println( io, fNameWL );
	end
	
	# @infiltrate
	
	return fName;
end

function loops_MC_methods_WL2d( divNum = 64; dosIncrInit = 1, dosIncrMin = 0.001, cAreaInit = 0, nDim = 2, isFileNameOnly = false, fMainOutside = "", wlHistDosType = WLHistDos2DFull, wlHistDosArgs = (), wlResetInterval = 1000 )	
	flipChecker = WL2dHistStructFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = wlResetInterval, wlHistDosType = wlHistDosType, wlHistDosArgs = wlHistDosArgs );
	flipProposer = OneFlipProposer();
	initializer = genMeanFieldInitializer( cAreaInit );
	updaterType = SingleUpdater;
	auxDataType = WangLandauAuxData;
	
	itController = WangLandauItController( dosIncrMin, flipChecker );
	fName = loops_MC_methods_Base( divNum; updaterType = updaterType, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, auxDataType = auxDataType, itController = itController, nDim = nDim, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	
	itNum = itController.itRef[];
	
	println( "f = ", flipChecker.dosIncrRef[] );
	
	# save( fName, "histArr", flipChecker.histArr, "dosArr", flipChecker.dosArr );
	save( fName, "histArr", getHistArr( flipChecker ), "dosArr", getDosArr( flipChecker ) );
	
	# open( dirLog * fNameFileLstWL, "w" ) do io
		# println( io, fNameWL );
	# end
	
	return fName;
end

function loops_MC_methods_WL2dZoned( divNum = 64; dosIncrInit = 1, dosIncrMin = 0.001, cAreaInit = 0, nDim = 2, isFileNameOnly = false, fMainOutside = "", EMinRatio = -2.0, EMaxRatio = 0.2, EOverlap = 0.06, numZones = 8 )
	updaterType = SingleUpdater;
	auxDataType = WangLandauAuxData;
	flipProposer = OneFlipProposer();
	if !isFileNameOnly
		EWidth = ( (EMaxRatio - EMinRatio) + EOverlap * (numZones - 1) ) / numZones;
		
		EIntervalLst = [ Tuple( EMinRatio .+ [ (iE-1)*EWidth, iE * EWidth ] .- EOverlap .* (iE-1) ) for iE = 1 : numZones ];
		EMidLst = [ ( EIntervalLst[iE][1] + EIntervalLst[iE][2] ) / 2 for iE = 1 : numZones ];
		probLst = EMidLst ./ 4 .+ 1/2;
		
		configLst = loops_MC_WL_ScanConfig( divNum; nDim = nDim );
		
		flipCheckerLst = [ WL2dHistStructFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = 1000, wlHistDosType = WLHistDos2DZoned, wlHistDosArgs = EIntervalLst[iE] ) for iE = 1 : numZones ];
		initializerLst = [ PresetInitializer( configLst[ flipCheckerLst[iE].histDos.idMin ] ) for iE = 1 : numZones ];
		
		itControllerLst = [ WangLandauItController( dosIncrMin, flipCheckerLst[iE] ) for iE = 1 : numZones ];
		auxDataLst = Vector{auxDataType}(undef,numZones);
		fMainLst = Vector{String}(undef,numZones);
		
		println("Loops_MC zoning starts:")
		Threads.@threads for iE = 1 : numZones
			fMainLst[iE], auxDataLst[iE] = loops_MC_methods_Base_detailed_output( divNum; updaterType = updaterType, flipChecker = flipCheckerLst[iE], flipProposer = flipProposer, initializer = initializerLst[iE], auxDataType = auxDataType, itController = itControllerLst[iE], nDim = nDim, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
		end
		
		dosArrLst = [ (auxDataLst[iE].dataSampleLst)[2][end] for iE = 1 : numZones ];
		
		idMinLst = [ flipCheckerLst[iE].histDos.idMin for iE = 1 : numZones ];
		idMaxLst = [ flipCheckerLst[iE].histDos.idMax for iE = 1 : numZones ];
		grdNum = flipCheckerLst[1].grdNum;
		dosArrFull = zeros(grdNum-1);
		
		dosArrFull[idMinLst[1]:idMaxLst[1]] .= dosArrLst[1];
		@views for iE = 2 : numZones
			lnOverlap = idMaxLst[iE-1] - idMinLst[iE] + 1;
			shOverlap = mean( dosArrLst[iE][1:lnOverlap] - dosArrLst[iE-1][end-lnOverlap+1:end] );
			dosArrLst[iE] .-= shOverlap;
			dosArrFull[idMinLst[iE]:idMaxLst[iE-1]] .= ( dosArrFull[idMinLst[iE]:idMaxLst[iE-1]] .+ dosArrLst[iE][1:lnOverlap] ) ./ 2;
			dosArrFull[(idMaxLst[iE-1]+1) : idMaxLst[iE]] .= dosArrLst[iE][lnOverlap+1:end];
		end
		
		for id = idMaxLst[end]+1 : length(dosArrFull)
			if dosArrFull[id] == 0
				dosArrFull[id] = dosArrFull[length(dosArrFull)+1-id];
			end
		end
		
		itNum = itControllerLst[1].itRef[];
	end
	
	fMainGlued = "loops_WL2dZonesGluedThrd";
	
	flipCheckerGluedDummy = WL2dHistStructFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = 1000, wlHistDosType = WLHistDos2DFull );
	initializerDummy = PresetInitializer( [zeros(Bool,divNum, divNum)] );
	itControllerDummy = WangLandauItController( dosIncrMin, flipCheckerGluedDummy )
	
	fNameGlued = loops_MC_methods_Base( divNum; updaterType = updaterType, flipChecker = flipCheckerGluedDummy, flipProposer = flipProposer, initializer = initializerDummy, auxDataType = auxDataType, itController = itControllerDummy, nDim = nDim, isFileNameOnly = true, fMainOutside = fMainGlued );
	
	if !isFileNameOnly
		save( fNameGlued, "dosArrGlued", dosArrFull );
		println( "f = ", flipCheckerLst[end].dosIncrRef[] );
	end
	
	if isFileNameOnly
		return fNameGlued
	else
		return fMainLst;
	end
end

function genWLZone_idMinMaxLst( histType::Type{<:WLHistDosZonedInE}, divNum::Int64, EMinRatio::Real, EMaxRatio::Real, EOverlapRatio::Real, numZones::Int64 )
	EWidth = (EMaxRatio - EMinRatio) / ( numZones - (numZones-1) * EOverlapRatio );
	EOverlap = EWidth * EOverlapRatio;
	
	EIntervalLst = [ Tuple( min.( EMinRatio .+ [ (iE-1)*EWidth, iE * EWidth ] .- EOverlap .* (iE-1), 2 ) ) for iE = 1 : numZones ];
	
	idMinLst = ones(Int64, numZones);
	idMaxLst = similar(idMinLst);
	
	for iZone = 1 : numZones
		idMinLst[iZone], idMaxLst[iZone] = getIdMinMaxWLHist( histType, divNum, EIntervalLst[iZone]... );
	end
	
	return idMinLst, idMaxLst;
end

function getAttrValLstLttcWLReplica( divNum::Int64, nDim::Int64; itController::ItController, flipCheckerArr::Array{<:FlipChecker}, numZones::Int64, numWalksEach::Int64 )
	attrLstLttc, valLstLttc = genAttrValLstLttc( divNum, nDim );
	attrLstItCtrl, valLstItCtrl = getAttrValLstItController( itController );
	attrLstFlip, valLstFlip = getFlipCheckerAttrValLst( flipCheckerArr[1] );
	# attrLstInit, valLstInit = getAttrValInitializer( initializer );
	attrLstWLReplica = ["nZones", "nWalks"];
	valLstWLReplica = [numZones, numWalksEach];
	
	attrLstLst = [attrLstLttc, attrLstItCtrl, attrLstWLReplica, attrLstFlip];
	valLstLst = [valLstLttc, valLstItCtrl, valLstWLReplica, valLstFlip];
	
	attrLst = append!( attrLstLst... );
	valLst = append!( valLstLst... );
	
	return attrLst, valLst;
end

function loops_MC_methods_WL2dReplica( divNum = 64; dosIncrInit = 1, dosIncrMin = 0.001, cAreaInit = 0, nDim = 2, isFileNameOnly = false, fMainOutside = "", EMinRatio = -2.0, EMaxRatio = 2.0, EOverlapRatio = 0.7, numZones = 8, numWalksEach = 3, itExchange = 1000, wlResetInterval = 1000, histCutoffThres = 0.5, D_hist = 1, isCheckNan = false, fMod = "" )
	updaterType = SingleUpdater;
	auxDataType = WangLandauAuxData;
	flipProposer = OneFlipProposer();
	
	fModCheckNan = "checkNan";
	if isCheckNan
		fMod = strAppendWith_( fMod, fModCheckNan );
	end	
	
	# params = ParamsLoops( divNum, nDim );
	
	histDosType = WLHistDosZonedInE{nDim,D_hist};
	histDosFullType = WLHistDosFull{nDim,D_hist};
	flipCheckerType = WL2dHistStructFlipChecker{histDosType};
	
	idMinLst, idMaxLst = genWLZone_idMinMaxLst( histDosType, divNum, EMinRatio, EMaxRatio, EOverlapRatio, numZones );
	
	fMainGlued = "loops_WL2dZonesGluedReplica";
	if !isFileNameOnly
		# configLst = loops_MC_WL_ScanConfig( divNum, idMinLst, idMaxLst; nDim = nDim );
		configLst = loops_MC_WL_ScanConfig_Zoned( divNum, idMinLst, idMaxLst; nDim = nDim );
		fMainOutside = fMainGlued;
	else
		configLst = nothing;
	end
	
	fNameOut, histDosFinalLst, fLast = loops_MC_methods_WL2dReplicaBase( divNum; dosIncrInit = dosIncrInit, dosIncrMin = dosIncrMin, cAreaInit = cAreaInit, nDim = nDim, idMinLst = idMinLst, idMaxLst = idMaxLst, configZoneLst = configLst, numWalksEach = numWalksEach, itExchange = itExchange, wlResetInterval = wlResetInterval, histCutoffThres = histCutoffThres, D_hist = D_hist, isCheckNan = isCheckNan, fMainOutside = fMainOutside, isFileNameOnly = isFileNameOnly, fMod = fMod );
		
	if !isFileNameOnly
		histDosFull = genWLHistDos2DZonedFull( divNum );
		dosArrFull = histDosFull.dosArr;
		
		# iMatchLst = findGluePtsDosArrReplica( dosArrLst, idMinLst, idMaxLst );
		# glueDosArrReplicaFromGluePts!( dosArrFull, dosArrLst, idMinLst, idMaxLst, iMatchLst )
		iMatchLst, dosShLst = findGluePtsDosArrReplica(histDosFinalLst);
		dosArrFull = glueDosArrReplicaFromGluePts( histDosFinalLst, iMatchLst, dosShLst );
		
		# itNum = itController.itRef[];
		save( fNameOut, "dosArrGlued", dosArrFull, "dosArrZonedLst", getDosArr.(histDosFinalLst) );
		println( "f = ", fLast );
	end
	
	return fNameOut;
end

function loops_MC_methods_WL2dReplicaBase( divNum = 64; dosIncrInit = 1, dosIncrMin = 0.001, cAreaInit = 0, nDim = 2, idMinLst::Vector{Int64}, idMaxLst::Vector{Int64}, configZoneLst, numWalksEach = 3, itExchange = 1000, wlResetInterval = 1000, histCutoffThres = 0.5, D_hist = 1, isCheckNan = false, isFileNameOnly = false, fMainOutside::String = "", fMod = "" )
	numZones = length(idMinLst);
	if !isFileNameOnly && ( numZones != length(idMaxLst) || numZones != length(configZoneLst) )
		error( "idMin, idMax, or configZoneLst length mismatched" );
	end
	updaterType = SingleUpdater;
	auxDataType = WangLandauAuxData;
	flipProposer = OneFlipProposer();	
	
	params = ParamsLoops( divNum, nDim );
	
	histDosType = WLHistDosZonedInE{nDim,D_hist};
	histDosFullType = WLHistDosFull{nDim,D_hist};
	flipCheckerType = WL2dHistStructFlipChecker{histDosType};

	updaterLst = [ updaterType(params) for iE = 1 : numZones, iW = 1 : numWalksEach ];
	flipCheckerLst = [ WLFriendsFlipChecker{flipCheckerType}( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = wlResetInterval, histCutoffThres = histCutoffThres, wlHistDosArgs = ( idMinLst[iE], idMaxLst[iE] ) ) for iE = 1 : numZones, iW = 1 : numWalksEach ];
	
	histDosFinalLst = getHistDos.( @views( flipCheckerLst[:,1] ) );	
	
	if isnothing(configZoneLst)
		initializerLst = [ ConstantInitializer() for iE = 1 : numZones ];
	else
		initializerLst = [ PresetInitializer( configZoneLst[ iE ] ) for iE = 1 : numZones ];
	end
	
	auxDataLst = [ auxDataType( flipCheckerLst[iZ,iW] ) for iZ = 1 : numZones, iW = 1 : numWalksEach ];
	
	@views for iE = 1 : numZones, iW = 1 : numWalksEach
		setFriendCheckerLst!( flipCheckerLst[iE, iW], flipCheckerLst[iE,:]... );
	end
	
	fMainLst = Vector{String}(undef,numZones);
	
	if isCheckNan
		errCheckItController = NanCheckWLItController( flipCheckerLst );
		itControllerWL = WangLandauReplicaFineSampleItController( flipCheckerLst, dosIncrMin, itExchange );
		itController = JointItController( [errCheckItController, itControllerWL] );
		# fModOut = join( [fModOut, "nancheck"], "_" );
	else
		itController = WangLandauReplicaItController( flipCheckerLst, dosIncrMin, itExchange );
	end
	
	println("Loops_MC zoning starts:")
	fName = loops_MC_NoPrefabHelper_Replica( numZones, numWalksEach; params = params, updaterLst = updaterLst, flipCheckerLst = flipCheckerLst, flipProposer = flipProposer, initializerLst = initializerLst, auxDataLst = auxDataLst, itController = itController, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside, fMod = fMod );
	
	itNum = itController.itRef[];
	
	dosIncrLast = getDosIncr( flipCheckerLst[end] );
	dosIncrLst = getDosIncr.( flipCheckerLst );
	
	paramsToSaveLst = [itNum, dosIncrLast, dosIncrLst];
	paramsToSaveNameLst = ["itNumFinal","dosIncrLast" ,"dosIncrLst"];
	
	if isCheckNan
		isNanChecked = !testItNotDone( errCheckItController );
		push!(paramsToSaveLst, isNanChecked);
		push!( paramsToSaveNameLst, "isNanDetected" );
	end
	
	fMainMiscDataWL = "loops_WL2dReplica_MiscData";
	fNameMiscDataWL = loops_MC_NoPrefabHelper_Replica( numZones, numWalksEach; params = params, updaterLst = updaterLst, flipCheckerLst = flipCheckerLst, flipProposer = flipProposer, initializerLst = initializerLst, auxDataLst = auxDataLst, itController = itController, isFileNameOnly = true, fMainOutside = fMainMiscDataWL, fMod = fMod );
	Utils.saveJld2Lst( fNameMiscDataWL, paramsToSaveNameLst, paramsToSaveLst );
	
	return fName, histDosFinalLst, dosIncrLast;
end

function loops_MC_NoPrefabHelper_Replica( numZones, numWalksEach; params::ParamsLoops, updaterLst::Array{<:LoopsUpdater}, flipCheckerLst::Array{<:FlipChecker}, flipProposer::FlipProposer = OneFlipProposer, auxDataLst::Array{<:AuxData}, initializerLst::Vector{<:BLinkInitializer}, itController::Union{WangLandauReplicaItController,FindConfigReplicaItController,JointItController{<:Tuple{ErrCheckItController,AbstractWLReplicaItController}}}, fMod = "", isFileNameOnly::Bool = false, fMainOutside::Union{String, Nothing}= "" )
	# fModOut = copy( fMod );
	# if !isnothing(errCheckItController)
		# itController = JointItController( [errCheckItController, itController] );
		# fModOut = join( [fModOut, "nancheck"], "_" );
	# end
	
	attrLstReplica, valLstReplica = getAttrValLstLttcWLReplica( params.divNum, params.nDim; itController = itController, numZones = numZones, numWalksEach = numWalksEach, flipCheckerArr = flipCheckerLst );
	
	fNameOutside = fNameFunc( fMainOutside, attrLstReplica, valLstReplica, jld2Type; fMod = fMod );
	
	if !isFileNameOnly
		
		fNameLst = similar( flipCheckerLst, String );
		fModOutLst = similar( fNameLst );
		attrLstLst = similar(fNameLst, Vector{String});
		valLstLst = similar(fNameLst, Vector);
		
		for iZ = 1 : numZones, iW = 1 : numWalksEach
			updater = updaterLst[iZ, iW];
			flipChecker = flipCheckerLst[iZ, iW];
			initializer = initializerLst[iZ];
			fModOut = getFModLoopsMC( fMod, typeof( updater ); flipChecker = flipChecker, flipProposer = flipProposer );
			
			fMain = fMainLoopsMC;
			attrLst, valLst = genAttrLstLttcFullUpdater( params.divNum, params.nDim, flipChecker, initializer; itController = itController );
			fNameLst[iZ, iW] = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
			fModOutLst[iZ, iW] = fModOut;
			attrLstLst[iZ,iW] = attrLst;
			valLstLst[iZ,iW] = valLst;
		end
		
		# attrLst, valLst = genAttrLstLttcFullUpdater( params.divNum, params.nDim, flipChecker, initializer; itController = itController );
		
		# if isFileNameOnly
			# fNameOutside = fNameFunc( fMainOutside, attrLst, valLst, jld2Type; fMod = fModOut );
			# return fNameOutside;
		# end
		
		bLinkDataLst = [ genAuxData( BLinkAuxData, params, itController ) for iZ = 1 : numZones, iW = 1 : numWalksEach ];
		# BfieldLst, linkLst, linkFerroLst = bLinkData.dataLst;
		
		for iZ = 1 : numZones, iW = 1 : numWalksEach
			BfieldLst, linkLst, linkFerroLst = bLinkDataLst[iZ,iW].dataLst;
			initializeBL( initializerLst[iZ], getBfieldLst( bLinkDataLst[iZ,iW] ), params );
			updateLinkFrom0ByBAllDims( getBfieldLst( bLinkDataLst[iZ,iW] ), getLinkLst( bLinkDataLst[iZ,iW] ), getLinkFerroLst( bLinkDataLst[iZ,iW] ), params );
			calcAuxData!( auxDataLst[iZ,iW], params, BfieldLst, linkLst, linkFerroLst );
		end
		walksIdLst = [1:numWalksEach;];
		flipCheckerFriendsLst = [ @view( flipCheckerLst[iZ,:] ) for iZ = 1 : numZones ];
		flipCheckerMasterLst = @view( flipCheckerLst[:,1] );
		lnMaster = length(flipCheckerMasterLst);
		idMidMaster = Int64(floor(lnMaster / 2));
		idDisplayMasterLst = pushfirst!( Int64.( floor.( [1:4;]./4 * lnMaster ) ), 1 );
		# idDisplayMasterLst = [1,idMidMaster,lnMaster];
		
		resetItControl( itController );
		
		itNum, itNumSample, itNumStartSample = getItNumLst( itController );
		
		# itExchange = itController.itExchange;
		itExchange = getItExchange( itController );
		
		itSample = 1;
		itStepOutput = 10000;
		itStepGC = 10000 * itExchange;
		@time while( testItNotDone( itController ) )
			# @time begin
			it = itController.itRef[];
			if mod(it,itStepOutput) == 0 || mod(it,itStepOutput) == 1
				# print( "it = ", it, ", fmax = ", maximum( getDosIncr.( flipCheckerLst ) ), ", fmin = ", minimum( getDosIncr.( flipCheckerLst ) ), ", ", roundKeepInt.( getDosIncr.( flipCheckerMasterLst ); digits=5 ), " \r" )
				print( "it = ", it, ", fmax = ", maximum( getDosIncr.( flipCheckerLst ) ), ", fmin = ", minimum( getDosIncr.( flipCheckerLst ) ) );
				print(", ", roundKeepInt.( getDosIncr.( flipCheckerMasterLst[idDisplayMasterLst] ); digits=5 ) );
				print( " \r" );
				# ", itWL = " , flipCheckerLst[1,1].selfChecker.wlCounter[]
			end
			
			# @timeit tmr "sample" 
			if testItDoSample( itController )
				for iZ = 1 : numZones, iW = 1 : numWalksEach
					storeAuxDataSample( bLinkDataLst[iZ,iW], itSample, it );
					storeAuxDataSample( auxDataLst[iZ,iW], itSample, it );
				end
				itSample += 1;
			end
			
			# @timeit tmr "startSample" 
			if testItDoStartSample( itController )
				for iZ = 1 : numZones, iW = 1 : numWalksEach
					storeAuxDataStartSample( bLinkDataLst[iZ,iW], it );
					storeAuxDataStartSample( auxData[iZ,iW], it );
				end
			end
			
			# @timeit tmr "exchange" 
			if testItDoExchange( itController )
				for iZ = 1 : numZones - 1
					shuffle!( walksIdLst );
					for iW = 1 : numWalksEach
						iWPair = walksIdLst[iW];
						flipCheckExchange( flipCheckerLst[iZ,iW], flipCheckerLst[iZ+1,iWPair], params, bLinkDataLst[iZ,iW], bLinkDataLst[iZ+1,iWPair] );
					end
				end
			end
			
			# @time Threads.@threads 
			# @time 
			# for iSeg = 1 : length(bLinkDataLst)
				# storeAuxDataNum( bLinkDataLst[iSeg], it );
				# storeAuxDataNum( auxDataLst[iSeg], it );
			# end
			
			# println("updateLoops")
			# @time begin
			# Threads.@threads 
			
			# @time 
			GC.enable(false)
			# @timeit tmr "updateLoop" 
			Threads.@threads for iSeg = 1 : length( bLinkDataLst )
				# GC.enable(false)
				BfieldLst, linkLst, linkFerroLst = bLinkDataLst[iSeg].dataLst;
				# println("old wlCunter = ", flipCheckerLst[iSeg].selfChecker.wlCounter[])
				for itInt = 1 : itExchange
					updateLoops( updaterLst[iSeg], flipCheckerLst[iSeg], flipProposer, auxDataLst[iSeg], BfieldLst, linkLst, linkFerroLst, params );
				end
				# println("new wlCunter = ", flipCheckerLst[iSeg].selfChecker.wlCounter[])
				if flipCheckerLst[iSeg].selfChecker.wlCounter[] == 1
					flipCheckerLst[iSeg].selfChecker.wlCounter[] = itExchange;
				else
					flipCheckerLst[iSeg].selfChecker.wlCounter[] += itExchange;
				end
				# GC.enable(true)
			end
			# Threads.@threads for iSeg = 1 : length( bLinkDataLst )
				GC.enable(true)
			# end
			# GC.gc()
			
			# @time Threads.@threads 
			# @time 
			# @timeit tmr "wlSynchronize" 
			Threads.@threads for iZ = 1 : numZones
				wangLandauHistResetSynchronized( flipCheckerFriendsLst[iZ] );
			end
			if flipCheckerLst[1].selfChecker.wlCounter[] % flipCheckerLst[1].selfChecker.wlResetInterval == 0			
				# @infiltrate
			end
			
			# @timeit tmr "advanceIt" 
			for itInt = 1 : itExchange
				advanceItControl( itController );
			end
			
			if mod( it, itStepGC ) == 0 || mod( it, itStepGC ) == 1
				# Threads.@threads for iSeg = 1 : length( bLinkDataLst )
					# GC.enable(true)
					GC.gc();
					# GC.enable(false)
				# end
			end
			# @infiltrate mod( it, 50000 ) == 0
			# end
			# @infiltrate
			# if it > 1e6
				# break;
			# end
		end
		
		# if testItDoSample( itController )
		it = itController.itRef[];
		for iZ = 1 : numZones, iW = 1 : numWalksEach
			storeAuxDataSample( bLinkDataLst[iZ,iW], itSample, it );
			storeAuxDataSample( auxDataLst[iZ,iW], itSample, it );
		end
		itSample += 1;
		
		for iZ = 1 : numZones
			wangLandauHistResetSynchronized( flipCheckerFriendsLst[iZ] );
		end
		
		wlBundledAuxData = BundledArrAuxData( auxDataLst );
		bLinkBundledAuxData = BundledArrAuxData( bLinkDataLst );
		
		saveAuxDataAll( wlBundledAuxData, attrLstReplica, valLstReplica; fMod = fMod );
		saveAuxDataAll( bLinkBundledAuxData, attrLstReplica, valLstReplica; fMod = fMod );
		
		println("it total: ", itController.itRef[]);
	
		# end
	end
	
	# for iZ = 1 : numZones, iW = 1 : numWalksEach
		# attrLst = attrLstLst[iZ,iW];
		# valLst = valLstLst[iZ,iW];
		# fModOutZones = join( [fModOutLst[iZ, iW], string(iZ), string(iW) ] );
		# saveAuxDataAll( bLinkDataLst[iZ,iW], attrLst, valLst; fMod = fModOutZones );
		# saveAuxDataAll( auxDataLst[iZ,iW], attrLst, valLst; fMod = fModOutZones );
	# end
	
	# if isFileNameOnly
		# return fNameOutside;
	# end
	
	# return fNameLst;
	return fNameOutside
end

function loops_MC_methods_WL2dReplica_MPI( divNum = 64; dosIncrInit = 1, dosIncrMin = 0.001, cAreaInit = 0, nDim = 2, isFileNameOnly = false, fMainOutside = "", EMinRatio = -2.0, EMaxRatio = 2.0, EOverlapRatio = 0.7, numZones = 8, numWalksEach = 3, itExchange = 1000, wlResetInterval = 1000, histCutoffThres = 0.5, D_hist = 1 )
	updaterType = SingleUpdater;
	auxDataType = WangLandauAuxData;
	flipProposer = OneFlipProposer(); 
	
	
	params = ParamsLoops( divNum, nDim );
	
	# histDosType = WLHistDos2DZoned;
	# histDosType = WLHistDosZonedInE1dDos;
	histDosType = WLHistDosZonedInE{D_hist};
	histDosFullType = WLHistDosFull{D_hist};
	flipCheckerType = WL2dHistStructFlipChecker{histDosType};
	
	if !isFileNameOnly
		EWidth = (EMaxRatio - EMinRatio) / ( numZones - (numZones-1) * EOverlapRatio );
		EOverlap = EWidth * EOverlapRatio;
		
		EIntervalLst = [ Tuple( EMinRatio .+ [ (iE-1)*EWidth, iE * EWidth ] .- EOverlap .* (iE-1) ) for iE = 1 : numZones ];
		EMidLst = [ ( EIntervalLst[iE][1] + EIntervalLst[iE][2] ) / 2 for iE = 1 : numZones ];
		probLst = EMidLst ./ 4 .+ 1/2;
		
		configLst = loops_MC_WL_ScanConfig( divNum; nDim = nDim );
		
		updaterLst = [ updaterType(params) for iE = 1 : numZones, iW = 1 : numWalksEach ]
		flipCheckerLst = [ WLFriendsFlipChecker{flipCheckerType}( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = wlResetInterval, histCutoffThres = histCutoffThres, wlHistDosArgs = EIntervalLst[iE] ) for iE = 1 : numZones, iW = 1 : numWalksEach ];
		initializerLst = [ PresetInitializer( configLst[ getHistDos( flipCheckerLst[iE,1] ).idMin ] ) for iE = 1 : numZones ];
		auxDataLst = [ auxDataType( flipCheckerLst[iZ,iW] ) for iZ = 1 : numZones, iW = 1 : numWalksEach ];
		
		@views for iE = 1 : numZones, iW = 1 : numWalksEach
			setFriendCheckerLst!( flipCheckerLst[iE, iW], flipCheckerLst[iE,:]... );
		end
		
		itController = WangLandauReplicaItController( flipCheckerLst, dosIncrMin, itExchange );
		fMainLst = Vector{String}(undef,numZones);
		
		println("Loops_MC zoning starts:")
		fName = loops_MC_NoPrefabHelper_Replica( numZones, numWalksEach; params = params, updaterLst = updaterLst, flipCheckerLst = flipCheckerLst, flipProposer = flipProposer, initializerLst = initializerLst, auxDataLst = auxDataLst, itController = itController, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
		
		histDosFinalLst = getHistDos.( @views( flipCheckerLst[:,1] ) );
		dosArrLst = [ (auxDataLst[iE,1].dataSampleLst)[2][end] for iE = 1 : numZones ];
		
		idMinLst = [ getHistDos( flipCheckerLst[iE,1] ).idMin for iE = 1 : numZones ];
		idMaxLst = [ getHistDos( flipCheckerLst[iE,1] ).idMax for iE = 1 : numZones ];
		histDosFull = genWLHistDos2DZonedFull( divNum );
		dosArrFull = histDosFull.dosArr;
		
		# iMatchLst = findGluePtsDosArrReplica( dosArrLst, idMinLst, idMaxLst );
		# glueDosArrReplicaFromGluePts!( dosArrFull, dosArrLst, idMinLst, idMaxLst, iMatchLst )
		iMatchLst, dosShLst = findGluePtsDosArrReplica(histDosFinalLst);
		dosArrFull = glueDosArrReplicaFromGluePts( histDosFinalLst, iMatchLst, dosShLst );
		
		itNum = itController.itRef[];
	end
	
	fMainGlued = "loops_WL2dZonesGluedReplica";
	
	flipCheckerGluedDummy = WL2dHistStructFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = wlResetInterval, wlHistDosType = histDosFullType );
	initializerDummy = PresetInitializer( [zeros(Bool,divNum, divNum)] );
	attrLstMod = ["numZones", "numWalks"];
	valLstMod = [numZones, numWalksEach];
	itController = WangLandauItController( dosIncrMin, flipCheckerGluedDummy );
	
	fNameGlued = loops_MC_methods_Base( divNum; updaterType = updaterType, flipChecker = flipCheckerGluedDummy, flipProposer = flipProposer, initializer = initializerDummy, auxDataType = auxDataType, itController = itController, nDim = nDim, isFileNameOnly = true, fMainOutside = fMainGlued, attrLstMod = attrLstMod, valLstMod = valLstMod );
	
	if !isFileNameOnly
		save( fNameGlued, "dosArrGlued", dosArrFull, "dosArrZonedLst", dosArrLst );
		println( "f = ", getDosIncr( flipCheckerLst[end] ) );
	end
	
	if isFileNameOnly
		return fNameGlued
	else
		return fName;
	end
end

function loops_MC_NoPrefabHelper_Replica_MPI( numZones, numWalksEach; params::ParamsLoops, updaterLst::Array{<:LoopsUpdater}, flipCheckerLst::Array{<:FlipChecker}, flipProposer::FlipProposer = OneFlipProposer, auxDataLst::Array{<:AuxData}, initializerLst::Vector{<:BLinkInitializer}, itController::WangLandauReplicaItController, fMod = "", isFileNameOnly::Bool = false, fMainOutside::Union{String, Nothing}= "" )
	fNameLst = similar( flipCheckerLst, String );
	fModOutLst = similar( fNameLst );
	attrLstLst = similar(fNameLst, Vector{String});
	valLstLst = similar(fNameLst, Vector);
	for iZ = 1 : numZones, iW = 1 : numWalksEach
		updater = updaterLst[iZ, iW];
		flipChecker = flipCheckerLst[iZ, iW];
		initializer = initializerLst[iZ];
		fModOut = getFModLoopsMC( fMod, typeof( updater ); flipChecker = flipChecker, flipProposer = flipProposer );
		
		fMain = fMainLoopsMC;
		attrLst, valLst = genAttrLstLttcFullUpdater( params.divNum, params.nDim, flipChecker, initializer; itController = itController );
		fNameLst[iZ, iW] = fNameFunc( fMain, attrLst, valLst, jld2Type; fMod = fModOut );
		fModOutLst[iZ, iW] = fModOut;
		attrLstLst[iZ,iW] = attrLst;
		valLstLst[iZ,iW] = valLst;
	end
	
	# attrLst, valLst = genAttrLstLttcFullUpdater( params.divNum, params.nDim, flipChecker, initializer; itController = itController );
	
	# if isFileNameOnly
		# fNameOutside = fNameFunc( fMainOutside, attrLst, valLst, jld2Type; fMod = fModOut );
		# return fNameOutside;
	# end
	
	bLinkDataLst = [ genAuxData( BLinkAuxData, params, itController ) for iZ = 1 : numZones, iW = 1 : numWalksEach ];
	# BfieldLst, linkLst, linkFerroLst = bLinkData.dataLst;
	
	for iZ = 1 : numZones, iW = 1 : numWalksEach
		BfieldLst, linkLst, linkFerroLst = bLinkDataLst[iZ,iW].dataLst;
		initializeBL( initializerLst[iZ], getBfieldLst( bLinkDataLst[iZ,iW] ), params );
		updateLinkFrom0ByBAllDims( getBfieldLst( bLinkDataLst[iZ,iW] ), getLinkLst( bLinkDataLst[iZ,iW] ), getLinkFerroLst( bLinkDataLst[iZ,iW] ), params );
		calcAuxData!( auxDataLst[iZ,iW], params, BfieldLst, linkLst, linkFerroLst );
	end
	walksIdLst = [1:numWalksEach;];
	flipCheckerFriendsLst = [ @view( flipCheckerLst[iZ,:] ) for iZ = 1 : numZones ];
	flipCheckerMasterLst = @view( flipCheckerLst[:,1] );
	lenMaster = length(flipCheckerMasterLst);
	idMidMaster = Int64(floor(lenMaster / 2));
	idDisplayMasterLst = [1,idMidMaster,lnMaster];
	
	resetItControl( itController );
	
	itNum, itNumSample, itNumStartSample = getItNumLst( itController );
	
	itSample = 1;
	@time while( testItNotDone( itController ) )
		# @time begin
		it = itController.itRef[];
		if mod(it,1000) == 0
			# print( "it = ", it, ", fmax = ", maximum( getDosIncr.( flipCheckerLst ) ), ", fmin = ", minimum( getDosIncr.( flipCheckerLst ) ), ", ", roundKeepInt.( getDosIncr.( flipCheckerMasterLst ); digits=5 ), " \r" )
			print( "it = ", it, ", fmax = ", maximum( getDosIncr.( flipCheckerLst ) ), ", fmin = ", minimum( getDosIncr.( flipCheckerLst ) ) );
			print(", ", roundKeepInt.( getDosIncr( flipCheckerMasterLst[idDisplayMasterLst] ); digits=5 ) );
			print( " \r" );
		end
		
		if testItDoSample( itController )
			for iZ = 1 : numZones, iW = 1 : numWalksEach
				storeAuxDataSample( bLinkDataLst[iZ,iW], itSample, it );
				storeAuxDataSample( auxDataLst[iZ,iW], itSample, it );
			end
			itSample += 1;
		end
		
		if testItDoStartSample( itController )
			for iZ = 1 : numZones, iW = 1 : numWalksEach
				storeAuxDataStartSample( bLinkDataLst[iZ,iW], it );
				storeAuxDataStartSample( auxData[iZ,iW], it );
			end
		end
		
		if testItDoExchange( itController )
			for iZ = 1 : numZones - 1
				shuffle!( walksIdLst );
				for iW = 1 : numWalksEach
					iWPair = walksIdLst[iW];
					flipCheckExchange( flipCheckerLst[iZ,iW], flipCheckerLst[iZ+1,iWPair], params, bLinkDataLst[iZ,iW], bLinkDataLst[iZ+1,iWPair] );
				end
			end
		end
		
		# @time Threads.@threads 
		# @time 
		for iSeg = 1 : length(bLinkDataLst)
			storeAuxDataNum( bLinkDataLst[iSeg], it );
			storeAuxDataNum( auxDataLst[iSeg], it );
		end
		
		# println("updateLoops")
		# @time begin
		# Threads.@threads 
		# @time 
		Threads.@threads for iSeg = 1 : length( bLinkDataLst )
			BfieldLst, linkLst, linkFerroLst = bLinkDataLst[iSeg].dataLst;
			updateLoops( updaterLst[iSeg], flipCheckerLst[iSeg], flipProposer, auxDataLst[iSeg], BfieldLst, linkLst, linkFerroLst, params );
		end
		
		# @time Threads.@threads 
		# @time 
		for iZ = 1 : numZones
			wangLandauHistResetSynchronized( flipCheckerFriendsLst[iZ] );
		end
		
		advanceItControl( itController );
		@infiltrate mod( it, 50000 ) == 0
		# end
	end
	
	# if testItDoSample( itController )
	it = itController.itRef[];
	for iZ = 1 : numZones, iW = 1 : numWalksEach
		storeAuxDataSample( bLinkDataLst[iZ,iW], itSample, it );
		storeAuxDataSample( auxDataLst[iZ,iW], itSample, it );
	end
	itSample += 1;
	
	for iZ = 1 : numZones
		wangLandauHistResetSynchronized( flipCheckerFriendsLst[iZ] );
	end
	
	# end
	
	# for iZ = 1 : numZones, iW = 1 : numWalksEach
		# attrLst = attrLstLst[iZ,iW];
		# valLst = valLstLst[iZ,iW];
		# fModOutZones = join( [fModOutLst[iZ, iW], string(iZ), string(iW) ] );
		# saveAuxDataAll( bLinkDataLst[iZ,iW], attrLst, valLst; fMod = fModOutZones );
		# saveAuxDataAll( auxDataLst[iZ,iW], attrLst, valLst; fMod = fModOutZones );
	# end
	
	# if isFileNameOnly
		# return fNameOutside;
	# end
	
	return fNameLst;
end

# function findGluePtsDosArrReplica( histDosLst::AbstractVector{WLHistDosZonedInE1dDos} )
	# diffLst1 = zeros(0);
	# diffLst2 = zeros(0);
	# diffDiff = zeros(0);
	# iMatchLst = zeros(Int64, length( histDosLst ) - 1);
	# dosShLst = zeros(length(iMatchLst));
	# for iZ = 1 : length( dosShLst )
		# lnOverlap = histDosLst[iZ].idMax - histDosLst[iZ+1].idMin + 1;
		# resize!.( [diffLst1, diffLst2, diffDiff], lnOverlap );
		
		# iDiff = 1;
		# for iOver = 1 : lnOverlap
			# diffLst1[iOver] = arrDiffCalc( getDosArr(histDosLst[iZ]), length(getDosArr(histDosLst[iZ])) -
			# lnOverlap + iOver );
			# diffLst2[iOver] = arrDiffCalc( getDosArr(histDosLst[iZ+1]), iOver );
		# end
		# diffDiff .= abs.( diffLst2 .- diffLst1 );
		# _, iMatchLst[iZ] = findmin( diffDiff );
		# iMatchLst[iZ] += histDosLst[iZ+1].idMin - 1;
	# end
	
	# return iMatchLst;
# end

function findGluePtsDosArrReplica( histDosLst::AbstractVector{<:WLHistDosZonedInE{nDim,1} where {nDim}} )
	return findGluePtsDosArrReplica( getDosArr.( histDosLst ), (h->h.idMin).(histDosLst), (h->h.idMax).(histDosLst) )
end

function findGluePtsDosArrReplica( histDosLst::AbstractVector{<:WLHistDosZonedInE{nDim,2} where{nDim}} )
	dosArr2dLst = getDosArr.( histDosLst );
	dosMinLst = minimum.(a -> a !=0 ? a : Inf, dosArr2dLst);
	for iHist = 1 : length( dosArr2dLst )
		for i2d = 1 : length( dosArr2dLst[iHist] )
			if dosArr2dLst[iHist][i2d] != 0
				dosArr2dLst[iHist][i2d] -= dosMinLst[iHist];
				dosArr2dLst[iHist][i2d] += log(2);
			end
		end
	end
	
	dosMaxLst = maximum.(dosArr2dLst);

	dosArrExpShMaxLst = [ similar( dosArr2dLst[ii] ) for ii = 1 : length(histDosLst) ];
	for iHist = 1 : length( dosArr2dLst )
		for i2d = 1 : length( dosArr2dLst[iHist] )
			if dosArr2dLst[iHist][i2d] == 0
				dosArrExpShMaxLst[iHist][i2d] = 0;
			else
				dosArrExpShMaxLst[iHist][i2d] = exp(dosArr2dLst[iHist][i2d] - dosMaxLst[iHist]);
			end
		end
	end

	# dosArr1dLst = ( d -> dropdims( log.( sum( exp.( d ); dims = 1 ) ); dims = 1 ) ).(dosArr2dLst);
	dosArr1dLst = ( (d,mx) -> dropdims( log.( sum( d; dims = 1 ) ); dims = 1 ) .+ mx ).(dosArrExpShMaxLst, dosMaxLst);
	return findGluePtsDosArrReplica( dosArr1dLst, (h->h.idMin).(histDosLst), (h->h.idMax).(histDosLst) )
end

# function findGluePtsDosArrReplica( histDosLst::AbstractVector{<:WLHistDosZonedInE{nDim,2} where{nDim}} )
	# dosArr2dLst = getDosArr.( histDosLst );
	# dosMinLst = minimum.(a -> a !=0 ? a : Inf, dosArr2dLst);
	# # dosArr2dLst = ( (x,y)->((a,b)->(a != 0 ? a -= b - log(2) : nothing)).(x,y) ).(dosArr2dLst, dosMinLst);
	# for iHist = 1 : length( dosArr2dLst )
		# for i2d = 1 : length( dosArr2dLst[iHist] )
			# if dosArr2dLst[iHist][i2d] != 0
				# dosArr2dLst[iHist][i2d] -= dosMinLst[iHist];
				# dosArr2dLst[iHist][i2d] += log(2);
			# end
		# end
	# end
	
	# dosArr1dLst = ( d -> dropdims( log.( sum( exp.( d ); dims = 1 ) ); dims = 1 ) ).(dosArr2dLst);
	# return findGluePtsDosArrReplica( dosArr1dLst, (h->h.idMin).(histDosLst), (h->h.idMax).(histDosLst) )
# end

function findGluePtsDosArrReplica( dosArrLst, idMinLst, idMaxLst )
	diffLst1 = zeros(0);
	diffLst2 = zeros(0);
	diffDiff = zeros(0);
	iMatchLst = zeros(Int64, length(dosArrLst) - 1);
	dosShLst = zeros(length(dosArrLst) - 1);
	for iZ = 1 : length( dosArrLst ) - 1
		lnOverlap = idMaxLst[iZ] - idMinLst[iZ+1] + 1;
		resize!.( [diffLst1, diffLst2, diffDiff], lnOverlap );
		
		iDiff = 1;
		for iOver = 1 : lnOverlap
			diffLst1[iOver] = arrDiffCalc( dosArrLst[iZ], length(dosArrLst[iZ]) - lnOverlap + iOver );
			diffLst2[iOver] = arrDiffCalc( dosArrLst[iZ+1], iOver );
		end
		diffDiff .= abs.( diffLst2 .- diffLst1 );
		_, iMatchLst[iZ] = findmin( diffDiff );
		iMatchLst[iZ] += idMinLst[iZ+1]-1;
	end
	
	for iZ = 1 : length( dosArrLst ) - 1
		idNxt = iMatchLst[iZ] - idMinLst[iZ+1] + 1;
		idThis = iMatchLst[iZ] - idMinLst[iZ] + 1;
		dosShLst[iZ] = dosArrLst[iZ+1][idNxt] - dosArrLst[iZ][idThis];
	end
	
	return iMatchLst, dosShLst;
end

function glueDosArrReplicaFromGluePts( histDosLst::AbstractVector{<:WLHistDosZonedInE{nDim,1} where {nDim}}, iMatchLst::Vector{Int64}, dosShLst::Vector{Float64} )
	dosArrFull = zeros( histDosLst[1].numLnkVal - 1 );
	dosArrLst = getDosArr.( histDosLst );
	
	iMatchLstTailHead = copy( iMatchLst );
	pushfirst!(iMatchLstTailHead, histDosLst[1].idMin-1 )
	push!(iMatchLstTailHead, histDosLst[end].idMax);
	
	for iZ = 1 : length(dosArrLst)
		if iZ > 1
			ii = iMatchLstTailHead[iZ];
			iThis = getZonedIdFromFull( histDosLst[iZ], ii );
			getDosArr( histDosLst[iZ] ) .-= getDosArr( histDosLst[iZ] )[iThis] - dosArrFull[ii];
		end
		for ii = iMatchLstTailHead[iZ]+1 : iMatchLstTailHead[iZ+1]
			# iThis = ii - histDosLst[iZ].idMin + 1
			iThis = getZonedIdFromFull( histDosLst[iZ], ii );
			dosArrFull[ii] = getDosArr( histDosLst[iZ] )[iThis];
		end
	end
	
	for ii = iMatchLstTailHead[end] + 1 : length(dosArrFull)
		if dosArrFull[ii] == 0
			dosArrFull[ii] = dosArrFull[ length(dosArrFull)-ii+1 ];
		end
	end
	
	return dosArrFull;
end

function glueDosArrReplicaFromGluePts( histDosLst::AbstractVector{<: (WLHistDosZonedInE{nDim,2} where {nDim})}, iMatchLst::Vector{Int64}, dosShLst::Vector{Float64} )
	dosArrFull = zeros( size( histDosLst[1].histArr, 1 ), histDosLst[1].numLnkVal - 1 );
	dosArrLst = getDosArr.( histDosLst );
	
	dosShLstAcc = accumulate(+,dosShLst);
	
	iMatchLstTailHead = copy( iMatchLst );
	pushfirst!(iMatchLstTailHead, histDosLst[1].idMin-1 )
	push!(iMatchLstTailHead, histDosLst[end].idMax);
	
	@views for iZ = 1 : length(dosArrLst)
		if iZ > 1
			ii = iMatchLstTailHead[iZ];
			iThis = getZonedIdFromFull( histDosLst[iZ], ii );
			# shDos = mean( getDosArr( histDosLst[iZ] )[:,iThis] - dosArrFull[:,ii] );
			for iDos = 1 : length( getDosArr(histDosLst[iZ]) )
				if getDosArr(histDosLst[iZ])[iDos] != 0
					getDosArr(histDosLst[iZ])[iDos] -= dosShLstAcc[iZ-1];
				end
			end
			(a -> a!=0 ? a -= dosShLstAcc[iZ-1] : nothing).( getDosArr( histDosLst[iZ] ) );
			# (a -> a!=0 ? a -= shDos : nothing),getDosArr( histDosLst[iZ] ) .-= shDos;
		end
		for ii = iMatchLstTailHead[iZ]+1 : iMatchLstTailHead[iZ+1]
			iThis = getZonedIdFromFull( histDosLst[iZ], ii );
			dosArrFull[:,ii] = getDosArr( histDosLst[iZ] )[:,iThis];
		end
	end
	
	lnLinks = size(dosArrFull)[end];
	@views for ii = iMatchLstTailHead[end] + 1 : lnLinks
		# if dosArrFull[ii] == 0
			dosArrFull[:,ii] = dosArrFull[ :, lnLinks-ii+1 ];
		# end
	end
	
	return dosArrFull;
end

function glueDosArrReplicaFromGluePts!( dosArrFull, dosArrLst, idMinLst, idMaxLst, iMatchLst )
	iMatchLstTailHead = copy( iMatchLst );
	pushfirst!(iMatchLstTailHead, idMinLst[1]-1 )
	push!(iMatchLstTailHead, idMaxLst[end]);
	# @infiltrate
	# dosArrFull[1] = dosArrLst[1][1];
	for iZ = 1 : length(dosArrLst)
		if iZ > 1
			ii = iMatchLstTailHead[iZ];
			iThis = ii - idMinLst[iZ] + 1
			dosArrLst[iZ] .-= dosArrLst[iZ][iThis] - dosArrFull[ii];
		end
		for ii = iMatchLstTailHead[iZ]+1 : iMatchLstTailHead[iZ+1]
			iThis = ii - idMinLst[iZ] + 1
			dosArrFull[ii] = dosArrLst[iZ][iThis];
		end
	end
	
	for ii = iMatchLstTailHead[end] + 1 : length(dosArrFull)
		if dosArrFull[ii] == 0
			dosArrFull[ii] = dosArrFull[ length(dosArrFull)-ii+1 ];
		end
	end
end

function loops_MC_WL_ScanConfig( divNum = 64; nDim = 2, isFileNameOnly = false, fMainOutside = "" )
	flipChecker = WL2dHistStructFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = 1000, wlHistDosType = WLHistDos2DFull );
	params = ParamsLoops( divNum, nDim );
	
	flipProposer = OneFlipProposer();
	initializer = ConstantInitializer();
	updater = SingleUpdater(params);
	configAuxData = BConfigOfLnkValAuxData{params.nDim}( flipChecker );
	
	itController = FindConfigItController( configAuxData );
		
	fName = loops_MC_NoPrefabHelper_Base( ; params = params, updater = updater, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, auxData = configAuxData, itController = itController, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	itNum = itController.itRef[];
	
	return configAuxData.configArr;
end

function loops_MC_WL_ScanConfig( divNum::Int64, histDosLst::AbstractVector{<:WLHistDosZonedInE}; nDim = 2 )
	idMinLst = (h->h.idMin).(histDosLst);
	idMaxLst = (h->h.idMax).(histDosLst);
	
	loops_MC_WL_ScanConfig( divNum, idMinLst, idMaxLst; nDim = nDim );
end

function loops_MC_WL_ScanConfig( divNum::Int64, idMinLst::Vector{Int64}, idMaxLst::Vector{Int64}; isFileNameOnly = false, fMainOutside = "", nDim = 2 )
	flipChecker = WL2dHistStructFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = 1000, wlHistDosType = WLHistDosFull{nDim,1} );
	params = ParamsLoops( divNum, nDim );
	
	flipProposer = OneFlipProposer();
	initializer = ConstantInitializer();
	updater = SingleUpdater(params);
	configAuxData = BConfigOfLnkValZonedAuxData{params.nDim}( idMinLst, idMaxLst, 0, 0, 0 );
	
	itController = FindConfigItController( configAuxData );
		
	fName = loops_MC_NoPrefabHelper_Base( ; params = params, updater = updater, flipChecker = flipChecker, flipProposer = flipProposer, initializer = initializer, auxData = configAuxData, itController = itController, isFileNameOnly = isFileNameOnly, fMainOutside = fMainOutside );
	itNum = itController.itRef[];
	
	return configAuxData.configArr;
end

function loops_MC_WL_ScanConfig_Zoned( divNum::Int64, idMinLst::Vector{Int64}, idMaxLst::Vector{Int64}; isFileNameOnly = false, fMainOutside = "", nDim = 2, numWalksEach = 3, numZonesPre = 4, overlapRatioPre = 0.5 )
	flipChecker = WL2dHistStructFlipChecker( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = 1000, wlHistDosType = WLHistDosFull{nDim,1} );
	params = ParamsLoops( divNum, nDim );
	
	flipProposer = OneFlipProposer();
	initializer = ConstantInitializer();
	updater = SingleUpdater(params);
	
	idMinBnd = idMinLst[1];
	idMaxBnd = idMaxLst[end];
	idWidth = (idMaxBnd - idMinBnd) / ( numZonesPre - (numZonesPre-1) * overlapRatioPre );
	idOverlap = idWidth * overlapRatioPre;
	
	idIntervalLst = [ Int64( floor( min( idMinBnd + ( iZ + iZSh ) * idWidth - idOverlap * (iZ-1), idMaxBnd ) ) ) for iZ = 1 : numZonesPre, iZSh = -1 : 0 ];
	
	idMinLstPre = idIntervalLst[:,1];
	idMaxLstPre = idIntervalLst[:,2];
	idMinLstPre[1] = idMinLst[1];
	idMaxLstPre[end] = idMaxLst[end];
	
	# EWidth = (EMaxRatio - EMinRatio) / ( numZones - (numZones-1) * EOverlapRatio );
	# EOverlap = EWidth * EOverlapRatio;
	
	# EIntervalLst = [ Tuple( min.( EMinRatio .+ [ (iE-1)*EWidth, iE * EWidth ] .- EOverlap .* (iE-1), 2 ) ) for iE = 1 : numZones ];
	
	preConfigZoneLst = loops_MC_WL_ScanConfig( divNum, idMinLstPre, idMaxLstPre; nDim = nDim );
		
	configArrFinal = loops_MC_methods_WL2dReplicaScanBase( divNum; nDim = nDim, idMinLstPre = idMinLstPre, idMaxLstPre = idMaxLstPre, preConfigZoneLst = preConfigZoneLst, idMinLst = idMinLst, idMaxLst = idMaxLst, numWalksEach = numWalksEach )
	
	return configArrFinal;
end

function loops_MC_methods_WL2dReplicaScanBase( divNum = 64; dosIncrInit = 1, dosIncrMin = 0.001, cAreaInit = 0, nDim = 2, idMinLstPre::Vector{Int64}, idMaxLstPre::Vector{Int64}, idMinLst::Vector{Int64}, idMaxLst::Vector{Int64}, preConfigZoneLst, numWalksEach = 3, itExchange = 1000, wlResetInterval = 1000, histCutoffThres = 0.5 )
	D_hist = 1;
	numZonesPre = length(idMinLstPre);
	if numZonesPre != length(idMaxLstPre) || numZonesPre != length(preConfigZoneLst)
		error( "idMin, idMax, or configZoneLst length mismatched" );
	end
	updaterType = SingleUpdater;
	configAuxData = BConfigOfLnkValZonedAuxData{nDim}( idMinLst, idMaxLst, 0, 0, 0 );
	flipProposer = OneFlipProposer();	
	
	params = ParamsLoops( divNum, nDim );
	
	histDosType = WLHistDosZonedInE{nDim,D_hist};
	histDosFullType = WLHistDosFull{nDim,D_hist};
	flipCheckerType = WL2dHistStructFlipChecker{histDosType};

	updaterLst = [ updaterType(params) for iE = 1 : numZonesPre, iW = 1 : numWalksEach ];
	flipCheckerLst = [ WLFriendsFlipChecker{flipCheckerType}( divNum, nDim; histStdRatioThres = 0.15, wlResetInterval = wlResetInterval, histCutoffThres = histCutoffThres, wlHistDosArgs = ( idMinLstPre[iE], idMaxLstPre[iE] ) ) for iE = 1 : numZonesPre, iW = 1 : numWalksEach ];
	
	histDosFinalLst = getHistDos.( @views( flipCheckerLst[:,1] ) );	
	
	initializerLst = [ PresetInitializer( preConfigZoneLst[ iE ] ) for iE = 1 : numZonesPre ];
	auxDataLst = [ configAuxData for iZ = 1 : numZonesPre, iW = 1 : numWalksEach ];
	
	@views for iE = 1 : numZonesPre, iW = 1 : numWalksEach
		setFriendCheckerLst!( flipCheckerLst[iE, iW], flipCheckerLst[iE,:]... );
	end
	
	itController = FindConfigReplicaItController( configAuxData, flipCheckerLst, itExchange );
	fMainLst = Vector{String}(undef,numZonesPre);
	
	println("Loops_MC zoning starts:")
	fName = loops_MC_NoPrefabHelper_Replica( numZonesPre, numWalksEach; params = params, updaterLst = updaterLst, flipCheckerLst = flipCheckerLst, flipProposer = flipProposer, initializerLst = initializerLst, auxDataLst = auxDataLst, itController = itController );
	
	itNum = itController.itRef[];
	
	return configAuxData.configArr;
end
