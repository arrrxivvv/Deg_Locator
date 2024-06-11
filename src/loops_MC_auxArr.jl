#Loops_MC module

function getAuxDataSummaryName( noAuxDataType::Type{NoAuxData} )
	return "noAuxData";
end

function getAuxDataNameLst( noAuxData::NoAuxData )
	return String[];
end

renewAuxJldVarLst!( noAuxData::NoAuxData ) = nothing;

function flipAuxData!( noAuxData::NoAuxData, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	nothing;
end

function calcAuxData!( noAuxData::NoAuxData, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	nothing;
end

# storeAuxDataSample( noAuxData::NoAuxData, itSample::Int64, it::Int64 ) = nothing;
# storeAuxDataSampleDataOnly( noAuxData::NoAuxData, itSample::Int64, it::Int64 ) = nothing;
storeAuxDataSampleNoBndCheck( noAuxData::NoAuxData, itSample::Int64 ) = nothing;
checkDataSampleExtend!( noAuxData::NoAuxData, itSample::Int64 ) = nothing;

storeAuxDataStartSampleNoBndCheck( noAuxData::NoAuxData, itStartSample::Int64 ) = nothing;
checkDataStartSampleExtend!( noAuxData::NoAuxData, itStartSample::Int64 ) = nothing;

storeAuxDataNumNoBndCheck( noAuxData::NoAuxData, it::Int64 ) = nothing;
checkDataNumExtend!( noAuxData::NoAuxData, itStartSample::Int64 ) = nothing;

getJldVarSampleLst( noAuxData::NoAuxData ) = nothing;

getJldVarStartSampleLst( noAuxData::NoAuxData ) = nothing;

getJldVarNumLst( noAuxData::NoAuxData ) = nothing;

saveAuxDataAll( auxData::AuxData, attrLst::Vector{String}, valLst::Vector, fMod::String ) = nothing;




ZakArrAuxData( params::ParamsLoops, flipChecker::FlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) = ZakArrAuxData( params, itNum, itNumSample, itNumStartSample );

genAuxData( zakAuxDataType::Type{ZakArrAuxData}, params::ParamsLoops, flipChecker::FlipChecker, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) = genAuxData( zakAuxDataType, params, itNum, itNumSample, itNumStartSample );
genAuxData( zakAuxDataType::Type{ZakArrAuxData}, params::ParamsLoops, itNum::Int64, itNumSample::Int64, itNumStartSample::Int64 ) = zakAuxDataType( params, itNum, itNumSample, itNumStartSample );

function getAuxDataSummaryName( zakAuxDataType::Type{ZakArrAuxData} )
	return "zakAux";
end

function getAuxDataNameLst( zakAuxDataType::Type{<:ZakArrAuxData} )
	return ["zakArr"];
end

function renewAuxDataOutLst!( zakAuxData::ZakArrAuxData )
	# @infiltrate
	zakAuxData.dataSampleOutLst[1] = cat( zakAuxData.dataSampleLst[1]...; dims = 4 );
	zakAuxData.dataStartSampleOutLst[1] = cat( zakAuxData.dataStartSampleLst[1]...; dims = 4 );
	zakAuxData.dataNumOutLst[1] = cat( zakAuxData.dataNumLst[1]...; dims = 2 );
	
	GC.gc();
end

function flipAuxData!( zakAuxData::ZakArrAuxData, flipProposer::FlipProposer, params::ParamsLoops, dim::Int64, pos::CartesianIndex{D}, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	# calcAuxData!( zakAuxData, params, BfieldLst, linkLst, linkFerroLst );
	nothing;
end

function calcAuxData!( zakAuxData::ZakArrAuxData, params::ParamsLoops, BfieldLst::Vector{Array{Bool,D}}, linkLst::Vector{Array{Bool,D}}, linkFerroLst::Matrix{Array{Bool,D}} ) where {D}
	zakAuxData.zakArr .= false;
	@views for iDim = 1 : params.nDim
		xorIdFun = (a, b) -> xor( a, BfieldLst[iDim][b] ) ;
		for iD = 1 : params.divLst[iDim]
			zakAuxData.zakArr[:,:,iDim] .= xorIdFun.( zakAuxData.zakArr[:,:,iDim], params.posSlcLst[iDim][iD] );
		end
	end
end

function storeAuxDataSampleNoBndCheck( zakAuxData::ZakArrAuxData, itSample::Int64 )
	zakAuxData.zakArrSampleLst[itSample] .= zakAuxData.zakArr;
end

# function checkDataSampleExtend!( zakAuxData::ZakArrAuxData, itSample::Int64 )
	# if itSample > length( zakAuxData.zakArrSampleLst )
		# push!( zakAuxData.zakArrSampleLst, similar(zakAuxData.zakArr) );
	# end
# end

function storeAuxDataStartSampleNoBndCheck( zakAuxData::ZakArrAuxData, itStartSample::Int64 )
	zakAuxData.zakArrStartSampleLst[itStartSample] .= zakAuxData.zakArr;
end

# function checkDataStartSampleExtend!( zakAuxData.ZakArrAuxData )
	# if itStartSample > length( zakAuxData.zakArrSampleLst )
		# push!( zakAuxData.zakArrStartSampleLst, similar(zakAuxData.zakArr) );
	# end
# end

function storeAuxDataNumNoBndCheck( zakAuxData::ZakArrAuxData, it::Int64 )
	@views for dim = 1 : zakAuxData.nDim
		zakAuxData.zakMeanLst[it][dim] = mean( zakAuxData.zakArr[:,:,dim] );
	end
end

# function checkDataNumExtend!( zakAuxData::ZakArrAuxData, it::Int64 )
	# if it > length(zakAuxData.zakMeanLst)
		# push!( zakAuxData.zakMeanLst, zeros(zakAuxData.nDim) );
	# end
# end




function getAuxDataSummaryName( auxDataType::Type{<:BLinkAuxData} )
	return "BLink";
end

function getAuxDataNameLst( auxDataType::Type{<:BLinkAuxData} )
	return ["Bfield","link","linkFerro"]
end

getBfieldLst( auxData::BLinkAuxData ) = auxData.dataLst[1];
getLinkLst( auxData::BLinkAuxData ) = auxData.dataLst[2];
getLinkFerroLst( auxData::BLinkAuxData ) = auxData.dataLst[3];

function calcAuxData!( auxData::BLinkAuxData )
	nothing;
end

function renewAuxDataOutLst!( auxData::BLinkAuxData )
	posLst = auxData.params.posLst;
	colonLst = ntuple(ii->Colon(),auxData.nDim);
	divTup = auxData.params.divTup;
	iDimLstLst = [ CartesianIndices( auxData.dataLst[iD] ) for iD = 1 : length(auxData.dataLst) ];
	dimLstLst = [ size(iDimLstLst[iD]) for iD = 1 : length(auxData.dataLst) ];
	itNum = length(auxData.dataNumLst[1]);
	itSampleNum = length(auxData.dataSampleLst[1]);
	itStartSampleNum = length(auxData.dataStartSampleLst[1]);
	
	for iD = 1 : length(auxData.dataLst)
		auxData.dataSampleOutLst[iD] = zeros( Bool, divTup..., dimLstLst[iD]..., itSampleNum );
		auxData.dataStartSampleOutLst[iD] = zeros( Bool, divTup..., dimLstLst[iD]..., itStartSampleNum );
		auxData.dataNumOutLst[iD] = zeros( Bool, itNum, dimLstLst[iD]... );
		Threads.@threads for it = 1 : itSampleNum 
			for iDim in iDimLstLst[iD]
				auxData.dataSampleOutLst[iD][colonLst..., iDim, it] = auxData.dataSampleLst[iD][it][iDim];
			end
		end
		Threads.@threads for it = 1 : itStartSampleNum
			for iDim in iDimLstLst[iD]
				auxData.dataStartSampleOutLst[iD][colonLst..., iDim, it] .= auxData.dataStartSampleLst[iD][it][iDim] ;
			end
		end
		Threads.@threads for it = 1 : itNum
			for iDim in iDimLstLst[iD]
				auxData.dataNumOutLst[iD][it, iDim] = auxData.dataNumLst[iD][it][iDim];
			end
		end
		# @infiltrate
	end
	# @infiltrate
	
	GC.gc();
end

function storeAuxDataSampleNoBndCheck( auxData::BLinkAuxData, itSample::Int64 )
	for iD = 1 : length(auxData.dataLst)
		( (a,b) -> (a .= b) ).(auxData.dataSampleLst[iD][itSample], auxData.dataLst[iD]);
	end
end

function storeAuxDataStartSampleNoBndCheck( auxData::BLinkAuxData, itStartSample::Int64 )
	# @infiltrate
	for iD = 1 : length(auxData.dataLst)
		( (a,b) -> (a .= b) ).(auxData.dataStartSampleLst[iD][itStartSample], auxData.dataLst[iD]);
	end
end

function storeAuxDataNumNoBndCheck( auxData::BLinkAuxData, it::Int64 )
	# @infiltrate
	for iDim = 1 : auxData.nDimB
		auxData.dataNumLst[1][it][iDim] = sum( auxData.dataLst[1][iDim] );
	end
	for iDim = 1 : auxData.nDim
		auxData.dataNumLst[2][it][iDim] = sum( auxData.dataLst[2][iDim] );
	end
	for iDimB = 1 : auxData.nDimB, iDimLayer = 1 : auxData.nDimLayer
		auxData.dataNumLst[3][it][iDimLayer,iDimB] = sum( auxData.dataLst[3][iDimLayer,iDimB] );
	end
end
