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

storeAuxDataSample( noAuxData::NoAuxData, itSample::Int64 ) = nothing;

storeAuxDataStartSample( noAuxData::NoAuxData, itStartSample::Int64 ) = nothing;

storeAuxDataNum( noAuxData::NoAuxData, it::Int64 ) = nothing;

getJldVarSampleLst( noAuxData::NoAuxData ) = nothing;

getJldVarStartSampleLst( noAuxData::NoAuxData ) = nothing;

getJldVarNumLst( noAuxData::NoAuxData ) = nothing;

saveAuxDataAll( auxData::AuxData, attrLst::Vector{String}, valLst::Vector, fMod::String ) = nothing;




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

function storeAuxDataSample( zakAuxData::ZakArrAuxData, itSample::Int64 )
	zakAuxData.zakArrSampleLst[itSample] .= zakAuxData.zakArr;
end

function storeAuxDataStartSample( zakAuxData::ZakArrAuxData, itStartSample::Int64 )
	zakAuxData.zakArrStartSampleLst[itStartSample] .= zakAuxData.zakArr;
end

function storeAuxDataNum( zakAuxData::ZakArrAuxData, it::Int64 )
	@views for dim = 1 : zakAuxData.nDim
		zakMeanLst[it][dim] = mean( zakAuxData.zakArr[:,:,dim] );
	end
end
