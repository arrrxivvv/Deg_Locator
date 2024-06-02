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

storeAuxDataSample( noAuxData::NoAuxData, itSample::Int64, it::Int64 ) = nothing;
storeAuxDataSampleDataOnly( noAuxData::NoAuxData, itSample::Int64, it::Int64 ) = nothing;

storeAuxDataStartSample( noAuxData::NoAuxData, itStartSample::Int64 ) = nothing;

storeAuxDataNum( noAuxData::NoAuxData, it::Int64 ) = nothing;

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

function storeAuxDataSampleDataOnly( zakAuxData::ZakArrAuxData, itSample::Int64 )
	if itSample > length( zakAuxData.zakArrSampleLst )
		push!( zakAuxData.zakArrSampleLst, similar(zakAuxData.zakArr) );
	end
	zakAuxData.zakArrSampleLst[itSample] .= zakAuxData.zakArr;
end

function storeAuxDataStartSample( zakAuxData::ZakArrAuxData, itStartSample::Int64 )
	if itStartSample > length( zakAuxData.zakArrSampleLst )
		push!( zakAuxData.zakArrStartSampleLst, similar(zakAuxData.zakArr) );
	end
	zakAuxData.zakArrStartSampleLst[itStartSample] .= zakAuxData.zakArr;
end

function storeAuxDataNum( zakAuxData::ZakArrAuxData, it::Int64 )
	if it > length(zakAuxData.zakMeanLst)
		push!( zakAuxData.zakMeanLst, zeros(zakAuxData.nDim) );
	end
	@views for dim = 1 : zakAuxData.nDim
		zakMeanLst[it][dim] = mean( zakAuxData.zakArr[:,:,dim] );
	end
end
