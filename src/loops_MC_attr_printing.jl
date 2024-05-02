#Loops_MC funcs
#printing attribute values

function summarizeArrAttr( arr::Array )
	if isempty(arr)
		throw( ArgumentError( "array is empty" ) );
	elseif length(arr) == 1
		return arr[1];
	elseif length(arr) > 1
		return arr[[1,end]];
	end
end

function roundKeepInt( num::Number; digits )
	if isinteger(num)
		return Int64(num);
	else 
		return round(num; digits = digits);
	end
end

function summarizeArrAttrRound( arr::Array{<:Number}; digs = 3 )
	arrAttr = summarizeArrAttr( arr );
	if !(eltype(arr) <: Int)
		arrAttr = roundKeepInt.( arrAttr; digits = digs );
	end
	return arrAttr;
end

function summarizeArrBoolAttr( arr::Array{<:Union{Bool,Int}} )
	if isempty(arr)
		throw( ArgumentError( "array is empty" ) );
	elseif length(arr) > 2 
		throw( ArgumentError( "array should be true/false or +1/-1 sign" ) );
	elseif length(arr) == 1
		return arr[1];
	else
		return "both";
	end
end
