function funAxisDropDims( arr, fun; dims )
	return dropdims( fun( arr; dims = dims ); dims = Tuple(dims) );
end

function maxDropDims( arr; dims )
	return funAxisDropDims( arr, maximum; dims = dims );
end

function minDropDims( arr; dims )
	return funAxisDropDims( arr, minimum; dims = dims );
end

function sumDropDims( arr; dims )
	return funAxisDropDims( arr, sum; dims = dims );
end
