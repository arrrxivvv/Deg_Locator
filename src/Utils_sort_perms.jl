function permute1d!( lst, lstTmp, ixLst )
	lstTmp .= lst;
	for id = 1 : length(lst)
		lst[id] = lstTmp[ixLst[id]];
	end
end

function permuteCol2d!( arr, arrTmp, ixLst )
	arrTmp .= arr;
	ln1 = size( arr, 1 );
	ln2 = size( arr, 2 );
	# id = 1;
	# @infiltrate
	# selectdim( arr, pDim, id );
	# selectdim( arr, pDim, 1 ) .= selectdim( arrTmp, pDim, 1 );
	for i2 = 1 : ln2
		for i1 = 1 : ln1
			arr[i1,i2] = arrTmp[i1,ixLst[i2]];
		end
	end
end

function permuteArrNd!( arr, arrTmp, pDim, ixLst )
	arrTmp .= arr;
	lnDim = size( arr, pDim );
	id = 1;
	# @infiltrate
	# selectdim( arr, pDim, id );
	selectdim( arr, pDim, 1 ) .= selectdim( arrTmp, pDim, 1 );
	for id = 1 : lnDim
		selectdim( arr, pDim, id ) .= selectdim( arrTmp, pDim, ixLst[id] );
		# selectdim( arr, pDim, id ) .= selectdim( arrTmp, pDim, id );
		# selectdim( arr, pDim, id ) .= getindex.( Ref( selectdim( arrTmp, pDim, id ) ), ixLst );
	end
end

function SortABbyA( key, object )
	indSort = sortperm( key );
	key .= key[indSort];
	object .= object[indSort,..];
	return key, object;
end

function EigenSort!( valLSt, vecLst )
	indSort = sortperm( valLSt );
	valLSt = valLSt[indSort];
	vecLst = vecLst[..,indSort];
	# return valLSt, vecLst
end
# function EigenSort!( valLSt, vecLst )
	# indSort = sortperm( valLSt );
	# valLSt .= valLSt[indSort];
	# vecLst .= vecLst[..,indSort];
	# # return valLSt, vecLst
# end