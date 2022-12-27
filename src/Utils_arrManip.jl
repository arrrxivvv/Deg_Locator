function arrSlcLst( arr, iDim, sz )
	return [ selectdim( arr, iDim, iZ )
		for iZ = 1 : sz ];
end

function arrShAllLst( arr, nDim )
	return [ ShiftedArrays.circshift( arr, 
		ntuple( ( x -> (x == iDim ? -1 : 0) ), nDim ) ) 
		for iDim = 1 : nDim];
end 

function assignArrOfArrs!( arrLst1, arrLst2 )
	( (a,b) -> (a .= b) ).( arrLst1, arrLst2 );
end
