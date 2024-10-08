module ImgProcessing

using ImageFiltering
using OffsetArrays
using Utils

struct PrefabbedFilts
	sobelFiltLst::Vector{OffsetMatrix{Int64, Matrix{Int64}}};
	gaussFiltLstNoOffset::Vector{Matrix{Float64}};
	gaussFiltLst::Vector{OffsetMatrix{Float64, Matrix{Float64}}};

	function PrefabbedFilts()
		sobelFiltLstRaw = [ zeros(Int64, 3,3) for ii = 1 : 2 ];
		sobelFiltLstBase = [1, 2, 1];
		
		idSobelLst = [1,3];
		
		for ii = 1 : 2
			sobelFiltLstRaw[1][:,idSobelLst[ii]] .= sobelFiltLstBase .* (-1)^(ii-1);
			sobelFiltLstRaw[2,][idSobelLst[ii],:] .= sobelFiltLstBase .* (-1)^(ii-1);
		end
		
		sobelFiltLst = reflect.( centered.( sobelFiltLstRaw ) );
		
		gaussFiltLstRaw = Vector{Matrix{Float64}}(undef,4);
		
		gaussFiltLstRaw[1] = ones(1,1);
		gaussFiltLstRaw[2] = [1.0 2 1; 2 4 2; 1 2 1];
		gaussFiltLstRaw[2] .= gaussFiltLstRaw[2] ./ 16;
		gaussFiltLstRaw[3] = [ 1 4 7 4 1; 4 16 26 16 4; 7 26 41 26 7; 4 16 26 16 4; 1 4 7 4 1 ];
		gaussFiltLstRaw[3] .= gaussFiltLstRaw[3] ./ 273;
		gaussFiltLstRaw[4] = [ 0 0 1 2 1 0 0; 0 3 13 22 13 3 0; 1 13 59 97 59 13 1; 2 22 97 159 97 22 2; 1 13 59 97 59 13 1; 0 3 13 22 13 3 0; 0 0 1 2 1 0 0; ];
		gaussFiltLstRaw[4] .= gaussFiltLstRaw[4] ./ 1003;
		
		gaussFiltLst = centered.( gaussFiltLstRaw );
		
		new( sobelFiltLst, gaussFiltLstRaw, gaussFiltLst );
	end
end

const prefabbedFilts = PrefabbedFilts();

const sobelFiltLst = prefabbedFilts.sobelFiltLst;
const gaussFiltLst = prefabbedFilts.gaussFiltLst;
const gaussFiltLstNoOffset = prefabbedFilts.gaussFiltLstNoOffset;

function genBoxFit( sz::Int64 )
	len = 2*sz-1;
	val = 1 / len^2;
	return centered( fill( val, ( len, len ) ) );
end

const boxFiltLst = genBoxFit.( [1:4;] );


function convolute!( resultArr::AbstractArray{<:Number}, imgArr::AbstractArray{<:Number}, filtReflArr::AbstractArray{<:Number} )
	imfilter!(resultArr, imgArr, filtReflArr, "circular");
end

function nonMaxSuppress!( maxSuppressedArr::AbstractArray{<:Number}, imgArr::AbstractArray{<:Number}; ln::Int64 = 1 )
	maxVal = 0;
	sz1, sz2 = size(imgArr);
	for jj = 1 : sz2
		for ii = 1 : sz1
			maxVal = imgArr[ii,jj];
			for j2 = jj - ln : jj + ln, i2 = ii - ln : ii + ln
				jFilt = wrapIntInd( j2, sz2 );
				iFilt = wrapIntInd( i2, sz1 );
				
				maxVal = max( maxVal, imgArr[iFilt, jFilt] );
			end
			if maxVal != imgArr[ii,jj]
				maxSuppressedArr[ii,jj] = 0;
			else
				maxSuppressedArr[ii,jj] = imgArr[ii,jj];
			end
		end
	end
end

end
