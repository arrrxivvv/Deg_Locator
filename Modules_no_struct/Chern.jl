using LinearAlgebra
using EllipsisNotation

function SortABbyA( key, object )
	indSort = sortperm( key );
	key .= key[indSort];
	object .= object[indSort,..];
	return key, object
end

function EigenSort!( valLSt, vecLst )
	indSort = sortperm( valLSt );
	valLSt .= valLSt[indSort];
	vecLst .= vecLst[..,indSort];
	# return valLSt, vecLst
end

function movedPosFun( arr, bndShft )
	bndShftDir = sign(bndShft);
	if bndShftDir > 0
		movedPos = 1;
	else 
		movedPos = length(arr);
	end
end

function RollE_Vec( eLst, vecLst, bndShft )
	bndShftDir = sign(bndShft);
	movedPos = movedPosFun( eLst, bndShft );
	
	for rep = 1 : abs(bndShft)
		eLstSh = circshift( eLst, bndShftDir );
		eLst[movedPos] = eLst[movedPos] -bndShftDir * 2 * pi;
		vecLstSh = circshift( vecLst, bndShftDir );
	end
	return eLstSh, vecLstSh
end

function ShftOrNot( ShftFunOpt::Function, quasiELst, quasiELstTmp, ind1, ind2, bndShft )
	movedPos = movedPosFun( eLst, bndShft );
	indsThis = [ind1 ind2;]
	if  (ind2 != 0)
		prevDir = 1;
	else
		prevDir = 0;
	end
	indsPrev = deepcopy( indsThis );
	indsPrev[prevDir] -= 1;
	
	quasiELstPrev = quasiELst[indsPrev[1], indsPrev[2], :];

	return ShftFunOpt( quasiELstPrev, quasiELstTmp, movedPos, bndShft )
end

function ShftOrNotSum( quasiELstPrev, quasiELstTmp, bndShft )
	movedPos = movedPosFun( eLst, bndShft )
	quasiELstShft, _ = RollE_Vec( quasiELstTmp, quasiELstTmp, bndShft );
	
	diffVal = abs( sum( quasiELstTmp - quasiELstPrev ) );
	diff_shft = abs( sum( quasiELstShft- quasiELstPrev ) );
	
	return ( diff_shft < diffVal )
end

function Chern( N, U_mat_fun::Function, ang_divide, isFlq::Bool = true, isBandMatch::Bool = true )
	k_num = N;
	ang_ln = ang_divide + 1;
	ang_step = ang2pi / ang_divide;
	ang_lst = [0 : ang_step : 2*pi];
	
	U_mat_temp = zeros( Complex{Float64}, (k_num,k_num) );
	quasiELst = zeros( Complex{Float64}, (ang_ln, ang_ln, N) );
	eVecLst = zeros( Complex{Float64}, (ang_ln, ang_ln, N, N) );
	
	bndShftMem = 0;
	
	for ind1 = 1 : ang_ln
		ang1 = ang_lst[ind1];
		for ind1 = 1 : ang_ln
			ang2 = ang_lst[ind2];
			U_mat_temp = U_mat_fun( ang1, ang2 );
			
			eValLstTmp,eVecLstTmp = eigen(U_mat_temp);
			eVecLstTmp = transpose(eVecLstTmp);
			if isFlq
				quasiELstTmp = log.(eValLstTmp) / (-1im);
			else
				quasiELstTmp = eValLstTmp;
			end
			quasiELstImTmp = imag.(quasiELstTmp);
			quasiELstTmp = real.(quasiELstTmp);
			#quasiELstTmp = eValLstTmp
			
			# (quasiELst[ind1][ind2], eVecLst[ind1][ind2]) = SortABbyA( quasiELstTmp, eVecLstTmp )			
			SortABbyA( quasiELstTmp, eVecLstTmp );

			quasiELstTmp, eVecLstTmp = RollE_Vec( quasiELstTmp, eVecLstTmp, bndShftMem );
			
			if isFlq && isBandMatch && N > 2 && !(ind1==0 && ind2==0)
				for bndShft = [1 -1;]
					if ShftOrNot( ShftOrNotSum, quasiELst, quasiELstTmp, ind1, ind2, bndShft )
						quasiELstTmp, eVecLstTmp = RollE_Vec( quasiELstTmp, eVecLstTmp, bndShft );
						bndShftMem += bndShft;
						break
					end
				end
			end
			
			quasiELst[ind1,ind2,:] .= quasiELstTmp;
			eVecLst[ind1,ind2,:,:] .= eVecLstTmp;
		end
	end

	link1Lst = zeros( Complex(Float64), (ang_ln, ang_ln-1) );
	link2Lst = zeros( Complex(Float64), (ang_ln-1, ang_ln) );	
			
	chernLst = zeros( Complex(Float64),(N) );
	uProdLst = zeros( Complex(Float64), (ang_ln,ang_ln,N) );
	f12Lst = zeros( Complex(Float64), (ang_ln,ang_ln,N) );
	for indN = 1:N
		f12Sum = 0;
		for ind1 = 1:ang_ln
			for ind2 = 1:ang_ln-1
				dotProd = dot( eVecLst[ind1,ind2,indN], eVecLst[ind1,ind2+1,indN] );
				link1Lst[ind1,ind2] = dotProd / abs( dotProd )	;
			end
		end
		for ind1 = 1:ang_ln-1
			for ind2 = 1:ang_ln
				dotProd = dot( eVecLst[ind1,ind2,indN], eVecLst[ind1+1,ind2,indN] );
				link2Lst[ind1,ind2] = dotProd / abs( dotProd )	;
			end
		end
		
		for ind1 = 1:ang_ln-1
			for ind2 = 1:ang_ln-1
				f12 = log( link1Lst[ind1,ind2] / link1Lst[ind1+1,ind2] * link2Lst[ind1,ind2+1] / link2Lst[ind1,ind2] );
				f12Sum = f12Sum + f12;
				
				uProdLst[ind1,ind2,indN] = link1Lst[ind1,ind2] / link1Lst[ind1+1,ind2] * link2Lst[ind1,ind2+1] / link2Lst[ind1,ind2];
				f12Lst[ind1,ind2,indN] = f12;
			end
		end
		chernLst[indN] = f12Sum / (2 * pi * 1im );
	end

	if isFlq
		chernLst += 1 / N;
	end
	return (chernLst, quasiELst, eVecLst)
end
