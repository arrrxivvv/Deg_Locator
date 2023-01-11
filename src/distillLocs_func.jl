using DegLocatorDiv
using Utils
using JLD
using JLD2
using UMat
# using Infiltrator
using Statistics

Mlst = [10];

function locLstPurify_detailedOutput( posLocLstRaw, negLocLstRaw, posNlst, param_divide; isReversed = false )
	itNum = size(posNlst, 1);
	N = size(posNlst, 2);
	nLst = ones(Int64, itNum, N);
	param_divide_lim = param_divide[1];
	
	locLstRawPol = [ posLocLstRaw, negLocLstRaw ];
	locLstPols = [ [ [ sortslices( locLstRawPol[iPol][it][n], dims=1 ) for n = 1:N ] for it = 1:itNum ] for iPol = 1 : 2 ];
	# @infiltrate
	
	pureNlst = zeros( Int64, itNum, N-1 );
	pureNlst[:,1] = posNlst[:,1];
	pureNlst[:,N-1] = posNlst[:,N];
	for n = 2 : N-2
		pureNlst[:,n] = posNlst[:,n] - pureNlst[:,n-1];
	end
	
	locLstPolsDirty = [ [ [ [ locLstPols[pol][it][n][id,:] for id = 1 : size(locLstPols[pol][it][n],1) ] for n = 1 : N ] for it = 1:itNum ] for pol = 1 : 2 ];
	locLstPolsPure = [ [ [ Vector{Vector{Int64}}(undef,0) for n = 1 : N-1 ] for it = 1:itNum ] for iPol = 1 : 2 ];
	if isReversed
		for pol = 1:2, it = 1 : itNum
			reverse!( locLstPolsDirty[pol][it] );
		end
	end
	
	polInvLst = [2,1];
	locDiffAbs = [0,0,0];
	diffLimMax = 3;
	diffLimExtended = 4;
	matchLst = [mX, mSameN, mPropagate];
	lnMatch = length(matchLst);
	isExtendedArr = [false, true];
	diffLimMaxArr = [[diffLimMax, diffLimExtended] for iM = 1:lnMatch];
	diffLimMaxArr[end] .+= [3,6];
	diffLimStartArr = [[0,2] for iM = 1:lnMatch];
	unMatchedLst = [];
	unMatchedNumPrev = 0;
	unMatchedClosest = [];
	isFinished = false;
	isExtended = false;
	
	lnDirtyPrevLst = zeros(Int64, 2);
	
	# Threads.@threads for it = 1 : itNum
	for it = 1 : itNum
		print("\r", it, ", ");
		isExtended = false;
		for iLim = 1 : lnMatch, iExt = 1:2
			isExtended = isExtendedArr[iExt];
			mMatch = matchLst[iLim];
			diffLimEnd = diffLimMaxArr[iLim][iExt];
			diffLimStart = diffLimStartArr[iLim][iExt];
			for diffLim = diffLimStart : diffLimEnd
				lnDirtyPrevLst .= 0;
				for n = 2 : N
					for pol = 1 : 2
						polInv = polInvLst[pol];
						
						locLstPrev = locLstPolsDirty[polInv][it][n-1];
						if mMatch == mPropagate
							nThisLst = [n,n-1,collect(n+1:N)...];
						elseif mMatch == mSameN
							nThisLst = [n-1];
						else
							nThisLst = [n]
						end
						for nThis in nThisLst
							locLstThis = locLstPolsDirty[pol][it][nThis];
							if mMatch == mSameN && pol == 2
								break;
							end
							idDirtyStart = 1;
							
							if mMatch == mExtended && lnDirtyPrevLst[pol] == 0 && length(locLstPrev) == length(locLstThis)
								# @infiltrate
								for id = 1 : length(locLstPrev)
									push!( locLstPolsPure[polInv][it][n-1], locLstPrev[id] );
								end
								empty!(locLstPrev);
								empty!(locLstThis);
							end
							
							idPure = 1;
							while idPure <= length( locLstPrev )
								locPure = locLstPrev[idPure];
								
								isMatched = false;
								idDirty = idDirtyStart;
								idDirtyRev = idRetard( idDirty, locLstThis );
								
								isMatched, idDirty, idDirtyRev, idDirtyStart = scanMatchFull!( locPure, locLstThis, idDirty, idDirtyRev, diffLim, param_divide, isExtended ); 
								
								if isMatched
									if mMatch != mSameN
										for nPush = (n-1) : (nThis-1)
											push!( locLstPolsPure[polInv][it][nPush], locPure );
										end
									end
									deleteat!( locLstPrev, idPure );
								else
									idPure += 1;
								end
							end
						end
					end
					for pol = 1 : 2
						lnDirtyPrevLst[pol] = length(locLstPolsDirty[pol][it][n-1]);
					end
				end
				isFinished = all( isempty.(locLstPolsDirty[1][it]) ) && all( isempty.( locLstPolsDirty[2][it] ) );
				if isFinished
					break;
				end
				# @infiltrate cond = it >= 2
			end
			if isFinished
				break;
			end
			# isExtended = true;
		end
	end
	# @infiltrate	
	
	pureNlstPolsFinal = [ zeros( Int64, itNum, N-1 ) for pol = 1 : 2 ];
	pureNlstPolsComb = [ zeros( Int64, itNum, N-1 ) for pol = 1 : 2 ];
	locLstPolsComb = deepcopy(locLstPolsPure);
	for pol = 1 : 2
		for it = 1 : itNum
			for n = 1 : N-1
				for id = 1 : length( locLstPolsDirty[pol][it][n] )
					locDirty = locLstPolsDirty[pol][it][n][id];
					for nPush = n:N-1
						push!( locLstPolsComb[pol][it][nPush],locDirty  );
					end
				end
			end
		end
	end
	# @infiltrate
	for pol = 1 : 2
		for it = 1 : itNum
			for n = 1 : N-1
				pureNlstPolsFinal[pol][it,n] = length( locLstPolsPure[pol][it][n] );
				pureNlstPolsComb[pol][it,n] = length( locLstPolsComb[pol][it][n] );
			end
		end
	end
	
	if isReversed
		for pol = 1:2, it = 1 : itNum
			reverse!( locLstPolsDirty[pol][it] );
			reverse!( locLstPolsPure[pol][it] );
		end
	end
	
	# @infiltrate
	locLstPolsPureFinal = [ [ [ isempty(locLstPolsPure[pol][it][n]) ? [] : copy( hcat( locLstPolsPure[pol][it][n]... )' ) for n = 1 : N-1 ] for it = 1 : itNum] for pol = 1 : 2];
	locLstPolsCombFinal = [ [ [ isempty(locLstPolsPure[pol][it][n]) ? [] : copy( hcat( locLstPolsComb[pol][it][n]... )' ) for n = 1 : N-1 ] for it = 1 : itNum] for pol = 1 : 2];
	# @infiltrate
	
	return locLstPolsPureFinal, pureNlstPolsFinal, locLstPolsPure, locLstPolsDirty, locLstPolsComb, pureNlstPolsComb, locLstPolsCombFinal;
end

function locLstPurify( posLocLstRaw, negLocLstRaw, posNlst, param_divide )
	locLstPolsPureFinal, pureNlstPolsFinal, locLstPolsPure, locLstPolsDirty, locLstPolsComb, pureNlstPolsComb, locLstPolsCombFinal = locLstPurify_detailedOutput( posLocLstRaw, negLocLstRaw, posNlst, param_divide );
	return locLstPolsCombFinal, pureNlstPolsComb;
end

function locLstPurify_20211130( posLocLstRaw, negLocLstRaw, posNlst, param_divide )
	itNum = size(posNlst, 1);
	N = size(posNlst, 2);
	nLst = ones(Int64, itNum, N);
	param_divide_lim = param_divide[1];
	
	locLstRawPol = [ posLocLstRaw, negLocLstRaw ];
	locLstPols = [ [ [ sortslices( locLstRawPol[iPol][it][n], dims=1 ) for n = 1:N ] for it = 1:itNum ] for iPol = 1 : 2 ];
	
	pureNlst = zeros( Int64, itNum, N-1 );
	pureNlst[:,1] = posNlst[:,1];
	pureNlst[:,N-1] = posNlst[:,N];
	for n = 2 : N-2
		pureNlst[:,n] = posNlst[:,n] - pureNlst[:,n-1];
	end
	
	# locLstPolsPure = [ [ [ zeros( Int64, pureNlst[it,n], 3 ) for n = 1:N-1 ] for it = 1:itNum ] for iPol = 1 : 2 ];
	locLstPolsPure = [ [ [ Vector{Vector{Int64}}(undef,0) for n = 1 : N-1 ] for it = 1:itNum ] for iPol = 1 : 2 ];
	locLstPolsPureFinal = [ [ [ Vector{Vector{Int64}}(undef,0) for n = 1 : N-1 ] for it = 1:itNum ] for iPol = 1 : 2 ];
	diffBaskets = [ [] for ix = 1 : param_divide_lim ];
	
	polInvLst = [2,1];
	locDiffAbs = [0,0,0]
	diffLimMax = 2;
	diffLimExtended = 5;
	unMatchedLst = [];
	unMatchedNumPrev = 0;
	unMatchedClosest = [];
	
	# Threads.@threads for it = 1 : itNum
	for it = 1 : itNum
		print("\r", it, ", ");
		for pol = 1 : 2
			locLstPolsPure[pol][it][1] = [ locLstPols[pol][it][1][id,:] for id = 1 : posNlst[it,1] ];
			locLstPolsPure[pol][it][N-1] = [ locLstPols[pol][it][N][id,:] for id = 1 : posNlst[it,N] ];;
		end
		for n = 2 : N - 2
			# print(n, ",");
			for pol = 1 : 2
				polInv = polInvLst[pol];
				
				lnDirty = posNlst[it,n];
				lnPure = pureNlst[it,n-1];
				idDirtyStart = 1;
				
				# @infiltrate cond = it >= 811
				pureLstNext = [ locLstPols[pol][it][n][id,:] for id = 1 : size(locLstPols[pol][it][n],1) ];
				locLstPolsPure[pol][it][n] = pureLstNext;
				empty!( unMatchedLst );
				empty!( unMatchedClosest );
				# @infiltrate
				
				for idPure = 1 : lnPure
					# print(idPure, ",")
					locPure = locLstPolsPure[polInv][it][n-1][idPure];
					# @infiltrate cond = n >= 4
					
					empty!.( diffBaskets );
					
					diffLim = 0;
					isMatched = false;
					isMatchedExtended = false;
					idDirty = idDirtyStart;
					idDirtyRev = idRetard( idDirty, pureLstNext );
					# @infiltrate cond = n>=4 && idPure >= 19
					for diffLimEnd = [diffLimMax, diffLimExtended]
						# @infiltrate cond = n>=4 && idPure >= 19
						while diffLim <= diffLimEnd && !isMatched
							isMatched, idDirty, idDirtyRev, idDirtyStart, locDiffAbs = scanMatchFull!( locPure, pureLstNext, idDirty, idDirtyRev, diffLim, param_divide_lim, diffBaskets, isMatchedExtended );
							diffLim += 1;
							# @infiltrate cond = n>=4 && idPure >= 19 
						end
						# @infiltrate cond = n>=4 && idPure >= 19 
						if isMatched
							break;
						else
							isMatchedExtended = true;
							diffLim = 1;
						end
					end
					
					# @infiltrate cond = n >= 4
					# @infiltrate cond = !isMatched;
					if isMatchedExtended
						push!( unMatchedLst, idPure );
						push!( unMatchedClosest, sort( locDiffAbs, rev=true ) );
					end
				end
				# @infiltrate cond = n>=5 && pol == 2
				popNum = min( length( locLstPolsPure[polInv][it][n-1] ) - pureNlst[it,n-1], length(unMatchedLst) );
				popIdLst = sortperm( unMatchedClosest, by = x -> [maximum(x[2:3]),minimum(x[2:3]),x[1]], rev=true );
				deleteat!( locLstPolsPure[polInv][it][n-1], unMatchedLst[1:popNum] );
				# @infiltrate
			end
		end
		# @infiltrate
	end
	
	locLstPolsPureFinal = [ [ [ copy( hcat( locLstPolsPure[pol][it][n]... )' ) for n = 1 : N-1 ] for it = 1 : itNum] for pol = 1 : 2];
	pureNlstPolsFinal = [ zeros( Int64, itNum, N-1 ) for pol = 1 : 2 ];
	for pol = 1 : 2
		for it = 1 : itNum
			for n = 1 : N-1
				pureNlstPolsFinal[pol][it,n] = length( locLstPolsPure[pol][it][n] );
			end
		end
	end
	# @infiltrate
	
	return locLstPolsPureFinal, pureNlstPolsFinal;
end

function idRetard( id, arr )
	idChecked = 1;
	if !isempty(arr)
		idChecked = mod( id - 1 - 1, length(arr) ) + 1;
	end
	return idChecked;
end

function idRecheck( id, arr )
	idChecked = 1;
	if !isempty(arr)
		idChecked = mod( id - 1, length(arr) ) + 1;
	end
	return idChecked;
end

function idAdvance( id, arr )
	idChecked = 1;
	if !isempty(arr)
		idChecked = mod( id - 1 + 1, length(arr) ) + 1;
	end
	return idChecked;
end

function diffWrap( x1, x2, param_divide_lim )
	sh = div(param_divide_lim - 1, 2);
	return mod( x1 - x2 + sh, param_divide_lim ) - sh;
end

function scanMatchFull!( locPure, pureLstNext, idDirty, idDirtyRev, diffLim, param_divide, isExtended )
	isMatched = false;
	
	isAdvArr = [true, false];
	idDirtyArr = [idDirty, idDirtyRev];
	idDirtyStart = idDirty;
	for iRev = 1 : 2
		# @infiltrate cond = n >= 1 && idPure >= 10
		if !isMatched
			isMatched, idDirtyArr[iRev], idDirtyStart, locDiff = scanMatch!( locPure, pureLstNext, idDirtyArr[iRev], diffLim, param_divide, isExtended; advOrRet = isAdvArr[iRev] );
		end
	end
	
	return isMatched, idDirtyArr[1], idDirtyArr[2], idDirtyStart;
end

function scanMatchFull!( locPure, pureLstNext, idDirty, idDirtyRev, diffLim, param_divide_lim, diffBaskets, isMatchedExtended )
	iDiff = 1;
	isMatched = false;
	while iDiff <= length( diffBaskets[diffLim+1] ) && !isMatched
		# @infiltrate cond = n >= 4 && idPure >= 6
		# @infiltrate cond = locPure == [53,13,5];
		idPureDiff = diffBaskets[diffLim+1][iDiff];
		isMatched, overflowed, idDirtyStart, locDiff = matchAndPop( locPure, pureLstNext, idPureDiff, diffLim, param_divide_lim, diffBaskets, isMatchedExtended; isFromBasket = true );
		iDiff += 1;
	end
	isAdv = true;
	if !isMatched
		isMatched, idDirty, idDirtyStart, locDiff = scanMatch!( locPure, pureLstNext, idDirty, diffLim, param_divide_lim, diffBaskets, isMatchedExtended; advOrRet = true );
	end
	# @infiltrate cond = n >= 1 && idPure >= 10

	isAdv = false;
	if !isMatched
		isMatched, idDirtyRev, idDirtyStart, locDiff = scanMatch!( locPure, pureLstNext, idDirtyRev, diffLim, param_divide_lim, diffBaskets, isMatchedExtended; advOrRet = false );
	end
	
	return isMatched, idDirty, idDirtyRev, idDirtyStart, locDiff;
end

function scanMatch!( locPure, pureLstNext, idDirty, diffLim, param_divide, isExtended; advOrRet )
	overflowed =false;
	isMatched = false;
	idDirtyStart = idDirty;
	idAdvRet = advOrRet ? idAdvance : idRetard;
	
	locDiff = zeros(Int64, 3);
	idCount = 0;
	while !isMatched && !overflowed
		isMatched, overflowed, idDirtyStart, locDiff = matchAndPop( locPure, pureLstNext, idDirty, diffLim, param_divide, isExtended; advOrRet = advOrRet );
		idCount += 1;
		# @infiltrate cond = locPure == [30,6,6] && diffLim >= 1
		if overflowed || idCount > length( pureLstNext )
			break;
		end
		idDirty = idAdvRet(idDirty, pureLstNext);
	end	
	return isMatched, idDirty, idDirtyStart, locDiff;
end

function scanMatch!( locPure, pureLstNext, idDirty, diffLim, param_divide_lim, diffBaskets, isMatchedExtended; advOrRet )
	overflowed =false;
	isMatched = false;
	idDirtyStart = idDirty;
	idAdvRet = advOrRet ? idAdvance : idRetard;
	
	locDiff = zeros(Int64, 3);
	while !isMatched && !overflowed
		isMatched, overflowed, idDirtyStart, locDiff = matchAndPop( locPure, pureLstNext, idDirty, diffLim, param_divide_lim, diffBaskets, isMatchedExtended;  advOrRet = advOrRet );
		# @infiltrate cond = locPure == [53,13,5];
		if overflowed
			break;
		end
		idDirty = idAdvRet(idDirty, pureLstNext);
	end	
	return isMatched, idDirty, idDirtyStart, locDiff;
end

function matchAndPop( locPure, pureLstNext, idDirty, diffLim, param_divide, isExtended; advOrRet=true )
	isMatched = false;
	overflowed = true;
	idDirtyStart = 1;
	locDiffAbs = zeros( Int64, 3 );
	if !isempty( pureLstNext )
		locDirty = pureLstNext[idDirty];
		idDirtyStart = idDirty;
		isMatched, dx, dMax, locDiffAbs = dirtyPureMatch( locPure, locDirty, diffLim, param_divide, isExtended );
		overflowed = false;
		if isMatched
			deleteat!( pureLstNext, idDirty );
			idDirtyStart = idRecheck( idDirty, pureLstNext );
		else
			overflowed = advOrRet ? (dx < -diffLim) : (dx > diffLim);
		end
	end
	return isMatched, overflowed, idDirtyStart, locDiffAbs;
end

function matchAndPop( locPure, pureLstNext, idDirty, diffLim, param_divide_lim, diffBaskets, isMatchedExtended; advOrRet=true, isFromBasket = false )
	locDirty = pureLstNext[idDirty];
	idDirtyStart = idDirty;
	isMatched, dx, dMax, locDiffAbs = dirtyPureMatch( locPure, locDirty, diffLim, param_divide_lim, isMatchedExtended );
	overflowed = false;
	if isMatched
		if !isMatchedExtended
			deleteat!( pureLstNext, idDirty );
		end
		idDirtyStart = idRecheck( idDirty, pureLstNext );
	else
		overflowed = advOrRet ? (dx < -diffLim) : (dx > diffLim);
		if !overflowed && !isFromBasket
			push!( diffBaskets[dMax+1], idDirty );
		end
	end
	return isMatched, overflowed, idDirtyStart, locDiffAbs;
end

function locLstPurify_oldBaskets( posLocLstRaw, negLocLstRaw, posNlst, param_divide )
	itNum = size(posNlst, 1);
	N = size(posNlst, 2);
	nLst = ones(Int64, itNum, N);
	param_divide_lim = param_divide[1];
	
	locLstRawPol = [ posLocLstRaw, negLocLstRaw ];
	locLstPols = [ [ [ sortslices( locLstRawPol[iPol][it][n], dims=1 ) for n = 1:N ] for it = 1:itNum ] for iPol = 1 : 2 ];
	
	pureNlst = zeros( Int64, itNum, N-1 );
	pureNlst[:,1] = posNlst[:,1];
	pureNlst[:,N-1] = posNlst[:,N];
	for n = 2 : N-2
		pureNlst[:,n] = posNlst[:,n] - pureNlst[:,n-1];
	end
	
	locLstPolsPure = [ [ [ zeros( Int64, pureNlst[it,n], 3 ) for n = 1:N-1 ] for it = 1:itNum ] for iPol = 1 : 2 ];
	diffBaskets = [ [] for ix = 1 : param_divide_lim ];
	
	polInvLst = [2,1];
	
	# Threads.@threads for it = 1 : itNum
	for it = 1 : itNum
		print("\r", it, ", ");
		for pol = 1 : 2
			locLstPolsPure[pol][it][N-1] .= locLstPols[pol][it][N];
		end
		for n = 2 : N - 2
			print(n, ",");
			for pol = 1 : 2
				polInv = polInvLst[pol];
				
				lnDirty = posNlst[it,n];
				lnPure = pureNlst[it,n-1];
				idDirty = 1;
				idPure2 = 1;
				idPureTrace = idPure2;
				
				pureLstNext = [];
				
				for idPure = 1 : lnPure
					print(idPure, ",")
					locPure = locLstPolsPure[polInv][it][n-1][idPure,:];
					# @infiltrate cond = n >= 4
					
					for iDiff = 1 : param_divide_lim
						diffBaskets[iDiff] = [];
					end
					
					idPureTmp = 1;
					diffLim = 1;
					diffLimMax = 3;
					isMatched = false;
					while diffLim <= diffLimMax && !isMatched
						idPureTrace = idPure2 - 1;
						iDiff = 1;
						while iDiff <= length( diffBaskets[diffLim+1] ) && !isMatched
							# @infiltrate cond = n >= 4 && idPure >= 6
							idPureDiff = diffBaskets[diffLim+1][iDiff];
							locDirty = pureLstNext[idPureDiff];
							isMatched, dx = dirtyPureMatch( locPure, locDirty, diffLim );
							if isMatched
								idPure2 = popNextLst!( pureLstNext, idPureDiff, idPure2 );
							end
							iDiff += 1;
						end
						overflowed =false;
						while idDirty <= lnDirty && !isMatched && !overflowed
							locDirty = locLstPols[pol][it][n][idDirty,:];
							isMatched, dx, dMax = dirtyPureMatch( locPure, locDirty, diffLim );
							if !isMatched
								overflowed = dx < -diffLim;
								if overflowed
									break;
								end
								push!( diffBaskets[dMax+1], idPure2 );
								idPure2 = pushNextLst!( pureLstNext, locDirty, idPure2 );
							end
							idDirty += 1;
						end
						
						overflowed = false;
						while idPureTrace >= 1 && !isMatched && !overflowed
							locDirty = pureLstNext[idPureTrace];
							isMatched, dx, dMax = dirtyPureMatch( locPure, locDirty, diffLim );
							if isMatched
								idPure2 = popNextLst!( pureLstNext, idPureTrace, idPure2 );
							else
								overflowed = dx > diffLim;
								push!( diffBaskets[dMax+1], idPureTrace );
							end
							idPureTrace -= 1;
						end
						diffLim += 1;
					end
					# @infiltrate cond = n >= 4
				end
				while idDirty <= lnDirty
					locDirty = locLstPols[pol][it][n][idDirty,:];
					push!( pureLstNext, locDirty );
					idDirty += 1;
				end
				
				for iPure = 1 : length( pureLstNext )
					locLstPolsPure[pol][it][n][iPure,:] .= pureLstNext[iPure];
				end
				# @infiltrate
			end
		end
	end
	
	return locLstPolsPure[1], locLstPolsPure[2], pureNlst;
end

function dirtyPureMatch( locPure, locDirty, diffLim, param_divide, isExtended )
	locDiff = ( (x1,x2,divide) -> diffWrap( x1, x2, divide ) ).(locPure, locDirty, param_divide);
	locDiffAbs = abs.(locDiff);
	diffLimY = min( diffLim, 1 );
	if !isExtended
		diffVec = [diffLim, diffLimY, diffLimY];
	else
		diffVec = [diffLim, diffLim, diffLim];
	end
	isMatched = all( locDiffAbs .<= diffVec );
	dx = locDiff[1];
	dMax = maximum( locDiffAbs );
	# dMax = locDiffAbs[1];
	return isMatched, dx, dMax, locDiffAbs;
end

# function dirtyPureMatch( locPure, locDirty, diffLim, param_divide_lim, isExtended )
	# locDiff = ( (x1,x2) -> diffWrap( x1, x2, param_divide_lim ) ).(locPure, locDirty);
	# locDiffAbs = abs.(locDiff);
	# diffLimY = min( diffLim, 1 );
	# if !isExtended
		# diffVec = [diffLim, diffLimY, diffLimY];
	# else
		# diffVec = [diffLim, diffLim, diffLim];
	# end
	# isMatched = all( locDiffAbs .<= diffVec );
	# dx = locDiff[1];
	# dMax = maximum( locDiffAbs );
	# # dMax = locDiffAbs[1];
	# return isMatched, dx, dMax, locDiffAbs;
# end

function pushNextLst!( pureLstNext, locDirty, idPure2 )
	push!( pureLstNext, locDirty );
	return idPure2 + 1;
end

function popNextLst!( pureLstNext, idDel, idPure2 )
	deleteat!( pureLstNext, idDel );
	return idPure2 - 1;
end

function locLstPurify_old1( posLocLstRaw, negLocLstRaw, posNlst )
	itNum = size(posNlst, 1);
	N = size(posNlst, 2);
	nLst = ones(Int64, itNum, N);
	
	locLstRawPol = [ posLocLstRaw, negLocLstRaw ];
	locLstPols = [ [ [ sortslices( locLstRawPol[iPol][it][n], dims=1 ) for n = 1:N ] for it = 1:itNum ] for iPol = 1 : 2 ];
	
	pureNlst = zeros( Int64, itNum, N-1 );
	pureNlst[:,1] = posNlst[:,1];
	pureNlst[:,N-1] = posNlst[:,N];
	for n = 2 : N-2
		pureNlst[:,n] = posNlst[:,n] - pureNlst[:,n-1];
	end
	
	locLstPolsPure = [ [ [ zeros( Int64, pureNlst[it,n], 3 ) for n = 1:N-1 ] for it = 1:itNum ] for iPol = 1 : 2 ];
	# negLocLstPure = [ [ zeros( Int64, pureNlst[it,n], 3 ) for n = 1:N-1 ] for it = 1:itNum ];
	
	# locLstPols = [posLocLst, negLocLst];
	# locLstPolsPure = [posLocLstPure, negLocLstPure];
	polInvLst = [2,1];
	# @infiltrate
	
	# Threads.@threads for it = 1 : itNum
	for it = 1 : itNum
		print("\r", it, ", ");
		for pol = 1 : 2
			locLstPolsPure[pol][it][1,:] = locLstPols[pol][it][1,:];
			locLstPolsPure[pol][it][N-1,:] = locLstPols[pol][it][N,:];
		end
		for n = 2 : N - 2
			print(n, ",");
			for pol = 1 : 2
				polInv = polInvLst[pol];
				
				backLogLocs = [ [] for it = 1 : 3 ];
				backLogsNext = [ [] for it = 1 : 3 ];
				backLogsExpired = [];
				
				lnDirty = posNlst[it,n];
				lnPure = pureNlst[it,n-1];
				idDirty = 1;
				idPure2 = 1;
				# locDirty = locLstPols[pol][it][n][idDirty,:];
				idPure = 1;
				# locPure = locLstPolsPure[polInv][it][n-1][idPure,:];
				# xPure = locPure[1];
				xPure = 1;
				xPurePrev = 1;
				for idPure = 1 : lnPure
					print(idPure, ",")
					locPure = locLstPolsPure[polInv][it][n-1][idPure,:];
					xPure = locPure[1];
					# @infiltrate cond = n >= 4
					if idPure > 1
						if xPure - xPurePrev == 1
							backLogLocs[1] = backLogsNext[2];
							backLogLocs[2] = backLogsNext[3];
							backLogLocs[3] = [];
							expiredEnd = 1;
						elseif xPure - xPurePrev == 2
							backLogLocs[1] = backLogsNext[3];
							backLogLocs[2] = [];
							backLogLocs[3] = [];
							expiredEnd = 2;
						elseif xPure == xPurePrev
							backLogLocs .= backLogsNext;
							expiredEnd = 0;
						else
							backLogLocs = [ [] for iBack = 1 : 3 ];
							expiredEnd = 3;
						end
						idPure2 = flushBackLog!( backLogsNext, backLogsExpired, expiredEnd, locLstPolsPure[pol][it][n], idPure2 );
						backLogsNext = [ [] for iBack = 1 : 3 ];
					end
					isMatched = false;
					for iBack = 1 : 3
						for jj = 1 : length( backLogLocs[iBack] )
							locDirty = backLogLocs[iBack][jj];
							if isMatched
								push!( backLogsNext[iBack], locDirty );
							else
								isMatched, idPure2 = dirtyPureMatch!( locPure, locDirty, backLogsNext, locLstPolsPure[pol][it][n], idPure2 );
							end
						end
					end
					# @infiltrate cond = n >= 4
					
					while idDirty <= lnDirty && !isMatched
						locDirty = locLstPols[pol][it][n][idDirty,:];
						isMatched, idPure2 = dirtyPureMatch!( locPure, locDirty, backLogsNext, locLstPolsPure[pol][it][n], idPure2 );
						idDirty += 1;
					end
					
					xPurePrev = xPure;
					# @infiltrate cond = n >= 4
				end
				expiredEnd = 3;
				idPure2 = flushBackLog!( backLogsNext, backLogsExpired, expiredEnd, locLstPolsPure[pol][it][n], idPure2 );
				# @infiltrate cond = n >= 4
				idPure2 = flushToPure!( [ locLstPols[pol][it][n][jj,:] for jj = idDirty : lnDirty ], locLstPolsPure[pol][it][n], idPure2 );
				# @infiltrate cond = n >= 4
			end
		end
	end
	
	return locLstPolsPure[1], locLstPolsPure[2], pureNlst;
end

function dirtyPureMatch_old!( locPure, locDirty, backLogsNext, pureLstNext, idPure2 )
	isMatched = false;
	# xPure = locPure[1];
	logId = 1;
	if all( abs.( locDirty - locPure ) .<= 1 )
		isMatched = true;
	else
		isMatched = false;
		dx = locPure[1] - locDirty[1];
		if abs(dx) >= 2
			pureLstNext[idPure2,:] = locDirty;
			idPure2 += 1;
		else
			if dx == 1
				logId = 1;
			elseif dx == 0
				logId = 2;
			elseif dx == -1
				logId = 3;
			end
			push!(backLogsNext[logId], locDirty);
		end
	end
	return isMatched, idPure2;
end

function flushBackLog!( backLogsNext, backLogsExpired, expiredEnd, lstPure, idPure2 )
	backLogsExpired = collect( Iterators.flatten( backLogsNext[1:expiredEnd] ) );
	idPure2 = flushToPure!( backLogsExpired, lstPure, idPure2 );
	
	return idPure2;
end

function flushToPure!( flushedLst, lstPure, idPure2 )
	# @infiltrate cond = idPure2 >= 28
	for iFl = 1 : length(flushedLst)
		lstPure[idPure2,:] = flushedLst[iFl];
		idPure2 += 1;
	end
	# @infiltrate cond = idPure2 >= 28
	return idPure2
end

function distillLocsFromFile( avgNum, N, param_divide, seed; fMain = "deg", fMod = "", fExt = jld2Type, dim = 3, pow = 2, avgNumLoc = avgNum, scale = degOptDefaultLst[1], ratio = degOptDefaultLst[2], alpha = degOptDefaultLst[3], enumSaveMem = memNone, thresNM = rtFndDefaultLst[1], thresEDeg = rtFndDefaultLst[2], isNewFlux = false )
	attrLst, valLst = fAttrOptLstFunc( N, param_divide, avgNum, seed; dim = dim, scale = scale, ratio = ratio, alpha = alpha, enumSaveMem = enumSaveMem, thresNM = thresNM, thresEDeg = thresEDeg );
	varFileName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	@info(varFileName);
	if isNewFlux
		locLstPols = load( varFileName, "locLstPol" );
		locNPol = load( varFileName, "NLstPol" );
		posLocLst = locLstPols[1];
		negLocLst = locLstPols[2];
		posNlst = locNPol[1];
		negNlst = locNPol[2];
	else
		posLocLst = load(varFileName, "posLocLst");
		negLocLst = load(varFileName, "negLocLst");
		posNlst = load(varFileName, "posNlst");
		negNlst = load(varFileName, "negNlst");
		locLstPols = [posLocLst, negLocLst];
		locNPol = [posNlst, negNlst];
	end
	
	# @infiltrate
	locNCumPol = distillLocN( avgNumLoc, N, locNPol );
	
	dim = 3;
	
	pureLocLstPols, pureNlstPols = locLstPurify( posLocLst, negLocLst, posNlst, param_divide );
	
	pureNlstAvg = vcat( mean.( pureNlstPols, dims=1 )... );
	
	# @infiltrate
	# locLstPurePol = [posLocPureLst, negLocPureLst];
	outFileName = fNameFunc( fMain * "_locDistilled", attrLst, valLst, fExt; fMod = fMod );
	
	# outFileName = string( fileNameMain, "locDistilled", "_", fileNameAttr, jld2Type );
	# jldsave( outFileName; locDistilledLst = pureLocLstPols, NDistilledLst = pureNlstPols );
	save( outFileName, "locDistilledLst", pureLocLstPols, "NDistilledLst", pureNlstPols );
	
	pureNlstName = fNameFunc( "gDistilled", attrLst, valLst, fExt; fMod = fMod );
	# pureNlstName = string( "gDistilled", "_", fileNameAttr, npyType );
	npzwrite( pureNlstName, pureNlstAvg );
	
	return pureLocLstPols, pureNlstPols;
end

function distillLocsFromFile_byEnergy( avgNum, N, param_divide, seed; fileNameMain = "deg", fileNameMod = "", dim = 3, pow = 2 )
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum, seed; dim = dim );
	fileNameAttr = string( fileNameAttr, fileNameMod );
	degFileName = string( fileNameMain, "_", fileNameAttr );
	varFileName = string( degFileName, jldType );
	@info(varFileName);
	posLocLst = load(varFileName, "posLocLst");
	negLocLst = load(varFileName, "negLocLst");
	posNlst = load(varFileName, "posNlst");
	negNlst = load(varFileName, "negNlst");
	H_GUE_lst = load(varFileName, "H_GUE_lst");
	
	locLstPol = [posLocLst, negLocLst];
	locNPol = [posNlst, negNlst];
	locNCumPol = distillLocN( avgNum, N, locNPol );
	
	dim = 3;
	locLstWhichPol = [[ [ zeros(Bool,locNPol[id][iA,n]) for n = 2 : N-2 ] for iA = 1 : avgNum ] for id = 1:2];
	
	param_dim = dim;
	param_min = fill(0.0,param_dim);
	param_max = fill(2*pi,param_dim);
	param_step = ( param_max .- param_min ) ./ param_divide;
	param_grids = [ collect( range( param_min[iDim], param_max[iDim] - param_step[iDim], length = param_divide[iDim] ) ) for iDim = 1:param_dim ];
	posLst = CartesianIndices(Tuple(param_divide));
	param_mesh = [ [ param_grids[j][ind[j]] for j = 1:param_dim ] for ind in posLst ];
	
	for idPol = 1 : 2
		locLst = locLstPol[idPol];
		println(avgNum);
		for iA = 1 : avgNum
			for n = 2 : N-2
				# itCum = 1;
				for it = 1 : size( locLst[iA][n],1 )
					print( '\r', iA, " ", n, " ", it );
					loc = locLst[iA][n][it,:];
					param = param_mesh[loc...] + 1/2 * param_step;
					Htest = Hmat_3comb( param , H_GUE_lst[iA] );
					sol = eigen(Htest);
					vals = sol.values;
					isThis =  vals[n+1]-vals[n] < vals[n] - vals[n-1];
					locLstWhichPol[idPol][iA][n-1][it] = isThis;
					# print(it," ");
				end
			end
		end
	end
	
	locLstCumPol = [[ [ 
		if n == 1
			locLstPol[id][iA][1];
		elseif n == N-1
			locLstPol[id][iA][N];
		else 
			indices = findall( locLstWhichPol[id][iA][n-1] );
			locLstPol[id][iA][n][indices,:];
		end 
		for n = 1 : N-1 ] 
		for iA = 1 : avgNum ] 
	for id = 1:2 ];
	
	outFileName = string( "locDistilled", "_", fileNameAttr, jldType );
	save( outFileName, "locDistilledLst", locLstCumPol );
	
	return locLstCumPol;
end

function distillLocsDiff1( avgNum, N, param_divide, seed; fileNameMOd = "", dim = 3, pow = 2 )
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed; dim = dim );
	fileNameAttr = string( fileNameAttr, fileNameMod );
	degFileName = string( "deg", "_", fileNameAttr );
	varFileName = string( degFileName, jldType );
	@info(varFileName);
	posLocLst = load(varFileName, "posLocLst");
	negLocLst = load(varFileName, "negLocLst");
	posNlst = load(varFileName, "posNlst");
	negNlst = load(varFileName, "negNlst");
	
	
end

function distillLocsFromFile_old( avgNum, N, param_divide, seed; fileNameMod = "", dim = 3, pow = 2 )
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed; dim = dim );
	fileNameAttr = string( fileNameAttr, fileNameMod );
	degFileName = string( "deg", "_", fileNameAttr );
	varFileName = string( degFileName, jldType );
	@info(varFileName);
	posLocLst = load(varFileName, "posLocLst");
	negLocLst = load(varFileName, "negLocLst");
	posNlst = load(varFileName, "posNlst");
	negNlst = load(varFileName, "negNlst");
	H_GUE_lst = load(varFileName, "H_GUE_lst");
	
	locLstPol = [posLocLst, negLocLst];
	locNPol = [posNlst, negNlst];
	locNCumPol = distillLocN( avgNum, N, locNPol );
	
	dim = 3;
	locLstCumPol = [[ [ zeros(Int64, locNCumPol[1][iA,n],dim) for n = 1 : N-1 ] for iA = 1 : avgNum ] for id = 1:2	];
	
	param_dim = dim;
	param_min = fill(0.0,param_dim);
	param_max = fill(2*pi,param_dim);
	param_step = ( param_max .- param_min ) ./ param_divide;
	param_grids = [ collect( range( param_min[iDim], param_max[iDim] - param_step[iDim], length = param_divide[iDim] ) ) for iDim = 1:param_dim ];
	posLst = CartesianIndices(Tuple(param_divide));
	param_mesh = [ [ param_grids[j][ind[j]] for j = 1:param_dim ] for ind in posLst ];
	
	for idPol = 1 : 2
		locLst = locLstPol[idPol];
		for iA = 1 : avgNum
			for n = 2 : N
				itCum = 1;
				for it = 1 : size( locLst[iA][n],1 )
					loc = locLst[iA][n][it,:];
					param = param_mesh[loc...] + 1/2 * param_step;
					Htest = Hmat_3comb( param , H_GUE_lst[iA] );
					sol = eigen(Htest);
					vals = sol.values;
					if vals[n+1]-vals[n] < vals[n] - vals[n-1]
						locLstCumPol[idPol][iA][n][itCum,:] = loc;
						itCum += 1;
					end
				end
			end
		end
	end
	
	return locLstCumPol;
end

function whichLocsFromFile( avgNum, N, param_divide, seed; fileNameMod = "", dim = 3, pow = 2 )
	fileNameAttr = fileNameAttrFunc( N, param_divide, avgNum,  seed; dim = dim );
	fileNameAttr = string( fileNameAttr, fileNameMod );
	degFileName = string( "deg", "_", fileNameAttr );
	varFileName = string( degFileName, jldType );
	@info(varFileName);
	posLocLst = load(varFileName, "posLocLst");
	negLocLst = load(varFileName, "negLocLst");
	posNlst = load(varFileName, "posNlst");
	negNlst = load(varFileName, "negNlst");
	H_GUE_lst = load(varFileName, "H_GUE_lst");
	
	locLstPol = [posLocLst, negLocLst];
	locNPol = [posNlst, negNlst];
	locNCumPol = distillLocN( avgNum, N, locNPol );
	
	dim = 3;
	locLstWhichPol = [[ [ zeros(Bool,locNPol[id][iA,n]) for n = 2 : N-2 ] for iA = 1 : avgNum ] for id = 1:2];
	
	param_dim = dim;
	param_min = fill(0.0,param_dim);
	param_max = fill(2*pi,param_dim);
	param_step = ( param_max .- param_min ) ./ param_divide;
	param_grids = [ collect( range( param_min[iDim], param_max[iDim] - param_step[iDim], length = param_divide[iDim] ) ) for iDim = 1:param_dim ];
	posLst = CartesianIndices(Tuple(param_divide));
	param_mesh = [ [ param_grids[j][ind[j]] for j = 1:param_dim ] for ind in posLst ];
	
	for idPol = 1 : 2
		locLst = locLstPol[idPol];
		println(avgNum);
		for iA = 1 : avgNum
			for n = 2 : N-2
				# itCum = 1;
				for it = 1 : size( locLst[iA][n],1 )
					print( '\r', iA, " ", n, " ", it );
					loc = locLst[iA][n][it,:];
					param = param_mesh[loc...] + 1/2 * param_step;
					Htest = Hmat_3comb( param , H_GUE_lst[iA] );
					sol = eigen(Htest);
					vals = sol.values;
					isThis =  vals[n+1]-vals[n] < vals[n] - vals[n-1];
					locLstWhichPol[idPol][iA][n-1][it] = isThis;
					# print(it," ");
				end
			end
		end
	end
	
	return locLstWhichPol;
end

function distillLocN( avgNum, N, locNPol )
	locNCumPol = [zeros( Int64, avgNum, N-1 ) for ip = 1 : 2];
	for iPol = 1 : 2
		locNCumPol[iPol][:,1] = locNPol[iPol][:,1];
		for n = 2 : N-1
			locNCumPol[iPol][:,n] = locNPol[iPol][:,n] - locNCumPol[iPol][:,n-1];
		end
	end
	return locNCumPol;
end

function distillLocsFromWhich( locLstPol, locLstWhichPol, avgNum, N )
	locLstCumPol = [[ [ 
		if n == 1
			locLstPol[id][iA][1];
		elseif n == N-1
			locLstPol[id][iA][N];
		else 
			indices = findall( locLstWhichPol[id][iA][n-1] );
			locLstPol[id][iA][n][indices,:];
		end 
		for n = 1 : N-1 ] 
		for iA = 1 : avgNum ] 
	for id = 1:2 ];
	
	return locLstCumPol;
end

function loopsClassifyGOE( posLocLst, param_divide, N )
	posLocLstPackedArr = [ circshift( hcat( (dim-1)*ones(size(posLocLst[dim,id3][n],1)), posLocLst[dim,id3][n].- 0.5 ), [0,-(dim-1)] ) for dim = 1 : 3, id3 = 1:param_divide, n = 1 : N ];
	posLocLstPackedLst = reshape( posLocLstPackedArr, (dim*param_divide, N) );
	posLocLstPackedLstSorted = reshape( posLocLstPackedArr, (dim*param_divide, N) );
end

function locNvarFromFile( N, param_divide, itNum, seed; dim = nothing, alpha = nothing, fMod = "", isDistilled = true, fExt=jld2Type, opt = "distill2" )
	if isDistilled && opt == "distill2"
		fNameMain = fLocDistilled;
	elseif isDistilled && opt == "distill"
		fNameMain = fLocDistilledOld;
	else
		fNameMain = fDeg;
	end
	attrLst, valLst = fNameAttrLstFunc( N, param_divide, itNum, seed; dim = dim, alpha = alpha );
	# attrLst = ["N", "param_divide", "instanceNum", "seed"];
	# valLst = [N, param_divide, itNum, seed];
	# if !isnothing(dim)
		# insert!(attrLst,1,"dim");
		# insert!(valLst,1,dim);
	# end
	# if !isnothing(alpha)
		# append!(attrLst, ["alpha"] );
		# append!(valLst,alpha);
	# end
	fName = fNameFunc( fNameMain, attrLst, valLst, fExt; fMod = fMod );
	
	locNlst = load( fName, varNDistilled );
	
	posNlst = locNlst[1];
	
	gVarLst = std( posNlst; dims = 1 );
	
	fGvarNpy = fNameFunc( fGvar, attrLst, valLst, npyType; fMod = fMod );
	npzwrite( fGvarNpy, gVarLst );
end

function parityGOE_resave_fromFile( N, param_divide, itNum, seed; fMod = "", dim = nothing, fExt = jld2Type, alpha = 0 )
	divNum = param_divide[1];
	fMain = "deg_GOE_3d_full";
	attrLst, valLst = fAttrValFunc( N, param_divide, itNum, seed, dim; alpha = alpha );
	# attrLst = deepcopy(attrLstBase);
	# valLst = [N, param_divide, itNum, seed];
	# if !isnothing(dim)
		# attrLst = insert!( attrLst, 1, "dim" );
		# valLst = insert!( valLst, 1, dim );
	# end
	# @infiltrate
	fName = fNameFunc( fMain, attrLst, valLst, fExt; fMod );
	@info(fName);
	tFull = @timed begin
		parity2dLst = load( fName, "parity2dLst" ) ;
	end
	@info("read file: ");
	@info( timeMemStr( tFull.time, tFull.bytes ) );
	# @infiltrate
	tFull = @timed begin
		parity2dArr = [ real( parity2dLst[it][d][x,y][n] ) for it=1:itNum, d=1:dim, x=1:divNum, y=1:divNum, n=1:N ];
	end
	@info("lst to arr: ");
	@info( timeMemStr( tFull.time, tFull.bytes ) )
	oFmain = fNameZacArr;
	oFname = fNameFunc(oFmain, attrLst, valLst, jld2Type; fMod);
	save( oFname, "parity2dArr", parity2dArr );
end

function parityAvg_fromFile( N, param_divide, itNum, seed; dim = nothing, fMod = "", fExt = jld2Type, alpha = 0 )
	attrLst, valLst = fAttrValFunc( N, param_divide, itNum, seed, dim; alpha = alpha );
	fNameMain = fNameFunc( fNameZacArr, attrLst, valLst, fExt; fMod = fMod );
	parityArr = load( fNameMain, varNameParityArr );
	sumDimLst = (3,4);
	zacAvgArr = dropdims( mean( parityArr; dims = sumDimLst ); dims = sumDimLst );
	oFmain = fNameZacAvg;
	oFname = fNameFunc( oFmain, attrLst, valLst, jld2Type; fMod = fMod );
	save( oFname, varNameZacAvg, zacAvgArr );
	# @infiltrate
end

function collisionPurifyFromFile( mSz, param_divide, itNum, seed; dim = nothing, fMod = "", fExt = jld2Type, alpha = 0, divColl = 10, paramMax = 2*pi, bndThrRat = 0.1, fidlThr = 0.99999 )
	attrLst, valLst = fAttrValFunc( mSz, param_divide, itNum, seed, dim; alpha = alpha );
	
	fMain = fDeg;
	fExt = jld2Type;
	
	fVarName = fNameFunc( fMain, attrLst, valLst, fExt; fMod = fMod );
	
	posLocLst = load( fVarName, varNamePosLoc );
	negLocLst = load( fVarName, varNameNegLoc );
	posNlst = load( fVarName, varNamePosN );
	negNlst = load( fVarName, varNameNegN );
	H_GUE_lst = load( fVarName, "H_GUE_lst" );
	
	# @infiltrate
	nLev = size(posNlst,2);
	if isnothing(dim)
		dim = 3;
	end
	ratio = ones(dim);
	if !isa( divColl, Array )
		divLstHash = fill(divColl,dim);
	else
		divLstHash = divColl;
	end
	if !isa( paramMax, Array )
		paramMaxLstHash = fill(paramMax, dim);
	else
		paramMaxLstHash = paramMax;
	end
	stepLstHash = paramMaxLstHash ./ divLstHash;
	
	idHash = Array{Vector{Int64},dim}(undef,divLstHash...);
	for ii in eachindex(idHash)
		idHash[ii] = zeros(Int64,0);
	end
	
	locLstPol = [posLocLst, negLocLst];
	NLstPol = [posNlst,negNlst];
	NLstPolOrg = deepcopy(NLstPol);
	
	locLstPolOut = [ [ Vector{Array{Float64,2}}(undef,nLev) for it = 1:itNum ] for iPol = 1:2 ];
	vLstPolOut = [ [ Vector{Array{ComplexF64,3}}(undef,nLev) for it = 1:itNum ] for iPol = 1:2 ];
	fidlLstPol = [ [Vector{Vector{Float64}}(undef,nLev) for it = 1:itNum] for iPol = 1:2 ];
	
	HmatThr = threaded_zeros( ComplexF64, mSz, mSz );
	vLstThr = threaded_zeros( ComplexF64, mSz, mSz );
	eLstThr = threaded_zeros( Float64, mSz );
	eigWorks = eigenPreworkStruct!( getThrInst(HmatThr), getThrInst(eLstThr), getThrInst(vLstThr) );
	eigWorksThr = thrStructCpyTheRest( eigWorks );
	for iPol = 1:2
		locLstLst = locLstPol[iPol];
		NLstLst = NLstPol[iPol];
		for it = 1 : itNum, m = 1 : nLev
			HmatFun = ( xLst -> Hmat_3comb_ratio!( getThrInst( HmatThr ), xLst, H_GUE_lst[it], ratio ) );
			vLstTmp = zeros( ComplexF64, mSz, 2, NLstLst[it,m] );
			Threads.@threads for iLoc = 1 : NLstLst[it,m]
				@views HmatFun( locLstLst[it][m][iLoc,:] );
				eigenZheevrStruct!( getThrInst(HmatThr), getThrInst(eLstThr), getThrInst(vLstThr), getThrInst(eigWorksThr) );
				@views vLstTmp[:,:,iLoc] = vLstThr[:,m:m+1];
			end
			locLstPolOut[iPol][it][m], vLstPolOut[iPol][it][m], fidlLstPol[iPol][it][m] = collisionPurifyPerLevel!( idHash, locLstLst[it][m], vLstTmp, paramMaxLstHash, divLstHash, stepLstHash; bndThrRat = bndThrRat, fidlThr = fidlThr, isRecordFidl = true );
			NLstLst[it,m] = size( locLstPolOut[iPol][it][m], 1 )
		end
		# vLstLst = [ [ zeros(ComplexF64,N,2,NLstLst[it,m]) for m = 1 : length[locLstLst[it]] ] for it = 1 : itNum ];
	end
	fidlLstAllPol = [ vcat.(fidlLstPol[iPol]...) for iPol = 1:2 ];
	
	oFmain = fMain * "_collisionPurified";
	attrLstOut = copy(attrLst);
	valLstOut = copy(valLst);
	push!(attrLstOut, "divColl", "bndThrRat");
	push!(valLstOut, divColl, bndThrRat);
	oFname = fNameFunc( oFmain, attrLstOut, valLstOut, fExt; fMod = fMod );
	save( oFname, varNamePosLoc, locLstPolOut[1], varNameNegLoc, locLstPolOut[2], varNamePosN, NLstPol[1], varNameNegN, NLstPol[2], "fidlLstPol", fidlLstPol, "fidlLstAllPol", fidlLstAllPol );
	# @infiltrate
	return locLstPolOut, NLstPol, fidlLstPol;
end

function collisionPurify!( idHash, locLstLst, vLstLst, mSz, paramMaxLst, divLst, stepLst; bndThrRat = 0.01, fidlThr = 0.99999 )
	locLstLstPure = Vector{Matrix{Float64}}(undef,mSz-1);
	vLstLstPure = Vector{Array{ComplexF64,3}}(undef,mSz-1);
	
	for m = 1 : mSz-1
		locLstLstPure[m], vLstLstPure[m] = collisionPUrifyPerLevel!( idHash, locLstLst[m], vLstLst[m], paramMaxLst, divLst, stepLst; bndThrRat = bndThrRat, fidlThr = fidlThr );
	end
	return locLstLstPure, vLstLstPure;
end

function collisionPurifyPerLevel!( idHash, locLst, vLst, paramMaxLst, divLst, stepLst; bndThrRat = 0.01, fidlThr = 0.99999, isRecordFidl = false )
	for idLst in idHash
		empty!(idLst);
	end
	
	if isempty(locLst)
		return
	end
	
	nDim = length(divLst);
	bndThrAbsMin = bndThrRat * stepLst;
	bndThrAbsMax = (1 - bndThrRat) * stepLst;
	
	# if isRecordFidl
	fidlLst = zeros(0);
	# end
	
	idTmp = zeros(Int64,nDim);
	idTmpSh = similar(idTmp);
	idBndSh = zeros(Int64,nDim);
	idRemain = zeros(nDim);
	idShLst = zeros(Int64,nDim,2^nDim);
	iShThis = 1;
	# iShNext = 1;
	isBndLst = fill(false,nDim);
	idLin = 1;
	
	idNonColl = zeros(Int64,0);
	for ii = 1 : size(locLst,1)
		# @infiltrate
		idTmp .= floor.( mod.( @view(locLst[ii,:]), paramMaxLst ) ./ stepLst ) .+ 1;
		idRemain .= @view(locLst[ii,:]) .- (idTmp .- 1) .* stepLst;
		idLin = getLinInd( idTmp, divLst, nDim );
		
		isCollide = false;
		for iHsh = 1 : length( idHash[idLin] )
			# @infiltrate
			fidl = fidelity2by2( vLst[:,:,ii], vLst[:,:,idHash[idLin][iHsh]] );
			if isRecordFidl
				push!(fidlLst,fidl);
			end
			@views if fidl > fidlThr
				isCollide = true;
				break;
			end
		end
		if isCollide
			continue;
		end
		push!(idNonColl,ii);
		
		idBndSh .= 0;
		for iDim = 1 : nDim
			if idRemain[iDim] < bndThrAbsMin[iDim]
				idBndSh[iDim] = -1;
			elseif idRemain[iDim] > bndThrAbsMax[iDim]
				idBndSh[iDim] = 1;
			end
		end
		# isBndLst .= idRemain .> bndThrAbs;
		@views idShLst[:,1] .= 0;
		iShThis = 1;
		for iDim = 1:nDim
			if idBndSh[iDim] != 0
				for iShBnd = 1 : iShThis
					@views idShLst[:,iShBnd+iShThis] .= idShLst[:,iShBnd];
					idShLst[iDim,iShBnd+iShThis] = idBndSh[iDim];
				end
				iShThis *= 2;
			end
		end
		
		for iSh = 1 : iShThis
			idTmpSh .= mod.( @view(idShLst[:,iSh]) .+ idTmp .- 1, divLst ) .+ 1;
			idLin = getLinInd( idTmpSh, divLst, nDim );
			push!(idHash[idLin],ii);
		end
	end
	
	# @infiltrate
	return locLst[idNonColl,:], vLst[:,:,idNonColl], fidlLst;
end

function fidelity2by2( v1lst, v2lst; numVec = 2 )
	H22 = zeros(ComplexF64,numVec,numVec);
	for i1 = 1:numVec
		for i2 = 1:numVec
			for jj = 1:numVec
				@views H22[i1,i2] += 1/4 * dot( v1lst[:,i1], v2lst[:,jj] ) * dot( v2lst[:,jj], v1lst[:,i2] );
			end
		end
	end
	return real( tr(sqrt(H22)) );
end
