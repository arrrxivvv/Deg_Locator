function deltaN( N, param_divide, posLocLst, negLocLst )
	pow = 2;
	return deltaNpower( N, param_divide, posLocLst, negLocLst, pow );
end

function deltaNpower( N, param_divide, posLocLst, negLocLst, pow )
	delNlst = [zeros( N, param_divide ) for idPol=1:2];
	
	locLst = [posLocLst, negLocLst];
	
	for n = 1:N
		for idPol = 1:2
			locX = sort( locLst[idPol][n][:,1] );
			lastPos = 1;
			numDeg = 0;
			for loc in locX
				delNlst[idPol][n,lastPos:loc-1] .= numDeg;
				numDeg += 1;
				delNlst[idPol][n,loc] = numDeg;
				lastPos = loc+1;
			end
			delNlst[idPol][n,lastPos:param_divide] .= numDeg;
		end
	end
	delNlstTotal = delNlst[1]-delNlst[2];
	
	return delNlstTotal.^pow;
end

# function deltaN( N, param_divide, posLocLst, negLocLst )
	# delNlst = [zeros( N, param_divide ) for idPol=1:2];
	
	# locLst = [posLocLst, negLocLst];
	
	# for n = 1:N
		# for idPol = 1:2
			# locX = sort( locLst[idPol][n][:,1] );
			# lastPos = 1;
			# numDeg = 0;
			# for loc in locX
				# delNlst[idPol][n,lastPos:loc-1] .= numDeg;
				# numDeg += 1;
				# delNlst[idPol][n,loc] = numDeg;
				# lastPos = loc+1;
			# end
			# delNlst[idPol][n,lastPos:param_divide] .= numDeg;
		# end
	# end
	# delNlstTotal = delNlst[1]-delNlst[2];
	
	# return delNlstTotal.^2;
# end
