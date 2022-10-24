indlst = zeros(Int64,3);
for indlst[1] = 1:10, indlst[2] = 1:100, indlst[3] = 1:100
	x = indlst[1];
	y = indlst[2];
	print("x, y, z = ", string(x), " ", string(y), " ", string(indlst[3]), "\r");
end	