
Threads.@threads for it = 1 : 1000
	for x = 1 : 1000
		if mod(x,100) == 0 && mod(it,10) == 0
			print("it = ", string(it), ", x = ", string(x), "\r");
		end
	end
end
