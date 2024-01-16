function timeMemStr( t, mem )
	"   time: " * string(t) * " seconds, memory: " * Base.format_bytes(mem)
end

macro timeInfo( ex )
	exTimed = 
		quote
			tFull = @timed $(esc(ex));
			@info( timeMemStr( tFull.time, tFull.bytes ) );
		end
	return ( exTimed );
end
