function sim!(P, C, dt)
    for p in P
        integrate!(p, p.param, Float32(dt))
        record!(p)
    end
    for c in C
        forward!(c, c.param)
        record!(c)
    end
end

function sim!(P, C; dt = 0.1ms, duration = 10ms)
<<<<<<< HEAD
	sized = duration/dt
	for p in P
		if hasproperty(p, :spike_raster)
			p.spike_raster::Vector{Int32} = zeros(trunc(Int, sized))
		end
		for t = 0ms:dt:duration#(duration - dt)
			integrate!(p, p.param, Float32(dt))
			record!(p)
		end
	end
=======
    for t = 0ms:dt:(duration - dt)
        sim!(P, C, dt)
    end
>>>>>>> 8135f6c... remove trailing sim!/train!
end

	#println(size(spike_raster))
	#println(P.sized)

    #if hasproperty(P, :spike_raster)
		#P.spike_raster::Vector{SNNInt} = zeros(trunc(Int, size))
		#spike_raster::VFT = fill(0, N)

		#P.spike_raster::Vector{Int32} = zeros(trunc(Int, size))
		#println(P.spike_raster)
    #end
    #for t = 0ms:dt:(duration - dt)
    #    sim!(P, C, dt)
    #end

function train!(P, C, dt, t = 0)
    for p in P
        integrate!(p, p.param, Float32(dt))
        record!(p)
    end
    for c in C
        forward!(c, c.param)
        plasticity!(c, c.param, Float32(dt), Float32(t))
        record!(c)
    end
end

function train!(P, C; dt = 0.1ms, duration = 10ms)
    for t = 0ms:dt:(duration - dt)
        train!(P, C, Float32(dt), Float32(t))
    end
end
