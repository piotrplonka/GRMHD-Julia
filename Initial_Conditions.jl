using BenchmarkTools
using Base.Threads
using StaticArrays

#Size of grid
N1 = 128*2
N2 = 64*2
N3 = 32*2

#Grid Limits (in gravitational Radius)
N1_min = 10
N1_max = 100

N2_min = 0.1
N2_max = 2

N3_min = 0
N3_max = 2*3.1415

#Number of primitive variables
NP = 8

# Grid Initialization
grid = zeros(N1,N2,N3,NP)

#Coordinates
N1_grid = collect(range(N1_min, N1_max, length=N1))
N2_grid = collect(range(N2_min, N2_max, length=N2))
N3_grid = collect(range(N3_min, N3_max, length=N3))

# Primitive variables (density, internal energy, velocity, magnetic fields)
@threads for i in 1:N1
    for j in 1:N2
        for k in 1:N3
            # Wzrost gęstości i energii w przestrzeni
            grid[i, j, k, 1] = 0.5 + 0.5 * (i / N1) + 0.2 * rand()  # Gęstość (rho)
            grid[i, j, k, 2] = 0.5 + 0.5 * (j / N2) + 0.2 * rand()  # Energia wewnętrzna (u)
            grid[i, j, k, 3] = 0.05 * rand()  # Prędkość U1
            grid[i, j, k, 4] = 0.05 * rand()  # Prędkość U2
            grid[i, j, k, 5] = 0.05 * rand()  # Prędkość U3
            grid[i, j, k, 6] = 0.5 * rand()  # Pole magnetyczne B1
            grid[i, j, k, 7] = 0.5 * rand()  # Pole magnetyczne B2
            grid[i, j, k, 8] = 0.5 * rand()  # Pole magnetyczne B3
        end
    end
end

for i in 1:N1
	for j in 1:N2
		for k in 1:N3
			#println("x: ", grid[i, j, k, :])
		end
	end
end
