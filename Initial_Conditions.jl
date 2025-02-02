using BenchmarkTools
using Base.Threads
using StaticArrays

#Size of grid
N1 = 50
N2 = 50
N3 = 50

#Grid Limits (in gravitational Radius)
N1_min = 0
N1_max = 10

N2_min = 0
N2_max = 10

N3_min = 0
N3_max = 10

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
            grid[i, j, k, 1] =  2+ 0.1*sin(i+j+k)    # Gęstość (rho)
            grid[i, j, k, 2] =  2+ 0.2*cos(i+j+k)    # Energia wewnętrzna (u)
            grid[i, j, k, 3] = 0.3- 0.1*sin(i+j+k)   # Prędkość U1
            grid[i, j, k, 4] = 0.3-0.1*sin(i+j+k)    # Prędkość U2
            grid[i, j, k, 5] = 0.3- 0.1*sin(i+j+k)   # Prędkość U3
            grid[i, j, k, 6] = 0.3+ 0.01*sin(i+j+k)  # Pole magnetyczne B1
            grid[i, j, k, 7] = 0.3+ 0.01*sin(i+j+k)  # Pole magnetyczne B2
            grid[i, j, k, 8] = 0.3+ 0.01*sin(i+j+k)  # Pole magnetyczne B3
        end
    end
end



#Iteracja i=62, j=30, k=8
#[2.144299779477927, 2.389710151903291, 0.061200220522073474, 0.061200220522073474, 0.061200220522073474, 0.005429977947792654, 0.005429977947792654, 0.005429977947792654]
#Zbiegło!!! w 87 iteracji

