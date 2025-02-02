using Plots
gr()  
using BenchmarkTools
using Base.Threads
using StaticArrays


N1 = 100 
N2 = 100  
N3 = 100  

N1_min = -10.0
N1_max =  10.0

N2_min = -10.0
N2_max =  10.0

N3_min = -10.0
N3_max =  10.0


NP = 8


grid = zeros(N1, N2, N3, NP)  


N1_grid = collect(range(N1_min, N1_max, length=N1))
N2_grid = collect(range(N2_min, N2_max, length=N2))
N3_grid = collect(range(N3_min, N3_max, length=N3))


star_radius = 6.0
star_center = (0.0, 0.0, 0.0)


for i in 1:N1
    for j in 1:N2
        for k in 1:N3
            x = N1_grid[i]
            y = N2_grid[j]
            z = N3_grid[k]
            
            r = sqrt((x - star_center[1])^2 + (y - star_center[2])^2 + (z - star_center[3])^2)
            
            if r <= star_radius
                density = 2.0 * exp(-r/2)
                internal_energy = 2.0 
            else
                density = 1.0e-1
                internal_energy = 1.0e-1
            end
            grid[i, j, k, 1] = density
            grid[i, j, k, 2] = internal_energy

            if r > 0 && r <= star_radius
                v0 = 0.1  
                
                vx = -v0 * (x / r) 
                vy = -v0 * (y / r) 
                vz = -v0 * (z / r) 
            else
                vx = 0.0
                vy = 0.0
                vz = 0.0
            end
            grid[i, j, k, 3] = vx
            grid[i, j, k, 4] = vy
            grid[i, j, k, 5] = vz

            Bx = 0.0 
            By = 0.0 
            Bz = 0.0
            grid[i, j, k, 6] = Bx
            grid[i, j, k, 7] = By
            grid[i, j, k, 8] = Bz
        end
    end
end

using Plots
gr()  

mid_y = div(N2, 2)
density_slice = grid[:, mid_y, :, 1]  # Zmienna 1 to gęstość


p = heatmap(N1_grid, N3_grid, density_slice',
    xlabel = "x (w jednostkach promienia grawitacyjnego)",
    ylabel = "z (w jednostkach promienia grawitacyjnego)",
    title = "Rozkład gęstości (przekrój w płaszczyźnie x-z)",
    colorbar_title = "Gęstość")
savefig(p, "density_plot.png") 
