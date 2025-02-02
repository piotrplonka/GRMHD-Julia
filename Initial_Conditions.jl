using Plots
gr()  # Ustawienie backendu GR
using BenchmarkTools
using Base.Threads
using StaticArrays

# Size of grid
N1 = 100  # Oś x
N2 = 100  # Oś y
N3 = 100  # Oś z

# Grid Limits (in gravitational Radius)
N1_min = -10.0
N1_max =  10.0

N2_min = -10.0
N2_max =  10.0

N3_min = -10.0
N3_max =  10.0

# Number of primitive variables
NP = 8

# Grid Initialization
grid = zeros(N1, N2, N3, NP)  # Zmieniono nazwę 'grid' na 'grid'

# Coordinates (kartezjański układ)
N1_grid = collect(range(N1_min, N1_max, length=N1))
N2_grid = collect(range(N2_min, N2_max, length=N2))
N3_grid = collect(range(N3_min, N3_max, length=N3))

# Ustalmy promień gwiazdy i środek
star_radius = 6.0
star_center = (0.0, 0.0, 0.0)

# Primitive variables (density, internal energy, velocities, magnetic fields)
# Oznaczenia:
#   1: ρ (density)
#   2: u (internal energy)
#   3: u¹ (velocity in x-direction)
#   4: u² (velocity in y-direction)
#   5: u³ (velocity in z-direction)
#   6: B¹ (magnetic field in x-direction)
#   7: B² (magnetic field in y-direction)
#   8: B³ (magnetic field in z-direction)
for i in 1:N1
    for j in 1:N2
        for k in 1:N3
            # Pobieramy kartezjańskie współrzędne
            x = N1_grid[i]
            y = N2_grid[j]
            z = N3_grid[k]
            
            # Obliczamy odległość od środka gwiazdy
            r = sqrt((x - star_center[1])^2 + (y - star_center[2])^2 + (z - star_center[3])^2)
            
            # Ustalanie gęstości i energii wewnętrznej:
            if r <= star_radius
                # Wewnątrz gwiazdy – rozkład sferyczny (np. eksponencjalny spadek gęstości)
                density = 2.0 * exp(-r/2)
                internal_energy = 2.0 + 0.2 * cos(r)  # można dodać perturbacje zależne od r
            else
                # Poza gwiazdą – niska gęstość otoczenia (atmosfera)
                density = 1.0e-1
                internal_energy = 1.0e-1
            end
            grid[i, j, k, 1] = density
            grid[i, j, k, 2] = internal_energy

            # Ustalanie prędkości: chcemy, by wewnątrz gwiazdy prędkość była skierowana do środka (zapadanie)
            if r > 0 && r <= star_radius
                v0 = 0.2  # amplituda prędkości zapadania
                # Składowe prędkości skierowane do środka
                vx = -v0 * (x / r) + 0.01 * sin(x + y + z)
                vy = -v0 * (y / r) + 0.01 * cos(x + y + z)
                vz = -v0 * (z / r) + 0.01 * sin(x - y + z)
            else
                # Poza gwiazdą przyjmujemy zerową prędkość
                vx = 0.0
                vy = 0.0
                vz = 0.0
            end
            grid[i, j, k, 3] = vx
            grid[i, j, k, 4] = vy
            grid[i, j, k, 5] = vz

            # Ustalanie pola magnetycznego: przyjmujemy stałe pole pionowe (wzdłuż osi z) z małymi perturbacjami
            Bx = 0.0 + 0.001 * sin(x + y + z)
            By = 0.0 + 0.001 * cos(x + y + z)
            Bz = 0.3 + 0.001 * sin(x - y + z)
            grid[i, j, k, 6] = Bx
            grid[i, j, k, 7] = By
            grid[i, j, k, 8] = Bz
        end
    end
end

using Plots
gr()  # Ustawienie backendu GR

# Teraz rysujemy wykres w płaszczyźnie x-z dla gęstości:
mid_y = div(N2, 2)
density_slice = grid[:, mid_y, :, 1]  # Zmienna 1 to gęstość

# Rysowanie mapy ciepła:
p = heatmap(N1_grid, N3_grid, density_slice',
    xlabel = "x (w jednostkach promienia grawitacyjnego)",
    ylabel = "z (w jednostkach promienia grawitacyjnego)",
    title = "Rozkład gęstości (przekrój w płaszczyźnie x-z)",
    colorbar_title = "Gęstość")

# Zapisz wykres do pliku (np. w formacie PNG):
savefig(p, "density_plot.png")  # Zapisuje wykres w bieżącym katalogu
