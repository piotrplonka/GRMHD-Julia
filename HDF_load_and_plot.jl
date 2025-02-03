using HDF5
using Plots
using Interpolations

# Ustawienie wysokiej rozdzielczości dla backendu (opcjonalnie)
default(dpi=300)

function load_dump_HDF(num::Int64)
    filename = joinpath("dumps", "dump$(num).h5")
    # Otwarcie pliku HDF5 w trybie odczytu
    h5open(filename, "r") do file
        # Wczytanie danych z pliku (zakładamy, że dane są w zestawie 'grid_data')
        dataset = file["grid_data"]
        grid_data = read(dataset)
        return grid_data
    end
end

# Funkcja do interpolacji danych – zwiększenie rozdzielczości o zadany współczynnik (factor)
function interpolate_slice(data_slice, factor::Int=4)
    # Pobranie rozmiarów oryginalnej macierzy
    nx, nz = size(data_slice)
    # Oryginalne osie
    x = 1:nx
    z = 1:nz
    # Utworzenie interpolatora liniowego
    itp = interpolate(data_slice, BSpline(Linear()))
    # Nowe osie (zwiększona liczba punktów)
    x_new = range(first(x), last(x), length=nx*factor)
    z_new = range(first(z), last(z), length=nz*factor)
    # Interpolacja – iterujemy po nowej siatce (uwaga: w pętli najpierw iterujemy po z_new)
    data_interp = [itp(x_val, z_val) for z_val in z_new, x_val in x_new]
    return x_new, z_new, data_interp
end

for i in 1:100
    grid_data_loaded = load_dump_HDF(i)
    println("Wczytano dane z dump_$i.h5")

    # Przyjmujemy, że w grid_data_loaded wymiary są:
    #   N1 – oś x, N2 – oś y, N3 – oś z, a ostatni wymiar to kanały danych.
    N1 = size(grid_data_loaded, 1)  # liczba punktów w osi x
    N2 = size(grid_data_loaded, 2)  # liczba punktów w osi y
    N3 = size(grid_data_loaded, 3)  # liczba punktów w osi z

    mid_y = div(N2, 2)  # Przekrój w połowie osi y

    # Wyodrębniamy poszczególne kanały:
    # Kanał 1: gęstość (rho)
    density_slice = grid_data_loaded[:, mid_y, :, 1]
    # Kanał 2: energia wewnętrzna
    energy_slice  = grid_data_loaded[:, mid_y, :, 2]
    # Kanały 3,4,5: składowe prędkości (ux, uy, uz)
    ux_slice = grid_data_loaded[:, mid_y, :, 3]
    uy_slice = grid_data_loaded[:, mid_y, :, 4]
    uz_slice = grid_data_loaded[:, mid_y, :, 5]
    # Obliczamy prędkość całkowitą
    velocity_slice = sqrt.(ux_slice.^2 .+ uy_slice.^2 .+ uz_slice.^2)

    # Interpolacja – zwiększamy rozdzielczość (tutaj przykład dla współczynnika 10)
    factor = 10
    x_new, z_new, density_interp = interpolate_slice(density_slice, factor)
    _,      _, energy_interp   = interpolate_slice(energy_slice, factor)
    _,      _, velocity_interp = interpolate_slice(velocity_slice, factor)

    # Aby uniknąć problemów z logarytmowaniem zer (lub bardzo małych wartości)
    eps = 1e-12
    density_log = log10.(density_interp .+ eps)
    energy_log  = log10.(energy_interp  .+ eps)
    # Dla prędkości możesz pozostawić skalę liniową – lub też przeskalować, jeśli potrzeba

    # Tworzymy wykresy z układem 1x2 dla density i energy
    p1 = heatmap(x_new, z_new, density_log,
        xlabel = "x (jednostki promienia grawitacyjnego)",
        ylabel = "z (jednostki promienia grawitacyjnego)",
        title = "Gęstość (log₁₀ skala)",
        colorbar_title = "log₁₀(ρ)",
        aspect_ratio = :equal,
        c = :viridis,
        framestyle = :box,
        grid = false)
        
    p2 = heatmap(x_new, z_new, energy_log,
        xlabel = "x (jednostki promienia grawitacyjnego)",
        ylabel = "z (jednostki promienia grawitacyjnego)",
        title = "Energia wewnętrzna (log₁₀ skala)",
        colorbar_title = "log₁₀(E)",
        aspect_ratio = :equal,
        c = :inferno,
        framestyle = :box,
        grid = false)

    # Łączymy oba wykresy w jeden rysunek (layout 1 x 2)
    p_combined = plot(p1, p2, layout = (1, 2), size=(1200,500))
    # Zapisujemy rysunek do pliku PNG (wyższa rozdzielczość dzięki default dpi)
    savefig(p_combined, "density_energy_plot_$i.png")

    # Dodatkowy wykres dla prędkości całkowitej (skala liniowa)
    p_vel = heatmap(x_new, z_new, velocity_interp,
        xlabel = "x (jednostki promienia grawitacyjnego)",
        ylabel = "z (jednostki promienia grawitacyjnego)",
        title = "Prędkość całkowita",
        colorbar_title = "|v|",
        aspect_ratio = :equal,
        c = :plasma,
        framestyle = :box,
        grid = false,
        size=(600,500))
    savefig(p_vel, "velocity_plot_$i.png")
end

