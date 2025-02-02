using HDF5
using Plots

function load_dump_HDF(num::Int64)
    filename = joinpath("dumps", "dump$num.h5")
    
    # Otwarcie pliku HDF5 w trybie odczytu
    h5open(filename, "r") do file
        # Wczytanie danych z pliku (w tym przypadku dane są w zestawie 'grid_data')
        dataset = file["grid_data"]
        grid_data = read(dataset)
        return grid_data
    end
end

	
for i in 1:32
	grid_data_loaded = load_dump_HDF(i)
	println("Wczytano dane z dump_$i.h5")

	# Przykładowe rozmiary siatki
	N1 = size(grid_data_loaded, 1)  # Rozmiar w osi x
	N2 = size(grid_data_loaded, 2)  # Rozmiar w osi y
	N3 = size(grid_data_loaded, 3)  # Rozmiar w osi z

	# Zakładając, że mamy jedną jednostkę długości w każdej osi
	N1_grid = 1:N1  # Przykładowa siatka x
	N3_grid = 1:N3  # Przykładowa siatka z

	mid_y = div(N2, 2)  # Przekrój w połowie osi y
	density_slice = grid_data_loaded[:, mid_y, :, 1]  # Wybór przekroju wzdłuż osi x-z

	# Tworzenie wykresu
	p = heatmap(N1_grid, N3_grid, density_slice',
	    xlabel = "x (w jednostkach promienia grawitacyjnego)",
	    ylabel = "z (w jednostkach promienia grawitacyjnego)",
	    title = "Rozkład gęstości (przekrój w płaszczyźnie x-z)",
	    colorbar_title = "Gęstość")

	# Zapisanie wykresu do pliku PNG
	savefig(p, "density_plot_$i.png")
end
