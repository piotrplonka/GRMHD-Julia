using BenchmarkTools
using Base.Threads
using Profile
using LinearAlgebra  
using HDF5

function save_dump_HDF(num::Int64, grid::AbstractArray)
    folder = "dumps"           
    mkpath(folder)              
    filename = joinpath(folder, "dump$num.h5")  
    h5open(filename, "w") do file
        dataset = create_dataset(file, "grid_data", Float64, size(grid))
        write(dataset, grid)
    end
    println("Dane zostały zapisane do pliku: $filename")
end


