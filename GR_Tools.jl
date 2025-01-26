using BenchmarkTools
using Base.Threads
using StaticArrays

include("Matrix_operations.jl")
include("Christoffel_Symbols.jl")
include("eos.jl")

function Metric_Calculation(r::AbstractVector, θ::AbstractVector, φ::AbstractVector, a::Float64)
    
    # Dimensions of metric
    Nx = length(r)
    Ny = length(θ)
    Nz = length(φ)
    
    # Zero values for metric
    gcov = zeros(Nx, Ny, Nz, 4, 4)
    println("Shape of gcov: ", size(gcov))
    
    # Calculations
    for i in 1:Nx
        for j in 1:Ny
         @threads   for k in 1:Nz
                gcov[i, j, k, :, :] = Kerr_Schild_metric_cov(r[i], θ[j], φ[k], a)
            end
        end
    end
    
    return gcov
end

function Christoffel_Symbols_Calculation(r::AbstractVector, θ::AbstractVector, φ::AbstractVector, a::Float64)
    
    # Dimensions of metric
    Nx = length(r)
    Ny = length(θ)
    Nz = length(φ)
    
    # Zero values for Christoffel Symbols
    CSymbols = zeros(Nx, Ny, Nz, 4, 4, 4)
    println("Shape of Christoffel Symbols: ", size(CSymbols))
    
    # Calculations
    for i in 1:Nx
        for j in 1:Ny
            for k in 1:Nz
            	for K1 in 1:4
            		for K2 in 1:4
            			@threads for K3 in 1:4
                			CSymbols[i, j, k, K1, K2, K3] = Kerr_Schild_Christoffel_Symbols(K1,K2,K3,r[i], θ[j], φ[k], a)
                		end
                	end
                end
            end
        end
    end
    
    return CSymbols
end
