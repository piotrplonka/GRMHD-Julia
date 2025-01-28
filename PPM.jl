using BenchmarkTools
using Base.Threads
using Profile
using LinearAlgebra  

include("Matrix_operations.jl")
include("Christoffel_Symbols.jl")
include("structs.jl")
include("Initial_Conditions.jl")
include("fluxlimiter.jl")
include("eos.jl")


function PPM_N1(Grid_Of_Values::Array{Float64, 4}, i::Int64, j::Int64, k::Int64)
	
	#PPM
	κ::Float64 = 1/3
	
	#δu calculations 
	δu_plus_1_2  ::Vector = Grid_Of_Values[i+1, j, k, :] .- Grid_Of_Values[i,   j, k, :]
	δu_sub_1_2   ::Vector = Grid_Of_Values[i,   j, k, :] .- Grid_Of_Values[i-1, j, k, :]
	δu_plus_3_2  ::Vector = Grid_Of_Values[i+2, j, k, :] .- Grid_Of_Values[i+1, j, k, :]
	δu_sub_3_2   ::Vector = Grid_Of_Values[i-1, j, k, :] .- Grid_Of_Values[i-2, j, k, :]                

	#r calculations 
	r            ::Vector = (Grid_Of_Values[i,   j, k, :] .- Grid_Of_Values[i-1, j, k, :]) ./ (Grid_Of_Values[i+1, j, k, :] .- Grid_Of_Values[i,   j, k, :])
	r_plus_1     ::Vector = (Grid_Of_Values[i+1, j, k, :] .- Grid_Of_Values[i,   j, k, :]) ./ (Grid_Of_Values[i+2, j, k, :] .- Grid_Of_Values[i+1, j, k, :])
	r_sub_1      ::Vector = (Grid_Of_Values[i-1, j, k, :] .- Grid_Of_Values[i-2, j, k, :]) ./ (Grid_Of_Values[i,   j, k, :] .- Grid_Of_Values[i-1, j, k, :])

	#Final Interpolated Values
	u_plus_1_2_L ::Vector = Grid_Of_Values[i,   j, k, :] .+ minmod(r       )/4 .* ((1-κ)*  δu_sub_1_2  .+ (1+κ) * δu_plus_1_2)
	u_plus_1_2_R ::Vector = Grid_Of_Values[i+1, j, k, :] .- minmod(r_plus_1 )/4 .* ((1-κ)* δu_plus_3_2 .+ (1+κ) * δu_plus_1_2)                
	u_sub_1_2_L  ::Vector = Grid_Of_Values[i-1, j, k, :] .+ minmod(r_sub_1  )/4 .* ((1-κ)* δu_sub_3_2  .+ (1+κ) * δu_sub_1_2 )
	u_sub_1_2_R  ::Vector = Grid_Of_Values[i,   j, k, :] .- minmod(r        )/4 .* ((1-κ)* δu_plus_1_2 .+ (1+κ) * δu_sub_1_2 )
	

return u_plus_1_2_L, u_plus_1_2_R, u_sub_1_2_L, u_sub_1_2_R

end



function PPM_N2(Grid_Of_Values::Array{Float64, 4}, i::Int64, j::Int64, k::Int64)
	
	#PPM
	κ::Float64 = 1/3
	
        δu_plus_1_2  ::Vector  = Grid_Of_Values[i, j+1, k, :] .- Grid_Of_Values[i, j,   k, :]
        δu_sub_1_2   ::Vector  = Grid_Of_Values[i, j,   k, :] .- Grid_Of_Values[i, j-1, k, :]
        δu_plus_3_2  ::Vector  = Grid_Of_Values[i, j+2, k, :] .- Grid_Of_Values[i, j+1, k, :]
        δu_sub_3_2   ::Vector  = Grid_Of_Values[i, j-1, k, :] .- Grid_Of_Values[i, j-2, k, :]

        r            ::Vector  = (Grid_Of_Values[i, j,   k, :] .- Grid_Of_Values[i, j-1, k, :]) ./ (Grid_Of_Values[i, j+1, k, :] .- Grid_Of_Values[i, j,   k, :])
        r_plus_1     ::Vector  = (Grid_Of_Values[i, j+1, k, :] .- Grid_Of_Values[i, j,   k, :]) ./ (Grid_Of_Values[i, j+2, k, :] .- Grid_Of_Values[i, j+1, k, :])
        r_sub_1      ::Vector  = (Grid_Of_Values[i, j-1, k, :] .- Grid_Of_Values[i, j-2, k, :]) ./ (Grid_Of_Values[i, j,   k, :] .- Grid_Of_Values[i, j-1, k, :])

        u_plus_1_2_L ::Vector  = Grid_Of_Values[i, j,   k, :] .+ minmod(r        )/4 .* ((1-κ)* δu_sub_1_2  .+ (1+κ) * δu_plus_1_2)
        u_plus_1_2_R ::Vector  = Grid_Of_Values[i, j+1, k, :] .- minmod(r_plus_1 )/4 .* ((1-κ)* δu_plus_3_2 .+ (1+κ) * δu_plus_1_2)
        u_sub_1_2_L  ::Vector  = Grid_Of_Values[i, j-1, k, :] .+ minmod(r_sub_1  )/4 .* ((1-κ)* δu_sub_3_2  .+ (1+κ) * δu_sub_1_2 )
        u_sub_1_2_R  ::Vector  = Grid_Of_Values[i, j,   k, :] .- minmod(r        )/4 .* ((1-κ)* δu_plus_1_2 .+ (1+κ) * δu_sub_1_2 )
	

return u_plus_1_2_L, u_plus_1_2_R, u_sub_1_2_L, u_sub_1_2_R

end



function PPM_N3(Grid_Of_Values::Array{Float64, 4}, i::Int64, j::Int64, k::Int64)
	
	#PPM
	κ::Float64 = 1/3
	
	δu_plus_1_2  ::Vector = Grid_Of_Values[i, j, k+1, :] .- Grid_Of_Values[i, j, k,   :]
	δu_sub_1_2   ::Vector = Grid_Of_Values[i, j, k,   :] .- Grid_Of_Values[i, j, k-1, :]
	δu_plus_3_2  ::Vector = Grid_Of_Values[i, j, k+2, :] .- Grid_Of_Values[i, j, k+1, :]
	δu_sub_3_2   ::Vector = Grid_Of_Values[i, j, k-1, :] .- Grid_Of_Values[i, j, k-2, :]

	r            ::Vector = (Grid_Of_Values[i, j, k,   :] .- Grid_Of_Values[i, j, k-1, :]) ./ (Grid_Of_Values[i, j, k+1, :] .- Grid_Of_Values[i, j, k,   :])
	r_plus_1     ::Vector = (Grid_Of_Values[i, j, k+1, :] .- Grid_Of_Values[i, j, k,   :]) ./ (Grid_Of_Values[i, j, k+2, :] .- Grid_Of_Values[i, j, k+1, :])
	r_sub_1      ::Vector = (Grid_Of_Values[i, j, k-1, :] .- Grid_Of_Values[i, j, k-2, :]) ./ (Grid_Of_Values[i, j, k,   :] .- Grid_Of_Values[i, j, k-1, :])

	u_plus_1_2_L ::Vector = Grid_Of_Values[i, j, k,   :] .+ minmod(r        )/4 .* ((1-κ)* δu_sub_1_2  .+ (1+κ) * δu_plus_1_2)
	u_plus_1_2_R ::Vector = Grid_Of_Values[i, j, k+1, :] .- minmod(r_plus_1 )/4 .* ((1-κ)* δu_plus_3_2 .+ (1+κ) * δu_plus_1_2)
	u_sub_1_2_L  ::Vector = Grid_Of_Values[i, j, k-1, :] .+ minmod(r_sub_1  )/4 .* ((1-κ)* δu_sub_3_2  .+ (1+κ) * δu_sub_1_2 )
	u_sub_1_2_R  ::Vector = Grid_Of_Values[i, j, k,   :] .- minmod(r        )/4 .* ((1-κ)* δu_plus_1_2 .+ (1+κ) * δu_sub_1_2 )

return u_plus_1_2_L, u_plus_1_2_R, u_sub_1_2_L, u_sub_1_2_R

end
