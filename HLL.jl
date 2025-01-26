using BenchmarkTools
using Base.Threads
using Profile

include("Matrix_operations.jl")
include("Christoffel_Symbols.jl")
include("structs.jl")
include("Initial_Conditions.jl")
include("fluxlimiter.jl")

function HLL()
a = 0.8
buffer_U = zeros(NP)
buffer_U_true =zeros(N1,N2,N3,NP)

@views @threads for i in 3:(N1-2)
   for j in 3:(N2-2)
          for k in 3:(N3-2)
                # Call PtoU with the grid values and store the results directly in buffer_U
                PtoU(grid[i, j, k, :], buffer_U, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a), eos)       
                buffer_U_true[i, j, k, :] = buffer_U
                
                #PPM METHOD
                κ = 1/3
                
                #δu calculations in N1 - direction
                δu_i_plus_1_2_N1 = grid[i+1, j, k,:] .- grid[i, j, k,  :]
                δu_i_sub_1_2_N1  = grid[i, j, k,  :] .- grid[i-1, j, k,:]
                δu_i_plus_3_2_N1 = grid[i+2, j, k,:] .- grid[i+1, j, k,:]
                δu_i_sub_3_2_N1  = grid[i-1, j, k,:] .- grid[i-2, j, k,:]                
                
                #r calculations in N1 - direction
		r_i_N1        = (grid[i, j, k,  :] .- grid[i-1, j, k,:]) ./ (grid[i+1, j, k,:] .- grid[i, j, k,  :])
                r_i_plus_1_N1 = (grid[i+1, j, k,:] .- grid[i, j, k,  :]) ./ (grid[i+2, j, k,:] .- grid[i+1, j, k,:])
                r_i_sub_1_N1  = (grid[i-1, j, k,:] .- grid[i-2, j, k,:]) ./ (grid[i, j, k,  :] .- grid[i-1, j, k,:])
                
                #Final N1 interpolation
                u_i_plus_1_2_L_N1 = grid[i, j, k,  :] .+ minmod(r_i_N1       )/4 .* ((1-κ)* δu_i_sub_1_2_N1  .+ (1+κ) * δu_i_plus_1_2_N1)
                u_i_plus_1_2_R_N1 = grid[i+1, j, k,:] .- minmod(r_i_plus_1_N1)/4 .* ((1-κ)* δu_i_plus_3_2_N1 .+ (1+κ) * δu_i_plus_1_2_N1)                
                u_i_sub_1_2_L_N1  = grid[i-1, j, k,:] .+ minmod(r_i_sub_1_N1 )/4 .* ((1-κ)* δu_i_sub_3_2_N1  .+ (1+κ) * δu_i_sub_1_2_N1 )
                u_i_sub_1_2_R_N1  = grid[i, j, k,  :] .- minmod(r_i_N1       )/4 .* ((1-κ)* δu_i_plus_1_2_N1 .+ (1+κ) * δu_i_sub_1_2_N1 )                   

                #δu calculations in N2 - direction
                δu_i_plus_1_2_N2 = grid[i, j+1, k,:] .- grid[i, j, k,  :]
                δu_i_sub_1_2_N2  = grid[i, j, k,  :] .- grid[i, j-1, k,:]
                δu_i_plus_3_2_N2 = grid[i, j+2, k,:] .- grid[i, j+1, k,:]
                δu_i_sub_3_2_N2  = grid[i, j-1, k,:] .- grid[i, j-2, k,:]                
                
                #r calculations in N2 - direction
		r_i_N2        = (grid[i, j, k,  :] .- grid[i, j-1, k,:]) ./ (grid[i, j+1, k,:] .- grid[i, j, k,  :])
                r_i_plus_1_N2 = (grid[i, j+1, k,:] .- grid[i, j, k,  :]) ./ (grid[i, j+2, k,:] .- grid[i, j+1, k,:])
                r_i_sub_1_N2  = (grid[i, j-1, k,:] .- grid[i, j-2, k,:]) ./ (grid[i, j, k,  :] .- grid[i, j-1, k,:])
                
                #Final N2 interpolation
                u_i_plus_1_2_L_N2 = grid[i, j, k,  :] .+ minmod(r_i_N2       )/4 .* ((1-κ)* δu_i_sub_1_2_N2  .+ (1+κ) * δu_i_plus_1_2_N2)
                u_i_plus_1_2_R_N2 = grid[i, j+1, k,:] .- minmod(r_i_plus_1_N2)/4 .* ((1-κ)* δu_i_plus_3_2_N2 .+ (1+κ) * δu_i_plus_1_2_N2)                
                u_i_sub_1_2_L_N2  = grid[i, j-1, k,:] .+ minmod(r_i_sub_1_N2 )/4 .* ((1-κ)* δu_i_sub_3_2_N2  .+ (1+κ) * δu_i_sub_1_2_N2 )
                u_i_sub_1_2_R_N2  = grid[i, j, k,  :] .- minmod(r_i_N2       )/4 .* ((1-κ)* δu_i_plus_1_2_N2 .+ (1+κ) * δu_i_sub_1_2_N2 )    

                #δu calculations in N3 - direction
                δu_i_plus_1_2_N3 = grid[i, j, k+1,:] .- grid[i, j, k,  :]
                δu_i_sub_1_2_N3  = grid[i, j, k,  :] .- grid[i, j, k-2,:]
                δu_i_plus_3_2_N3 = grid[i, j, k+2,:] .- grid[i, j, k+1,:]
                δu_i_sub_3_2_N3  = grid[i, j, k-1,:] .- grid[i, j, k-2,:]                
                
                #r calculations in N3 - direction
		r_i_N3        = (grid[i, j, k,  :] .- grid[i, j, k-1,:]) ./ (grid[i, j, k+1,:] .- grid[i, j, k,  :])
                r_i_plus_1_N3 = (grid[i, j, k+1,:] .- grid[i, j, k,  :]) ./ (grid[i, j, k+2,:] .- grid[i, j, k+1,:])
                r_i_sub_1_N3  = (grid[i, j, k-1,:] .- grid[i, j, k-2,:]) ./ (grid[i, j, k,  :] .- grid[i, j, k-1,:])
                
                #Final N3 interpolation
                u_i_plus_1_2_L_N3 = grid[i, j, k,  :] .+ minmod(r_i_N3       )/4 .* ((1-κ)* δu_i_sub_1_2_N3  .+ (1+κ) * δu_i_plus_1_2_N3)
                u_i_plus_1_2_R_N3 = grid[i, j, k+1,:] .- minmod(r_i_plus_1_N3)/4 .* ((1-κ)* δu_i_plus_3_2_N3 .+ (1+κ) * δu_i_plus_1_2_N3)                
                u_i_sub_1_2_L_N3  = grid[i, j, k-1,:] .+ minmod(r_i_sub_1_N3 )/4 .* ((1-κ)* δu_i_sub_3_2_N3  .+ (1+κ) * δu_i_sub_1_2_N3 )
                u_i_sub_1_2_R_N3  = grid[i, j, k,  :] .- minmod(r_i_N3       )/4 .* ((1-κ)* δu_i_plus_1_2_N3 .+ (1+κ) * δu_i_sub_1_2_N3 )    
                
            end
        end
    end
end

@time HLL()
				
