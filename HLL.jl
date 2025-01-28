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
include("PPM.jl")

function HLL()
a = 0.8
buffer_U = zeros(NP)
buffer_U_true = zeros(N1,N2,N3,NP)

buffer_Fx_LL=zeros(NP)
buffer_Fx_LR=zeros(NP)
buffer_Fx_RL=zeros(NP)
buffer_Fx_RR=zeros(NP)

buffer_Fy_LL=zeros(NP)
buffer_Fy_LR=zeros(NP)
buffer_Fy_RL=zeros(NP)
buffer_Fy_RR=zeros(NP)

buffer_Fz_LL=zeros(NP)
buffer_Fz_LR=zeros(NP)
buffer_Fz_RL=zeros(NP)
buffer_Fz_RR=zeros(NP)



ULL_N1 = zeros(8)
ULR_N1 = zeros(8)
URL_N1 = zeros(8)
URR_N1 = zeros(8)

ULL_N2 = zeros(8)
ULR_N2 = zeros(8)
URL_N2 = zeros(8)
URR_N2 = zeros(8)               

ULL_N3 = zeros(8)
ULR_N3 = zeros(8)
URL_N3 = zeros(8)
URR_N3 = zeros(8)    


#Call PtoU with the grid values and store the results directly in buffer_U
#PtoU(grid[i, j, k, :], buffer_U, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a), eos)       
#buffer_U_true[i, j, k, :] = buffer_U
#if norm(UtoP(buffer_U, grid[i, j, k, :] .*(1.1), Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a), eos) .-grid[i, j, k, :]) >1
#println("Źle zbiega!")
#end


@views  for i in 3:(N1-2)
        	for j in 3:(N2-2)
       @threads 	for k in 3:(N3-2)
                             
                #N1 interpolation
                u_i_plus_1_2_L_N1, u_i_plus_1_2_R_N1, u_i_sub_1_2_L_N1, u_i_sub_1_2_R_N1 = PPM_N1(grid, i, j, k)

                #N2 interpolation
                u_i_plus_1_2_L_N2, u_i_plus_1_2_R_N2, u_i_sub_1_2_L_N2, u_i_sub_1_2_R_N2 = PPM_N2(grid, i, j, k)

                #N3 interpolation
                u_i_plus_1_2_L_N3, u_i_plus_1_2_R_N3, u_i_sub_1_2_L_N3, u_i_sub_1_2_R_N3 = PPM_N3(grid, i, j, k)
                
                
                #Calculating Fluxes using interpolated values
                PtoFx(u_i_sub_1_2_L_N1,  buffer_Fx_LL, Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i-1])/2, N2_grid[j], N3_grid[k], a), eos)
                PtoFx(u_i_sub_1_2_R_N1,  buffer_Fx_LR, Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i-1])/2, N2_grid[j], N3_grid[k], a), eos)
                PtoFx(u_i_plus_1_2_L_N1, buffer_Fx_RL, Kerr_Schild_metric_cov((N1_grid[i+1] + N1_grid[i])/2  , N2_grid[j], N3_grid[k], a), eos)
                PtoFx(u_i_plus_1_2_R_N1, buffer_Fx_RR, Kerr_Schild_metric_cov((N1_grid[i+1] + N1_grid[i])/2  , N2_grid[j], N3_grid[k], a), eos)
                
                PtoFy(u_i_sub_1_2_L_N2,  buffer_Fy_LL, Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j]   + N2_grid[j-1])/2, N3_grid[k], a), eos)
                PtoFy(u_i_sub_1_2_R_N2,  buffer_Fy_LR, Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j]   + N2_grid[j-1])/2, N3_grid[k], a), eos)
                PtoFy(u_i_plus_1_2_L_N2, buffer_Fy_RL, Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j+1] + N2_grid[j])/2,   N3_grid[k], a), eos)
                PtoFy(u_i_plus_1_2_R_N2, buffer_Fy_RR, Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j+1] + N2_grid[j])/2,   N3_grid[k], a), eos)
                
                PtoFz(u_i_sub_1_2_L_N3,  buffer_Fz_LL, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k]   + N3_grid[k-1])/2, a), eos)
                PtoFz(u_i_sub_1_2_R_N3,  buffer_Fz_LR, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k]   + N3_grid[k-1])/2, a), eos)
                PtoFz(u_i_plus_1_2_L_N3, buffer_Fz_RL, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k+1] + N3_grid[k])/2,   a), eos)
                PtoFz(u_i_plus_1_2_R_N3, buffer_Fz_RR, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k+1] + N3_grid[k])/2,   a), eos)
                
                
                #Velocities of the Magnetosonic Waves
                V_LL_N1 = MagnetosonicSpeed(u_i_sub_1_2_L_N1,  Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i-1])/2, N2_grid[j], N3_grid[k], a), eos)
                V_LR_N1 = MagnetosonicSpeed(u_i_sub_1_2_R_N1,  Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i-1])/2, N2_grid[j], N3_grid[k], a), eos)
                V_RL_N1 = MagnetosonicSpeed(u_i_plus_1_2_L_N1, Kerr_Schild_metric_cov((N1_grid[i+1] + N1_grid[i])/2  , N2_grid[j], N3_grid[k], a), eos)
                V_RR_N1 = MagnetosonicSpeed(u_i_plus_1_2_R_N1, Kerr_Schild_metric_cov((N1_grid[i+1] + N1_grid[i])/2  , N2_grid[j], N3_grid[k], a), eos)

                V_LL_N2 = MagnetosonicSpeed(u_i_sub_1_2_L_N2,  Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j]   + N2_grid[j-1])/2, N3_grid[k], a), eos)
                V_LR_N2 = MagnetosonicSpeed(u_i_sub_1_2_R_N2,  Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j]   + N2_grid[j-1])/2, N3_grid[k], a), eos)
                V_RL_N2 = MagnetosonicSpeed(u_i_plus_1_2_L_N2, Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j+1] + N2_grid[j])/2,   N3_grid[k], a), eos)
                V_RR_N2 = MagnetosonicSpeed(u_i_plus_1_2_R_N2, Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j+1] + N2_grid[j])/2,   N3_grid[k], a), eos)

                V_LL_N3 = MagnetosonicSpeed(u_i_sub_1_2_L_N3,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k]   + N3_grid[k-1])/2, a), eos)
                V_LR_N3 = MagnetosonicSpeed(u_i_sub_1_2_R_N3,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k]   + N3_grid[k-1])/2, a), eos)
                V_RL_N3 = MagnetosonicSpeed(u_i_plus_1_2_L_N3, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k+1] + N3_grid[k])/2,   a), eos)
                V_RR_N3 = MagnetosonicSpeed(u_i_plus_1_2_R_N3, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k+1] + N3_grid[k])/2,   a), eos)


                #UR and UL calculations
                PtoU(u_i_sub_1_2_L_N1,  ULL_N1,  Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i-1])/2, N2_grid[j], N3_grid[k], a), eos)
                PtoU(u_i_sub_1_2_R_N1,  ULR_N1,  Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i-1])/2, N2_grid[j], N3_grid[k], a), eos)
                PtoU(u_i_plus_1_2_L_N1, URL_N1,  Kerr_Schild_metric_cov((N1_grid[i+1] + N1_grid[i])/2  , N2_grid[j], N3_grid[k], a), eos)
                PtoU(u_i_plus_1_2_R_N1, URR_N1,  Kerr_Schild_metric_cov((N1_grid[i+1] + N1_grid[i])/2  , N2_grid[j], N3_grid[k], a), eos)
                
                PtoU(u_i_sub_1_2_L_N2,  ULL_N2,  Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j]   + N2_grid[j-1])/2, N3_grid[k], a), eos)
                PtoU(u_i_sub_1_2_R_N2,  ULR_N2,  Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j]   + N2_grid[j-1])/2, N3_grid[k], a), eos)
                PtoU(u_i_plus_1_2_L_N2, URL_N2,  Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j+1] + N2_grid[j])/2,   N3_grid[k], a), eos)
                PtoU(u_i_plus_1_2_R_N2, URR_N2,  Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j+1] + N2_grid[j])/2,   N3_grid[k], a), eos)
                
                PtoU(u_i_sub_1_2_L_N3,  ULL_N3,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k]   + N3_grid[k-1])/2, a), eos)
                PtoU(u_i_sub_1_2_R_N3,  ULR_N3,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k]   + N3_grid[k-1])/2, a), eos)
                PtoU(u_i_plus_1_2_L_N3, URL_N3,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k+1] + N3_grid[k])/2,   a), eos)
                PtoU(u_i_plus_1_2_R_N3, URR_N3,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k+1] + N3_grid[k])/2,   a), eos)
     
            end
        end
    end
end


@time HLL()

				
