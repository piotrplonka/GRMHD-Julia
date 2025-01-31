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
include("velocity_composition.jl")

function HLL()
a = 0.8
buffer_P = zeros(N1,N2,N3,8)
buffer_U = zeros(N1,N2,N3,8)
buffer_flux = zeros(N1,N2,N3,8,3)

buffer_Fx_R = zeros(NP)
buffer_Fx_L = zeros(NP)

buffer_Fy_R = zeros(NP)
buffer_Fy_L = zeros(NP)

buffer_Fz_R = zeros(NP)
buffer_Fz_L = zeros(NP)

UL_N1 = zeros(8)
UR_N1 = zeros(8)
UL_N2 = zeros(8)
UR_N2 = zeros(8)
UL_N3 = zeros(8)
UR_N3 = zeros(8)

#Stupid!!
stupid_U = zeros(8)
stupid_P = zeros(8)

for i in 1:N1
	for j in 1:N2
		@threads for k in 1:N3
			PtoU(grid[i, j, k, :],  stupid_U,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a), eos)
			buffer_U[i, j, k, :] .= stupid_U
		end
	end
end



for i in 1:N1
	for j in 1:N2
		@threads for k in 1:N3
			buffer_P[i, j, k, :] .= grid[i, j, k, :]
		end
	end
end


for i in 3:(N1-2)
	for j in 3:(N2-2)
@threads 	for k in 3:(N3-2)
                             
                #N1 interpolation
                q_i_plus_1_2_N1, q_i_sub_1_2_N1 = WENOZ_vec(grid[i-2, j, k, :], grid[i-1, j, k, :], grid[i, j, k, :], grid[i+1, j, k, :], grid[i+2, j, k, :])
                
                #N2 interpolation
                q_i_plus_1_2_N2, q_i_sub_1_2_N2 = WENOZ_vec(grid[i, j-2, k, :], grid[i, j-1, k, :], grid[i, j, k, :], grid[i, j+1, k, :], grid[i, j+2, k, :])

                #N3 interpolation
                q_i_plus_1_2_N3, q_i_sub_1_2_N3 = WENOZ_vec(grid[i, j, k-2, :], grid[i, j, k-1, :], grid[i, j, k, :], grid[i, j, k+1, :], grid[i, j, k+2, :])
                
                
                #Calculating Fluxes using interpolated values
                PtoFx(q_i_plus_1_2_N1, buffer_Fx_R, Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i+1])/2, N2_grid[j], N3_grid[k], a), eos)
                PtoFx(q_i_sub_1_2_N1,  buffer_Fx_L, Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i-1])/2, N2_grid[j], N3_grid[k], a), eos)
                
                PtoFy(q_i_plus_1_2_N2, buffer_Fy_R, Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j]   + N2_grid[j+1])/2, N3_grid[k], a), eos)
                PtoFy(q_i_sub_1_2_N2,  buffer_Fy_L, Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j]   + N2_grid[j-1])/2, N3_grid[k], a), eos)

                PtoFz(q_i_plus_1_2_N3, buffer_Fz_R, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k]   + N3_grid[k+1])/2, a), eos)
                PtoFz(q_i_sub_1_2_N3,  buffer_Fz_L, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k]   + N3_grid[k-1])/2, a), eos)

                
                #Velocities of the Magnetosonic Waves
                V_L_N1 = MagnetosonicSpeed(q_i_sub_1_2_N1,   Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i-1])/2, N2_grid[j], N3_grid[k], a), eos)
                V_R_N1 = MagnetosonicSpeed(q_i_plus_1_2_N1,  Kerr_Schild_metric_cov((N1_grid[i+1] + N1_grid[i])/2  , N2_grid[j], N3_grid[k], a), eos)

                V_L_N2 = MagnetosonicSpeed(q_i_sub_1_2_N2,   Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j]   + N2_grid[j-1])/2, N3_grid[k], a), eos)
                V_R_N2 = MagnetosonicSpeed(q_i_plus_1_2_N2,  Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j+1] + N2_grid[j])/2,   N3_grid[k], a), eos)

                V_L_N3 = MagnetosonicSpeed(q_i_sub_1_2_N3,   Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k]   + N3_grid[k-1])/2, a), eos)
                V_R_N3 = MagnetosonicSpeed(q_i_plus_1_2_N3,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k+1] + N3_grid[k])/2,   a), eos)


                #UR and UL calculations
                PtoU(q_i_sub_1_2_N1,  UL_N1,  Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i-1])/2, N2_grid[j], N3_grid[k], a), eos)
                PtoU(q_i_plus_1_2_N1, UR_N1,  Kerr_Schild_metric_cov((N1_grid[i+1] + N1_grid[i])/2  , N2_grid[j], N3_grid[k], a), eos)
                
                PtoU(q_i_sub_1_2_N2,  UL_N2,  Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j]   + N2_grid[j-1])/2, N3_grid[k], a), eos)
                PtoU(q_i_plus_1_2_N2, UR_N2,  Kerr_Schild_metric_cov(N1_grid[i], (N2_grid[j+1] + N2_grid[j])/2,   N3_grid[k], a), eos)
                
                PtoU(q_i_sub_1_2_N3,  UL_N3,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k]   + N3_grid[k-1])/2, a), eos)
                PtoU(q_i_plus_1_2_N3, UR_N3,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k+1] + N3_grid[k])/2,   a), eos)
                
                #Velocities for N1, N2, N3
                C_max_N1, C_min_N1 = cmin_cmax_N1(q_i_plus_1_2_N1, q_i_sub_1_2_N1, V_L_N1, V_R_N1, i, j, k, a, eos::Polytrope)
                C_max_N2, C_min_N2 = cmin_cmax_N2(q_i_plus_1_2_N2, q_i_sub_1_2_N2, V_L_N2, V_R_N2, i, j, k, a, eos::Polytrope)
                C_max_N3, C_min_N3 = cmin_cmax_N3(q_i_plus_1_2_N3, q_i_sub_1_2_N3, V_L_N3, V_R_N3, i, j, k, a, eos::Polytrope)
                
                buffer_flux[i,j,k,:,1] .= (C_min_N1.*buffer_Fx_R + C_max_N1.*buffer_Fx_L - C_max_N1*C_min_N1*(UR_N1 .- UL_N1))/(C_max_N1 + C_min_N1)
                buffer_flux[i,j,k,:,2] .= (C_min_N2.*buffer_Fy_R + C_max_N2.*buffer_Fy_L - C_max_N2*C_min_N2*(UR_N2 .- UL_N2))/(C_max_N2 + C_min_N2) 
                buffer_flux[i,j,k,:,3] .= (C_min_N3.*buffer_Fz_R + C_max_N3.*buffer_Fz_L - C_max_N3*C_min_N3*(UR_N3 .- UL_N3))/(C_max_N3 + C_min_N3)
                
            end
        end
    end


for i in 1:N1
	for j in 1:N2
	@threads for k in 1:N3
			 buffer_P[i, j, k, :] .= UtoP(buffer_U[i,j,k,:], buffer_P[i,j,k,:].*1.3,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a), eos::Polytrope)
		end
	end
end


end




@time HLL()

				
