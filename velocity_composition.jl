using LinearAlgebra

function cmin_cmax_N1(q_i_sub_1_2::AbstractVector,q_i_plus_1_2::AbstractVector,V_L::Float64,V_R::Float64, i::Int64, j::Int64, k::Int64, a::Float64, eos::Polytrope)

        lorL::Float64 = Lorentz_factor(q_i_sub_1_2,  Kerr_Schild_metric_cov((N1_grid[i]   + N1_grid[i-1])/2, N2_grid[j], N3_grid[k], a), eos)
       	lorR::Float64 = Lorentz_factor(q_i_plus_1_2, Kerr_Schild_metric_cov((N1_grid[i+1] + N1_grid[i])/2  , N2_grid[j], N3_grid[k], a), eos)
       	        
        vL::Float64 = q_i_sub_1_2[3]  / lorL
        vR::Float64 = q_i_plus_1_2[3] / lorL
            	
        sigma_S_L::Float64 = V_L^2 / ( lorL^2 * (1-V_L^2))
        sigma_S_R::Float64 = V_R^2 / ( lorR^2 * (1-V_R^2))
        	
        C_max::Float64 =  max( (vL + sqrt(sigma_S_L * (1-vL^2 + sigma_S_L)) ) / (1 + sigma_S_L), (vR + sqrt(sigma_S_R * (1-vR^2 + sigma_S_R)) ) / (1 + sigma_S_R)) 
        C_min::Float64 = -min( (vL - sqrt(sigma_S_L * (1-vL^2 + sigma_S_L)) ) / (1 + sigma_S_L), (vR - sqrt(sigma_S_R * (1-vR^2 + sigma_S_R)) ) / (1 + sigma_S_R)) 
        
	return C_max, C_min
end

function cmin_cmax_N2(q_i_sub_1_2::AbstractVector,q_i_plus_1_2::AbstractVector,V_L::Float64,V_R::Float64, i::Int64, j::Int64, k::Int64, a::Float64, eos::Polytrope)

        lorL::Float64 = Lorentz_factor(q_i_sub_1_2,  Kerr_Schild_metric_cov(N1_grid[i],  (N2_grid[j]   + N2_grid[j-1])/2, N3_grid[k], a), eos)
       	lorR::Float64 = Lorentz_factor(q_i_plus_1_2, Kerr_Schild_metric_cov(N1_grid[i] , (N2_grid[j+1] + N2_grid[j])/2, N3_grid[k], a), eos)
       	        
        vL::Float64 = q_i_sub_1_2[3]  / lorL
        vR::Float64 = q_i_plus_1_2[3] / lorL
            	
        sigma_S_L::Float64 = V_L^2 / ( lorL^2 * (1-V_L^2))
        sigma_S_R::Float64 = V_R^2 / ( lorR^2 * (1-V_R^2))
        	
        C_max::Float64 =  max( (vL + sqrt(sigma_S_L * (1-vL^2 + sigma_S_L)) ) / (1 + sigma_S_L), (vR + sqrt(sigma_S_R * (1-vR^2 + sigma_S_R)) ) / (1 + sigma_S_R)) 
        C_min::Float64 = -min( (vL - sqrt(sigma_S_L * (1-vL^2 + sigma_S_L)) ) / (1 + sigma_S_L), (vR - sqrt(sigma_S_R * (1-vR^2 + sigma_S_R)) ) / (1 + sigma_S_R)) 
        
	return C_max, C_min
end

function cmin_cmax_N3(q_i_sub_1_2::AbstractVector,q_i_plus_1_2::AbstractVector,V_L::Float64,V_R::Float64, i::Int64, j::Int64, k::Int64, a::Float64, eos::Polytrope)

        lorL::Float64 = Lorentz_factor(q_i_sub_1_2,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k-1]+N3_grid[k])/2, a), eos)
       	lorR::Float64 = Lorentz_factor(q_i_plus_1_2, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], (N3_grid[k+1]+N3_grid[k])/2, a), eos)
       	        
        vL::Float64 = q_i_sub_1_2[3]  / lorL
        vR::Float64 = q_i_plus_1_2[3] / lorL
            	
        sigma_S_L::Float64 = V_L^2 / ( lorL^2 * (1-V_L^2))
        sigma_S_R::Float64 = V_R^2 / ( lorR^2 * (1-V_R^2))
        	
        C_max::Float64 =  max( (vL + sqrt(sigma_S_L * (1-vL^2 + sigma_S_L)) ) / (1 + sigma_S_L), (vR + sqrt(sigma_S_R * (1-vR^2 + sigma_S_R)) ) / (1 + sigma_S_R)) 
        C_min::Float64 = -min( (vL - sqrt(sigma_S_L * (1-vL^2 + sigma_S_L)) ) / (1 + sigma_S_L), (vR - sqrt(sigma_S_R * (1-vR^2 + sigma_S_R)) ) / (1 + sigma_S_R)) 
        
	return C_max, C_min
end
