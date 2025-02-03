using Base: mod1

# Funkcja opakowująca indeksy – zakładamy, że siatka ma długość n
periodic(i, n) = mod1(i, n)

function cmin_cmax_N1(q_i_plus_1_2::AbstractVector,
                      q_i_sub_1_2::AbstractVector,
                      CL::Float64,
                      CR::Float64,
                      i::Int64,
                      j::Int64,
                      k::Int64,
                      a::Float64,
                      eos::Polytrope)
    
    # Zakładamy, że N1_grid to tablica z wartościami w kierunku 1.
    n1 = length(N1_grid)
    
    # Używamy opakowania indeksów przy wyliczaniu średniej współrzędnej:
    xL = (N1_grid[periodic(i, n1)] + N1_grid[periodic(i-1, n1)])/2
    xR = (N1_grid[periodic(i+1, n1)] + N1_grid[periodic(i, n1)])/2

    lorL = Lorentz_factor(q_i_sub_1_2, Kerr_Schild_metric_cov(xL, N2_grid[j], N3_grid[k], a), eos)
    lorR = Lorentz_factor(q_i_plus_1_2, Kerr_Schild_metric_cov(xR, N2_grid[j], N3_grid[k], a), eos)
    
    vL = q_i_sub_1_2[3] / lorL
    vR = q_i_plus_1_2[3] / lorR
    
    sigma_S_L = CL^2 / (lorL^2 * (1-CL^2))
    sigma_S_R = CR^2 / (lorR^2 * (1-CR^2))
    
    C_max = max( (vL + sqrt(sigma_S_L*(1-vL^2+sigma_S_L))) / (1+sigma_S_L),
                 (vR + sqrt(sigma_S_R*(1-vR^2+sigma_S_R))) / (1+sigma_S_R) )
    C_min = -min( (vL - sqrt(sigma_S_L*(1-vL^2+sigma_S_L))) / (1+sigma_S_L),
                  (vR - sqrt(sigma_S_R*(1-vR^2+sigma_S_R))) / (1+sigma_S_R) )
    
    return C_max, C_min
end

function cmin_cmax_N2(q_i_plus_1_2::AbstractVector,
                      q_i_sub_1_2::AbstractVector,
                      CL::Float64,
                      CR::Float64,
                      i::Int64,
                      j::Int64,
                      k::Int64,
                      a::Float64,
                      eos::Polytrope)
    
    n2 = length(N2_grid)
    
    yL = (N2_grid[periodic(j, n2)] + N2_grid[periodic(j-1, n2)])/2
    yR = (N2_grid[periodic(j+1, n2)] + N2_grid[periodic(j, n2)])/2

    lorL = Lorentz_factor(q_i_sub_1_2, Kerr_Schild_metric_cov(N1_grid[i], yL, N3_grid[k], a), eos)
    lorR = Lorentz_factor(q_i_plus_1_2, Kerr_Schild_metric_cov(N1_grid[i], yR, N3_grid[k], a), eos)
    
    vL = q_i_sub_1_2[4] / lorL
    vR = q_i_plus_1_2[4] / lorR
    
    sigma_S_L = CL^2 / (lorL^2 * (1-CL^2))
    sigma_S_R = CR^2 / (lorR^2 * (1-CR^2))
    
    C_max = max( (vL + sqrt(sigma_S_L*(1-vL^2+sigma_S_L))) / (1+sigma_S_L),
                 (vR + sqrt(sigma_S_R*(1-vR^2+sigma_S_R))) / (1+sigma_S_R) )
    C_min = -min( (vL - sqrt(sigma_S_L*(1-vL^2+sigma_S_L))) / (1+sigma_S_L),
                  (vR - sqrt(sigma_S_R*(1-vR^2+sigma_S_R))) / (1+sigma_S_R) )
    
    return C_max, C_min
end

function cmin_cmax_N3(q_i_plus_1_2::AbstractVector,
                      q_i_sub_1_2::AbstractVector,
                      CL::Float64,
                      CR::Float64,
                      i::Int64,
                      j::Int64,
                      k::Int64,
                      a::Float64,
                      eos::Polytrope)
    
    n3 = length(N3_grid)
    
    zL = (N3_grid[periodic(k-1, n3)] + N3_grid[k])/2  # tu można również opakować k, jeśli potrzeba
    zR = (N3_grid[periodic(k+1, n3)] + N3_grid[k])/2

    lorL = Lorentz_factor(q_i_sub_1_2, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], zL, a), eos)
    lorR = Lorentz_factor(q_i_plus_1_2, Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], zR, a), eos)
    
    vL = q_i_sub_1_2[5] / lorL
    vR = q_i_plus_1_2[5] / lorR
    
    sigma_S_L = CL^2 / (lorL^2 * (1-CL^2))
    sigma_S_R = CR^2 / (lorR^2 * (1-CR^2))
    
    C_max = max( (vL + sqrt(sigma_S_L*(1-vL^2+sigma_S_L))) / (1+sigma_S_L),
                 (vR + sqrt(sigma_S_R*(1-vR^2+sigma_S_R))) / (1+sigma_S_R) )
    C_min = -min( (vL - sqrt(sigma_S_L*(1-vL^2+sigma_S_L))) / (1+sigma_S_L),
                  (vR - sqrt(sigma_S_R*(1-vR^2+sigma_S_R))) / (1+sigma_S_R) )
    
    return C_max, C_min
end

