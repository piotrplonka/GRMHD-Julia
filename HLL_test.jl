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
include("GR_Flux.jl")
include("HDF.jl")

# Funkcja pomocnicza do opakowywania indeksu (periodyczność)
periodic(i, n) = mod1(i, n)  # mod1 zapewnia, że wynik należy do 1:n

function HLL()
    a = 0.8

    # Prealokacja buforów – zakładamy, że zmienne globalne N1,N2,N3, NP, grid, N1_grid, N2_grid, N3_grid, eos są już zdefiniowane
    buffer_P = zeros(N1, N2, N3, 8)
    buffer_U = zeros(N1, N2, N3, 8)
    buffer_flux = zeros(N1, N2, N3, 8, 3)
    buffer_flux_GR = zeros(N1, N2, N3, 5)
    buffer_flux_magnetic = zeros(N1, N2, N3, 3)

    buffer_Fx_R = zeros(NP)
    buffer_Fx_L = zeros(NP)
    buffer_Fy_R = zeros(NP)
    buffer_Fy_L = zeros(NP)
    buffer_Fz_R = zeros(NP)
    buffer_Fz_L = zeros(NP)

    # Przygotowanie pomocniczych wektorów – tworzymy je lokalnie wewnątrz pętli gdy to możliwe
    stupid_U = zeros(8)
    stupid_P = zeros(8)

    dx = (N1_grid[end] - N1_grid[1]) / N1
    dy = (N2_grid[end] - N2_grid[1]) / N2
    dz = (N3_grid[end] - N3_grid[1]) / N3
    c = 1.0
    tot = 0.0
    dt = 0.1 * dx / c

    T = 100 * dt
    println("Wykonuje przerzucenie z siatki primitive variable.")

    @threads for i in 1:N1
        for j in 1:N2
            for k in 1:N3
                @inbounds begin
                    buffer_P[i, j, k, :] .= grid[i, j, k, :]
                    if buffer_P[i, j, k, 1] < 0
                        println("Ujemna wartość density na pozycji ($i,$j,$k): ", buffer_P[i, j, k, 1])
                        buffer_P[i, j, k, 1] = 1e-5
                    end
                    if buffer_P[i, j, k, 2] < 0
                        println("Ujemna wartość pressure na pozycji ($i,$j,$k): ", buffer_P[i, j, k, 2])
                        buffer_P[i, j, k, 2] = 1e-5
                    end
                end
            end
        end
    end

    println("Wykonuje obliczenie wartości zachowanych.")
    @threads for i in 1:N1
        for j in 1:N2
            for k in 1:N3
                @inbounds begin
                    local_local_U = similar(stupid_U)
                    # Używamy bieżących współrzędnych – zakładamy, że N1_grid, N2_grid, N3_grid są tablicami
                    PtoU(buffer_P[i, j, k, :], local_local_U,
                         Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a), eos)
                    buffer_U[i, j, k, :] .= local_local_U
                end
            end
        end
    end

    number_of_iteration::Int64 = 0

    while tot < T
        println("Krok:", tot/T)
        save_dump_HDF(number_of_iteration, buffer_P)
        println("Wykonuje obliczenie fluxów dla kroku: ", number_of_iteration)
        # Zmieniamy pętle, aby przebiegały po całej siatce – dzięki periodycznym warunkom brzegowym
        @threads for i in 1:N1
            for j in 1:N2
                for k in 1:N3
                    @inbounds begin

                        # Przy definiowaniu rekonstrukcji korzystamy z opakowanych indeksów
                        q_i_plus_1_2_N1_R = zeros(8)
                        q_i_sub_1_2_N1_L  = zeros(8)
                        q_i_plus_1_2_N2_R = zeros(8)
                        q_i_sub_1_2_N2_L  = zeros(8)
                        q_i_plus_1_2_N3_R = zeros(8)
                        q_i_sub_1_2_N3_L  = zeros(8)

                        # Dla kierunku 1 (x)
                        _, q_i_sub_1_2_N1_L = WENOZ_vec( 
                            buffer_P[periodic(i-2, N1), j, k, :],
                            buffer_P[periodic(i-1, N1), j, k, :],
                            buffer_P[i, j, k, :],
                            buffer_P[periodic(i+1, N1), j, k, :],
                            buffer_P[periodic(i+2, N1), j, k, :]
                        )
                        q_i_sub_1_2_N1_L[1:2] = max.(q_i_sub_1_2_N1_L[1:2], 1e-7)

                        q_i_plus_1_2_N1_R, _ = WENOZ_vec(
                            buffer_P[periodic(i-1, N1), j, k, :],
                            buffer_P[i, j, k, :],
                            buffer_P[periodic(i+1, N1), j, k, :],
                            buffer_P[periodic(i+2, N1), j, k, :],
                            buffer_P[periodic(i+3, N1), j, k, :]
                        )
                        q_i_plus_1_2_N1_R[1:2] = max.(q_i_plus_1_2_N1_R[1:2], 1e-7)

                        # Dla kierunku 2 (y)
                        _, q_i_sub_1_2_N2_L = WENOZ_vec(
                            buffer_P[i, periodic(j-2, N2), k, :],
                            buffer_P[i, periodic(j-1, N2), k, :],
                            buffer_P[i, j, k, :],
                            buffer_P[i, periodic(j+1, N2), k, :],
                            buffer_P[i, periodic(j+2, N2), k, :]
                        )
                        q_i_sub_1_2_N2_L[1:2] = max.(q_i_sub_1_2_N2_L[1:2], 1e-7)

                        q_i_plus_1_2_N2_R, _ = WENOZ_vec(
                            buffer_P[i, periodic(j-1, N2), k, :],
                            buffer_P[i, j, k, :],
                            buffer_P[i, periodic(j+1, N2), k, :],
                            buffer_P[i, periodic(j+2, N2), k, :],
                            buffer_P[i, periodic(j+3, N2), k, :]
                        )
                        q_i_plus_1_2_N2_R[1:2] = max.(q_i_plus_1_2_N2_R[1:2], 1e-7)

                        # Dla kierunku 3 (z)
                        _, q_i_sub_1_2_N3_L = WENOZ_vec(
                            buffer_P[i, j, periodic(k-2, N3), :],
                            buffer_P[i, j, periodic(k-1, N3), :],
                            buffer_P[i, j, k, :],
                            buffer_P[i, j, periodic(k+1, N3), :],
                            buffer_P[i, j, periodic(k+2, N3), :]
                        )
                        q_i_sub_1_2_N3_L[1:2] = max.(q_i_sub_1_2_N3_L[1:2], 1e-7)

                        q_i_plus_1_2_N3_R, _ = WENOZ_vec(
                            buffer_P[i, j, periodic(k-1, N3), :],
                            buffer_P[i, j, k, :],
                            buffer_P[i, j, periodic(k+1, N3), :],
                            buffer_P[i, j, periodic(k+2, N3), :],
                            buffer_P[i, j, periodic(k+3, N3), :]
                        )
                        q_i_plus_1_2_N3_R[1:2] = max.(q_i_plus_1_2_N3_R[1:2], 1e-7)

                        # Przygotowanie lokalnych wektorów flux (dla uproszczenia zakładamy, że lokalne bufor fluxów mają rozmiar NP)
                        local_Fx_R = similar(buffer_Fx_R)
                        local_Fx_L = similar(buffer_Fx_L)
                        local_Fy_R = similar(buffer_Fy_R)
                        local_Fy_L = similar(buffer_Fy_L)
                        local_Fz_R = similar(buffer_Fz_R)
                        local_Fz_L = similar(buffer_Fz_L)

                        # Obliczenia fluxów w kierunku N1 (x)
                        PtoFx(q_i_plus_1_2_N1_R, local_Fx_R,
                              Kerr_Schild_metric_cov((N1_grid[periodic(i, N1)] + N1_grid[periodic(i+1, N1)])/2,
                                                     N2_grid[j], N3_grid[k], a), eos)
                        PtoFx(q_i_sub_1_2_N1_L, local_Fx_L,
                              Kerr_Schild_metric_cov((N1_grid[periodic(i, N1)] + N1_grid[periodic(i-1, N1)])/2,
                                                     N2_grid[j], N3_grid[k], a), eos)

                        # Obliczenia fluxów w kierunku N2 (y)
                        PtoFy(q_i_plus_1_2_N2_R, local_Fy_R,
                              Kerr_Schild_metric_cov(N1_grid[i],
                                                     (N2_grid[periodic(j, N2)] + N2_grid[periodic(j+1, N2)])/2,
                                                     N3_grid[k], a), eos)
                        PtoFy(q_i_sub_1_2_N2_L, local_Fy_L,
                              Kerr_Schild_metric_cov(N1_grid[i],
                                                     (N2_grid[periodic(j, N2)] + N2_grid[periodic(j-1, N2)])/2,
                                                     N3_grid[k], a), eos)

                        # Obliczenia fluxów w kierunku N3 (z)
                        PtoFz(q_i_plus_1_2_N3_R, local_Fz_R,
                              Kerr_Schild_metric_cov(N1_grid[i],
                                                     N2_grid[j],
                                                     (N3_grid[periodic(k, N3)] + N3_grid[periodic(k+1, N3)])/2,
                                                     a), eos)
                        PtoFz(q_i_sub_1_2_N3_L, local_Fz_L,
                              Kerr_Schild_metric_cov(N1_grid[i],
                                                     N2_grid[j],
                                                     (N3_grid[periodic(k, N3)] + N3_grid[periodic(k-1, N3)])/2,
                                                     a), eos)

                        # Obliczenie prędkości dźwięku – analogicznie z wykorzystaniem opakowywania indeksów
                        V_L_N1 = SoundSpeed(q_i_sub_1_2_N1_L,
                            Kerr_Schild_metric_cov((N1_grid[periodic(i, N1)] + N1_grid[periodic(i-1, N1)])/2,
                                                   N2_grid[j], N3_grid[k], a), eos)
                        V_R_N1 = SoundSpeed(q_i_plus_1_2_N1_R,
                            Kerr_Schild_metric_cov((N1_grid[periodic(i+1, N1)] + N1_grid[periodic(i, N1)])/2,
                                                   N2_grid[j], N3_grid[k], a), eos)

                        V_L_N2 = SoundSpeed(q_i_sub_1_2_N2_L,
                            Kerr_Schild_metric_cov(N1_grid[i],
                                                   (N2_grid[periodic(j, N2)] + N2_grid[periodic(j-1, N2)])/2,
                                                   N3_grid[k], a), eos)
                        V_R_N2 = SoundSpeed(q_i_plus_1_2_N2_R,
                            Kerr_Schild_metric_cov(N1_grid[i],
                                                   (N2_grid[periodic(j+1, N2)] + N2_grid[periodic(j, N2)])/2,
                                                   N3_grid[k], a), eos)

                        V_L_N3 = SoundSpeed(q_i_sub_1_2_N3_L,
                            Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j],
                                                   (N3_grid[periodic(k, N3)] + N3_grid[periodic(k-1, N3)])/2,
                                                   a), eos)
                        V_R_N3 = SoundSpeed(q_i_plus_1_2_N3_R,
                            Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j],
                                                   (N3_grid[periodic(k+1, N3)] + N3_grid[periodic(k, N3)])/2,
                                                   a), eos)

                        UL_N1 = zeros(8)
                        UR_N1 = zeros(8)
                        UL_N2 = zeros(8)
                        UR_N2 = zeros(8)
                        UL_N3 = zeros(8)
                        UR_N3 = zeros(8)

                        PtoU(q_i_sub_1_2_N1_L, UL_N1,
                              Kerr_Schild_metric_cov((N1_grid[periodic(i, N1)] + N1_grid[periodic(i-1, N1)])/2,
                                                     N2_grid[j], N3_grid[k], a), eos)
                        PtoU(q_i_plus_1_2_N1_R, UR_N1,
                              Kerr_Schild_metric_cov((N1_grid[periodic(i+1, N1)] + N1_grid[periodic(i, N1)])/2,
                                                     N2_grid[j], N3_grid[k], a), eos)

                        PtoU(q_i_sub_1_2_N2_L, UL_N2,
                              Kerr_Schild_metric_cov(N1_grid[i],
                                                     (N2_grid[periodic(j, N2)] + N2_grid[periodic(j-1, N2)])/2,
                                                     N3_grid[k], a), eos)
                        PtoU(q_i_plus_1_2_N2_R, UR_N2,
                              Kerr_Schild_metric_cov(N1_grid[i],
                                                     (N2_grid[periodic(j+1, N2)] + N2_grid[periodic(j, N2)])/2,
                                                     N3_grid[k], a), eos)

                        PtoU(q_i_sub_1_2_N3_L, UL_N3,
                              Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j],
                                                     (N3_grid[periodic(k, N3)] + N3_grid[periodic(k-1, N3)])/2,
                                                     a), eos)
                        PtoU(q_i_plus_1_2_N3_R, UR_N3,
                              Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j],
                                                     (N3_grid[periodic(k+1, N3)] + N3_grid[periodic(k, N3)])/2,
                                                     a), eos)

                        # Obliczamy minimalne i maksymalne prędkości fal – analogicznie
                        C_max_N1, C_min_N1 = cmin_cmax_N1(q_i_plus_1_2_N1_R, q_i_sub_1_2_N1_L, V_L_N1, V_R_N1, i, j, k, a, eos::Polytrope)
                        C_max_N2, C_min_N2 = cmin_cmax_N2(q_i_plus_1_2_N2_R, q_i_sub_1_2_N2_L, V_L_N2, V_R_N2, i, j, k, a, eos::Polytrope)
                        C_max_N3, C_min_N3 = cmin_cmax_N3(q_i_plus_1_2_N3_R, q_i_sub_1_2_N3_L, V_L_N3, V_R_N3, i, j, k, a, eos::Polytrope)

                        # Obliczenie fluxu metodą HLL – analogicznie jak wcześniej
                        if C_max_N1 < 0
                            buffer_flux[i, j, k, 1:5, 1] .= local_Fx_R[1:5]
                        elseif C_min_N1 < 0
                            buffer_flux[i, j, k, 1:5, 1] .= local_Fx_L[1:5]
                        else
                            buffer_flux[i, j, k, 1:5, 1] .= (C_min_N1 .* local_Fx_R[1:5] +
                                                             C_max_N1 .* local_Fx_L[1:5] -
                                                             C_max_N1 * C_min_N1 .* (UR_N1[1:5] .- UL_N1[1:5])) /
                                                            (C_max_N1 + C_min_N1)
                        end

                        if C_max_N2 < 0
                            buffer_flux[i, j, k, 1:5, 2] .= local_Fy_R[1:5]
                        elseif C_min_N2 < 0
                            buffer_flux[i, j, k, 1:5, 2] .= local_Fy_L[1:5]
                        else
                            buffer_flux[i, j, k, 1:5, 2] .= (C_min_N2 .* local_Fy_R[1:5] +
                                                             C_max_N2 .* local_Fy_L[1:5] -
                                                             C_max_N2 * C_min_N2 .* (UR_N2[1:5] .- UL_N2[1:5])) /
                                                            (C_max_N2 + C_min_N2)
                        end

                        if C_max_N3 < 0
                            buffer_flux[i, j, k, 1:5, 3] .= local_Fz_R[1:5]
                        elseif C_min_N3 < 0
                            buffer_flux[i, j, k, 1:5, 3] .= local_Fz_L[1:5]
                        else
                            buffer_flux[i, j, k, 1:5, 3] .= (C_min_N3 .* local_Fz_R[1:5] +
                                                             C_max_N3 .* local_Fz_L[1:5] -
                                                             C_max_N3 * C_min_N3 .* (UR_N3[1:5] .- UL_N3[1:5])) /
                                                            (C_max_N3 + C_min_N3)
                        end

                        buffer_flux_magnetic[i, j, k, 1] = PtoF_Bx(q_i_plus_1_2_N1_R,
                            Kerr_Schild_metric_cov((N1_grid[periodic(i+1, N1)] + N1_grid[periodic(i, N1)])/2,
                                                   N2_grid[j], N3_grid[k], a), eos::Polytrope) -
                            PtoF_Bx(q_i_sub_1_2_N1_L,
                            Kerr_Schild_metric_cov((N1_grid[periodic(i, N1)] + N1_grid[periodic(i-1, N1)])/2,
                                                   N2_grid[j], N3_grid[k], a), eos::Polytrope)

                        buffer_flux_magnetic[i, j, k, 2] = PtoF_By(q_i_plus_1_2_N2_R,
                            Kerr_Schild_metric_cov(N1_grid[i],
                                                   (N2_grid[periodic(j+1, N2)] + N2_grid[periodic(j, N2)])/2,
                                                   N3_grid[k], a), eos::Polytrope) -
                            PtoF_By(q_i_sub_1_2_N2_L,
                            Kerr_Schild_metric_cov(N1_grid[i],
                                                   (N2_grid[periodic(j, N2)] + N2_grid[periodic(j-1, N2)])/2,
                                                   N3_grid[k], a), eos::Polytrope)

                        buffer_flux_magnetic[i, j, k, 3] = PtoF_Bz(q_i_plus_1_2_N3_R,
                            Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j],
                                                   (N3_grid[periodic(k+1, N3)] + N3_grid[periodic(k, N3)])/2,
                                                   a), eos::Polytrope) -
                            PtoF_Bz(q_i_sub_1_2_N3_L,
                            Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j],
                                                   (N3_grid[periodic(k, N3)] + N3_grid[periodic(k-1, N3)])/2,
                                                   a), eos::Polytrope)
                    end
                end
            end
        end

        println("Wykonuje przeliczenie nowych wartości zachowanych dla kroku: ", number_of_iteration)

        # Aktualizacja wartości U – również z periodycznym opakowaniem indeksów przy pobieraniu fluxów
        @threads for i in 1:N1
            for j in 1:N2
                for k in 1:N3
                    @inbounds begin
                        # Dla różnic centralnych stosujemy opakowywanie indeksów
                        buffer_U[i, j, k, 1:5] .= buffer_U[i, j, k, 1:5] -
                            (dt/dx) * (buffer_flux[i, j, k, 1:5, 1] - buffer_flux[periodic(i-1, N1), j, k, 1:5, 1]) -
                            (dt/dy) * (buffer_flux[i, j, k, 1:5, 2] - buffer_flux[i, periodic(j-1, N2), k, 1:5, 2]) -
                            (dt/dz) * (buffer_flux[i, j, k, 1:5, 3] - buffer_flux[i, j, periodic(k-1, N3), 1:5, 3])
                            
                        buffer_U[i, j, k, 6] += (dt/dx) * (buffer_flux_magnetic[i, j, k, 1] -
                                                           buffer_flux_magnetic[periodic(i-1, N1), j, k, 1])
                        buffer_U[i, j, k, 7] += (dt/dy) * (buffer_flux_magnetic[i, j, k, 2] -
                                                           buffer_flux_magnetic[i, periodic(j-1, N2), k, 2])
                        buffer_U[i, j, k, 8] += (dt/dz) * (buffer_flux_magnetic[i, j, k, 3] -
                                                           buffer_flux_magnetic[i, j, periodic(k-1, N3), 3])
                    end
                end
            end
        end

        println("Wykonuje obliczenie wartości primitive variable z wartości zachowanych: ", number_of_iteration)

        @threads for i in 1:N1
            for j in 1:N2
                for k in 1:N3
                    @inbounds begin
                        buffer_P[i, j, k, :] .= UtoP_LU(buffer_U[i, j, k, :],
                                                        buffer_P[i, j, k, :],
                                                        Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a),
                                                        eos::Polytrope)
                        buffer_P[i, j, k, 1:2] .= max.(buffer_P[i, j, k, 1:2],1e-7)
                    end
                end
            end
        end

        number_of_iteration += 1
        tot += dt
    end

end

@time HLL()

