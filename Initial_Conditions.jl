using Random
using LinearAlgebra   # sqrt, acos, itd.
using BenchmarkTools, Base.Threads, StaticArrays

# ===============================
# Parametry symulacji
# ===============================
floor    = 1e-8
outer    = 1e-4

# Rozmiary przestrzeni – przyjmujemy, że:
# oś X i Z odpowiadają kierunkom poprzecznym,
# oś Y – kierunek dżetu.
box_X  = 50.0   # zakres X: od -50 do 50
box_Y  = 100.0  # zakres Y: np. od -10 do 90 (możesz modyfikować)
box_Z  = 50.0   # zakres Z: od -50 do 50

R_max    = 50.0   # maksymalny promień, po którym gęstość przechodzi do wartości outer
R_eng    = 10.0   # promień rdzenia dżetu (wewnątrz którego gęstość jest wyższa)
Temp     = 0.0    # parametr temperatury (tu nie używany)
Rho0     = 5e-3
Rho_core = 0.2
U0       = 5e-2
angle_jet = 0.1   # kąt rozwarcia stożka dżetu (w radianach)
vmax     = 0.3

# ===============================
# Parametry siatki 3D
# ===============================
N1, N2, N3 = 100, 100, 100   # N1: oś X, N2: oś Y, N3: oś Z
NP = 8                     # liczba zmiennych (1: gęstość, 2: energia, 3-5: prędkości, 6-8: pole B)

# Zakresy fizyczne – możesz je modyfikować
N1_min, N1_max = -box_X, box_X
N2_min, N2_max = -10.0, box_Y-10.0   # przykładowy zakres dla osi Y (np. przesunięty, by bazę dżetu umieścić poniżej 0)
N3_min, N3_max = -box_Z, box_Z

# Tworzymy tablicę na dane
grid = zeros(Float64, N1, N2, N3, NP)

# Tworzymy wektory współrzędnych dla każdej osi
N1_grid = collect(range(N1_min, N1_max, length=N1))
N2_grid = collect(range(N2_min, N2_max, length=N2))
N3_grid = collect(range(N3_min, N3_max, length=N3))

# ===============================
# Inicjalizacja siatki – definiujemy stożek dżetowy
# ===============================
# Założenie: dżet skierowany jest w stronę dodatniej osi Y.
# W danym punkcie (X, Y, Z) obliczamy kąt theta między wektorem pozycji a osią Y.
# Jeśli theta < angle_jet, to punkt należy do wnętrza stożka dżetu.
for i in 1:N1
    X = N1_grid[i]
    for j in 1:N2
        Y = N2_grid[j]
        for k in 1:N3
            Z = N3_grid[k]
            
            # Obliczamy długość wektora pozycji – dodajemy mały offset, aby uniknąć dzielenia przez zero
            R_tot = sqrt(X^2 + Y^2 + Z^2) + 1e-10

            # Obliczamy kąt theta między wektorem (X, Y, Z) a osią Y
            theta = acos(clamp(Y / R_tot, -1.0, 1.0))
            
            # Jeśli punkt znajduje się wewnątrz stożka (i kierunek dżetu to Y>0)
            if theta < angle_jet && Y > 0
                # Ustalamy gęstość – tutaj przyjmujemy, że w obrębie rdzenia dżetu (np. w obrębie pewnego promienia) gęstość jest wyższa
                # Możesz przyjąć inny profil zależny od odległości od osi dżetu
                # W poniższym przykładzie wykorzystujemy odległość poprzeczną od osi dżetu:
                r_cyl = sqrt(X^2 + Z^2)  # odległość od osi Y
                if r_cyl < R_eng
                    density = Rho_core
                elseif r_cyl < R_max
                    density = Rho0 * min((R_max / (r_cyl + 1e-10))^2, Rho0 * 10)
                else
                    density = outer
                end

                # Energia – przyjmujemy wartość U0 z niewielkim szumem
                energy = U0 * (1 + randn() * 3e-3)

                # Prędkość – nadajemy prędkość tylko wewnątrz stożka.
                # W tym przykładzie skaluje się ona ze współczynnikiem zależnym od pozycji wzdłuż osi Y
                # (możesz przyjąć inny profil w zależności od potrzeb)
                v = vmax * (Y / (R_eng + 1e-10))
                v = min(v, vmax)  # ograniczenie maksymalnej prędkości
                gamma = 1 / sqrt(1 - v^2)  # czynnik Lorentza (dla symulacji relatywistycznych)
                # Ustalamy kierunek prędkości – dla stożka dżetowego przyjmujemy, że skierowana jest od początku układu w stronę punktu,
                # ale tylko jeśli punkt należy do stożka; dzięki temu prędkość ma komponenty zarówno poprzeczne, jak i równoległe do Y.
                vx = gamma * v * (X / R_tot)
                vy = gamma * v * (Y / R_tot)
                vz = gamma * v * (Z / R_tot)
            else
                # Poza stożkiem dżetu – przyjmujemy warunki tła
                density = outer
                energy = U0 * (1 + randn()*3e-3)
                vx = 0.0
                vy = 0.0
                vz = 0.0
            end

            # Wypełniamy tablicę grid zgodnie z ustalonym porządkiem:
            grid[i, j, k, 1] = density * (1 + randn() * 3e-3)
            grid[i, j, k, 2] = energy
            grid[i, j, k, 3] = vx   # składowa vx
            grid[i, j, k, 4] = vy   # składowa vy
            grid[i, j, k, 5] = vz   # składowa vz
            grid[i, j, k, 6] = 0.0  # Bx
            grid[i, j, k, 7] = 0.0  # By
            grid[i, j, k, 8] = 0.0  # Bz
        end
    end
end

