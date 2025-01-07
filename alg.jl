using BenchmarkTools
using Base.Threads
using Profile

include("Matrix_operations.jl")
include("Christoffel_Symbols.jl")
include("structs.jl")


r = range(10.0, 1000.0, length = 256)  
theta = range(0.1, 1.0, length = 128)  
phi = range(0.1, 1.0, length = 64)   

x = [0.7, 0.2, 0.3, 0.4, 0.5, 2.0, 3.0, 4.0]
x_guess = [0.7, 0.2, 0.3, 0.4, 0.5, 2.0, 3.0, 4.0] * 0.865

local_buffer = zeros(8)

@time @threads for i in r
    for j in theta
        for k in phi
            metric = Kerr_Schild_metric_cov(i, j, k, 0.95)
            PtoFx(x, local_buffer, metric,eos::Polytrope)
            PtoFy(x, local_buffer, metric,eos::Polytrope)
            PtoFz(x, local_buffer, metric,eos::Polytrope)
            PtoU(x, local_buffer, metric,eos::Polytrope)
            UtoP(local_buffer, x_guess, metric,eos::Polytrope)
            for z1 in 1:4
                for z2 in 1:4
                    for z3 in 1:4
                        Kerr_Schild_Christoffel_Symbols(z1, z2, z3, i, j, k, 0.95)
                    end
                end
            end
        end
    end
end



