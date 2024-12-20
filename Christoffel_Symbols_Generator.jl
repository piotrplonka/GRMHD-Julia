using LinearAlgebra
using Base.Threads
using SymPy

function ChristoffelSymbol(g, xx)
    n = length(xx)
    ig = inv(g)  
    Γ = zeros(Sym, n, n, n)  
    for i in 1:n
        for j in 1:n
            for k in 1:n
                Γ[i, j, k] = (1 // 2) * sum(
                    ig[i, s] * (-diff(g[j, k], xx[s]) + diff(g[j, s], xx[k]) + diff(g[s, k], xx[j]))
                    for s in 1:n
                )
            end
        end
    end
    return simplify.(Γ)
end


function replace_print(expr)
    expr_str = string(expr)
    patterns = [
        ("sin(θ)", "sine_θ"),
        ("sin(2*θ)", "sine_2θ"),
        ("sin(3*θ)", "sine_3θ"),
        ("sin(4*θ)", "sine_4θ"),
        ("cos(θ)", "cos_θ"),
        ("cos(2*θ)", "cos_2θ"),
        ("cos(3*θ)", "cos_3θ"),
        ("cos(4*θ)", "cos_4θ"),
        ("tan(θ)", "tangent_θ")
    ]
    
    for (pattern, replacement) in patterns
        if count(x -> occursin(pattern, x), split(expr_str)) > 1
            expr_str = Base.replace(expr_str, pattern => replacement)
            println("\t\t$replacement = $pattern")
        end
    end
    
    for i in 2:19
        power_patterns = [
            ("r^$i", "r_$i"),
            ("sine_θ^$i", "sine_θ_$i"),
            ("cos_θ^$i", "cos_θ_$i"),
            ("a^$i", "a_$i")
        ]
        
        for (pattern, replacement) in power_patterns
            if count(x -> occursin(pattern, x), split(expr_str)) > 1
                expr_str = Base.replace(expr_str, pattern => replacement)
                println("\t\t$replacement = $pattern")
            end
        end
    end
    return expr_str
end

function replace_print_metric(expr)
    expr_str = string(expr)
    patterns = [
        ("sin(θ)", "sine_θ"),
        ("sin(2*θ)", "sine_2θ"),
        ("sin(3*θ)", "sine_3θ"),
        ("sin(4*θ)", "sine_4θ"),
        ("cos(θ)", "cos_θ"),
        ("cos(2*θ)", "cos_2θ"),
        ("cos(3*θ)", "cos_3θ"),
        ("cos(4*θ)", "cos_4θ"),
        ("tan(θ)", "tangent_θ")
    ]
    
    for (pattern, replacement) in patterns
        if count(x -> occursin(pattern, x), split(expr_str)) > 1
            expr_str = Base.replace(expr_str, pattern => replacement)
            println("\t$replacement = $pattern")
        end
    end
    
    for i in 2:19
        power_patterns = [
            ("r^$i", "r_$i"),
            ("sine_θ^$i", "sine_θ_$i"),
            ("cos_θ^$i", "cos_θ_$i"),
            ("a^$i", "a_$i")
        ]
        
        for (pattern, replacement) in power_patterns
            if count(x -> occursin(pattern, x), split(expr_str)) > 1
                expr_str = Base.replace(expr_str, pattern => replacement)
                println("\t$replacement = $pattern")
            end
        end
    end
    return expr_str
end

function replace(expr)
    
    expr_str = string(expr)
    
    patterns = [
        ("sin(θ)", "sine_θ"),
        ("sin(2*θ)", "sine_2θ"),
        ("sin(3*θ)", "sine_3θ"),
        ("sin(4*θ)", "sine_4θ"),
        ("cos(θ)", "cos_θ"),
        ("cos(2*θ)", "cos_2θ"),
        ("cos(3*θ)", "cos_3θ"),
        ("cos(4*θ)", "cos_4θ"),
        ("tan(θ)", "tangent_θ")
    ]


    for (pattern, replacement) in patterns
        if count(x -> occursin(pattern, x), split(expr_str)) > 1
            expr_str = Base.replace(expr_str, pattern => replacement)
        end
    end


    for i in 2:19
        power_patterns = [
            ("r^$i", "r_$i"),
            ("sine_θ^$i", "sine_θ_$i"),
            ("cos_θ^$i", "cos_θ_$i"),
            ("a^$i", "a_$i")
        ]

        for (pattern, replacement) in power_patterns
            if count(x -> occursin(pattern, x), split(expr_str)) > 1
                expr_str = Base.replace(expr_str, pattern => replacement)
            end
        end
    end

    return expr_str
end

function print_metric_cov(g,name_of_metric,parameters)
    parameter_string = join([p * "::Float64" for p in parameters], ", ")
    println("function $(name_of_metric)_metric_cov($parameter_string)")
    println("	g = zeros(4, 4)") 
    for i in 1:4
    	for j in 1:4
    		replace_print_metric(g[i, j])
    	end
    end
    for i in 1:4
        for j in 1:4
            if g[i, j] != 0
            	println("	g[$i, $j] = ", replace(g[i, j])) 
            end            
        end
    end
    println("	return g")
    println("end")
end

function print_christoffel_symbols(Γ,name_of_metric,parameters)
    parameter_string = join([p * "::Float64" for p in parameters], ", ")
    println("function $(name_of_metric)_Christoffel_Symbols(i::Int64, j::Int64, k::Int64, $parameter_string)")
    println("	Γ = zeros(4, 4, 4)")
    z=0
    for i in 1:4
        for j in 1:4
            for k in 1:4
            	if Γ[i, j, k] != 0
            		if z==0
                		println("	if i == $i && j == $j && k ==$k")
                		replace_print(Γ[i, j, k])
                		println("    		Γ[$i, $j, $k] = ", replace(Γ[i, j, k]))
                		z=1
                	else
                		println("	elseif i == $i && j == $j && k == $k")
                		replace_print(Γ[i, j, k])
                		println("    		Γ[$i, $j, $k] = ", replace(Γ[i, j, k]))	
                	end	
                end
            end
        end
    end
    println("	end")
    println("	return Γ[i, j, k]")
    println("end")
end

function print_christoffel_symbols_prevoius(Γ,name_of_metric,parameters)
    parameter_string = join([p * "::Float64" for p in parameters], ", ")
    println("function $(name_of_metric)_Christoffel_Symbols(i::Int64, j::Int64, k::Int64, $parameter_string)")
    println("	Γ = zeros(4, 4, 4)")
    z=0
    for i in 1:4
        for j in 1:4
            for k in 1:4
            	if Γ[i, j, k] != 0
            		if z==0
                		println("	if i == $i && j == $j && k == $k")
                		println("    		Γ[$i, $j, $k] = ", Γ[i, j, k])
                		z=1
                	end
                	println("	elseif i == $i && j == $j && k == $k")
                	println("    		Γ[$i, $j, $k] = ", Γ[i, j, k])	
                end
            end
        end
    end
    println("	end")
    println("	return Γ[i, j, k]")
    println("end")
end

function Kerr_BL()
    @syms t r θ φ m a
    m=1
    Σ = r^2 + a^2 * cos(θ)^2
    Δ_r = r^2 - 2*m*r + a^2
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] 
    g[1, 1] = 1 - 2*m*r / Σ
    g[1, 4] = 2*m*r*a*sin(θ)^2 / Σ
    g[2, 2] = -Σ / Δ_r
    g[3, 3] = -Σ
    g[4, 4] = -(r^2 + a^2 + 2*m*r*a^2*sin(θ)^2 / Σ) * sin(θ)^2
    g[4, 1] = g[1, 4]

    return simplify.(g)
end

function Kerr_BL_con()
    @syms t r θ φ m a
    m=1
    Σ = r^2 + a^2 * cos(θ)^2
    Δ_r = r^2 - 2*m*r + a^2
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] 
    g[1, 1] = 1 - 2*m*r / Σ
    g[1, 4] = 2*m*r*a*sin(θ)^2 / Σ
    g[2, 2] = -Σ / Δ_r
    g[3, 3] = -Σ
    g[4, 4] = -(r^2 + a^2 + 2*m*r*a^2*sin(θ)^2 / Σ) * sin(θ)^2
    g[4, 1] = g[1, 4]

    return simplify.(inv(g))
end

function Kerr_BL_concov()
    @syms t r θ φ m a
    m=1
    Σ = r^2 + a^2 * cos(θ)^2
    Δ_r = r^2 - 2*m*r + a^2
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] 
    g[1, 1] = 1 - 2*m*r / Σ
    g[1, 4] = 2*m*r*a*sin(θ)^2 / Σ
    g[2, 2] = -Σ / Δ_r
    g[3, 3] = -Σ
    g[4, 4] = -(r^2 + a^2 + 2*m*r*a^2*sin(θ)^2 / Σ) * sin(θ)^2
    g[4, 1] = g[1, 4]
    g_inv=inv(g)
    g_concov=Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] 
    for i in 1:4
    	for j in 1:4
    		g_concov[i,j]=g_inv[1,i]*g[1,j]+g_inv[2,i]*g[2,j]+g_inv[3,i]*g[3,j]+g_inv[4,i]*g[4,j]
    	end
    end	
    return simplify.(g_concov)
end

function Schwarzschild()
    @syms t r θ φ M G c  
    M=1
    c=1
    G=1
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    f = 1 - (2 * G * M / c^2) / r  
    g[1, 1] = simplify(-f * c^2)           # g_tt
    g[2, 2] = simplify(1 / f)              # g_rr
    g[3, 3] = r^2                          # g_θθ
    g[4, 4] = r^2 * sin(θ)^2               # g_φφ
    return simplify.(g)
end

function Schwarzschild_con()
    @syms t r θ φ M G c  
    M=1
    c=1
    G=1
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    f = 1 - (2 * G * M / c^2) / r  
    g[1, 1] = simplify(-f * c^2)           # g_tt
    g[2, 2] = simplify(1 / f)              # g_rr
    g[3, 3] = r^2                          # g_θθ
    g[4, 4] = r^2 * sin(θ)^2               # g_φφ
    return simplify.(inv(g))
end

function Schwarzschild_concov()
    @syms t r θ φ M G c  
    M=1
    c=1
    G=1
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    f = 1 - (2 * G * M / c^2) / r  
    g[1, 1] = simplify(-f * c^2)           # g_tt
    g[2, 2] = simplify(1 / f)              # g_rr
    g[3, 3] = r^2                          # g_θθ
    g[4, 4] = r^2 * sin(θ)^2               # g_φφ
    g_inv=inv(g)
    g_concov=Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] 
    for i in 1:4
    	for j in 1:4
    		g_concov[i,j]=g_inv[1,i]*g[1,j]+g_inv[2,i]*g[2,j]+g_inv[3,i]*g[3,j]+g_inv[4,i]*g[4,j]
    	end
    end	
    return simplify.(g_concov)
end


function Kerr_Schild()
    @syms t r θ φ a M c G
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    M=1
    c=1
    G=1
    g[1, 1] = -(1 - 2*r/(r^2 + a^2 * cos(θ)^2))
    g[1, 2] = 4*r/(r^2 + a^2 * cos(θ)^2)  
    g[2, 2] = (1 + 2*r/(r^2 + a^2 * cos(θ)^2))
    g[3, 3] = r^2 + a^2 * cos(θ)^2 
    g[4, 4] = sin(θ)^2 * ((r^2 + a^2 * cos(θ)^2) + a^2 * (1 + 2*r/(r^2 + a^2 * cos(θ)^2)) * sin(θ)^2)
    g[1, 4] = -4*a*r*sin(θ)^2/(r^2 + a^2 * cos(θ)^2)
    g[2, 4] = -2*a * (1 + 2*r/(r^2 + a^2 * cos(θ)^2)) * sin(θ)^2
    g[2, 1] = g[1, 2]
    g[4, 1] = g[1, 4]
    g[4, 2] = g[2, 4]
    g[3, 1] = g[1, 3] = 0
    g[3, 2] = g[2, 3] = 0
    g[3, 4] = g[4, 3] = 0
    return simplify.(g)
end

function Kerr_Schild_con()
    @syms t r θ φ a M c G
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    M=1
    c=1
    G=1
    g[1, 1] = -(1 - 2*r/(r^2 + a^2 * cos(θ)^2))
    g[1, 2] = 4*r/(r^2 + a^2 * cos(θ)^2)  
    g[2, 2] = (1 + 2*r/(r^2 + a^2 * cos(θ)^2))
    g[3, 3] = r^2 + a^2 * cos(θ)^2 
    g[4, 4] = sin(θ)^2 * ((r^2 + a^2 * cos(θ)^2) + a^2 * (1 + 2*r/(r^2 + a^2 * cos(θ)^2)) * sin(θ)^2)
    g[1, 4] = -4*a*r*sin(θ)^2/(r^2 + a^2 * cos(θ)^2)
    g[2, 4] = -2*a * (1 + 2*r/(r^2 + a^2 * cos(θ)^2)) * sin(θ)^2
    g[2, 1] = g[1, 2]
    g[4, 1] = g[1, 4]
    g[4, 2] = g[2, 4]
    g[3, 1] = g[1, 3] = 0
    g[3, 2] = g[2, 3] = 0
    g[3, 4] = g[4, 3] = 0
    return simplify.(inv(g))
end

function Kerr_Schild_concov()
    @syms t r θ φ a M c G
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    M=1
    c=1
    G=1
    g[1, 1] = -(1 - 2*r/(r^2 + a^2 * cos(θ)^2))
    g[1, 2] = 4*r/(r^2 + a^2 * cos(θ)^2)  
    g[2, 2] = (1 + 2*r/(r^2 + a^2 * cos(θ)^2))
    g[3, 3] = r^2 + a^2 * cos(θ)^2 
    g[4, 4] = sin(θ)^2 * ((r^2 + a^2 * cos(θ)^2) + a^2 * (1 + 2*r/(r^2 + a^2 * cos(θ)^2)) * sin(θ)^2)
    g[1, 4] = -4*a*r*sin(θ)^2/(r^2 + a^2 * cos(θ)^2)
    g[2, 4] = -2*a * (1 + 2*r/(r^2 + a^2 * cos(θ)^2)) * sin(θ)^2
    g[2, 1] = g[1, 2]
    g[4, 1] = g[1, 4]
    g[4, 2] = g[2, 4]
    g[3, 1] = g[1, 3] = 0
    g[3, 2] = g[2, 3] = 0
    g[3, 4] = g[4, 3] = 0
    g_inv=inv(g)
    g_concov=Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] 
    for i in 1:4
    	for j in 1:4
    		g_concov[i,j]=g_inv[1,i]*g[1,j]+g_inv[2,i]*g[2,j]+g_inv[3,i]*g[3,j]+g_inv[4,i]*g[4,j]
    	end
    end	
    return simplify.(g_concov)
end

function Kerr()
    @syms t r θ φ M G c J a
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    J=a*M*c
    M=1
    c=1
    G=1
    g[1, 1] = simplify((1 - (2 * G * M / c^2) * r / (r^2 + (J / (M * c))^2 * cos(θ)^2)) * c^2)
    g[1, 4] = simplify(r * (2 * G * M / c^2) * (J / (M * c)) * sin(θ)^2 * c / (r^2 + (J / (M * c))^2 * cos(θ)^2))
    g[2, 2] = simplify(-(r^2 + (J / (M * c))^2 * cos(θ)^2) / (r^2 - (2 * G * M / c^2) * r + (J / (M * c))^2))
    g[3, 3] = simplify(-(r^2 + (J / (M * c))^2 * cos(θ)^2))
    g[4, 4] = simplify(-sin(θ)^2 * (r^2 + (J / (M * c))^2 + sin(θ)^2 * (2 * G * M / c^2) * (J / (M * c))^2 * r / (r^2 + (J / (M * c))^2 * cos(θ)^2)))
    g[4, 1] = g[1, 4]
    g[4, 2] = 0
    g[4, 3] = 0
    return simplify.(g)
end

function Kerr_con()
    @syms t r θ φ M G c J a
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    J=a*M*c
    M=1
    c=1
    G=1
    g[1, 1] = simplify((1 - (2 * G * M / c^2) * r / (r^2 + (J / (M * c))^2 * cos(θ)^2)) * c^2)
    g[1, 4] = simplify(r * (2 * G * M / c^2) * (J / (M * c)) * sin(θ)^2 * c / (r^2 + (J / (M * c))^2 * cos(θ)^2))
    g[2, 2] = simplify(-(r^2 + (J / (M * c))^2 * cos(θ)^2) / (r^2 - (2 * G * M / c^2) * r + (J / (M * c))^2))
    g[3, 3] = simplify(-(r^2 + (J / (M * c))^2 * cos(θ)^2))
    g[4, 4] = simplify(-sin(θ)^2 * (r^2 + (J / (M * c))^2 + sin(θ)^2 * (2 * G * M / c^2) * (J / (M * c))^2 * r / (r^2 + (J / (M * c))^2 * cos(θ)^2)))
    g[4, 1] = g[1, 4]
    g[4, 2] = 0
    g[4, 3] = 0
    return simplify.(inv(g))
end

function Kerr_concov()
    @syms t r θ φ M G c J a
    g = Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    J=a*M*c
    M=1
    c=1
    G=1
    g[1, 1] = simplify((1 - (2 * G * M / c^2) * r / (r^2 + (J / (M * c))^2 * cos(θ)^2)) * c^2)
    g[1, 4] = simplify(r * (2 * G * M / c^2) * (J / (M * c)) * sin(θ)^2 * c / (r^2 + (J / (M * c))^2 * cos(θ)^2))
    g[2, 2] = simplify(-(r^2 + (J / (M * c))^2 * cos(θ)^2) / (r^2 - (2 * G * M / c^2) * r + (J / (M * c))^2))
    g[3, 3] = simplify(-(r^2 + (J / (M * c))^2 * cos(θ)^2))
    g[4, 4] = simplify(-sin(θ)^2 * (r^2 + (J / (M * c))^2 + sin(θ)^2 * (2 * G * M / c^2) * (J / (M * c))^2 * r / (r^2 + (J / (M * c))^2 * cos(θ)^2)))
    g[4, 1] = g[1, 4]
    g[4, 2] = 0
    g[4, 3] = 0
    g_inv=inv(g)
    g_concov=Sym[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] 
    for i in 1:4
    	for j in 1:4
    		g_concov[i,j]=g_inv[1,i]*g[1,j]+g_inv[2,i]*g[2,j]+g_inv[3,i]*g[3,j]+g_inv[4,i]*g[4,j]
    	end
    end	
    return simplify.(g_concov)
end


function print_metric_con(g,name_of_metric,parameters)
    parameter_string = join([p * "::Float64" for p in parameters], ", ")
    println("function $(name_of_metric)_metric_con($parameter_string)")
    println("	g = zeros(4, 4)") 
    for i in 1:4
        for j in 1:4
            if g[i, j] != 0
            	println("	g[$i, $j] = ", g[i, j]) 
            end            
        end
    end
    println("	return g")
    println("end")
end

function print_metric_concov(g,name_of_metric,parameters)
    parameter_string = join([p * "::Float64" for p in parameters], ", ")
    println("function $(name_of_metric)_metric_concov($parameter_string)")
    println("	g = zeros(4, 4)") 
    for i in 1:4
        for j in 1:4
            if g[i, j] != 0
            	println("	g[$i, $j] = ", g[i, j]) 
            end            
        end
    end
    println("	return g")
    println("end")
end

parameters = ["r", "θ", "φ"]
print_metric_cov(Schwarzschild(),"Schwarzschild",parameters)
println("	")


parameters = ["r", "θ", "φ", "a"]
print_metric_cov(Kerr_Schild(),"Kerr_Schild",parameters)
println("	")


parameters = ["r", "θ", "φ", "a"]
print_metric_cov(Kerr(),"Kerr",parameters)
println("	")


parameters = ["r", "θ", "φ", "a"]
print_metric_cov(Kerr_BL(),"Kerr_BL",parameters)
println("	")

parameters = ["r", "θ", "φ"]
print_metric_con(Schwarzschild_con(),"Schwarzschild",parameters)
println("	")


parameters = ["r", "θ", "φ", "a"]
print_metric_con(Kerr_Schild_con(),"Kerr_Schild",parameters)
println("	")


parameters = ["r", "θ", "φ", "a"]
print_metric_con(Kerr_con(),"Kerr",parameters)
println("	")


parameters = ["r", "θ", "φ", "a"]
print_metric_con(Kerr_BL_con(),"Kerr_BL",parameters)
println("	")

parameters = ["r", "θ", "φ"]
print_metric_concov(Schwarzschild_concov(),"Schwarzschild",parameters)
println("	")


parameters = ["r", "θ", "φ", "a"]
print_metric_concov(Kerr_Schild_concov(),"Kerr_Schild",parameters)
println("	")


parameters = ["r", "θ", "φ", "a"]
print_metric_concov(Kerr_concov(),"Kerr",parameters)
println("	")


parameters = ["r", "θ", "φ", "a"]
print_metric_concov(Kerr_BL_concov(),"Kerr_BL",parameters)
println("	")

g_kerr = Schwarzschild()
@syms t r θ φ a
Γ = ChristoffelSymbol(g_kerr, [t, r, θ, φ])
parameters = ["r", "θ", "φ"]#, "a","m"]
print_christoffel_symbols(Γ,"Schwarzschild",parameters)
println("	")



g_kerr = Kerr_Schild()
@syms t r θ φ a
Γ = ChristoffelSymbol(g_kerr, [t, r, θ, φ])
parameters = ["r", "θ", "φ", "a"]
print_christoffel_symbols(Γ,"Kerr_Schild",parameters)
println("	")



g_kerr = Kerr_BL()
@syms t r θ φ a
Γ = ChristoffelSymbol(g_kerr, [t, r, θ, φ])
parameters = ["r", "θ", "φ", "a"]
print_christoffel_symbols(Γ,"Kerr_BL",parameters)
println("	")



g_kerr = Kerr()
@syms t r θ φ a M
Γ = ChristoffelSymbol(g_kerr, [t, r, θ, φ])
parameters = ["r", "θ", "φ", "a"]
print_christoffel_symbols(Γ,"Kerr",parameters)
