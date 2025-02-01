using BenchmarkTools
using Base.Threads
using StaticArrays

include("Matrix_operations.jl")
include("Christoffel_Symbols.jl")
include("eos.jl")

function SourceGR(x::AbstractVector, gcov::Matrix{Float64}, radius::Float64, theta::Float64, phi::Float64, a_kerr::Float64, eos::Polytrope)
	
	#Parameters
	ρ::Float64  = x[1] #Density
	u::Float64  = x[2] #Internal Energy 
	u1::Float64 = x[3] #Contravariant Four-velocity in 1-direction
	u2::Float64 = x[4] #Contravariant Four-velocity in 2-direction
	u3::Float64 = x[5] #Contravariant Four-velocity in 3-direction   
	B1::Float64 = x[6] #Magnetic field in 1-direction
	B2::Float64 = x[7] #Magnetic field in 2-direction
	B3::Float64 = x[8] #Magnetic field in 3-direction

	#To find u0, we need to solve the quadratic equation g_ij*ui*uj = −1    
	a::Float64 = gcov[1,1]   
	b::Float64 = 2*gcov[1,2]*u1 + 2*gcov[1,3]*u2 + 2*gcov[1,4]*u3
	c::Float64 = gcov[2,2]*u1*u1 + 2*gcov[2,3]*u1*u2 + 2*gcov[2,4]*u1*u3 + gcov[3,3]*u2*u2 + gcov[4,4]*u3*u3 + 2*gcov[3,4]*u2*u3 + 1

	#Contravariant Four-velocity
	u0::Float64 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
       #u1::Float64 = x[3]
       #u2::Float64 = x[4]
       #u3::Float64 = x[5]

	#Covariant Four-velocity
	u_0::Float64 = u0*gcov[1,1] + u1*gcov[2,1] + u2*gcov[3,1] + u3*gcov[4,1]
	u_1::Float64 = u0*gcov[1,2] + u1*gcov[2,2] + u2*gcov[3,2] + u3*gcov[4,2]
	u_2::Float64 = u0*gcov[1,3] + u1*gcov[2,3] + u2*gcov[3,3] + u3*gcov[4,3]
	u_3::Float64 = u0*gcov[1,4] + u1*gcov[2,4] + u2*gcov[3,4] + u3*gcov[4,4]

	#Contravariant Four-magnetic field
	b0::Float64 = B1*u_1 + B2*u_2 + B3*u_3
	b1::Float64 = (B1 + b0*u1)/u0
	b2::Float64 = (B2 + b0*u2)/u0
	b3::Float64 = (B3 + b0*u3)/u0

	#Covariant Four-magnetic field    
	b_0::Float64 = b0*gcov[1,1] + b1*gcov[1,2] + b2*gcov[1,3] + b3*gcov[1,4]
	b_1::Float64 = b0*gcov[2,1] + b1*gcov[2,2] + b2*gcov[2,3] + b3*gcov[2,4]
	b_2::Float64 = b0*gcov[3,1] + b1*gcov[3,2] + b2*gcov[3,3] + b3*gcov[3,4]
	b_3::Float64 = b0*gcov[4,1] + b1*gcov[4,2] + b2*gcov[4,3] + b3*gcov[4,4]

	#Useful Values        
	bsq::Float64    = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3
	value::Float64  = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq # p + (½)*bsq
	sq_g::Float64   = sqrt_g(gcov) #square root of the determinant of the metric

	#Stress-Energy tensor Components
	SE_tensor::Matrix{Float64} = zeros(4,4)
	
	SE_tensor[1,1] = sq_g*(value*u0*u_0 + value2 - b1*b_0)
	SE_tensor[1,2] = sq_g*(value*u0*u_1 - b1*b_1)
	SE_tensor[1,3] = sq_g*(value*u0*u_2 - b1*b_2)
	SE_tensor[1,4] = sq_g*(value*u0*u_3 - b1*b_3)

	SE_tensor[2,1] = sq_g*(value*u1*u_0 - b1*b_0)
	SE_tensor[2,2] = sq_g*(value*u1*u_1 + value2 - b1*b_1)
	SE_tensor[2,3] = sq_g*(value*u1*u_2 - b1*b_2)
	SE_tensor[2,4] = sq_g*(value*u1*u_3 - b1*b_3)

	SE_tensor[3,1] = sq_g*(value*u2*u_0 - b2*b_0)
	SE_tensor[3,2] = sq_g*(value*u2*u_1 - b2*b_1)
	SE_tensor[3,3] = sq_g*(value*u2*u_2 + value2 - b2*b_2)
	SE_tensor[3,4] = sq_g*(value*u2*u_3 - b2*b_3)

	SE_tensor[4,1] = sq_g*(value*u3*u_0 - b3*b_0)
	SE_tensor[4,2] = sq_g*(value*u3*u_1 - b3*b_1)
	SE_tensor[4,3] = sq_g*(value*u3*u_2 - b3*b_2)
	SE_tensor[4,4] = sq_g*(value*u3*u_3 + value2 - b3*b_3)
	
	#Flux computation
	S_0::Float64 = 0
	S_1::Float64 = 0
	S_2::Float64 = 0
	S_3::Float64 = 0
	
	for ϰ in 1:4 
		for λ in 1:4
			S_0 = S_0 + SE_tensor[ϰ,λ]*Kerr_Schild_Christoffel_Symbols(λ, 1, ϰ,  radius, theta, phi, a_kerr)
		end
	end
	
	
	for ϰ in 1:4
		for λ in 1:4
			S_1 = S_1 + SE_tensor[ϰ,λ]*Kerr_Schild_Christoffel_Symbols(λ, 2, ϰ,  radius, theta, phi, a_kerr)
		end
	end 

	for ϰ in 1:4
		for λ in 1:4
			S_2 = S_2 + SE_tensor[ϰ,λ]*Kerr_Schild_Christoffel_Symbols(λ, 3, ϰ,  radius, theta, phi, a_kerr)
		end
	end
	


	for ϰ in 1:4
		for λ in 1:4
			S_3 = S_3 + SE_tensor[ϰ,λ]*Kerr_Schild_Christoffel_Symbols(λ, 4, ϰ,  radius, theta, phi, a_kerr)
		end
	end
	
	return S_0, S_1, S_2, S_3
end


