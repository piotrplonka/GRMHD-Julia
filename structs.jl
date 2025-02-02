using BenchmarkTools
using Base.Threads
using StaticArrays
using LinearAlgebra
include("Matrix_operations.jl")
include("Christoffel_Symbols.jl")
include("eos.jl")
include("LU.jl")

function PtoU(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
	#Parameters
	ρ ::Float64 = x[1] #Density
	u ::Float64 = x[2] #Internal Energy 
	u¹::Float64 = x[3] #Contravariant Four-velocity in 1-direction
	u²::Float64 = x[4] #Contravariant Four-velocity in 2-direction
	u³::Float64 = x[5] #Contravariant Four-velocity in 3-direction   
	B¹::Float64 = x[6] #Magnetic field in 1-direction
	B²::Float64 = x[7] #Magnetic field in 2-direction
	B³::Float64 = x[8] #Magnetic field in 3-direction

	#To find u⁰, we need to solve the quadratic equation g_ij*ui*uj = −1    
	a::Float64 = gcov[1,1]   
	b::Float64 = 2*gcov[1,2]*u¹ + 2*gcov[1,3]*u² + 2*gcov[1,4]*u³
	c::Float64 = gcov[2,2]*u¹*u¹ + 2*gcov[2,3]*u¹*u² + 2*gcov[2,4]*u¹*u³ + gcov[3,3]*u²*u² + gcov[4,4]*u³*u³ + 2*gcov[3,4]*u²*u³ + 1

	#Contravariant Four-velocity
	u⁰::Float64 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
       #u¹::Float64 = x[3]
       #u²::Float64 = x[4]
       #u³::Float64 = x[5]

	#Covariant Four-velocity
	u₀::Float64 = u⁰*gcov[1,1] + u¹*gcov[2,1] + u²*gcov[3,1] + u³*gcov[4,1]
	u₁::Float64 = u⁰*gcov[1,2] + u¹*gcov[2,2] + u²*gcov[3,2] + u³*gcov[4,2]
	u₂::Float64 = u⁰*gcov[1,3] + u¹*gcov[2,3] + u²*gcov[3,3] + u³*gcov[4,3]
	u₃::Float64 = u⁰*gcov[1,4] + u¹*gcov[2,4] + u²*gcov[3,4] + u³*gcov[4,4]

	#Contravariant Four-magnetic field
	b⁰::Float64 = B¹*u₁ + B²*u₂ + B³*u₃
	b¹::Float64 = (B¹ + b⁰*u¹)/u⁰
	b²::Float64 = (B² + b⁰*u²)/u⁰
	b³::Float64 = (B³ + b⁰*u³)/u⁰

	#Covariant Four-magnetic field    
	b₀::Float64 = b⁰*gcov[1,1] + b¹*gcov[1,2] + b²*gcov[1,3] + b³*gcov[1,4]
	b₁::Float64 = b⁰*gcov[2,1] + b¹*gcov[2,2] + b²*gcov[2,3] + b³*gcov[2,4]
	b₂::Float64 = b⁰*gcov[3,1] + b¹*gcov[3,2] + b²*gcov[3,3] + b³*gcov[3,4]
	b₃::Float64 = b⁰*gcov[4,1] + b¹*gcov[4,2] + b²*gcov[4,3] + b³*gcov[4,4]

	#Useful Values        
	bsq   ::Float64 = b⁰*b₀ + b¹*b₁ + b²*b₂ + b³*b₃
	value ::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq   # p + (1/2)*bsq
	sq_g  ::Float64 = sqrt_g(gcov)                    #square root of the determinant of the metric

	#Buffers
	buffer[1] = sq_g*ρ*u⁰ 
	buffer[2] = sq_g*((value)*u⁰*u₀ + value2 - b⁰*b₀)
	buffer[3] = sq_g*((value)*u⁰*u₁ - b⁰*b₁)
	buffer[4] = sq_g*((value)*u⁰*u₂ - b⁰*b₂)
	buffer[5] = sq_g*((value)*u⁰*u₃ - b⁰*b₃)
	buffer[6] = sq_g*B¹
	buffer[7] = sq_g*B²
	buffer[8] = sq_g*B³

end

function PtoFx(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
	#Parameters
	ρ ::Float64 = x[1] #Density
	u ::Float64 = x[2] #Internal Energy 
	u¹::Float64 = x[3] #Contravariant Four-velocity in 1-direction
	u²::Float64 = x[4] #Contravariant Four-velocity in 2-direction
	u³::Float64 = x[5] #Contravariant Four-velocity in 3-direction   
	B¹::Float64 = x[6] #Magnetic field in 1-direction
	B²::Float64 = x[7] #Magnetic field in 2-direction
	B³::Float64 = x[8] #Magnetic field in 3-direction

	#To find u⁰, we need to solve the quadratic equation g_ij*ui*uj = −1    
	a::Float64 = gcov[1,1]   
	b::Float64 = 2*gcov[1,2]*u¹ + 2*gcov[1,3]*u² + 2*gcov[1,4]*u³
	c::Float64 = gcov[2,2]*u¹*u¹ + 2*gcov[2,3]*u¹*u² + 2*gcov[2,4]*u¹*u³ + gcov[3,3]*u²*u² + gcov[4,4]*u³*u³ + 2*gcov[3,4]*u²*u³ + 1

	#Contravariant Four-velocity
	u⁰::Float64 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
       #u¹::Float64 = x[3]
       #u²::Float64 = x[4]
       #u³::Float64 = x[5]

	#Covariant Four-velocity
	u₀::Float64 = u⁰*gcov[1,1] + u¹*gcov[2,1] + u²*gcov[3,1] + u³*gcov[4,1]
	u₁::Float64 = u⁰*gcov[1,2] + u¹*gcov[2,2] + u²*gcov[3,2] + u³*gcov[4,2]
	u₂::Float64 = u⁰*gcov[1,3] + u¹*gcov[2,3] + u²*gcov[3,3] + u³*gcov[4,3]
	u₃::Float64 = u⁰*gcov[1,4] + u¹*gcov[2,4] + u²*gcov[3,4] + u³*gcov[4,4]

	#Contravariant Four-magnetic field
	b⁰::Float64 = B¹*u₁ + B²*u₂ + B³*u₃
	b¹::Float64 = (B¹ + b⁰*u¹)/u⁰
	b²::Float64 = (B² + b⁰*u²)/u⁰
	b³::Float64 = (B³ + b⁰*u³)/u⁰

	#Covariant Four-magnetic field    
	b₀::Float64 = b⁰*gcov[1,1] + b¹*gcov[1,2] + b²*gcov[1,3] + b³*gcov[1,4]
	b₁::Float64 = b⁰*gcov[2,1] + b¹*gcov[2,2] + b²*gcov[2,3] + b³*gcov[2,4]
	b₂::Float64 = b⁰*gcov[3,1] + b¹*gcov[3,2] + b²*gcov[3,3] + b³*gcov[3,4]
	b₃::Float64 = b⁰*gcov[4,1] + b¹*gcov[4,2] + b²*gcov[4,3] + b³*gcov[4,4]

	#Useful Values        
	bsq   ::Float64 = b⁰*b₀ + b¹*b₁ + b²*b₂ + b³*b₃
	value ::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq   # p + (1/2)*bsq
	sq_g  ::Float64 = sqrt_g(gcov)                    #square root of the determinant of the metric
	
	#Buffers    
	buffer[1] = sq_g*ρ*u¹
	buffer[2] = sq_g*(value*u¹*u₀ - b¹*b₀)
	buffer[3] = sq_g*(value*u¹*u₁ + value2 - b¹*b₁)
	buffer[4] = sq_g*(value*u¹*u₂ - b¹*b₂)
	buffer[5] = sq_g*(value*u¹*u₃ - b¹*b₃)

end

function PtoFy(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
	#Parameters
	ρ ::Float64 = x[1] #Density
	u ::Float64 = x[2] #Internal Energy 
	u¹::Float64 = x[3] #Contravariant Four-velocity in 1-direction
	u²::Float64 = x[4] #Contravariant Four-velocity in 2-direction
	u³::Float64 = x[5] #Contravariant Four-velocity in 3-direction   
	B¹::Float64 = x[6] #Magnetic field in 1-direction
	B²::Float64 = x[7] #Magnetic field in 2-direction
	B³::Float64 = x[8] #Magnetic field in 3-direction

	#To find u⁰, we need to solve the quadratic equation g_ij*ui*uj = −1    
	a::Float64 = gcov[1,1]   
	b::Float64 = 2*gcov[1,2]*u¹ + 2*gcov[1,3]*u² + 2*gcov[1,4]*u³
	c::Float64 = gcov[2,2]*u¹*u¹ + 2*gcov[2,3]*u¹*u² + 2*gcov[2,4]*u¹*u³ + gcov[3,3]*u²*u² + gcov[4,4]*u³*u³ + 2*gcov[3,4]*u²*u³ + 1

	#Contravariant Four-velocity
	u⁰::Float64 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
       #u¹::Float64 = x[3]
       #u²::Float64 = x[4]
       #u³::Float64 = x[5]

	#Covariant Four-velocity
	u₀::Float64 = u⁰*gcov[1,1] + u¹*gcov[2,1] + u²*gcov[3,1] + u³*gcov[4,1]
	u₁::Float64 = u⁰*gcov[1,2] + u¹*gcov[2,2] + u²*gcov[3,2] + u³*gcov[4,2]
	u₂::Float64 = u⁰*gcov[1,3] + u¹*gcov[2,3] + u²*gcov[3,3] + u³*gcov[4,3]
	u₃::Float64 = u⁰*gcov[1,4] + u¹*gcov[2,4] + u²*gcov[3,4] + u³*gcov[4,4]

	#Contravariant Four-magnetic field
	b⁰::Float64 = B¹*u₁ + B²*u₂ + B³*u₃
	b¹::Float64 = (B¹ + b⁰*u¹)/u⁰
	b²::Float64 = (B² + b⁰*u²)/u⁰
	b³::Float64 = (B³ + b⁰*u³)/u⁰

	#Covariant Four-magnetic field    
	b₀::Float64 = b⁰*gcov[1,1] + b¹*gcov[1,2] + b²*gcov[1,3] + b³*gcov[1,4]
	b₁::Float64 = b⁰*gcov[2,1] + b¹*gcov[2,2] + b²*gcov[2,3] + b³*gcov[2,4]
	b₂::Float64 = b⁰*gcov[3,1] + b¹*gcov[3,2] + b²*gcov[3,3] + b³*gcov[3,4]
	b₃::Float64 = b⁰*gcov[4,1] + b¹*gcov[4,2] + b²*gcov[4,3] + b³*gcov[4,4]

	#Useful Values        
	bsq   ::Float64 = b⁰*b₀ + b¹*b₁ + b²*b₂ + b³*b₃
	value ::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq   # p + (1/2)*bsq
	sq_g  ::Float64 = sqrt_g(gcov)                    #square root of the determinant of the metric
	
	#Buffers
	buffer[1] = sq_g*ρ*u²
	buffer[2] = sq_g*(value*u²*u₀ - b²*b₀)
	buffer[3] = sq_g*(value*u²*u₁ - b²*b₁)
	buffer[4] = sq_g*(value*u²*u₂ + value2 - b²*b₂)
	buffer[5] = sq_g*(value*u²*u₃ - b²*b₃)

end

function PtoFz(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
	#Parameters
	ρ ::Float64 = x[1] #Density
	u ::Float64 = x[2] #Internal Energy 
	u¹::Float64 = x[3] #Contravariant Four-velocity in 1-direction
	u²::Float64 = x[4] #Contravariant Four-velocity in 2-direction
	u³::Float64 = x[5] #Contravariant Four-velocity in 3-direction   
	B¹::Float64 = x[6] #Magnetic field in 1-direction
	B²::Float64 = x[7] #Magnetic field in 2-direction
	B³::Float64 = x[8] #Magnetic field in 3-direction

	#To find u⁰, we need to solve the quadratic equation g_ij*ui*uj = −1    
	a::Float64 = gcov[1,1]   
	b::Float64 = 2*gcov[1,2]*u¹ + 2*gcov[1,3]*u² + 2*gcov[1,4]*u³
	c::Float64 = gcov[2,2]*u¹*u¹ + 2*gcov[2,3]*u¹*u² + 2*gcov[2,4]*u¹*u³ + gcov[3,3]*u²*u² + gcov[4,4]*u³*u³ + 2*gcov[3,4]*u²*u³ + 1

	#Contravariant Four-velocity
	u⁰::Float64 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
       #u¹::Float64 = x[3]
       #u²::Float64 = x[4]
       #u³::Float64 = x[5]

	#Covariant Four-velocity
	u₀::Float64 = u⁰*gcov[1,1] + u¹*gcov[2,1] + u²*gcov[3,1] + u³*gcov[4,1]
	u₁::Float64 = u⁰*gcov[1,2] + u¹*gcov[2,2] + u²*gcov[3,2] + u³*gcov[4,2]
	u₂::Float64 = u⁰*gcov[1,3] + u¹*gcov[2,3] + u²*gcov[3,3] + u³*gcov[4,3]
	u₃::Float64 = u⁰*gcov[1,4] + u¹*gcov[2,4] + u²*gcov[3,4] + u³*gcov[4,4]

	#Contravariant Four-magnetic field
	b⁰::Float64 = B¹*u₁ + B²*u₂ + B³*u₃
	b¹::Float64 = (B¹ + b⁰*u¹)/u⁰
	b²::Float64 = (B² + b⁰*u²)/u⁰
	b³::Float64 = (B³ + b⁰*u³)/u⁰

	#Covariant Four-magnetic field    
	b₀::Float64 = b⁰*gcov[1,1] + b¹*gcov[1,2] + b²*gcov[1,3] + b³*gcov[1,4]
	b₁::Float64 = b⁰*gcov[2,1] + b¹*gcov[2,2] + b²*gcov[2,3] + b³*gcov[2,4]
	b₂::Float64 = b⁰*gcov[3,1] + b¹*gcov[3,2] + b²*gcov[3,3] + b³*gcov[3,4]
	b₃::Float64 = b⁰*gcov[4,1] + b¹*gcov[4,2] + b²*gcov[4,3] + b³*gcov[4,4]

	#Useful Values        
	bsq   ::Float64 = b⁰*b₀ + b¹*b₁ + b²*b₂ + b³*b₃
	value ::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq   # p + (1/2)*bsq
	sq_g  ::Float64 = sqrt_g(gcov)                    #square root of the determinant of the metric
	
	#Buffers
	buffer[1] = sq_g*ρ*u³
	buffer[2] = sq_g*(value*u³*u₀ - b³*b₀)
	buffer[3] = sq_g*(value*u³*u₁ - b³*b₁)
	buffer[4] = sq_g*(value*u³*u₂ - b³*b₂)
	buffer[5] = sq_g*(value*u³*u₃ + value2 - b³*b₃)

end

function PtoF_Bx_By_Bz(x::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
	#Parameters
	ρ ::Float64 = x[1] #Density
	u ::Float64 = x[2] #Internal Energy 
	u¹::Float64 = x[3] #Contravariant Four-velocity in 1-direction
	u²::Float64 = x[4] #Contravariant Four-velocity in 2-direction
	u³::Float64 = x[5] #Contravariant Four-velocity in 3-direction   
	B¹::Float64 = x[6] #Magnetic field in 1-direction
	B²::Float64 = x[7] #Magnetic field in 2-direction
	B³::Float64 = x[8] #Magnetic field in 3-direction

	#To find u⁰, we need to solve the quadratic equation g_ij*ui*uj = −1    
	a::Float64 = gcov[1,1]   
	b::Float64 = 2*gcov[1,2]*u¹ + 2*gcov[1,3]*u² + 2*gcov[1,4]*u³
	c::Float64 = gcov[2,2]*u¹*u¹ + 2*gcov[2,3]*u¹*u² + 2*gcov[2,4]*u¹*u³ + gcov[3,3]*u²*u² + gcov[4,4]*u³*u³ + 2*gcov[3,4]*u²*u³ + 1

	#Contravariant Four-velocity
	u⁰::Float64 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
       #u¹::Float64 = x[3]
       #u²::Float64 = x[4]
       #u³::Float64 = x[5]

	#Covariant Four-velocity
	u₀::Float64 = u⁰*gcov[1,1] + u¹*gcov[2,1] + u²*gcov[3,1] + u³*gcov[4,1]
	u₁::Float64 = u⁰*gcov[1,2] + u¹*gcov[2,2] + u²*gcov[3,2] + u³*gcov[4,2]
	u₂::Float64 = u⁰*gcov[1,3] + u¹*gcov[2,3] + u²*gcov[3,3] + u³*gcov[4,3]
	u₃::Float64 = u⁰*gcov[1,4] + u¹*gcov[2,4] + u²*gcov[3,4] + u³*gcov[4,4]

	#Contravariant Four-magnetic field
	b⁰::Float64 = B¹*u₁ + B²*u₂ + B³*u₃
	b¹::Float64 = (B¹ + b⁰*u¹)/u⁰
	b²::Float64 = (B² + b⁰*u²)/u⁰
	b³::Float64 = (B³ + b⁰*u³)/u⁰

	#Covariant Four-magnetic field    
	b₀::Float64 = b⁰*gcov[1,1] + b¹*gcov[1,2] + b²*gcov[1,3] + b³*gcov[1,4]
	b₁::Float64 = b⁰*gcov[2,1] + b¹*gcov[2,2] + b²*gcov[2,3] + b³*gcov[2,4]
	b₂::Float64 = b⁰*gcov[3,1] + b¹*gcov[3,2] + b²*gcov[3,3] + b³*gcov[3,4]
	b₃::Float64 = b⁰*gcov[4,1] + b¹*gcov[4,2] + b²*gcov[4,3] + b³*gcov[4,4]

	#Useful Values        
	bsq   ::Float64 = b⁰*b₀ + b¹*b₁ + b²*b₂ + b³*b₃
	value ::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq   # p + (1/2)*bsq
	sq_g  ::Float64 = sqrt_g(gcov)                    #square root of the determinant of the metric
	
	#Buffers
	buffer_1::Float64 = sq_g*(b⁰*u¹-b¹*u⁰ + b²*u¹-b¹*u² + b³*u¹-b¹*u³)
	buffer_2::Float64 = sq_g*(b⁰*u²-b²*u⁰ + b¹*u²-b²*u¹ + b³*u²-b²*u³)
	buffer_3::Float64 = sq_g*(b⁰*u³-b³*u⁰ + b¹*u³-b³*u¹ + b²*u³-b³*u²)
	return buffer_1, buffer_2, buffer_3
end

function Lorentz_factor(x::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
	#Parameters
	ρ::Float64  = x[1] #Density
	u::Float64  = x[2] #Internal Energy 
	u¹::Float64 = x[3] #Contravariant Four-velocity in 1-direction
	u²::Float64 = x[4] #Contravariant Four-velocity in 2-direction
	u³::Float64 = x[5] #Contravariant Four-velocity in 3-direction   
	B¹::Float64 = x[6] #Magnetic field in 1-direction
	B²::Float64 = x[7] #Magnetic field in 2-direction
	B³::Float64 = x[8] #Magnetic field in 3-direction

	#To find u⁰, we need to solve the quadratic equation g_ij*ui*uj = −1    
	a::Float64 = gcov[1,1]   
	b::Float64 = 2*gcov[1,2]*u¹ + 2*gcov[1,3]*u² + 2*gcov[1,4]*u³
	c::Float64 = gcov[2,2]*u¹*u¹ + 2*gcov[2,3]*u¹*u² + 2*gcov[2,4]*u¹*u³ + gcov[3,3]*u²*u² + gcov[4,4]*u³*u³ + 2*gcov[3,4]*u²*u³ + 1

	#For C=1 lorentz factor is:
	γ::Float64 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
	return γ
end


function Jacobian(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
	
	#Parameters
	ρ ::Float64 = x[1] #Density
	u ::Float64 = x[2] #Internal Energy 
	u¹::Float64 = x[3] #Contravariant Four-velocity in 1-direction
	u²::Float64 = x[4] #Contravariant Four-velocity in 2-direction
	u³::Float64 = x[5] #Contravariant Four-velocity in 3-direction   
	B¹::Float64 = x[6] #Magnetic field in 1-direction
	B²::Float64 = x[7] #Magnetic field in 2-direction
	B³::Float64 = x[8] #Magnetic field in 3-direction

	#To find u⁰, we need to solve the quadratic equation g_ij*ui*uj = −1    
	a::Float64 = gcov[1,1]   
	b::Float64 = 2*gcov[1,2]*u¹ + 2*gcov[1,3]*u² + 2*gcov[1,4]*u³
	c::Float64 = gcov[2,2]*u¹*u¹ + 2*gcov[2,3]*u¹*u² + 2*gcov[2,4]*u¹*u³ + gcov[3,3]*u²*u² + gcov[4,4]*u³*u³ + 2*gcov[3,4]*u²*u³ + 1

	#Contravariant Four-velocity
	u⁰::Float64 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
       #u¹::Float64 = x[3]
       #u²::Float64 = x[4]
       #u³::Float64 = x[5]

	#Covariant Four-velocity
	u₀::Float64 = u⁰*gcov[1,1] + u¹*gcov[2,1] + u²*gcov[3,1] + u³*gcov[4,1]
	u₁::Float64 = u⁰*gcov[1,2] + u¹*gcov[2,2] + u²*gcov[3,2] + u³*gcov[4,2]
	u₂::Float64 = u⁰*gcov[1,3] + u¹*gcov[2,3] + u²*gcov[3,3] + u³*gcov[4,3]
	u₃::Float64 = u⁰*gcov[1,4] + u¹*gcov[2,4] + u²*gcov[3,4] + u³*gcov[4,4]

	#Contravariant Four-magnetic field
	b⁰::Float64 = B¹*u₁ + B²*u₂ + B³*u₃
	b¹::Float64 = (B¹ + b⁰*u¹)/u⁰
	b²::Float64 = (B² + b⁰*u²)/u⁰
	b³::Float64 = (B³ + b⁰*u³)/u⁰

	#Covariant Four-magnetic field    
	b₀::Float64 = b⁰*gcov[1,1] + b¹*gcov[1,2] + b²*gcov[1,3] + b³*gcov[1,4]
	b₁::Float64 = b⁰*gcov[2,1] + b¹*gcov[2,2] + b²*gcov[2,3] + b³*gcov[2,4]
	b₂::Float64 = b⁰*gcov[3,1] + b¹*gcov[3,2] + b²*gcov[3,3] + b³*gcov[3,4]
	b₃::Float64 = b⁰*gcov[4,1] + b¹*gcov[4,2] + b²*gcov[4,3] + b³*gcov[4,4]

	#Useful Values        
	bsq   ::Float64 = b⁰*b₀ + b¹*b₁ + b²*b₂ + b³*b₃
	value ::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq   # p + (1/2)*bsq
	sq_g  ::Float64 = sqrt_g(gcov)                    #square root of the determinant of the metric
	
	#For u⁰
	value_sqrt::Float64 = sqrt(b^2 - 4*a*c)

	#For u⁰
	u⁰_x1::Float64 = 0
	u⁰_x2::Float64 = 0
	u⁰_x3::Float64 = (-2*gcov[1,2] - (-2*a*(2*gcov[2,2]*u¹ + 2*gcov[2,3]*u² + 2*gcov[2,4]*u³) + 2*gcov[1,2]*(b))/value_sqrt)/(2*a)
	u⁰_x4::Float64 = (-2*gcov[1,3] - (-2*a*(2*gcov[2,3]*u¹ + 2*gcov[3,3]*u² + 2*gcov[3,4]*u³) + 2*gcov[1,3]*(b))/value_sqrt)/(2*a)	
	u⁰_x5::Float64 = (-2*gcov[1,4] - (-2*a*(2*gcov[2,4]*u¹ + 2*gcov[3,4]*u² + 2*gcov[4,4]*u³) + 2*gcov[1,4]*(b))/value_sqrt)/(2*a)

	#For u¹
	u¹_x1::Float64 = 0
	u¹_x2::Float64 = 0
	u¹_x3::Float64 = 1
	u¹_x4::Float64 = 0
	u¹_x5::Float64 = 0	

	#For u²
	u²_x1::Float64 = 0
	u²_x2::Float64 = 0
	u²_x3::Float64 = 0
	u²_x4::Float64 = 1
	u²_x5::Float64 = 0	

	#For u³
	u³_x1::Float64 = 0
	u³_x2::Float64 = 0
	u³_x3::Float64 = 0
	u³_x4::Float64 = 0
	u³_x5::Float64 = 1

	#For u₀ 
	u₀_x1::Float64 = 0
	u₀_x2::Float64 = 0	
	u₀_x3::Float64 = u⁰_x3*gcov[1,1] + gcov[2,1]
	u₀_x4::Float64 = u⁰_x4*gcov[1,1] + gcov[3,1]
	u₀_x5::Float64 = u⁰_x5*gcov[1,1] + gcov[4,1]	
 
	#For u₁ 
	u₁_x1::Float64 = 0
	u₁_x2::Float64 = 0	
	u₁_x3::Float64 = u⁰_x3*gcov[1,2] + gcov[2,2]
	u₁_x4::Float64 = u⁰_x4*gcov[1,2] + gcov[3,2]
	u₁_x5::Float64 = u⁰_x5*gcov[1,2] + gcov[4,2]	

	#For u₂ 
	u₂_x1::Float64 = 0
	u₂_x2::Float64 = 0	
	u₂_x3::Float64 = u⁰_x3*gcov[1,3] + gcov[2,3]
	u₂_x4::Float64 = u⁰_x4*gcov[1,3] + gcov[3,3]
	u₂_x5::Float64 = u⁰_x5*gcov[1,3] + gcov[4,3]

	#For u₃ 
	u₃_x1::Float64 = 0
	u₃_x2::Float64 = 0	
	u₃_x3::Float64 = u⁰_x3*gcov[1,4] + gcov[2,4]
	u₃_x4::Float64 = u⁰_x4*gcov[1,4] + gcov[3,4]
	u₃_x5::Float64 = u⁰_x5*gcov[1,4] + gcov[4,4]
	
	#For b⁰ 
	b⁰_x1::Float64 = 0
	b⁰_x2::Float64 = 0
	b⁰_x3::Float64 = B¹*u₁_x3 + B²*u₂_x3 + B³*u₃_x3
	b⁰_x4::Float64 = B¹*u₁_x4 + B²*u₂_x4 + B³*u₃_x4	
	b⁰_x5::Float64 = B¹*u₁_x5 + B²*u₂_x5 + B³*u₃_x5	
	
	#For b¹	
	b¹_x1::Float64 = 0
	b¹_x2::Float64 = 0
	b¹_x3::Float64 = ((b⁰_x3*u¹ + u¹_x3*b⁰)*u⁰ - (B¹ + b⁰*u¹)*u⁰_x3)/(u⁰^2)
	b¹_x4::Float64 = ((b⁰_x4*u¹ + u¹_x4*b⁰)*u⁰ - (B¹ + b⁰*u¹)*u⁰_x4)/(u⁰^2)	
	b¹_x5::Float64 = ((b⁰_x5*u¹ + u¹_x5*b⁰)*u⁰ - (B¹ + b⁰*u¹)*u⁰_x5)/(u⁰^2)	

	#For b² 		
	b²_x1::Float64 = 0
	b²_x2::Float64 = 0
	b²_x3::Float64 = ((b⁰_x3*u² + u²_x3*b⁰)*u⁰ - (B² + b⁰*u²)*u⁰_x3)/(u⁰^2)
	b²_x4::Float64 = ((b⁰_x4*u² + u²_x4*b⁰)*u⁰ - (B² + b⁰*u²)*u⁰_x4)/(u⁰^2)	
	b²_x5::Float64 = ((b⁰_x5*u² + u²_x5*b⁰)*u⁰ - (B² + b⁰*u²)*u⁰_x5)/(u⁰^2)		

	#For b³ 		
	b³_x1::Float64 = 0
	b³_x2::Float64 = 0
	b³_x3::Float64 = ((b⁰_x3*u³ + u³_x3*b⁰)*u⁰ - (B³+b⁰*u³)*u⁰_x3)/(u⁰^2)
	b³_x4::Float64 = ((b⁰_x4*u³ + u³_x4*b⁰)*u⁰ - (B³+b⁰*u³)*u⁰_x4)/(u⁰^2)	
	b³_x5::Float64 = ((b⁰_x5*u³ + u³_x5*b⁰)*u⁰ - (B³+b⁰*u³)*u⁰_x5)/(u⁰^2)	

	#For b₀ 		
	b₀_x1::Float64 = 0
	b₀_x2::Float64 = 0
	b₀_x3::Float64 = b⁰_x3*gcov[1,1] + b¹_x3*gcov[1,2] + b²_x3*gcov[1,3] + b³_x3*gcov[1,4]
	b₀_x4::Float64 = b⁰_x4*gcov[1,1] + b¹_x4*gcov[1,2] + b²_x4*gcov[1,3] + b³_x4*gcov[1,4]	
	b₀_x5::Float64 = b⁰_x5*gcov[1,1] + b¹_x5*gcov[1,2] + b²_x5*gcov[1,3] + b³_x5*gcov[1,4]	
	
	#For b₁ 		
	b₁_x1::Float64 = 0
	b₁_x2::Float64 = 0
	b₁_x3::Float64 = b⁰_x3*gcov[2,1] + b¹_x3*gcov[2,2] + b²_x3*gcov[2,3] + b³_x3*gcov[2,4]
	b₁_x4::Float64 = b⁰_x4*gcov[2,1] + b¹_x4*gcov[2,2] + b²_x4*gcov[2,3] + b³_x4*gcov[2,4]
	b₁_x5::Float64 = b⁰_x5*gcov[2,1] + b¹_x5*gcov[2,2] + b²_x5*gcov[2,3] + b³_x5*gcov[2,4]		
	
	#For b₂ 		
	b₂_x1::Float64 = 0
	b₂_x2::Float64 = 0
	b₂_x3::Float64 = b⁰_x3*gcov[3,1] + b¹_x3*gcov[3,2] + b²_x3*gcov[3,3] + b³_x3*gcov[3,4]
	b₂_x4::Float64 = b⁰_x4*gcov[3,1] + b¹_x4*gcov[3,2] + b²_x4*gcov[3,3] + b³_x4*gcov[3,4]
	b₂_x5::Float64 = b⁰_x5*gcov[3,1] + b¹_x5*gcov[3,2] + b²_x5*gcov[3,3] + b³_x5*gcov[3,4]	
	
	#For b₃ 		
	b₃_x1::Float64 = 0
	b₃_x2::Float64 = 0
	b₃_x3::Float64 = b⁰_x3*gcov[4,1] + b¹_x3*gcov[4,2] + b²_x3*gcov[4,3] + b³_x3*gcov[4,4]
	b₃_x4::Float64 = b⁰_x4*gcov[4,1] + b¹_x4*gcov[4,2] + b²_x4*gcov[4,3] + b³_x4*gcov[4,4]
	b₃_x5::Float64 = b⁰_x5*gcov[4,1] + b¹_x5*gcov[4,2] + b²_x5*gcov[4,3] + b³_x5*gcov[4,4]
	
	#For bsq
	bsq_x1::Float64 = 0
	bsq_x2::Float64 = 0
	bsq_x3::Float64 = b⁰_x3*b₀ + b⁰*b₀_x3 + b¹_x3*b₁ + b¹*b₁_x3 + b²_x3*b₂ + b²*b₂_x3 + b³_x3*b₃ + b³*b₃_x3
	bsq_x4::Float64 = b⁰_x4*b₀ + b⁰*b₀_x4 + b¹_x4*b₁ + b¹*b₁_x4 + b²_x4*b₂ + b²*b₂_x4 + b³_x4*b₃ + b³*b₃_x4	
	bsq_x5::Float64 = b⁰_x5*b₀ + b⁰*b₀_x5 + b¹_x5*b₁ + b¹*b₁_x5 + b²_x5*b₂ + b²*b₂_x5 + b³_x5*b₃ + b³*b₃_x5	
	
	#For value
	value_x1::Float64 = 1
	value_x2::Float64 = eos.gamma	
	value_x3::Float64 = bsq_x3	
	value_x4::Float64 = bsq_x4
	value_x5::Float64 = bsq_x5
	
	#For value2
	value2_x1::Float64 = 0
	value2_x2::Float64 = eos.gamma - 1	
	value2_x3::Float64 = (0.5)*bsq_x3	
	value2_x4::Float64 = (0.5)*bsq_x4
	value2_x5::Float64 = (0.5)*bsq_x5	

	#Jacobian:
	buffer[1] = sq_g*u⁰
	buffer[6] = sq_g*(value_x1*u⁰*u₀ + value2_x1)
	buffer[11] = sq_g*(value_x1*u⁰*u₁)
	buffer[16] = sq_g*(value_x1*u⁰*u₂)
	buffer[21] = sq_g*(value_x1*u⁰*u₃)

	buffer[2] = 0
	buffer[7] = sq_g*(value_x2*u⁰*u₀ + value2_x2)
	buffer[12] = sq_g*(value_x2*u⁰*u₁)
	buffer[17] = sq_g*(value_x2*u⁰*u₂)
	buffer[22] = sq_g*(value_x2*u⁰*u₃)

	buffer[3] =  sq_g*ρ*u⁰_x3
	buffer[8] =  sq_g*((value_x3*(u⁰*u₀) + value*(u⁰_x3*u₀ + u₀_x3*u⁰)) + value2_x3 - (b⁰_x3*b₀ + b₀_x3*b⁰))
	buffer[13] = sq_g*((value_x3*(u⁰*u₁) + value*(u⁰_x3*u₁ + u₁_x3*u⁰)) - (b⁰_x3*b₁ + b₁_x3*b⁰))
	buffer[18] = sq_g*((value_x3*(u⁰*u₂) + value*(u⁰_x3*u₂ + u₂_x3*u⁰)) - (b⁰_x3*b₂ + b₂_x3*b⁰))
	buffer[23] = sq_g*((value_x3*(u⁰*u₃) + value*(u⁰_x3*u₃ + u₃_x3*u⁰)) - (b⁰_x3*b₃ + b₃_x3*b⁰))

	buffer[4] = sq_g*ρ*u⁰_x4
	buffer[9] = sq_g*((value_x4*(u⁰*u₀)  + value*(u⁰_x4*u₀ + u₀_x4*u⁰)) + value2_x4 - (b⁰_x4*b₀ + b₀_x4*b⁰))
	buffer[14] = sq_g*((value_x4*(u⁰*u₁) + value*(u⁰_x4*u₁ + u₁_x4*u⁰)) - (b⁰_x4*b₁ + b₁_x4*b⁰))
	buffer[19] = sq_g*((value_x4*(u⁰*u₂) + value*(u⁰_x4*u₂ + u₂_x4*u⁰)) - (b⁰_x4*b₂ + b₂_x4*b⁰))
	buffer[24] = sq_g*((value_x4*(u⁰*u₃) + value*(u⁰_x4*u₃ + u₃_x4*u⁰)) - (b⁰_x4*b₃ + b₃_x4*b⁰))

	buffer[5] = sq_g*ρ*u⁰_x5
	buffer[10] = sq_g*((value_x5*(u⁰*u₀) + value*(u⁰_x5*u₀ + u₀_x5*u⁰)) + value2_x5 - (b⁰_x5*b₀ + b₀_x5*b⁰))
	buffer[15] = sq_g*((value_x5*(u⁰*u₁) + value*(u⁰_x5*u₁ + u₁_x5*u⁰)) - (b⁰_x5*b₁ + b₁_x5*b⁰))
	buffer[20] = sq_g*((value_x5*(u⁰*u₂) + value*(u⁰_x5*u₂ + u₂_x5*u⁰)) - (b⁰_x5*b₂ + b₂_x5*b⁰))
	buffer[25] = sq_g*((value_x5*(u⁰*u₃) + value*(u⁰_x5*u₃ + u₃_x5*u⁰)) - (b⁰_x5*b₃ + b₃_x5*b⁰))

end

function invert_matrix_lu(A::AbstractMatrix)
    size(A, 1) == size(A, 2) || error("Macierz musi być kwadratowa!")
    lu_fact = lu(A)
    invA = inv(lu_fact)
    return invA
end


function UtoP(U::AbstractVector, initial_guess::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
    
    #Some values
    max_iter::Int64 = 10000 #Maximum amount of interations
    buff_jac = MVector{25, Float64}(zeros(25))
    buff_fun = MVector{8, Float64}(zeros(8))
    buff_jac_inv::Matrix{Float64} = zeros(5,5)
    buff_jac_reshape::Matrix{Float64} = zeros(5,5)
    
    sq_g::Float64 = sqrt_g(gcov) #square root of the determinant of the metric
    x = initial_guess #rename of x-vector
    
    #These three equations for the magnetic field can be solved analitically
    x[6] = U[6] / sq_g 
    x[7] = U[7] / sq_g 
    x[8] = U[8] / sq_g 
    
    #Newton-Rhapson algoritm for the non-linear system of equations 
    for i in 1:max_iter
        Jacobian(x, buff_jac, gcov, eos::Polytrope)
        buff_jac_reshape = reshape(buff_jac, 5, 5)'
        buff_jac_inv = invert_matrix_lu_manual(buff_jac_reshape)
        PtoU(x, buff_fun, gcov, eos::Polytrope)
        buff_fun[1:5] = buff_fun[1:5] - U[1:5]
        x[1:5] = x[1:5] - buff_jac_inv * buff_fun[1:5]
        if sqrt(buff_fun[1]^2 + buff_fun[2]^2 + buff_fun[3]^2 + buff_fun[4]^2 + buff_fun[5]^2) < 10^-8
        	#println("Zbiegło!!!") 
        	break
        end
    end
    if sqrt(buff_fun[1]^2 + buff_fun[2]^2 + buff_fun[3]^2 + buff_fun[4]^2 + buff_fun[5]^2) > 10^-8
    	println("Nie zbiegło!!!") 
    #	println(sqrt(buff_fun[1]^2 + buff_fun[2]^2 + buff_fun[3]^2 + buff_fun[4]^2 + buff_fun[5]^2))
    end
    
    return x
end
