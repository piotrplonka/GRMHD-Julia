using BenchmarkTools
using Base.Threads
using StaticArrays

include("Matrix_operations.jl")
include("Christoffel_Symbols.jl")
include("eos.jl")

function PtoU(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
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
	bsq::Float64 = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3
	value::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq # p + (1/2)*bsq
	sq_g::Float64 = sqrt_g(gcov) #square root of the determinant of the metric

	#Buffers
	buffer[1] = sq_g*ρ*u0 
	buffer[2] = sq_g*((value)*u0*u_0 + value2 - b0*b_0)
	buffer[3] = sq_g*((value)*u0*u_1 - b0*b_1)
	buffer[4] = sq_g*((value)*u0*u_2 - b0*b_2)
	buffer[5] = sq_g*((value)*u0*u_3 - b0*b_3)
	buffer[6] = sq_g*B1
	buffer[7] = sq_g*B2
	buffer[8] = sq_g*B3

end

function PtoFx(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
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
	bsq::Float64 = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3
	value::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq # p + (1/2)*bsq
	sq_g::Float64 = sqrt_g(gcov) #square root of the determinant of the metric
	
	#Buffers    
	buffer[1] = sq_g*ρ*u1
	buffer[2] = sq_g*(value*u1*u_0 - b1*b_0)
	buffer[3] = sq_g*(value*u1*u_1 + value2 - b1*b_1)
	buffer[4] = sq_g*(value*u1*u_2 - b1*b_2)
	buffer[5] = sq_g*(value*u1*u_3 - b1*b_3)

end

function PtoFy(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
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
	bsq::Float64 = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3
	value::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq # p + (1/2)*bsq
	sq_g::Float64 = sqrt_g(gcov) #square root of the determinant of the metric
	
	#Buffers
	buffer[1] = sq_g*ρ*u2
	buffer[2] = sq_g*(value*u2*u_0 - b2*b_0)
	buffer[3] = sq_g*(value*u2*u_1 - b2*b_1)
	buffer[4] = sq_g*(value*u2*u_2 + value2 - b2*b_2)
	buffer[5] = sq_g*(value*u2*u_3 - b2*b_3)

end

function PtoFz(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
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
	bsq::Float64 = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3
	value::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq # p + (1/2)*bsq
	sq_g::Float64 = sqrt_g(gcov) #square root of the determinant of the metric
	
	#Buffers
	buffer[1] = sq_g*ρ*u3
	buffer[2] = sq_g*(value*u3*u_0 - b3*b_0)
	buffer[3] = sq_g*(value*u3*u_1 - b3*b_1)
	buffer[4] = sq_g*(value*u3*u_2 - b3*b_2)
	buffer[5] = sq_g*(value*u3*u_3 + value2 - b3*b_3)

end

function PtoF_Bx_By_Bz(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
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
	
	#Useful Values        
	sq_g::Float64 = sqrt_g(gcov) #square root of the determinant of the metric
	
	#Buffers
	buffer[1] = sq_g*(b0*u1-b1*u0 + b2*u1-b1*u2 + b3*u1-b1*u3)
	buffer[2] = sq_g*(b0*u2-b2*u0 + b1*u2-b2*u1 + b3*u2-b2*u3)
	buffer[3] = sq_g*(b0*u3-b3*u0 + b1*u3-b3*u1 + b2*u3-b3*u2)
end

function Lorentz_factor(x::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
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

	#For C=1 lorentz factor is:
	γ::Float64 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
	return γ
end


function Jacobian(x::AbstractVector, buffer::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
	
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
	bsq::Float64 = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3
	value::Float64 = ρ + u + (eos.gamma - 1)*u + bsq # ρ + u + p 
	value2::Float64 = (eos.gamma - 1)*u + (1/2)*bsq # p + (1/2)*bsq
	sq_g::Float64 = sqrt_g(gcov) #square root of the determinant of the metric
	
	#For u0
	value_sqrt::Float64 = sqrt(b^2 - 4*a*c)

	#For u0
	u0_x1::Float64 = 0
	u0_x2::Float64 = 0
	u0_x3::Float64 = (-2*gcov[1,2] - (-2*a*(2*gcov[2,2]*u1 + 2*gcov[2,3]*u2 + 2*gcov[2,4]*u3) + 2*gcov[1,2]*(b))/value_sqrt)/(2*a)
	u0_x4::Float64 = (-2*gcov[1,3] - (-2*a*(2*gcov[2,3]*u1 + 2*gcov[3,3]*u2 + 2*gcov[3,4]*u3) + 2*gcov[1,3]*(b))/value_sqrt)/(2*a)	
	u0_x5::Float64 = (-2*gcov[1,4] - (-2*a*(2*gcov[2,4]*u1 + 2*gcov[3,4]*u2 + 2*gcov[4,4]*u3) + 2*gcov[1,4]*(b))/value_sqrt)/(2*a)

	#For u1
	u1_x1::Float64 = 0
	u1_x2::Float64 = 0
	u1_x3::Float64 = 1
	u1_x4::Float64 = 0
	u1_x5::Float64 = 0	

	#For u2
	u2_x1::Float64 = 0
	u2_x2::Float64 = 0
	u2_x3::Float64 = 0
	u2_x4::Float64 = 1
	u2_x5::Float64 = 0	

	#For u3
	u3_x1::Float64 = 0
	u3_x2::Float64 = 0
	u3_x3::Float64 = 0
	u3_x4::Float64 = 0
	u3_x5::Float64 = 1

	#For u_0 
	u_0_x1::Float64 = 0
	u_0_x2::Float64 = 0	
	u_0_x3::Float64 = u0_x3*gcov[1,1] + gcov[2,1]
	u_0_x4::Float64 = u0_x4*gcov[1,1] + gcov[3,1]
	u_0_x5::Float64 = u0_x5*gcov[1,1] + gcov[4,1]	
 
	#For u_1 
	u_1_x1::Float64 = 0
	u_1_x2::Float64 = 0	
	u_1_x3::Float64 = u0_x3*gcov[1,2] + gcov[2,2]
	u_1_x4::Float64 = u0_x4*gcov[1,2] + gcov[3,2]
	u_1_x5::Float64 = u0_x5*gcov[1,2] + gcov[4,2]	

	#For u_2 
	u_2_x1::Float64 = 0
	u_2_x2::Float64 = 0	
	u_2_x3::Float64 = u0_x3*gcov[1,3] + gcov[2,3]
	u_2_x4::Float64 = u0_x4*gcov[1,3] + gcov[3,3]
	u_2_x5::Float64 = u0_x5*gcov[1,3] + gcov[4,3]

	#For u_3 
	u_3_x1::Float64 = 0
	u_3_x2::Float64 = 0	
	u_3_x3::Float64 = u0_x3*gcov[1,4] + gcov[2,4]
	u_3_x4::Float64 = u0_x4*gcov[1,4] + gcov[3,4]
	u_3_x5::Float64 = u0_x5*gcov[1,4] + gcov[4,4]
	
	#For b0 
	b0_x1::Float64 = 0
	b0_x2::Float64 = 0
	b0_x3::Float64 = B1*u_1_x3 + B2*u_2_x3 + B3*u_3_x3
	b0_x4::Float64 = B1*u_1_x4 + B2*u_2_x4 + B3*u_3_x4	
	b0_x5::Float64 = B1*u_1_x5 + B2*u_2_x5 + B3*u_3_x5	
	
	#For b1	
	b1_x1::Float64 = 0
	b1_x2::Float64 = 0
	b1_x3::Float64 = ((b0_x3*u1 + u1_x3*b0)*u0 - (B1 + b0*u1)*u0_x3)/(u0^2)
	b1_x4::Float64 = ((b0_x4*u1 + u1_x4*b0)*u0 - (B1 + b0*u1)*u0_x4)/(u0^2)	
	b1_x5::Float64 = ((b0_x5*u1 + u1_x5*b0)*u0 - (B1 + b0*u1)*u0_x5)/(u0^2)	

	#For b2 		
	b2_x1::Float64 = 0
	b2_x2::Float64 = 0
	b2_x3::Float64 = ((b0_x3*u2 + u2_x3*b0)*u0 - (B2 + b0*u2)*u0_x3)/(u0^2)
	b2_x4::Float64 = ((b0_x4*u2 + u2_x4*b0)*u0 - (B2 + b0*u2)*u0_x4)/(u0^2)	
	b2_x5::Float64 = ((b0_x5*u2 + u2_x5*b0)*u0 - (B2 + b0*u2)*u0_x5)/(u0^2)		

	#For b3 		
	b3_x1::Float64 = 0
	b3_x2::Float64 = 0
	b3_x3::Float64 = ((b0_x3*u3 + u3_x3*b0)*u0 - (B3+b0*u3)*u0_x3)/(u0^2)
	b3_x4::Float64 = ((b0_x4*u3 + u3_x4*b0)*u0 - (B3+b0*u3)*u0_x4)/(u0^2)	
	b3_x5::Float64 = ((b0_x5*u3 + u3_x5*b0)*u0 - (B3+b0*u3)*u0_x5)/(u0^2)	

	#For b_0 		
	b_0_x1::Float64 = 0
	b_0_x2::Float64 = 0
	b_0_x3::Float64 = b0_x3*gcov[1,1] + b1_x3*gcov[1,2] + b2_x3*gcov[1,3] + b3_x3*gcov[1,4]
	b_0_x4::Float64 = b0_x4*gcov[1,1] + b1_x4*gcov[1,2] + b2_x4*gcov[1,3] + b3_x4*gcov[1,4]	
	b_0_x5::Float64 = b0_x5*gcov[1,1] + b1_x5*gcov[1,2] + b2_x5*gcov[1,3] + b3_x5*gcov[1,4]	
	
	#For b_1 		
	b_1_x1::Float64 = 0
	b_1_x2::Float64 = 0
	b_1_x3::Float64 = b0_x3*gcov[2,1] + b1_x3*gcov[2,2] + b2_x3*gcov[2,3] + b3_x3*gcov[2,4]
	b_1_x4::Float64 = b0_x4*gcov[2,1] + b1_x4*gcov[2,2] + b2_x4*gcov[2,3] + b3_x4*gcov[2,4]
	b_1_x5::Float64 = b0_x5*gcov[2,1] + b1_x5*gcov[2,2] + b2_x5*gcov[2,3] + b3_x5*gcov[2,4]		
	
	#For b_2 		
	b_2_x1::Float64 = 0
	b_2_x2::Float64 = 0
	b_2_x3::Float64 = b0_x3*gcov[3,1] + b1_x3*gcov[3,2] + b2_x3*gcov[3,3] + b3_x3*gcov[3,4]
	b_2_x4::Float64 = b0_x4*gcov[3,1] + b1_x4*gcov[3,2] + b2_x4*gcov[3,3] + b3_x4*gcov[3,4]
	b_2_x5::Float64 = b0_x5*gcov[3,1] + b1_x5*gcov[3,2] + b2_x5*gcov[3,3] + b3_x5*gcov[3,4]	
	
	#For b_3 		
	b_3_x1::Float64 = 0
	b_3_x2::Float64 = 0
	b_3_x3::Float64 = b0_x3*gcov[4,1] + b1_x3*gcov[4,2] + b2_x3*gcov[4,3] + b3_x3*gcov[4,4]
	b_3_x4::Float64 = b0_x4*gcov[4,1] + b1_x4*gcov[4,2] + b2_x4*gcov[4,3] + b3_x4*gcov[4,4]
	b_3_x5::Float64 = b0_x5*gcov[4,1] + b1_x5*gcov[4,2] + b2_x5*gcov[4,3] + b3_x5*gcov[4,4]
	
	#For bsq
	bsq_x1::Float64 = 0
	bsq_x2::Float64 = 0
	bsq_x3::Float64 = b0_x3*b_0 + b0*b_0_x3 + b1_x3*b_1 + b1*b_1_x3 + b2_x3*b_2 + b2*b_2_x3 + b3_x3*b_3 + b3*b_3_x3
	bsq_x4::Float64 = b0_x4*b_0 + b0*b_0_x4 + b1_x4*b_1 + b1*b_1_x4 + b2_x4*b_2 + b2*b_2_x4 + b3_x4*b_3 + b3*b_3_x4	
	bsq_x5::Float64 = b0_x5*b_0 + b0*b_0_x5 + b1_x5*b_1 + b1*b_1_x5 + b2_x5*b_2 + b2*b_2_x5 + b3_x5*b_3 + b3*b_3_x5	
	
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
	buffer[1] = sq_g*u0
	buffer[6] = sq_g*(value_x1*u0*u_0 + value2_x1)
	buffer[11] = sq_g*(value_x1*u0*u_1)
	buffer[16] = sq_g*(value_x1*u0*u_2)
	buffer[21] = sq_g*(value_x1*u0*u_3)

	buffer[2] = 0
	buffer[7] = sq_g*(value_x2*u0*u_0 + value2_x2)
	buffer[12] = sq_g*(value_x2*u0*u_1)
	buffer[17] = sq_g*(value_x2*u0*u_2)
	buffer[22] = sq_g*(value_x2*u0*u_3)

	buffer[3] =  sq_g*ρ*u0_x3
	buffer[8] =  sq_g*((value_x3*(u0*u_0) + value*(u0_x3*u_0 + u_0_x3*u0)) + value2_x3 - (b0_x3*b_0 + b_0_x3*b0))
	buffer[13] = sq_g*((value_x3*(u0*u_1) + value*(u0_x3*u_1 + u_1_x3*u0)) - (b0_x3*b_1 + b_1_x3*b0))
	buffer[18] = sq_g*((value_x3*(u0*u_2) + value*(u0_x3*u_2 + u_2_x3*u0)) - (b0_x3*b_2 + b_2_x3*b0))
	buffer[23] = sq_g*((value_x3*(u0*u_3) + value*(u0_x3*u_3 + u_3_x3*u0)) - (b0_x3*b_3 + b_3_x3*b0))

	buffer[4] = sq_g*ρ*u0_x4
	buffer[9] = sq_g*((value_x4*(u0*u_0) + value*(u0_x4*u_0 + u_0_x4*u0)) + value2_x4 - (b0_x4*b_0 + b_0_x4*b0))
	buffer[14] = sq_g*((value_x4*(u0*u_1) + value*(u0_x4*u_1 + u_1_x4*u0)) - (b0_x4*b_1 + b_1_x4*b0))
	buffer[19] = sq_g*((value_x4*(u0*u_2) + value*(u0_x4*u_2 + u_2_x4*u0)) - (b0_x4*b_2 + b_2_x4*b0))
	buffer[24] = sq_g*((value_x4*(u0*u_3) + value*(u0_x4*u_3 + u_3_x4*u0)) - (b0_x4*b_3 + b_3_x4*b0))

	buffer[5] = sq_g*ρ*u0_x5
	buffer[10] = sq_g*((value_x5*(u0*u_0) + value*(u0_x5*u_0 + u_0_x5*u0)) + value2_x5 - (b0_x5*b_0 + b_0_x5*b0))
	buffer[15] = sq_g*((value_x5*(u0*u_1) + value*(u0_x5*u_1 + u_1_x5*u0)) - (b0_x5*b_1 + b_1_x5*b0))
	buffer[20] = sq_g*((value_x5*(u0*u_2) + value*(u0_x5*u_2 + u_2_x5*u0)) - (b0_x5*b_2 + b_2_x5*b0))
	buffer[25] = sq_g*((value_x5*(u0*u_3) + value*(u0_x5*u_3 + u_3_x5*u0)) - (b0_x5*b_3 + b_3_x5*b0))

end

function UtoP(U::AbstractVector, initial_guess::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
    
    #Some values
    max_iter::Int64 = 100 #Maximum amount of interations
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
        buff_jac_inv = inverse_5x5(buff_jac_reshape)
        PtoU(x, buff_fun, gcov, eos::Polytrope)
        buff_fun[1:5]=buff_fun[1:5] - U[1:5]
        x[1:5] = x[1:5] - buff_jac_inv * buff_fun[1:5]
        if sqrt(buff_fun[1]^2 + buff_fun[2]^2 + buff_fun[3]^2 + buff_fun[4]^2 + buff_fun[5]^2) < 10^-10
        	break
        end
    end
    return x
end
