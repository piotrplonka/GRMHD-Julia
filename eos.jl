using LinearAlgebra
abstract type EOS end


struct Polytrope <:EOS
    gamma::Float64
end

function Pressure(u::Float64,eos::Polytrope)::Float64
    return (eos.gamma-1)*u
end

eos = Polytrope(5/3)

function SoundSpeed(x::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)
	
	#Parameters
	ρ ::Float64 = x[1] #Density
	u ::Float64 = x[2] #Internal Energy 
	u¹::Float64 = x[3] #Contravariant Four-velocity in 1-direction
	u²::Float64 = x[4] #Contravariant Four-velocity in 2-direction
	u³::Float64 = x[5] #Contravariant Four-velocity in 3-direction   
	B¹::Float64 = x[6] #Magnetic field in 1-direction
	B²::Float64 = x[7] #Magnetic field in 2-direction
	B³::Float64 = x[8] #Magnetic field in 3-direction
	
return sqrt((eos.gamma * (eos.gamma - 1) * u )/(ρ + eos.gamma * u))
end


function MagnetosonicSpeed(x::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)

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
	bsq     ::Float64 = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3
        EF      ::Float64 = ρ + eos.gamma*u
        Cs_mag2 ::Float64 = bsq/EF
        Cs2     ::Float64 = (eos.gamma * (eos.gamma - 1) * u )/(ρ + eos.gamma * u)
        
	return sqrt(Cs2 + Cs_mag2- Cs2*Cs_mag2) 
end
