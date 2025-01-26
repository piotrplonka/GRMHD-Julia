using LinearAlgebra
abstract type EOS end


struct Polytrope <:EOS
    gamma::Float64
end

function Pressure(u::Float64,eos::Polytrope)::Float64
    return (eos.gamma-1)*u
end

eos = Polytrope(5/3)

function SoundSpeed(rho::Float64,u::Float64,eos::Polytrope)::Float64
    return sqrt((eos.gamma * (eos.gamma - 1) * u )/(rho + eos.gamma * u))
end


function MagnetosonicSpeed(x::AbstractVector, gcov::Matrix{Float64}, eos::Polytrope)

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
        EF::Float64  = ρ + eos.gamma*u
        Cs_mag2::Float64 = bsq/EF
        Cs2::Float64 = (eos.gamma * (eos.gamma - 1) * u )/(ρ + eos.gamma * u)
        
	return sqrt(Cs2 + Cs_mag2 - Cs2*Cs_mag2)

	

end
