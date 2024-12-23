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
