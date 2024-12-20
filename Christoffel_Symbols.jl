function Schwarzschild_metric_cov(r::Float64, θ::Float64, φ::Float64)
	g = zeros(4, 4)
	g[1, 1] = (2.0 - r)/r
	g[2, 2] = r/(r - 2.0)
	g[3, 3] = r^2
	g[4, 4] = r^2*sin(θ)^2
	return g
end

function Kerr_Schild_metric_cov(r::Float64, θ::Float64, φ::Float64, a::Float64)
	g = zeros(4, 4)
	cos_θ = cos(θ)
	r_2 = r^2
	cos_θ_2 = cos_θ^2
	a_2 = a^2
	cos_θ = cos(θ)
	r_2 = r^2
	cos_θ_2 = cos_θ^2
	a_2 = a^2
	cos_θ = cos(θ)
	r_2 = r^2
	cos_θ_2 = cos_θ^2
	a_2 = a^2
	cos_θ = cos(θ)
	r_2 = r^2
	cos_θ_2 = cos_θ^2
	a_2 = a^2
	sine_θ = sin(θ)
	cos_θ = cos(θ)
	r_2 = r^2
	sine_θ_2 = sine_θ^2
	cos_θ_2 = cos_θ^2
	a_2 = a^2
	g[1, 1] = (-a_2*cos_θ_2 - r_2 + 2*r)/(a_2*cos_θ_2 + r_2)
	g[1, 2] = 4*r/(a^2*cos(θ)^2 + r^2)
	g[1, 4] = -4*a*r*sin(θ)^2/(a^2*cos(θ)^2 + r^2)
	g[2, 1] = 4*r/(a^2*cos(θ)^2 + r^2)
	g[2, 2] = (a_2*cos_θ_2 + r_2 + 2*r)/(a_2*cos_θ_2 + r_2)
	g[2, 4] = -2*a*(a_2*cos_θ_2 + r_2 + 2*r)*sin(θ)^2/(a_2*cos_θ_2 + r_2)
	g[3, 3] = a^2*cos(θ)^2 + r^2
	g[4, 1] = -4*a*r*sin(θ)^2/(a^2*cos(θ)^2 + r^2)
	g[4, 2] = -2*a*(a_2*cos_θ_2 + r_2 + 2*r)*sin(θ)^2/(a_2*cos_θ_2 + r_2)
	g[4, 4] = (a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2)*sine_θ_2/(a_2*cos_θ_2 + r_2)
	return g
end

	
function Kerr_metric_cov(r::Float64, θ::Float64, φ::Float64, a::Float64)
	g = zeros(4, 4)
	g[1, 1] = -2.0*r/(M^2*a^2*c^2*cos(θ)^2 + r^2) + 1
	g[1, 4] = 2.0*M*a*c*r*sin(θ)^2/(M^2*a^2*c^2*cos(θ)^2 + r^2)
	g[2, 2] = -(M^2*a^2*c^2*cos(θ)^2 + r^2)/(M^2*a^2*c^2 + r^2 - 2.0*r)
	g[3, 3] = -M^2*a^2*c^2*cos(θ)^2 - r^2
	g[4, 1] = 2.0*M*a*c*r*sin(θ)^2/(M^2*a^2*c^2*cos(θ)^2 + r^2)
	g[4, 4] = (-2.0*M^2*a^2*c^2*r*sin(θ)^2 + (-M^2*a^2*c^2 - r^2)*(M^2*a^2*c^2*cos(θ)^2 + r^2))*sin(θ)^2/(M^2*a^2*c^2*cos(θ)^2 + r^2)
	return g
end
	
function Kerr_BL_metric_cov(r::Float64, θ::Float64, φ::Float64, a::Float64)
	g = zeros(4, 4)
	g[1, 1] = (a^2*cos(θ)^2 + r^2 - 2*r)/(a^2*cos(θ)^2 + r^2)
	g[1, 4] = 2*a*r*sin(θ)^2/(a^2*cos(θ)^2 + r^2)
	g[2, 2] = -(a^2*cos(θ)^2 + r^2)/(a^2 + r^2 - 2*r)
	g[3, 3] = -a^2*cos(θ)^2 - r^2
	g[4, 1] = 2*a*r*sin(θ)^2/(a^2*cos(θ)^2 + r^2)
	g[4, 4] = (-2*a^2*r*sin(θ)^2 + (-a^2 - r^2)*(a^2*cos(θ)^2 + r^2))*sin(θ)^2/(a^2*cos(θ)^2 + r^2)
	return g
end
	
function Schwarzschild_metric_con(r::Float64, θ::Float64, φ::Float64)
	g = zeros(4, 4)
	g[1, 1] = -r/(r - 2.0)
	g[2, 2] = (r - 2.0)/r
	g[3, 3] = r^(-2)
	g[4, 4] = 1/(r^2*sin(θ)^2)
	return g
end
	
function Kerr_Schild_metric_con(r::Float64, θ::Float64, φ::Float64, a::Float64)
	g = zeros(4, 4)
	g[1, 1] = (-3*a^8*sin(θ)^2*cos(θ)^6 + a^8*cos(θ)^8 - 9*a^6*r^2*sin(θ)^2*cos(θ)^4 + 4*a^6*r^2*cos(θ)^6 - 12*a^6*r*sin(θ)^2*cos(θ)^4 + 2*a^6*r*cos(θ)^6 - 9*a^4*r^4*(1 - cos(4*θ))/8 + 6*a^4*r^4*cos(θ)^4 - 3*a^4*r^3*(1 - cos(4*θ)) + 6*a^4*r^3*cos(θ)^4 - 3*a^4*r^2*(1 - cos(4*θ))/2 - 3*a^2*r^6*sin(θ)^2 + 4*a^2*r^6*cos(θ)^2 - 12*a^2*r^5*sin(θ)^2 + 6*a^2*r^5*cos(θ)^2 - 12*a^2*r^4*sin(θ)^2 + r^8 + 2*r^7)/(3*a^8*sin(θ)^2*cos(θ)^6 - a^8*cos(θ)^8 + 9*a^6*r^2*sin(θ)^2*cos(θ)^4 - 4*a^6*r^2*cos(θ)^6 + 6*a^6*r*sin(θ)^2*cos(θ)^4 + 9*a^4*r^4*(1 - cos(4*θ))/8 - 6*a^4*r^4*cos(θ)^4 + 3*a^4*r^3*(1 - cos(4*θ))/2 + 5*a^4*r^2*(1 - cos(4*θ))/2 - 12*a^4*r^2*cos(θ)^4 + 3*a^2*r^6*sin(θ)^2 - 4*a^2*r^6*cos(θ)^2 + 6*a^2*r^5*sin(θ)^2 + 20*a^2*r^4*sin(θ)^2 - 24*a^2*r^4*cos(θ)^2 + 40*a^2*r^3*sin(θ)^2 - r^8 - 12*r^6)
	g[1, 2] = 4*r*(a^6*sin(θ)^2*cos(θ)^4 - a^6*cos(θ)^6 + a^4*r^2*(1 - cos(4*θ))/4 - 3*a^4*r^2*cos(θ)^4 + a^4*r*(1 - cos(4*θ))/4 + a^2*r^4*sin(θ)^2 - 3*a^2*r^4*cos(θ)^2 + 2*a^2*r^3*sin(θ)^2 - r^6)/(3*a^8*sin(θ)^2*cos(θ)^6 - a^8*cos(θ)^8 + 9*a^6*r^2*sin(θ)^2*cos(θ)^4 - 4*a^6*r^2*cos(θ)^6 + 6*a^6*r*sin(θ)^2*cos(θ)^4 + 9*a^4*r^4*(1 - cos(4*θ))/8 - 6*a^4*r^4*cos(θ)^4 + 3*a^4*r^3*(1 - cos(4*θ))/2 + 5*a^4*r^2*(1 - cos(4*θ))/2 - 12*a^4*r^2*cos(θ)^4 + 3*a^2*r^6*sin(θ)^2 - 4*a^2*r^6*cos(θ)^2 + 6*a^2*r^5*sin(θ)^2 + 20*a^2*r^4*sin(θ)^2 - 24*a^2*r^4*cos(θ)^2 + 40*a^2*r^3*sin(θ)^2 - r^8 - 12*r^6)
	g[1, 4] = 4*a*r*(a^4*cos(θ)^4 + 2*a^2*r^2*cos(θ)^2 + 2*a^2*r*cos(θ)^2 + r^4 + 2*r^3)/(-3*a^8*sin(θ)^2*cos(θ)^6 + a^8*cos(θ)^8 - 9*a^6*r^2*sin(θ)^2*cos(θ)^4 + 4*a^6*r^2*cos(θ)^6 - 6*a^6*r*sin(θ)^2*cos(θ)^4 - 9*a^4*r^4*(1 - cos(4*θ))/8 + 6*a^4*r^4*cos(θ)^4 - 3*a^4*r^3*(1 - cos(4*θ))/2 - 5*a^4*r^2*(1 - cos(4*θ))/2 + 12*a^4*r^2*cos(θ)^4 - 3*a^2*r^6*sin(θ)^2 + 4*a^2*r^6*cos(θ)^2 - 6*a^2*r^5*sin(θ)^2 - 20*a^2*r^4*sin(θ)^2 + 24*a^2*r^4*cos(θ)^2 - 40*a^2*r^3*sin(θ)^2 + r^8 + 12*r^6)
	g[2, 1] = 4*r*(a^6*sin(θ)^2*cos(θ)^4 - a^6*cos(θ)^6 + a^4*r^2*(1 - cos(4*θ))/4 - 3*a^4*r^2*cos(θ)^4 + a^4*r*(1 - cos(4*θ))/4 + a^2*r^4*sin(θ)^2 - 3*a^2*r^4*cos(θ)^2 + 2*a^2*r^3*sin(θ)^2 - r^6)/(3*a^8*sin(θ)^2*cos(θ)^6 - a^8*cos(θ)^8 + 9*a^6*r^2*sin(θ)^2*cos(θ)^4 - 4*a^6*r^2*cos(θ)^6 + 6*a^6*r*sin(θ)^2*cos(θ)^4 + 9*a^4*r^4*(1 - cos(4*θ))/8 - 6*a^4*r^4*cos(θ)^4 + 3*a^4*r^3*(1 - cos(4*θ))/2 + 5*a^4*r^2*(1 - cos(4*θ))/2 - 12*a^4*r^2*cos(θ)^4 + 3*a^2*r^6*sin(θ)^2 - 4*a^2*r^6*cos(θ)^2 + 6*a^2*r^5*sin(θ)^2 + 20*a^2*r^4*sin(θ)^2 - 24*a^2*r^4*cos(θ)^2 + 40*a^2*r^3*sin(θ)^2 - r^8 - 12*r^6)
	g[2, 2] = (-a^8*sin(θ)^2*cos(θ)^6 - a^8*cos(θ)^8 - 3*a^6*r^2*sin(θ)^2*cos(θ)^4 - 4*a^6*r^2*cos(θ)^6 + 2*a^6*r*cos(θ)^6 - 3*a^4*r^4*(1 - cos(4*θ))/8 - 6*a^4*r^4*cos(θ)^4 + 6*a^4*r^3*cos(θ)^4 - 3*a^4*r^2*(1 - cos(4*θ))/2 - a^2*r^6*sin(θ)^2 - 4*a^2*r^6*cos(θ)^2 + 6*a^2*r^5*cos(θ)^2 - 12*a^2*r^4*sin(θ)^2 - r^8 + 2*r^7)/(3*a^8*sin(θ)^2*cos(θ)^6 - a^8*cos(θ)^8 + 9*a^6*r^2*sin(θ)^2*cos(θ)^4 - 4*a^6*r^2*cos(θ)^6 + 6*a^6*r*sin(θ)^2*cos(θ)^4 + 9*a^4*r^4*(1 - cos(4*θ))/8 - 6*a^4*r^4*cos(θ)^4 + 3*a^4*r^3*(1 - cos(4*θ))/2 + 5*a^4*r^2*(1 - cos(4*θ))/2 - 12*a^4*r^2*cos(θ)^4 + 3*a^2*r^6*sin(θ)^2 - 4*a^2*r^6*cos(θ)^2 + 6*a^2*r^5*sin(θ)^2 + 20*a^2*r^4*sin(θ)^2 - 24*a^2*r^4*cos(θ)^2 + 40*a^2*r^3*sin(θ)^2 - r^8 - 12*r^6)
	g[2, 4] = 2*a*(a^6*cos(θ)^6 + 3*a^4*r^2*cos(θ)^4 + 3*a^2*r^4*cos(θ)^2 + 4*a^2*r^2*cos(θ)^2 + r^6 + 4*r^4)/(-3*a^8*sin(θ)^2*cos(θ)^6 + a^8*cos(θ)^8 - 9*a^6*r^2*sin(θ)^2*cos(θ)^4 + 4*a^6*r^2*cos(θ)^6 - 6*a^6*r*sin(θ)^2*cos(θ)^4 - 9*a^4*r^4*(1 - cos(4*θ))/8 + 6*a^4*r^4*cos(θ)^4 - 3*a^4*r^3*(1 - cos(4*θ))/2 - 5*a^4*r^2*(1 - cos(4*θ))/2 + 12*a^4*r^2*cos(θ)^4 - 3*a^2*r^6*sin(θ)^2 + 4*a^2*r^6*cos(θ)^2 - 6*a^2*r^5*sin(θ)^2 - 20*a^2*r^4*sin(θ)^2 + 24*a^2*r^4*cos(θ)^2 - 40*a^2*r^3*sin(θ)^2 + r^8 + 12*r^6)
	g[3, 3] = 1/(a^2*cos(θ)^2 + r^2)
	g[4, 1] = 4*a*r*(a^4*cos(θ)^4 + 2*a^2*r^2*cos(θ)^2 + 2*a^2*r*cos(θ)^2 + r^4 + 2*r^3)/(-3*a^8*sin(θ)^2*cos(θ)^6 + a^8*cos(θ)^8 - 9*a^6*r^2*sin(θ)^2*cos(θ)^4 + 4*a^6*r^2*cos(θ)^6 - 6*a^6*r*sin(θ)^2*cos(θ)^4 - 9*a^4*r^4*(1 - cos(4*θ))/8 + 6*a^4*r^4*cos(θ)^4 - 3*a^4*r^3*(1 - cos(4*θ))/2 - 5*a^4*r^2*(1 - cos(4*θ))/2 + 12*a^4*r^2*cos(θ)^4 - 3*a^2*r^6*sin(θ)^2 + 4*a^2*r^6*cos(θ)^2 - 6*a^2*r^5*sin(θ)^2 - 20*a^2*r^4*sin(θ)^2 + 24*a^2*r^4*cos(θ)^2 - 40*a^2*r^3*sin(θ)^2 + r^8 + 12*r^6)
	g[4, 2] = 2*a*(a^6*cos(θ)^6 + 3*a^4*r^2*cos(θ)^4 + 3*a^2*r^4*cos(θ)^2 + 4*a^2*r^2*cos(θ)^2 + r^6 + 4*r^4)/(-3*a^8*sin(θ)^2*cos(θ)^6 + a^8*cos(θ)^8 - 9*a^6*r^2*sin(θ)^2*cos(θ)^4 + 4*a^6*r^2*cos(θ)^6 - 6*a^6*r*sin(θ)^2*cos(θ)^4 - 9*a^4*r^4*(1 - cos(4*θ))/8 + 6*a^4*r^4*cos(θ)^4 - 3*a^4*r^3*(1 - cos(4*θ))/2 - 5*a^4*r^2*(1 - cos(4*θ))/2 + 12*a^4*r^2*cos(θ)^4 - 3*a^2*r^6*sin(θ)^2 + 4*a^2*r^6*cos(θ)^2 - 6*a^2*r^5*sin(θ)^2 - 20*a^2*r^4*sin(θ)^2 + 24*a^2*r^4*cos(θ)^2 - 40*a^2*r^3*sin(θ)^2 + r^8 + 12*r^6)
	g[4, 4] = (a^6*cos(θ)^6 + 3*a^4*r^2*cos(θ)^4 + 3*a^2*r^4*cos(θ)^2 + 12*a^2*r^2*cos(θ)^2 + r^6 + 12*r^4)/((-3*a^8*sin(θ)^2*cos(θ)^6 + a^8*cos(θ)^8 - 9*a^6*r^2*sin(θ)^2*cos(θ)^4 + 4*a^6*r^2*cos(θ)^6 - 6*a^6*r*sin(θ)^2*cos(θ)^4 - 9*a^4*r^4*(1 - cos(4*θ))/8 + 6*a^4*r^4*cos(θ)^4 - 3*a^4*r^3*(1 - cos(4*θ))/2 - 5*a^4*r^2*(1 - cos(4*θ))/2 + 12*a^4*r^2*cos(θ)^4 - 3*a^2*r^6*sin(θ)^2 + 4*a^2*r^6*cos(θ)^2 - 6*a^2*r^5*sin(θ)^2 - 20*a^2*r^4*sin(θ)^2 + 24*a^2*r^4*cos(θ)^2 - 40*a^2*r^3*sin(θ)^2 + r^8 + 12*r^6)*sin(θ)^2)
	return g
end
	
function Kerr_metric_con(r::Float64, θ::Float64, φ::Float64, a::Float64)
	g = zeros(4, 4)
	g[1, 1] = (1.0*M^4*a^4*c^4*cos(θ)^2 + 1.0*M^2*a^2*c^2*r^2*cos(θ)^2 + 1.0*M^2*a^2*c^2*r^2 + 2.0*M^2*a^2*c^2*r*sin(θ)^2 + 1.0*r^4)/(1.0*M^4*a^4*c^4*cos(θ)^2 + 1.0*M^2*a^2*c^2*r^2*cos(θ)^2 + 1.0*M^2*a^2*c^2*r^2 - 2.0*M^2*a^2*c^2*r*cos(θ)^2 + 1.0*r^4 - 2.0*r^3)
	g[1, 4] = 2.0*M*a*c*r/(1.0*M^4*a^4*c^4*cos(θ)^2 + 1.0*M^2*a^2*c^2*r^2*cos(θ)^2 + 1.0*M^2*a^2*c^2*r^2 - 2.0*M^2*a^2*c^2*r*cos(θ)^2 + 1.0*r^4 - 2.0*r^3)
	g[2, 2] = (-1.0*M^2*a^2*c^2 - 1.0*r^2 + 2.0*r)/(M^2*a^2*c^2*cos(θ)^2 + r^2)
	g[3, 3] = -1/(M^2*a^2*c^2*cos(θ)^2 + r^2)
	g[4, 1] = 2.0*M*a*c*r/(1.0*M^4*a^4*c^4*cos(θ)^2 + 1.0*M^2*a^2*c^2*r^2*cos(θ)^2 + 1.0*M^2*a^2*c^2*r^2 - 2.0*M^2*a^2*c^2*r*cos(θ)^2 + 1.0*r^4 - 2.0*r^3)
	g[4, 4] = (-1.0*M^2*a^2*c^2*cos(θ)^2 - 1.0*r^2 + 2.0*r)/((1.0*M^4*a^4*c^4*cos(θ)^2 + 1.0*M^2*a^2*c^2*r^2*cos(θ)^2 + 1.0*M^2*a^2*c^2*r^2 - 2.0*M^2*a^2*c^2*r*cos(θ)^2 + 1.0*r^4 - 2.0*r^3)*sin(θ)^2)
	return g
end
	
function Kerr_BL_metric_con(r::Float64, θ::Float64, φ::Float64, a::Float64)
	g = zeros(4, 4)
	g[1, 1] = (a^4*cos(θ)^2 + a^2*r^2*cos(θ)^2 + a^2*r^2 + 2*a^2*r*sin(θ)^2 + r^4)/(a^4*cos(θ)^2 + a^2*r^2*cos(θ)^2 + a^2*r^2 - 2*a^2*r*cos(θ)^2 + r^4 - 2*r^3)
	g[1, 4] = 2*a*r/(a^4*cos(θ)^2 + a^2*r^2*cos(θ)^2 + a^2*r^2 - 2*a^2*r*cos(θ)^2 + r^4 - 2*r^3)
	g[2, 2] = (-a^2 - r^2 + 2*r)/(a^2*cos(θ)^2 + r^2)
	g[3, 3] = -1/(a^2*cos(θ)^2 + r^2)
	g[4, 1] = 2*a*r/(a^4*cos(θ)^2 + a^2*r^2*cos(θ)^2 + a^2*r^2 - 2*a^2*r*cos(θ)^2 + r^4 - 2*r^3)
	g[4, 4] = (-a^2*cos(θ)^2 - r^2 + 2*r)/((a^4*cos(θ)^2 + a^2*r^2*cos(θ)^2 + a^2*r^2 - 2*a^2*r*cos(θ)^2 + r^4 - 2*r^3)*sin(θ)^2)
	return g
end
	
function Schwarzschild_metric_concov(r::Float64, θ::Float64, φ::Float64)
	g = zeros(4, 4)
	g[1, 1] = 1
	g[2, 2] = 1
	g[3, 3] = 1
	g[4, 4] = 1
	return g
end
	
function Kerr_Schild_metric_concov(r::Float64, θ::Float64, φ::Float64, a::Float64)
	g = zeros(4, 4)
	g[1, 1] = 1
	g[2, 2] = 1
	g[3, 3] = 1
	g[4, 4] = 1
	return g
end
	
function Kerr_metric_concov(r::Float64, θ::Float64, φ::Float64, a::Float64)
	g = zeros(4, 4)
	g[1, 1] = 1.00000000000000
	g[2, 2] = 1.00000000000000
	g[3, 3] = 1
	g[4, 4] = 1.00000000000000
	return g
end
	
function Kerr_BL_metric_concov(r::Float64, θ::Float64, φ::Float64, a::Float64)
	g = zeros(4, 4)
	g[1, 1] = 1
	g[2, 2] = 1
	g[3, 3] = 1
	g[4, 4] = 1
	return g
end
	
function Schwarzschild_Christoffel_Symbols(i::Int64, j::Int64, k::Int64, r::Float64, θ::Float64, φ::Float64)
	Γ = zeros(4, 4, 4)
	if i == 1 && j == 1 && k ==2
    		Γ[1, 1, 2] = 1.0/(r*(r - 2.0))
	elseif i == 1 && j == 2 && k == 1
    		Γ[1, 2, 1] = 1.0/(r*(r - 2.0))
	elseif i == 2 && j == 1 && k == 1
    		Γ[2, 1, 1] = 1.0*(r - 2.0)/r^3
	elseif i == 2 && j == 2 && k == 2
    		Γ[2, 2, 2] = -0.5/(r*(0.5*r - 1.0))
	elseif i == 2 && j == 3 && k == 3
    		Γ[2, 3, 3] = 2.0 - r
	elseif i == 2 && j == 4 && k == 4
    		Γ[2, 4, 4] = (2.0 - r)*sin(θ)^2
	elseif i == 3 && j == 2 && k == 3
    		Γ[3, 2, 3] = 1/r
	elseif i == 3 && j == 3 && k == 2
    		Γ[3, 3, 2] = 1/r
	elseif i == 3 && j == 4 && k == 4
    		Γ[3, 4, 4] = -sin(2*θ)/2
	elseif i == 4 && j == 2 && k == 4
    		Γ[4, 2, 4] = 1/r
	elseif i == 4 && j == 3 && k == 4
    		Γ[4, 3, 4] = 1/tan(θ)
	elseif i == 4 && j == 4 && k == 2
    		Γ[4, 4, 2] = 1/r
	elseif i == 4 && j == 4 && k == 3
    		Γ[4, 4, 3] = 1/tan(θ)
	end
	return Γ[i, j, k]
end
	
function Kerr_Schild_Christoffel_Symbols(i::Int64, j::Int64, k::Int64, r::Float64, θ::Float64, φ::Float64, a::Float64)
	Γ = zeros(4, 4, 4)
	if i == 1 && j == 1 && k ==1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[1, 1, 1] = -4*r*(a_2*cos_θ_2 - r_2)*(-a_6*sine_θ_2*cos_θ_4 + a_6*cos_θ_6 + a_4*r_2*(cos_4θ - 1)/4 + 3*a_4*r_2*cos_θ_4 + a_4*r*(cos_4θ - 1)/4 - a_2*r_4*sine_θ_2 + 3*a_2*r_4*cos_θ_2 - 2*a_2*r_3*sine_θ_2 + r_6)/((a_2*cos_θ_2 + r_2)^2*(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(cos_4θ - 1)/8 + 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(cos_4θ - 1)/2 + 5*a_4*r_2*(cos_4θ - 1)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6))
	elseif i == 1 && j == 1 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_7 = r^7
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[1, 1, 2] = (3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ_8 + 3*a_6*r_2*sine_θ_2*cos_θ_4 - 2*a_6*r_2*cos_θ_6 + 4*a_6*r*sine_θ_2*cos_θ_4 - 2*a_6*r*cos_θ_6 - 3*a_4*r_4*(1 - cos_4θ)/8 - 2*a_4*r_3*cos_θ_4 - a_4*r_2*(1 - cos_4θ)/2 - 3*a_2*r_6*sine_θ_2 + 2*a_2*r_6*cos_θ_2 - 4*a_2*r_5*sine_θ_2 + 2*a_2*r_5*cos_θ_2 + 4*a_2*r_4*sine_θ_2 + r_8 + 2*r_7)/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r_7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 1 && j == 1 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		sine_θ_4 = sine_θ^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[1, 1, 3] = 2*a_2*r*(3*a_6*sine_θ^6 - 6*a_6*sine_θ_4 + 3*a_6*sine_θ_2 - a_6*cos_θ_6 - 9*a_4*r_2*sine_θ_4 + 12*a_4*r_2*sine_θ_2 - 3*a_4*r_2 + 2*a_4*r*sine_θ_4 - 2*a_4*r + 6*a_2*r_4*sine_θ_2 - 3*a_2*r_4 - 4*a_2*r_3 - 4*a_2*r_2*sine_θ_2 - 16*a_2*r_2 - r_6 - 2*r_5 - 16*r_4)*sine_θ*cos_θ/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r^7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 1 && j == 1 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[1, 1, 4] = 8*a*r*(a_2*cos_θ_2 - r_2)*(-a_6*sine_θ_2*cos_θ_4 + a_6*cos_θ_6 + a_4*r_2*(cos_4θ - 1)/4 + 3*a_4*r_2*cos_θ_4 + a_4*r*(cos_4θ - 1)/4 - a_2*r_4*sine_θ_2 + 3*a_2*r_4*cos_θ_2 - 2*a_2*r_3*sine_θ_2 + r_6)*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(cos_4θ - 1)/8 + 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(cos_4θ - 1)/2 + 5*a_4*r_2*(cos_4θ - 1)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6))
	elseif i == 1 && j == 2 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_7 = r^7
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[1, 2, 1] = (3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ_8 + 3*a_6*r_2*sine_θ_2*cos_θ_4 - 2*a_6*r_2*cos_θ_6 + 4*a_6*r*sine_θ_2*cos_θ_4 - 2*a_6*r*cos_θ_6 - 3*a_4*r_4*(1 - cos_4θ)/8 - 2*a_4*r_3*cos_θ_4 - a_4*r_2*(1 - cos_4θ)/2 - 3*a_2*r_6*sine_θ_2 + 2*a_2*r_6*cos_θ_2 - 4*a_2*r_5*sine_θ_2 + 2*a_2*r_5*cos_θ_2 + 4*a_2*r_4*sine_θ_2 + r_8 + 2*r_7)/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r_7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 1 && j == 2 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_7 = r^7
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[1, 2, 2] = 4*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ_8 + 3*a_6*r_2*sine_θ_2*cos_θ_4 - 2*a_6*r_2*cos_θ_6 + 7*a_6*r*sine_θ_2*cos_θ_4 - a_6*r*cos_θ_6 - 3*a_4*r_4*(1 - cos_4θ)/8 - a_4*r_3*cos_θ_4 + a_4*r_2*(1 - cos_4θ)/4 - 3*a_2*r_6*sine_θ_2 + 2*a_2*r_6*cos_θ_2 - 7*a_2*r_5*sine_θ_2 + a_2*r_5*cos_θ_2 - 2*a_2*r_4*sine_θ_2 + r_8 + r_7)/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r_7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 1 && j == 2 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[1, 2, 3] = 4*a_2*r*(-6*a_4*sine_θ^4 + 9*a_4*sine_θ_2 - 3*a_4 + 9*a_2*r_2*sine_θ_2 - 6*a_2*r_2 + 14*a_2*r*sine_θ_2 - 8*a_2*r - 3*r_4 - 8*r_3 - 8*r_2)*sine_θ*cos_θ/(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 3*a_4*r_3*(1 - cos_4θ)/2 - 5*a_4*r_2*(1 - cos_4θ)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6)
	elseif i == 1 && j == 2 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_7 = r^7
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[1, 2, 4] = 2*a*(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ_8 - 3*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 10*a_6*r*sine_θ_2*cos_θ_4 + 2*a_6*r*cos_θ_6 + 3*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 + 6*a_4*r_3*cos_θ_4 - a_4*r_2*(1 - cos_4θ) + 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 + 10*a_2*r_5*sine_θ_2 + 6*a_2*r_5*cos_θ_2 + 8*a_2*r_4*sine_θ_2 + r_8 + 2*r_7)*sine_θ_2/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r_7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 1 && j == 3 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		sine_θ_4 = sine_θ^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[1, 3, 1] = 2*a_2*r*(3*a_6*sine_θ^6 - 6*a_6*sine_θ_4 + 3*a_6*sine_θ_2 - a_6*cos_θ_6 - 9*a_4*r_2*sine_θ_4 + 12*a_4*r_2*sine_θ_2 - 3*a_4*r_2 + 2*a_4*r*sine_θ_4 - 2*a_4*r + 6*a_2*r_4*sine_θ_2 - 3*a_2*r_4 - 4*a_2*r_3 - 4*a_2*r_2*sine_θ_2 - 16*a_2*r_2 - r_6 - 2*r_5 - 16*r_4)*sine_θ*cos_θ/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r^7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 1 && j == 3 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[1, 3, 2] = 4*a_2*r*(-6*a_4*sine_θ^4 + 9*a_4*sine_θ_2 - 3*a_4 + 9*a_2*r_2*sine_θ_2 - 6*a_2*r_2 + 14*a_2*r*sine_θ_2 - 8*a_2*r - 3*r_4 - 8*r_3 - 8*r_2)*sine_θ*cos_θ/(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 3*a_4*r_3*(1 - cos_4θ)/2 - 5*a_4*r_2*(1 - cos_4θ)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6)
	elseif i == 1 && j == 3 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[1, 3, 3] = 4*r_2*(a_6*sine_θ_2*cos_θ_4 - a_6*cos_θ_6 + a_4*r_2*(1 - cos_4θ)/4 - 3*a_4*r_2*cos_θ_4 + a_4*r*(1 - cos_4θ)/4 + a_2*r_4*sine_θ_2 - 3*a_2*r_4*cos_θ_2 + 2*a_2*r_3*sine_θ_2 - r_6)/(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 3*a_4*r_3*(1 - cos_4θ)/2 - 5*a_4*r_2*(1 - cos_4θ)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6)
	elseif i == 1 && j == 3 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[1, 3, 4] = 4*a^3*r*(4*a_4*sine_θ^4 - 5*a_4*sine_θ_2 + a_4 - 5*a_2*r_2*sine_θ_2 + 2*a_2*r_2 - 6*a_2*r*sine_θ_2 + r_4 + 4*r_2)*sine_θ^3*cos_θ/(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 3*a_4*r_3*(1 - cos_4θ)/2 - 5*a_4*r_2*(1 - cos_4θ)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6)
	elseif i == 1 && j == 4 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[1, 4, 1] = 8*a*r*(a_2*cos_θ_2 - r_2)*(-a_6*sine_θ_2*cos_θ_4 + a_6*cos_θ_6 + a_4*r_2*(cos_4θ - 1)/4 + 3*a_4*r_2*cos_θ_4 + a_4*r*(cos_4θ - 1)/4 - a_2*r_4*sine_θ_2 + 3*a_2*r_4*cos_θ_2 - 2*a_2*r_3*sine_θ_2 + r_6)*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(cos_4θ - 1)/8 + 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(cos_4θ - 1)/2 + 5*a_4*r_2*(cos_4θ - 1)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6))
	elseif i == 1 && j == 4 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_7 = r^7
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[1, 4, 2] = 2*a*(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ_8 - 3*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 10*a_6*r*sine_θ_2*cos_θ_4 + 2*a_6*r*cos_θ_6 + 3*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 + 6*a_4*r_3*cos_θ_4 - a_4*r_2*(1 - cos_4θ) + 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 + 10*a_2*r_5*sine_θ_2 + 6*a_2*r_5*cos_θ_2 + 8*a_2*r_4*sine_θ_2 + r_8 + 2*r_7)*sine_θ_2/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r_7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 1 && j == 4 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[1, 4, 3] = 4*a^3*r*(4*a_4*sine_θ^4 - 5*a_4*sine_θ_2 + a_4 - 5*a_2*r_2*sine_θ_2 + 2*a_2*r_2 - 6*a_2*r*sine_θ_2 + r_4 + 4*r_2)*sine_θ^3*cos_θ/(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 3*a_4*r_3*(1 - cos_4θ)/2 - 5*a_4*r_2*(1 - cos_4θ)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6)
	elseif i == 1 && j == 4 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[1, 4, 4] = 4*r*(r*(a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2) - (a_2*cos_θ_2 + r_2)*(a_2*(r + 1)*sine_θ_2 + 2*r*(a_2*cos_θ_2 + r_2)))*(a_6*sine_θ_2*cos_θ_4 - a_6*cos_θ_6 + a_4*r_2*(1 - cos_4θ)/4 - 3*a_4*r_2*cos_θ_4 + a_4*r*(1 - cos_4θ)/4 + a_2*r_4*sine_θ_2 - 3*a_2*r_4*cos_θ_2 + 2*a_2*r_3*sine_θ_2 - r_6)*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 2 && j == 1 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
    		Γ[2, 1, 1] = (a_2*cos_θ_2 - r_2)*(a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ_8 + 3*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 2*a_6*r*cos_θ_6 + 3*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 6*a_4*r_3*cos_θ_4 + 3*a_4*r_2*(1 - cos_4θ)/2 + a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r_5*cos_θ_2 + 12*a_2*r_4*sine_θ_2 + r_8 - 2*r^7)/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ_8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r_5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r_8 - 12*r_6))
	elseif i == 2 && j == 1 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[2, 1, 2] = 4*(-a_2*(-a_2*cos_θ_2 + r_2)*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4)*sine_θ_2 + r*(a_2*cos_θ_2 - r_2)*(a_6*sine_θ_2*cos_θ_4 - a_6*cos_θ_6 + a_4*r_2*(1 - cos_4θ)/4 - 3*a_4*r_2*cos_θ_4 + a_4*r*(1 - cos_4θ)/4 + a_2*r_4*sine_θ_2 - 3*a_2*r_4*cos_θ_2 + 2*a_2*r_3*sine_θ_2 - r_6))/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 2 && j == 1 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[2, 1, 3] = 4*a_2*r*(a_4*sine_θ_2 - a_4 + a_2*r_2*sine_θ_2 - 2*a_2*r_2 - 2*a_2*r*sine_θ_2 - r_4 - 8*r_2)*sine_θ*cos_θ/(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 3*a_4*r_3*(1 - cos_4θ)/2 - 5*a_4*r_2*(1 - cos_4θ)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6)
	elseif i == 2 && j == 1 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
    		Γ[2, 1, 4] = 2*a*(-a_2*cos_θ_2 + r_2)*(a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ_8 + 3*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 2*a_6*r*cos_θ_6 + 3*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 6*a_4*r_3*cos_θ_4 + 3*a_4*r_2*(1 - cos_4θ)/2 + a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r_5*cos_θ_2 + 12*a_2*r_4*sine_θ_2 + r_8 - 2*r^7)*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ_8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r_5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r_8 - 12*r_6))
	elseif i == 2 && j == 2 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[2, 2, 1] = 4*(-a_2*(-a_2*cos_θ_2 + r_2)*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4)*sine_θ_2 + r*(a_2*cos_θ_2 - r_2)*(a_6*sine_θ_2*cos_θ_4 - a_6*cos_θ_6 + a_4*r_2*(1 - cos_4θ)/4 - 3*a_4*r_2*cos_θ_4 + a_4*r*(1 - cos_4θ)/4 + a_2*r_4*sine_θ_2 - 3*a_2*r_4*cos_θ_2 + 2*a_2*r_3*sine_θ_2 - r_6))/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 2 && j == 2 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_7 = r^7
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[2, 2, 2] = (-7*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ_8 - 7*a_6*r_2*sine_θ_2*cos_θ_4 + 2*a_6*r_2*cos_θ_6 - 16*a_6*r*sine_θ_2*cos_θ_4 + 14*a_6*r*cos_θ_6 + 7*a_4*r_4*(1 - cos_4θ)/8 + 14*a_4*r_3*cos_θ_4 - 13*a_4*r_2*(1 - cos_4θ)/2 + 7*a_2*r_6*sine_θ_2 - 2*a_2*r_6*cos_θ_2 + 16*a_2*r_5*sine_θ_2 - 14*a_2*r_5*cos_θ_2 + 52*a_2*r_4*sine_θ_2 - r_8 - 14*r_7)/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r_7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 2 && j == 2 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		sine_θ_4 = sine_θ^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_7 = r^7
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[2, 2, 3] = 2*a_2*(-2*a_8*cos_θ_8 - 8*a_6*r_2*cos_θ_6 - 3*a_6*r*sine_θ^6 + 6*a_6*r*sine_θ_4 - 3*a_6*r*sine_θ_2 - 3*a_6*r*cos_θ_6 - 12*a_4*r_4*sine_θ_4 + 24*a_4*r_4*sine_θ_2 - 12*a_4*r_4 - 3*a_4*r_3*sine_θ_4 + 12*a_4*r_3*sine_θ_2 - 9*a_4*r_3 + 6*a_4*r_2*sine_θ_4 - 4*a_4*r_2*sine_θ_2 - 2*a_4*r_2 + 8*a_2*r_6*sine_θ_2 - 8*a_2*r_6 + 6*a_2*r_5*sine_θ_2 - 9*a_2*r_5 - 4*a_2*r_4*sine_θ_2 - 4*a_2*r_4 - 4*a_2*r_3*sine_θ_2 - 16*a_2*r_3 - 2*r_8 - 3*r_7 - 2*r_6 - 16*r_5)*sine_θ*cos_θ/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r_7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 2 && j == 2 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[2, 2, 4] = 2*a*(4*r*(-a_2*cos_θ_2 + r_2)*(a_6*sine_θ_2*cos_θ_4 - a_6*cos_θ_6 + a_4*r_2*(1 - cos_4θ)/4 - 3*a_4*r_2*cos_θ_4 + a_4*r*(1 - cos_4θ)/4 + a_2*r_4*sine_θ_2 - 3*a_2*r_4*cos_θ_2 + 2*a_2*r_3*sine_θ_2 - r_6) + (r*(a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2) + (a_2*cos_θ_2 + r_2)*(a_2*(-r - 1)*sine_θ_2 - 2*r*(a_2*cos_θ_2 + r_2)))*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4))*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 2 && j == 3 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[2, 3, 1] = 4*a_2*r*(a_4*sine_θ_2 - a_4 + a_2*r_2*sine_θ_2 - 2*a_2*r_2 - 2*a_2*r*sine_θ_2 - r_4 - 8*r_2)*sine_θ*cos_θ/(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 3*a_4*r_3*(1 - cos_4θ)/2 - 5*a_4*r_2*(1 - cos_4θ)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6)
	elseif i == 2 && j == 3 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		sine_θ_4 = sine_θ^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_7 = r^7
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[2, 3, 2] = 2*a_2*(-2*a_8*cos_θ_8 - 8*a_6*r_2*cos_θ_6 - 3*a_6*r*sine_θ^6 + 6*a_6*r*sine_θ_4 - 3*a_6*r*sine_θ_2 - 3*a_6*r*cos_θ_6 - 12*a_4*r_4*sine_θ_4 + 24*a_4*r_4*sine_θ_2 - 12*a_4*r_4 - 3*a_4*r_3*sine_θ_4 + 12*a_4*r_3*sine_θ_2 - 9*a_4*r_3 + 6*a_4*r_2*sine_θ_4 - 4*a_4*r_2*sine_θ_2 - 2*a_4*r_2 + 8*a_2*r_6*sine_θ_2 - 8*a_2*r_6 + 6*a_2*r_5*sine_θ_2 - 9*a_2*r_5 - 4*a_2*r_4*sine_θ_2 - 4*a_2*r_4 - 4*a_2*r_3*sine_θ_2 - 16*a_2*r_3 - 2*r_8 - 3*r_7 - 2*r_6 - 16*r_5)*sine_θ*cos_θ/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r_7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 2 && j == 3 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
    		Γ[2, 3, 3] = r*(a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ_8 + 3*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 2*a_6*r*cos_θ_6 + 3*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 6*a_4*r_3*cos_θ_4 + 3*a_4*r_2*(1 - cos_4θ)/2 + a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r_5*cos_θ_2 + 12*a_2*r_4*sine_θ_2 + r_8 - 2*r^7)/(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ_8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r_5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r_8 - 12*r_6)
	elseif i == 2 && j == 3 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[2, 3, 4] = 8*a^3*r_2*(3*a_2*sine_θ_2 - a_2 - r_2 + 2*r)*sine_θ^3*cos_θ/(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 3*a_4*r_3*(1 - cos_4θ)/2 - 5*a_4*r_2*(1 - cos_4θ)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6)
	elseif i == 2 && j == 4 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
    		Γ[2, 4, 1] = 2*a*(-a_2*cos_θ_2 + r_2)*(a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ_8 + 3*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 2*a_6*r*cos_θ_6 + 3*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 6*a_4*r_3*cos_θ_4 + 3*a_4*r_2*(1 - cos_4θ)/2 + a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r_5*cos_θ_2 + 12*a_2*r_4*sine_θ_2 + r_8 - 2*r^7)*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ_8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r_5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r_8 - 12*r_6))
	elseif i == 2 && j == 4 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[2, 4, 2] = 2*a*(4*r*(-a_2*cos_θ_2 + r_2)*(a_6*sine_θ_2*cos_θ_4 - a_6*cos_θ_6 + a_4*r_2*(1 - cos_4θ)/4 - 3*a_4*r_2*cos_θ_4 + a_4*r*(1 - cos_4θ)/4 + a_2*r_4*sine_θ_2 - 3*a_2*r_4*cos_θ_2 + 2*a_2*r_3*sine_θ_2 - r_6) + (r*(a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2) + (a_2*cos_θ_2 + r_2)*(a_2*(-r - 1)*sine_θ_2 - 2*r*(a_2*cos_θ_2 + r_2)))*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4))*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 2 && j == 4 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[2, 4, 3] = 8*a^3*r_2*(3*a_2*sine_θ_2 - a_2 - r_2 + 2*r)*sine_θ^3*cos_θ/(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ)/8 + 6*a_4*r_4*cos_θ_4 - 3*a_4*r_3*(1 - cos_4θ)/2 - 5*a_4*r_2*(1 - cos_4θ)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6)
	elseif i == 2 && j == 4 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
    		Γ[2, 4, 4] = (r*(a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2) - (a_2*cos_θ_2 + r_2)*(a_2*(r + 1)*sine_θ_2 + 2*r*(a_2*cos_θ_2 + r_2)))*(a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ_8 + 3*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 2*a_6*r*cos_θ_6 - 3*a_4*r_4*(cos_4θ - 1)/8 + 6*a_4*r_4*cos_θ_4 - 6*a_4*r_3*cos_θ_4 - 3*a_4*r_2*(cos_4θ - 1)/2 + a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r_5*cos_θ_2 + 12*a_2*r_4*sine_θ_2 + r_8 - 2*r^7)*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ_8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(cos_4θ - 1)/8 + 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(cos_4θ - 1)/2 + 5*a_4*r_2*(cos_4θ - 1)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r_5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r_8 + 12*r_6))
	elseif i == 3 && j == 1 && k == 1
		a_2 = a^2
    		Γ[3, 1, 1] = -8*a_2*r*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)^3
	elseif i == 3 && j == 1 && k == 2
		a_2 = a^2
    		Γ[3, 1, 2] = -16*a_2*r*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)^3
	elseif i == 3 && j == 1 && k == 4
		r_2 = r^2
		a_2 = a^2
    		Γ[3, 1, 4] = 16*a*r*(a_2 + r_2)*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r_2)^3
	elseif i == 3 && j == 2 && k == 1
		a_2 = a^2
    		Γ[3, 2, 1] = -16*a_2*r*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)^3
	elseif i == 3 && j == 2 && k == 2
		a_2 = a^2
    		Γ[3, 2, 2] = -8*a_2*r*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)^3
	elseif i == 3 && j == 2 && k == 3
    		Γ[3, 2, 3] = r/(a^2*cos(θ)^2 + r^2)
	elseif i == 3 && j == 2 && k == 4
		cos_θ = cos(θ)
		r_2 = r^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		a_4 = a^4
    		Γ[3, 2, 4] = 2*a*(a_4*(1 - cos_θ_2)^2 + 2*a_4*cos_θ_2 - a_4 + 2*a_2*r_2*cos_θ_2 + 2*a_2*r + r^4 + 2*r^3)*sin(θ)*cos_θ/(a_2*cos_θ_2 + r_2)^3
	elseif i == 3 && j == 3 && k == 2
    		Γ[3, 3, 2] = r/(a^2*cos(θ)^2 + r^2)
	elseif i == 3 && j == 3 && k == 3
		a_2 = a^2
    		Γ[3, 3, 3] = -a_2*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)
	elseif i == 3 && j == 4 && k == 1
		r_2 = r^2
		a_2 = a^2
    		Γ[3, 4, 1] = 16*a*r*(a_2 + r_2)*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r_2)^3
	elseif i == 3 && j == 4 && k == 2
		cos_θ = cos(θ)
		r_2 = r^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		a_4 = a^4
    		Γ[3, 4, 2] = 2*a*(a_4*(1 - cos_θ_2)^2 + 2*a_4*cos_θ_2 - a_4 + 2*a_2*r_2*cos_θ_2 + 2*a_2*r + r^4 + 2*r^3)*sin(θ)*cos_θ/(a_2*cos_θ_2 + r_2)^3
	elseif i == 3 && j == 4 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		a_4 = a^4
    		Γ[3, 4, 4] = (-a_2*(a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2)*sine_θ_2 + (a_2*cos_θ_2 + r_2)*(-2*a_4*cos_θ_2 + a_4 - 2*a_2*r_2*cos_θ_2 + 4*a_2*r*cos_θ_2 - 4*a_2*r - r^4))*sine_θ*cos_θ/(a_2*cos_θ_2 + r_2)^3
	elseif i == 4 && j == 1 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 1, 1] = 2*a*(a_2*cos_θ_2 - r_2)*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4)/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 4 && j == 1 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 1, 2] = 2*a*(-2*r*(a_6*cos_θ_6 + a_4*r_2*cos_θ_4 + 2*a_4*r*cos_θ_4 - a_2*r_4*cos_θ_2 - r_6 - 2*r_5) - (-a_2*cos_θ_2 + r_2)*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 12*a_2*r_2*cos_θ_2 + r_6 + 12*r_4))/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r_5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 4 && j == 1 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 1, 3] = 32*a*r*(-2*a_4*sine_θ^4 + 3*a_4*sine_θ_2 - a_4 + 3*a_2*r_2*sine_θ_2 - 2*a_2*r_2 + 2*a_2*r*sine_θ_2 - r_4 - 12*r_2)/((-24*a_8*sine_θ_2*cos_θ_6 + 8*a_8*cos_θ^8 - 72*a_6*r_2*sine_θ_2*cos_θ_4 + 32*a_6*r_2*cos_θ_6 - 48*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ) + 48*a_4*r_4*cos_θ_4 - 12*a_4*r_3*(1 - cos_4θ) - 20*a_4*r_2*(1 - cos_4θ) + 96*a_4*r_2*cos_θ_4 - 24*a_2*r_6*sine_θ_2 + 32*a_2*r_6*cos_θ_2 - 48*a_2*r^5*sine_θ_2 - 160*a_2*r_4*sine_θ_2 + 192*a_2*r_4*cos_θ_2 - 320*a_2*r_3*sine_θ_2 + 8*r^8 + 96*r_6)*tan(θ))
	elseif i == 4 && j == 1 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 1, 4] = 4*a_2*(-a_2*cos_θ_2 + r_2)*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4)*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 4 && j == 2 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 2, 1] = 2*a*(-2*r*(a_6*cos_θ_6 + a_4*r_2*cos_θ_4 + 2*a_4*r*cos_θ_4 - a_2*r_4*cos_θ_2 - r_6 - 2*r_5) - (-a_2*cos_θ_2 + r_2)*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 12*a_2*r_2*cos_θ_2 + r_6 + 12*r_4))/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r_5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 4 && j == 2 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		r_8 = r^8
		cos_θ_8 = cos_θ^8
		a_8 = a^8
		a_10 = a^10
    		Γ[4, 2, 2] = 2*a*(-a_6*cos_θ_6 - a_4*r_2*cos_θ_4 + 8*a_4*r*cos_θ_4 + a_2*r_4*cos_θ_2 - 4*a_2*r_2*cos_θ_2 + r_6 - 8*r_5 + 4*r_4)/(-3*a_10*sine_θ_2*cos_θ_8 + a_10*cos_θ^10 - 12*a_8*r_2*sine_θ_2*cos_θ_6 + 5*a_8*r_2*cos_θ_8 - 6*a_8*r*sine_θ_2*cos_θ_6 - 18*a_6*r_4*sine_θ_2*cos_θ_4 + 10*a_6*r_4*cos_θ_6 - 18*a_6*r_3*sine_θ_2*cos_θ_4 - 20*a_6*r_2*sine_θ_2*cos_θ_4 + 12*a_6*r_2*cos_θ_6 - 3*a_4*r_6*(1 - cos_4θ)/2 + 10*a_4*r_6*cos_θ_4 - 9*a_4*r_5*(1 - cos_4θ)/4 - 5*a_4*r_4*(1 - cos_4θ) + 36*a_4*r_4*cos_θ_4 - 5*a_4*r_3*(1 - cos_4θ) - 3*a_2*r_8*sine_θ_2 + 5*a_2*r_8*cos_θ_2 - 6*a_2*r^7*sine_θ_2 - 20*a_2*r_6*sine_θ_2 + 36*a_2*r_6*cos_θ_2 - 40*a_2*r_5*sine_θ_2 + r^10 + 12*r_8)
	elseif i == 4 && j == 2 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 2, 3] = 16*a*(-a_6*cos_θ_6 - 3*a_4*r_2*cos_θ_4 - 2*a_4*r*cos_θ_4 - 3*a_2*r_4*cos_θ_2 - 4*a_2*r_3*cos_θ_2 + 8*a_2*r_2*sine_θ_2 - 12*a_2*r_2*cos_θ_2 - r_6 - 2*r_5 - 12*r_4 - 24*r_3)/((-24*a_8*sine_θ_2*cos_θ_6 + 8*a_8*cos_θ^8 - 72*a_6*r_2*sine_θ_2*cos_θ_4 + 32*a_6*r_2*cos_θ_6 - 48*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ) + 48*a_4*r_4*cos_θ_4 - 12*a_4*r_3*(1 - cos_4θ) - 20*a_4*r_2*(1 - cos_4θ) + 96*a_4*r_2*cos_θ_4 - 24*a_2*r_6*sine_θ_2 + 32*a_2*r_6*cos_θ_2 - 48*a_2*r_5*sine_θ_2 - 160*a_2*r_4*sine_θ_2 + 192*a_2*r_4*cos_θ_2 - 320*a_2*r_3*sine_θ_2 + 8*r^8 + 96*r_6)*tan(θ))
	elseif i == 4 && j == 2 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 2, 4] = (-8*a_2*r*(-a_2*cos_θ_2 + r_2)*(a_4*cos_θ_4 + 2*a_2*r_2*cos_θ_2 + 2*a_2*r*cos_θ_2 + r_4 + 2*r_3)*sine_θ_2 + (r*(a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2) - (a_2*cos_θ_2 + r_2)*(a_2*(r + 1)*sine_θ_2 + 2*r*(a_2*cos_θ_2 + r_2)))*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 12*a_2*r_2*cos_θ_2 + r_6 + 12*r_4))/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 4 && j == 3 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 3, 1] = 32*a*r*(-2*a_4*sine_θ^4 + 3*a_4*sine_θ_2 - a_4 + 3*a_2*r_2*sine_θ_2 - 2*a_2*r_2 + 2*a_2*r*sine_θ_2 - r_4 - 12*r_2)/((-24*a_8*sine_θ_2*cos_θ_6 + 8*a_8*cos_θ^8 - 72*a_6*r_2*sine_θ_2*cos_θ_4 + 32*a_6*r_2*cos_θ_6 - 48*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ) + 48*a_4*r_4*cos_θ_4 - 12*a_4*r_3*(1 - cos_4θ) - 20*a_4*r_2*(1 - cos_4θ) + 96*a_4*r_2*cos_θ_4 - 24*a_2*r_6*sine_θ_2 + 32*a_2*r_6*cos_θ_2 - 48*a_2*r^5*sine_θ_2 - 160*a_2*r_4*sine_θ_2 + 192*a_2*r_4*cos_θ_2 - 320*a_2*r_3*sine_θ_2 + 8*r^8 + 96*r_6)*tan(θ))
	elseif i == 4 && j == 3 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_5 = r^5
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 3, 2] = 16*a*(-a_6*cos_θ_6 - 3*a_4*r_2*cos_θ_4 - 2*a_4*r*cos_θ_4 - 3*a_2*r_4*cos_θ_2 - 4*a_2*r_3*cos_θ_2 + 8*a_2*r_2*sine_θ_2 - 12*a_2*r_2*cos_θ_2 - r_6 - 2*r_5 - 12*r_4 - 24*r_3)/((-24*a_8*sine_θ_2*cos_θ_6 + 8*a_8*cos_θ^8 - 72*a_6*r_2*sine_θ_2*cos_θ_4 + 32*a_6*r_2*cos_θ_6 - 48*a_6*r*sine_θ_2*cos_θ_4 - 9*a_4*r_4*(1 - cos_4θ) + 48*a_4*r_4*cos_θ_4 - 12*a_4*r_3*(1 - cos_4θ) - 20*a_4*r_2*(1 - cos_4θ) + 96*a_4*r_2*cos_θ_4 - 24*a_2*r_6*sine_θ_2 + 32*a_2*r_6*cos_θ_2 - 48*a_2*r_5*sine_θ_2 - 160*a_2*r_4*sine_θ_2 + 192*a_2*r_4*cos_θ_2 - 320*a_2*r_3*sine_θ_2 + 8*r^8 + 96*r_6)*tan(θ))
	elseif i == 4 && j == 3 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 3, 3] = 2*a*r*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4)/(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6)
	elseif i == 4 && j == 3 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 3, 4] = (32*a_2*(4*r_2*(a_2 + r_2)*(a_4*cos_θ_4 + 2*a_2*r_2*cos_θ_2 + 2*a_2*r*cos_θ_2 + r_4 + 2*r_3) + (a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4)*(a_4*(1 - cos_θ_2)^2 + 2*a_4*cos_θ_2 - a_4 + 2*a_2*r_2*cos_θ_2 + 2*a_2*r + r_4 + 2*r_3))*sine_θ_2 - 8*(a_2*(a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2)*sine_θ_2 + (a_2*cos_θ_2 + r_2)*(2*a_4*cos_θ_2 - a_4 + 2*a_2*r_2*cos_θ_2 - 4*a_2*r*cos_θ_2 + 4*a_2*r + r_4))*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 12*a_2*r_2*cos_θ_2 + r_6 + 12*r_4))/((a_2*cos_θ_2 + r_2)^2*(24*a_8*sine_θ_2*cos_θ_6 - 8*a_8*cos_θ^8 + 72*a_6*r_2*sine_θ_2*cos_θ_4 - 32*a_6*r_2*cos_θ_6 + 48*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ) - 48*a_4*r_4*cos_θ_4 + 12*a_4*r_3*(1 - cos_4θ) + 20*a_4*r_2*(1 - cos_4θ) - 96*a_4*r_2*cos_θ_4 + 24*a_2*r_6*sine_θ_2 - 32*a_2*r_6*cos_θ_2 + 48*a_2*r^5*sine_θ_2 + 160*a_2*r_4*sine_θ_2 - 192*a_2*r_4*cos_θ_2 + 320*a_2*r_3*sine_θ_2 - 8*r^8 - 96*r_6)*tan(θ))
	elseif i == 4 && j == 4 && k == 1
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 4, 1] = 4*a_2*(-a_2*cos_θ_2 + r_2)*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4)*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 4 && j == 4 && k == 2
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 4, 2] = (-8*a_2*r*(-a_2*cos_θ_2 + r_2)*(a_4*cos_θ_4 + 2*a_2*r_2*cos_θ_2 + 2*a_2*r*cos_θ_2 + r_4 + 2*r_3)*sine_θ_2 + (r*(a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2) - (a_2*cos_θ_2 + r_2)*(a_2*(r + 1)*sine_θ_2 + 2*r*(a_2*cos_θ_2 + r_2)))*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 12*a_2*r_2*cos_θ_2 + r_6 + 12*r_4))/((a_2*cos_θ_2 + r_2)^2*(3*a_8*sine_θ_2*cos_θ_6 - a_8*cos_θ^8 + 9*a_6*r_2*sine_θ_2*cos_θ_4 - 4*a_6*r_2*cos_θ_6 + 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ)/8 - 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(1 - cos_4θ)/2 + 5*a_4*r_2*(1 - cos_4θ)/2 - 12*a_4*r_2*cos_θ_4 + 3*a_2*r_6*sine_θ_2 - 4*a_2*r_6*cos_θ_2 + 6*a_2*r^5*sine_θ_2 + 20*a_2*r_4*sine_θ_2 - 24*a_2*r_4*cos_θ_2 + 40*a_2*r_3*sine_θ_2 - r^8 - 12*r_6))
	elseif i == 4 && j == 4 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 4, 3] = (32*a_2*(4*r_2*(a_2 + r_2)*(a_4*cos_θ_4 + 2*a_2*r_2*cos_θ_2 + 2*a_2*r*cos_θ_2 + r_4 + 2*r_3) + (a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4)*(a_4*(1 - cos_θ_2)^2 + 2*a_4*cos_θ_2 - a_4 + 2*a_2*r_2*cos_θ_2 + 2*a_2*r + r_4 + 2*r_3))*sine_θ_2 - 8*(a_2*(a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2)*sine_θ_2 + (a_2*cos_θ_2 + r_2)*(2*a_4*cos_θ_2 - a_4 + 2*a_2*r_2*cos_θ_2 - 4*a_2*r*cos_θ_2 + 4*a_2*r + r_4))*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 12*a_2*r_2*cos_θ_2 + r_6 + 12*r_4))/((a_2*cos_θ_2 + r_2)^2*(24*a_8*sine_θ_2*cos_θ_6 - 8*a_8*cos_θ^8 + 72*a_6*r_2*sine_θ_2*cos_θ_4 - 32*a_6*r_2*cos_θ_6 + 48*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(1 - cos_4θ) - 48*a_4*r_4*cos_θ_4 + 12*a_4*r_3*(1 - cos_4θ) + 20*a_4*r_2*(1 - cos_4θ) - 96*a_4*r_2*cos_θ_4 + 24*a_2*r_6*sine_θ_2 - 32*a_2*r_6*cos_θ_2 + 48*a_2*r^5*sine_θ_2 + 160*a_2*r_4*sine_θ_2 - 192*a_2*r_4*cos_θ_2 + 320*a_2*r_3*sine_θ_2 - 8*r^8 - 96*r_6)*tan(θ))
	elseif i == 4 && j == 4 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		cos_4θ = cos(4*θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_3 = r^3
		r_4 = r^4
		cos_θ_4 = cos_θ^4
		a_4 = a^4
		r_6 = r^6
		cos_θ_6 = cos_θ^6
		a_6 = a^6
		a_8 = a^8
    		Γ[4, 4, 4] = 2*a*(r*(a_2*(a_2*cos_θ_2 + r_2 + 2*r)*sine_θ_2 + (a_2*cos_θ_2 + r_2)^2) - (a_2*cos_θ_2 + r_2)*(a_2*(r + 1)*sine_θ_2 + 2*r*(a_2*cos_θ_2 + r_2)))*(a_6*cos_θ_6 + 3*a_4*r_2*cos_θ_4 + 3*a_2*r_4*cos_θ_2 + 4*a_2*r_2*cos_θ_2 + r_6 + 4*r_4)*sine_θ_2/((a_2*cos_θ_2 + r_2)^2*(-3*a_8*sine_θ_2*cos_θ_6 + a_8*cos_θ^8 - 9*a_6*r_2*sine_θ_2*cos_θ_4 + 4*a_6*r_2*cos_θ_6 - 6*a_6*r*sine_θ_2*cos_θ_4 + 9*a_4*r_4*(cos_4θ - 1)/8 + 6*a_4*r_4*cos_θ_4 + 3*a_4*r_3*(cos_4θ - 1)/2 + 5*a_4*r_2*(cos_4θ - 1)/2 + 12*a_4*r_2*cos_θ_4 - 3*a_2*r_6*sine_θ_2 + 4*a_2*r_6*cos_θ_2 - 6*a_2*r^5*sine_θ_2 - 20*a_2*r_4*sine_θ_2 + 24*a_2*r_4*cos_θ_2 - 40*a_2*r_3*sine_θ_2 + r^8 + 12*r_6))
	end
	return Γ[i, j, k]
end
	
function Kerr_BL_Christoffel_Symbols(i::Int64, j::Int64, k::Int64, r::Float64, θ::Float64, φ::Float64, a::Float64)
	Γ = zeros(4, 4, 4)
	if i == 1 && j == 1 && k ==2
		sine_θ = sin(θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		a_2 = a^2
		a_4 = a^4
    		Γ[1, 1, 2] = (a_4*sine_θ_2 - a_4 + a_2*r_2*sine_θ_2 + r^4)/((a_2*cos(θ)^2 + r_2)^2*(a_2 + r_2 - 2*r))
	elseif i == 1 && j == 1 && k == 3
		a_2 = a^2
    		Γ[1, 1, 3] = -4*a_2*r*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)^2
	elseif i == 1 && j == 2 && k == 1
		sine_θ = sin(θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		a_2 = a^2
		a_4 = a^4
    		Γ[1, 2, 1] = (a_4*sine_θ_2 - a_4 + a_2*r_2*sine_θ_2 + r^4)/((a_2*cos(θ)^2 + r_2)^2*(a_2 + r_2 - 2*r))
	elseif i == 1 && j == 2 && k == 4
		cos_θ = cos(θ)
		r_2 = r^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
    		Γ[1, 2, 4] = a*(a^4*cos_θ_2 - a_2*r_2*cos_θ_2 - a_2*r_2 - 3*r^4)*sin(θ)^2/((a_2*cos_θ_2 + r_2)^2*(a_2 + r_2 - 2*r))
	elseif i == 1 && j == 3 && k == 1
		a_2 = a^2
    		Γ[1, 3, 1] = -4*a_2*r*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)^2
	elseif i == 1 && j == 3 && k == 4
    		Γ[1, 3, 4] = 2*a^3*r*sin(θ)^3*cos(θ)/(a^2*cos(θ)^2 + r^2)^2
	elseif i == 1 && j == 4 && k == 2
		cos_θ = cos(θ)
		r_2 = r^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
    		Γ[1, 4, 2] = a*(a^4*cos_θ_2 - a_2*r_2*cos_θ_2 - a_2*r_2 - 3*r^4)*sin(θ)^2/((a_2*cos_θ_2 + r_2)^2*(a_2 + r_2 - 2*r))
	elseif i == 1 && j == 4 && k == 3
    		Γ[1, 4, 3] = 2*a^3*r*sin(θ)^3*cos(θ)/(a^2*cos(θ)^2 + r^2)^2
	elseif i == 2 && j == 1 && k == 1
		cos_θ = cos(θ)
		r_2 = r^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
    		Γ[2, 1, 1] = (a_2*cos_θ_2 - r_2)*(-a_2 - r_2 + 2*r)/(a_2*cos_θ_2 + r_2)^3
	elseif i == 2 && j == 1 && k == 4
		cos_θ = cos(θ)
		r_2 = r^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
    		Γ[2, 1, 4] = a*(a_2*cos_θ_2 - r_2)*(a_2 + r_2 - 2*r)*sin(θ)^2/(a_2*cos_θ_2 + r_2)^3
	elseif i == 2 && j == 2 && k == 2
		cos_θ = cos(θ)
		r_2 = r^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
    		Γ[2, 2, 2] = (r*(a_2 + r_2 - 2*r) + (1 - r)*(a_2*cos_θ_2 + r_2))/((a_2*cos_θ_2 + r_2)*(a_2 + r_2 - 2*r))
	elseif i == 2 && j == 2 && k == 3
		a_2 = a^2
    		Γ[2, 2, 3] = -a_2*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)
	elseif i == 2 && j == 3 && k == 2
		a_2 = a^2
    		Γ[2, 3, 2] = -a_2*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)
	elseif i == 2 && j == 3 && k == 3
		r_2 = r^2
		a_2 = a^2
    		Γ[2, 3, 3] = r*(-a_2 - r_2 + 2*r)/(a_2*cos(θ)^2 + r_2)
	elseif i == 2 && j == 4 && k == 1
		cos_θ = cos(θ)
		r_2 = r^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
    		Γ[2, 4, 1] = a*(a_2*cos_θ_2 - r_2)*(a_2 + r_2 - 2*r)*sin(θ)^2/(a_2*cos_θ_2 + r_2)^3
	elseif i == 2 && j == 4 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
    		Γ[2, 4, 4] = (r*(2*a_2*r*sine_θ_2 + (a_2 + r_2)*(a_2*cos_θ_2 + r_2)) + (a_2*cos_θ_2 + r_2)*(a_2*r*sine_θ_2 - 2*a_2*r - a_2*sine_θ_2 - 2*r^3))*(a_2 + r_2 - 2*r)*sine_θ_2/(a_2*cos_θ_2 + r_2)^3
	elseif i == 3 && j == 1 && k == 1
		a_2 = a^2
    		Γ[3, 1, 1] = -8*a_2*r*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)^3
	elseif i == 3 && j == 1 && k == 4
		r_2 = r^2
		a_2 = a^2
    		Γ[3, 1, 4] = 8*a*r*(a_2 + r_2)*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r_2)^3
	elseif i == 3 && j == 2 && k == 2
		r_2 = r^2
		a_2 = a^2
    		Γ[3, 2, 2] = a_2*sin(2*θ)/((a_2 + r_2 - 2*r)*(a_2*cos(2*θ) + a_2 + 2*r_2))
	elseif i == 3 && j == 2 && k == 3
    		Γ[3, 2, 3] = r/(a^2*cos(θ)^2 + r^2)
	elseif i == 3 && j == 3 && k == 2
    		Γ[3, 3, 2] = r/(a^2*cos(θ)^2 + r^2)
	elseif i == 3 && j == 3 && k == 3
		a_2 = a^2
    		Γ[3, 3, 3] = -a_2*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r^2)
	elseif i == 3 && j == 4 && k == 1
		r_2 = r^2
		a_2 = a^2
    		Γ[3, 4, 1] = 8*a*r*(a_2 + r_2)*sin(2*θ)/(a_2*cos(2*θ) + a_2 + 2*r_2)^3
	elseif i == 3 && j == 4 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		a_4 = a^4
    		Γ[3, 4, 4] = (-a_2*(2*a_2*r*sine_θ_2 + (a_2 + r_2)*(a_2*cos_θ_2 + r_2))*sine_θ_2 + (a_2*cos_θ_2 + r_2)*(-2*a_4*cos_θ_2 + a_4 - 2*a_2*r_2*cos_θ_2 + 4*a_2*r*cos_θ_2 - 4*a_2*r - r^4))*sine_θ*cos_θ/(a_2*cos_θ_2 + r_2)^3
	elseif i == 4 && j == 1 && k == 2
		cos_2θ = cos(2*θ)
		r_2 = r^2
		a_2 = a^2
    		Γ[4, 1, 2] = 2*a*(-a_2*cos_2θ - a_2 + 2*r_2)/((a_2 + r_2 - 2*r)*(a_2*cos_2θ + a_2 + 2*r_2)^2)
	elseif i == 4 && j == 1 && k == 3
    		Γ[4, 1, 3] = -2*a*r/((a^2*cos(θ)^2 + r^2)^2*tan(θ))
	elseif i == 4 && j == 2 && k == 1
		cos_2θ = cos(2*θ)
		r_2 = r^2
		a_2 = a^2
    		Γ[4, 2, 1] = 2*a*(-a_2*cos_2θ - a_2 + 2*r_2)/((a_2 + r_2 - 2*r)*(a_2*cos_2θ + a_2 + 2*r_2)^2)
	elseif i == 4 && j == 2 && k == 4
		cos_θ = cos(θ)
		r_2 = r^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		a_4 = a^4
    		Γ[4, 2, 4] = (a_4*r*(1 - cos_θ_2)^2 + 2*a_4*r*cos_θ_2 - a_4*r - a_4*(1 - cos_θ_2)^2 - a_4*cos_θ_2 + a_4 + 2*a_2*r^3*cos_θ_2 - a_2*r_2*cos_θ_2 - a_2*r_2 + r^5 - 2*r^4)/((a_2*cos_θ_2 + r_2)^2*(a_2 + r_2 - 2*r))
	elseif i == 4 && j == 3 && k == 1
    		Γ[4, 3, 1] = -2*a*r/((a^2*cos(θ)^2 + r^2)^2*tan(θ))
	elseif i == 4 && j == 3 && k == 4
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_4 = r^4
		a_4 = a^4
    		Γ[4, 3, 4] = (4*a_2*r_2*(a_2 + r_2)*sine_θ_2 + (a_2*(2*a_2*r*sine_θ_2 + (a_2 + r_2)*(a_2*cos_θ_2 + r_2))*sine_θ_2 + (a_2*cos_θ_2 + r_2)*(2*a_4*cos_θ_2 - a_4 + 2*a_2*r_2*cos_θ_2 - 4*a_2*r*cos_θ_2 + 4*a_2*r + r_4))*(a_2*cos_θ_2 + r_2 - 2*r))/((a_2*cos_θ_2 + r_2)^2*(a_4*cos_θ_2 + a_2*r_2*cos_θ_2 + a_2*r_2 - 2*a_2*r*cos_θ_2 + r_4 - 2*r^3)*tan(θ))
	elseif i == 4 && j == 4 && k == 2
		cos_θ = cos(θ)
		r_2 = r^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		a_4 = a^4
    		Γ[4, 4, 2] = (a_4*r*(1 - cos_θ_2)^2 + 2*a_4*r*cos_θ_2 - a_4*r - a_4*(1 - cos_θ_2)^2 - a_4*cos_θ_2 + a_4 + 2*a_2*r^3*cos_θ_2 - a_2*r_2*cos_θ_2 - a_2*r_2 + r^5 - 2*r^4)/((a_2*cos_θ_2 + r_2)^2*(a_2 + r_2 - 2*r))
	elseif i == 4 && j == 4 && k == 3
		sine_θ = sin(θ)
		cos_θ = cos(θ)
		r_2 = r^2
		sine_θ_2 = sine_θ^2
		cos_θ_2 = cos_θ^2
		a_2 = a^2
		r_4 = r^4
		a_4 = a^4
    		Γ[4, 4, 3] = (4*a_2*r_2*(a_2 + r_2)*sine_θ_2 + (a_2*(2*a_2*r*sine_θ_2 + (a_2 + r_2)*(a_2*cos_θ_2 + r_2))*sine_θ_2 + (a_2*cos_θ_2 + r_2)*(2*a_4*cos_θ_2 - a_4 + 2*a_2*r_2*cos_θ_2 - 4*a_2*r*cos_θ_2 + 4*a_2*r + r_4))*(a_2*cos_θ_2 + r_2 - 2*r))/((a_2*cos_θ_2 + r_2)^2*(a_4*cos_θ_2 + a_2*r_2*cos_θ_2 + a_2*r_2 - 2*a_2*r*cos_θ_2 + r_4 - 2*r^3)*tan(θ))
	end
	return Γ[i, j, k]
end
	
