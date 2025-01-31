function MINMOD(q_im1,q_i,q_ip1)
    sp = q_ip1 - q_i
    sm = q_i - q_im1
    ssp = sign(sp)
    ssm = sign(sm)
    asp = abs(sp)
    asm = abs(sm)
    dU = 0.25 * (ssp + ssm) * min(asp,asm)
    return q_i - dU, q_i + dU    
end

function SQR(x)
    return x^2
end

function WENOZ(q_im2::T,q_im1::T,q_i::T,q_ip1::T,q_ip2::T) where T
    beta_coeff = (T(13.) / T(12.), T(0.25))
    beta = (
        beta_coeff[1] * SQR(q_im2 + q_i - T(2) * q_im1) +
        beta_coeff[2] * SQR(q_im2 + T(3) * q_i - T(4) * q_im1),
        
        beta_coeff[1] * SQR(q_im1 + q_ip1 - T(2) * q_i) +
        beta_coeff[2] * SQR(q_im1 - q_ip1),
        
        beta_coeff[1] * SQR(q_ip2 + q_i - T(2) * q_ip1) +
        beta_coeff[2] * SQR(q_ip2 + T(3) * q_i - T(4) * q_ip1)
    )

    epsL = T(1.0e-42)
    tau_5 = abs(beta[1] - beta[3])
    
    indicator = (
        tau_5 / (beta[1] + epsL),
        tau_5 / (beta[2] + epsL),
        tau_5 / (beta[3] + epsL)
    )

    f = (
        T(2.) * q_im2 - T(7.) * q_im1 + T(11.) * q_i,
        -q_im1 + T(5.) * q_i + T(2.) * q_ip1,
        T(2.) * q_i + T(5.) * q_ip1 - q_ip2
    )

    alpha = (
        T(0.1) * (T(1.) + SQR(indicator[1])),
        T(0.6) * (T(1.) + SQR(indicator[2])),
        T(0.3) * (T(1.) + SQR(indicator[3]))
    )
    
    alpha_sum = T(6.) * sum(alpha)
    ql_ip1 = sum(f .* alpha) / alpha_sum

    f = (
        T(2.) * q_ip2 - T(7.) * q_ip1 + T(11.) * q_i,
        -q_ip1 + T(5.) * q_i + T(2.) * q_im1,
        T(2.) * q_i + T(5.) * q_im1 - q_im2
    )

    alpha = (
        T(0.1) * (T(1.) + SQR(indicator[3])),
        T(0.6) * (T(1.) + SQR(indicator[2])),
        T(0.3) * (T(1.) + SQR(indicator[1]))
    )
    
    alpha_sum = T(6.) * sum(alpha)
    qr_i = sum(f .* alpha) / alpha_sum

    return qr_i,ql_ip1
end


function WENOZ_vec(q_im2::AbstractVector{T}, q_im1::AbstractVector{T}, q_i::AbstractVector{T}, q_ip1::AbstractVector{T}, q_ip2::AbstractVector{T}) where T

    beta_coeff = (T(13) / T(12), T(0.25))
    SQR(x) = x .^ 2
    
    beta = (
        beta_coeff[1] .* SQR.(q_im2 .+ q_i .- T(2) .* q_im1) +
        beta_coeff[2] .* SQR.(q_im2 .+ T(3) .* q_i .- T(4) .* q_im1),
        
        beta_coeff[1] .* SQR.(q_im1 .+ q_ip1 .- T(2) .* q_i) +
        beta_coeff[2] .* SQR.(q_im1 .- q_ip1),
        
        beta_coeff[1] .* SQR.(q_ip2 .+ q_i .- T(2) .* q_ip1) +
        beta_coeff[2] .* SQR.(q_ip2 .+ T(3) .* q_i .- T(4) .* q_ip1)
    )
    
    epsL = T(1.0e-42)
    tau_5 = abs.(beta[1] .- beta[3])
    
    indicator = (
        tau_5 ./ (beta[1] .+ epsL),
        tau_5 ./ (beta[2] .+ epsL),
        tau_5 ./ (beta[3] .+ epsL)
    )
    
    f = (
        T(2) .* q_im2 .- T(7) .* q_im1 .+ T(11) .* q_i,
        -q_im1 .+ T(5) .* q_i .+ T(2) .* q_ip1,
        T(2) .* q_i .+ T(5) .* q_ip1 .- q_ip2
    )
    
    alpha = (
        T(0.1) .* (T(1) .+ SQR.(indicator[1])),
        T(0.6) .* (T(1) .+ SQR.(indicator[2])),
        T(0.3) .* (T(1) .+ SQR.(indicator[3]))
    )
    
    alpha_sum = T(6) .* (alpha[1] .+ alpha[2] .+ alpha[3])
    ql_ip1 = (f[1] .* alpha[1] .+ f[2] .* alpha[2] .+ f[3] .* alpha[3]) ./ alpha_sum
    
    f = (
        T(2) .* q_ip2 .- T(7) .* q_ip1 .+ T(11) .* q_i,
        -q_ip1 .+ T(5) .* q_i .+ T(2) .* q_im1,
        T(2) .* q_i .+ T(5) .* q_im1 .- q_im2
    )
    
    alpha = (
        T(0.1) .* (T(1) .+ SQR.(indicator[3])),
        T(0.6) .* (T(1) .+ SQR.(indicator[2])),
        T(0.3) .* (T(1) .+ SQR.(indicator[1]))
    )
    
    alpha_sum = T(6) .* (alpha[1] .+ alpha[2] .+ alpha[3])
    qr_i = (f[1] .* alpha[1] .+ f[2] .* alpha[2] .+ f[3] .* alpha[3]) ./ alpha_sum
    
    return qr_i, ql_ip1
end



