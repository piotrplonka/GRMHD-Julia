using BenchmarkTools
using Profile
using Base.Threads

function u_pow0(x1::Float64,x2::Float64,x3::Float64,gcov::Matrix{Float64})
    a=gcov[1,1]
    b=2*(gcov[1,2]*x1+gcov[1,3]*x2+gcov[1,4]*x3)
    c=(gcov[2,2]*x1*x1+2*(gcov[2,3]*x1*x2 + gcov[2,4]*x1*x3 + gcov[4,4]*x3*x3 + gcov[3,4]*x2*x3)+1)
    return (-(2*(gcov[1,2]*x1+gcov[1,3]*x2+gcov[1,4]*x3))+sqrt((2*(gcov[1,2]*x1+gcov[1,3]*x2+gcov[1,4]*x3))^2-4*gcov[1,1]*(gcov[2,2]*x1*x1+2*(gcov[2,3]*x1*x2 + gcov[2,4]*x1*x3 + gcov[4,4]*x3*x3 + gcov[3,4]*x2*x3)+1)))/(2*gcov[1,1])
end



