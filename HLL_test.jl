using BenchmarkTools
using Base.Threads
using Profile
using LinearAlgebra  

#include("Matrix_operations.jl")
#include("Christoffel_Symbols.jl")
include("structs.jl")
#include("Initial_Conditions.jl")
#include("fluxlimiter.jl")
#include("eos.jl")
#include("PPM.jl")
#include("velocity_composition.jl")
#include("GR_Flux.jl")

a=0.8

#for i in 3:(N1-2)
   # for j in 3:(N2-2)
        #for k in 3:(N3-2)
		#metric_kerrros =  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a)
		#println(metric_kerrros)
		#zeros_U =zeros(8)
		#PtoU(grid[i,j,k,:],  zeros_U,  metric_kerrros, eos::Polytrope)
		#println(zeros_U)
		#println("i:",j,"j:",j,"k:",k)
		#println(P_values)
		#println(zeros_U)
		#println("Tu zaczyna się iteracja")
		#println("Iteracja i=$i, j=$j, k=$k")
		#println(grid[i,j,k,:]*1.1)
		#println(zeros_U)		
		#UtoP(zeros_U, grid[i,j,k,:]*1.1,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a), eos::Polytrope)	
  #      end
 #   end
#end


#Iteracja i=62, j=30, k=11
#[2.2685287494586586, 2.0279092042248346, -0.06302874945865838, -0.06302874945865838, -0.06302874945865838, 0.01785287494586584, 0.01785287494586584, 0.01785287494586584]
#[114995.08678214469, -2.0541085949912039e6, 88160.89810849521, -1.5482378959918973e8, -1.0128640003282715e8, 123.93110213361336, 123.93110213361336, 123.93110213361336]
#Nie zbiegło!!!

#i=62
#j=30
#k=11
#metric_kerrros=Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a)
#zeros_U=zeros(8)
##println("TEST!")
#PtoU(grid[i,j,k,:],  zeros_U,  metric_kerrros, eos::Polytrope)
#println(grid[i,j,k,:])
#println(zeros_U)

#[2.062298863144235, 1.843553822022577, -0.057298863144234886, -0.057298863144234886, -0.057298863144234886, 0.01622988631442349, 0.01622988631442349, 0.01622988631442349]
#[114995.08678214469, -2.0541085949912039e6, 88160.89810849521, -1.5482378959918973e8, -1.0128640003282715e8, 123.93110213361336, 123.93110213361336, 123.93110213361336]

#P1 = [2.062298863144235, 1.843553822022577, -0.057298863144234886, -0.057298863144234886, -0.057298863144234886, 0.01622988631442349, 0.01622988631442349, 0.01622988631442349]
#U1= [114995.08678214469, -2.0541085949912039e6, 88160.89810849521, -1.5482378959918973e8, -1.0128640003282715e8, 123.93110213361336, 123.93110213361336, 123.93110213361336]

#UtoP_LU(U1, P1*1.1,  Kerr_Schild_metric_cov(N1_grid[i], N2_grid[j], N3_grid[k], a), eos::Polytrope)


