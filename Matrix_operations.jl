using LinearAlgebra
using Base.Threads
using Distributed

#Equations for inverting a 4x4 matrix were prepared using Wolfram Mathematica.

#A = {{a11, a12, a13, a14}, {a21, a22, a23, a24}, {a31, a32, a33, a34}, {a41, a42, a43, a44}}
#Inverse[A]


function inverse_4x4(matrix::Matrix{Float64})
    a11, a12, a13, a14 = matrix[1, 1], matrix[1, 2], matrix[1, 3], matrix[1, 4]
    a21, a22, a23, a24 = matrix[2, 1], matrix[2, 2], matrix[2, 3], matrix[2, 4]
    a31, a32, a33, a34 = matrix[3, 1], matrix[3, 2], matrix[3, 3], matrix[3, 4]
    a41, a42, a43, a44 = matrix[4, 1], matrix[4, 2], matrix[4, 3], matrix[4, 4]

    det = a14*a23*a32*a41 - a13*a24*a32*a41 - a14*a22*a33*a41 + a12*a24*a33*a41 + a13*a22*a34*a41 - a12*a23*a34*a41 -
          a14*a23*a31*a42 + a13*a24*a31*a42 + a14*a21*a33*a42 - a11*a24*a33*a42 - a13*a21*a34*a42 + a11*a23*a34*a42 +
          a14*a22*a31*a43 - a12*a24*a31*a43 - a14*a21*a32*a43 + a11*a24*a32*a43 + a12*a21*a34*a43 - a11*a22*a34*a43 -
          a13*a22*a31*a44 + a12*a23*a31*a44 + a13*a21*a32*a44 - a11*a23*a32*a44 - a12*a21*a33*a44 + a11*a22*a33*a44

    if det == 0
        error("The Matrix 4x4 is singular!")
    end

    @inbounds inv_matrix = [
        ( -a24*a33*a42 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 + a22*a33*a44 )/ det,
        ( a14*a33*a42 - a13*a34*a42 - a14*a32*a43 + a12*a34*a43 + a13*a32*a44 - a12*a33*a44 ) / det,
        ( -a14*a23*a42 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 + a12*a23*a44 )/ det,
        ( a14*a23*a32 - a13*a24*a32 - a14*a22*a33 + a12*a24*a33 + a13*a22*a34 - a12*a23*a34 ) / det,

        ( a24*a33*a41 - a23*a34*a41 - a24*a31*a43 + a21*a34*a43 + a23*a31*a44 - a21*a33*a44 ) / det,
        ( -a14*a33*a41 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 + a11*a33*a44 )/ det,
        ( a14*a23*a41 - a13*a24*a41 - a14*a21*a43 + a11*a24*a43 + a13*a21*a44 - a11*a23*a44 ) / det,
        ( -a14*a23*a31 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 + a11*a23*a34 )/ det,

        ( -a24*a32*a41 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 + a21*a32*a44 )/ det,
        ( a14*a32*a41 - a12*a34*a41 - a14*a31*a42 + a11*a34*a42 + a12*a31*a44 - a11*a32*a44 ) / det,
        ( -a14*a22*a41 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 + a11*a22*a44 )/ det,
        ( a14*a22*a31 - a12*a24*a31 - a14*a21*a32 + a11*a24*a32 + a12*a21*a34 - a11*a22*a34 ) / det,

        ( a23*a32*a41 - a22*a33*a41 - a23*a31*a42 + a21*a33*a42 + a22*a31*a43 - a21*a32*a43 ) / det,
        ( -a13*a32*a41 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 + a11*a32*a43 )/ det,
        ( a13*a22*a41 - a12*a23*a41 - a13*a21*a42 + a11*a23*a42 + a12*a21*a43 - a11*a22*a43 ) / det,
        ( -a13*a22*a31 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33 )/ det
    ]
    
    return reshape(inv_matrix, 4, 4)
end

# For calculating sqrt(-g), where g is the determinant of the metric

function sqrt_g(matrix::Matrix{Float64})
    a11, a12, a13, a14 = matrix[1, 1], matrix[1, 2], matrix[1, 3], matrix[1, 4]
    a21, a22, a23, a24 = matrix[2, 1], matrix[2, 2], matrix[2, 3], matrix[2, 4]
    a31, a32, a33, a34 = matrix[3, 1], matrix[3, 2], matrix[3, 3], matrix[3, 4]
    a41, a42, a43, a44 = matrix[4, 1], matrix[4, 2], matrix[4, 3], matrix[4, 4]
    
    det = a14*a23*a32*a41 - a13*a24*a32*a41 - a14*a22*a33*a41 + a12*a24*a33*a41 + a13*a22*a34*a41 - a12*a23*a34*a41 -
          a14*a23*a31*a42 + a13*a24*a31*a42 + a14*a21*a33*a42 - a11*a24*a33*a42 - a13*a21*a34*a42 + a11*a23*a34*a42 +
          a14*a22*a31*a43 - a12*a24*a31*a43 - a14*a21*a32*a43 + a11*a24*a32*a43 + a12*a21*a34*a43 - a11*a22*a34*a43 -
          a13*a22*a31*a44 + a12*a23*a31*a44 + a13*a21*a32*a44 - a11*a23*a32*a44 - a12*a21*a33*a44 + a11*a22*a33*a44

    if det >= 0
        error("The determinant is non-negative, cannot compute sqrt of a non-negative value.")
    end

    return sqrt(-det)
end
