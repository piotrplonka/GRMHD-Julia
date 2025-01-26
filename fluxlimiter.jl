function minmod(r::AbstractVector)
    return max.(0, min.(1, r))  
end
