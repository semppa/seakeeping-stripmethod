"""
This is the way to do multiline commenting



"""

# Tmp location
struct SectionData
    framelocation::String
    ys::Vector{Float64}
    zs::Vector{Float64}
    B::Float64
    D::Float64
    A::Float64
end

function conformalMapping(section::SectionData, N::Int)
    # Generate initial solution as lewis coefficients

    # Loop until error is enough small:
        # Calculate new point angles from This

        # Generate the matrices and solve them to find new coefficients

        # Calculate new points from new coefficients

        # Calculate the error between calculated and actual points

    # Return the coefficients
end

function lewisCoefficients(B::Float64, D::Float64, A::Float64)
    H₀ = (B/2) / D
    σₛ = A / (D*B)

    c₁ = 3 + 4*σₛ / π + (1 - 4*σₛ/π)((H₀ - 1)/(H₀ + 1))^2

    a₃ = (-c₁ + 3 + sqrt(9 - c₁)) / c₁
    a₁ = (H₀ - 1) / (H₀ + 1) * (a₃ + 1)
    
    return [a₃, a₁]
end

