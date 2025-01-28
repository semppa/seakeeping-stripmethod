"""
The so called "main" function to drive the program

Code is written in pseudo-format for now
"""

# Read the input file for
    # Geometry:
        # If body plan is given, then use it
        # If STL or IGES is supplied, then generate a body plan of cross-sections
    # Input parameters:
        # Find out what is required

# Generate Sections and calculate their respective values

# For each section (Utilize GPU if possible)
    # Calculate the conformal mapping coefficients
    # Calculate the added masses and damping
    # Solve the RAOs
    # Others? (Resistance/added wave resistance/etc)

# Postprocess the data

struct SectionData
    framelocation::String
    ys::Vector{Float64}
    zs::Vector{Float64}
    B::Float64
    D::Float64
    A::Float64
end

include("conformalmapping.jl")

xs = [0, 0.8, 1.2, 1.2, 1.4, 2.2, 3.5, 5.0, 6.5, 7.4] #[0, 1.5, 2.3, 3.0, 3.3, 4.3, 7.4] 
ys = [10, 9.7, 8.4, 6.1, 4.3, 2.1, 1.2, 0.55, 0.15, 0] # [10, 9.7, 8.4, 6.1, 4.3, 2.1, 0] 

tmp = SectionData("lol", xs, ys, 2*7.4, 10, 0.235*7.4*10*2)

#conformalMapping(tmp, 5)

for b in 4:100
    conformalMapping(tmp, b)
end
