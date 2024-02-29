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

include("conformalmapping.jl")
