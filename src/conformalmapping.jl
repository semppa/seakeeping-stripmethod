"""
This is the way to do multiline commenting

"""

# Tmp location

using Plots

function conformalMapping(section::SectionData, N::Int)
    # Generate initial solution as lewis coefficients

    # Initialize data
    I = length(section.ys)
    xᵢ = section.ys
    yᵢ = section.zs
    θᵢ = Vector{Float64}(undef, I)

    for i in 1:I
        θᵢ[i] = 0
    end

    x₀ᵢ = Vector{Float64}(undef, I)
    y₀ᵢ = Vector{Float64}(undef, I)

    #xᵢ[end] = secion.B / 2
    #yᵢ[1] = section.D
    θᵢ[1] = 0
    θᵢ[end] = π/2

    #cosφᵢ = Vector{Float64}(undef, I)
    #sinφᵢ = Vector{Float64}(undef, I)

    a₂ₙ₋₁ = Vector{Float64}(undef, N+1)

    for i in 1:N+1
        a₂ₙ₋₁[i] = 0
    end 

    a₂ₙ₋₁[1] = 1
    a₂ₙ₋₁[2:3] = lewisCoefficients(section.B, section.D, section.A)
    
    F = section.B / (2 * sum(a₂ₙ₋₁))

    #a₂ₙ₋₁[1:3] = [1, -0.1898, 0.2701]
    #F = 6.85

    #println(F, a₂ₙ₋₁)

    error = 1e10
    prev_error = 1.1e10 
    loopCount = 0
    # Loop until error is enough small:
    while abs(prev_error - error) > 1e-3 && loopCount < 100
        prev_error = error
        # Calculate new angles
        for i in 2:I-1
            cosφᵢ = ( xᵢ[i+1] - xᵢ[i-1]) / sqrt((xᵢ[i+1] - xᵢ[i-1])^2 + (yᵢ[i+1] - yᵢ[i-1])^2)
            sinφᵢ = (-yᵢ[i+1] + yᵢ[i-1]) / sqrt((xᵢ[i+1] - xᵢ[i-1])^2 + (yᵢ[i+1] - yᵢ[i-1])^2)

            f = θ -> xᵢ[i]*cosφᵢ + F*cosφᵢ*sum([(-1)^n * a₂ₙ₋₁[n+1] * sin((2n-1)θ) for n in 0:N]) -
                     yᵢ[i]*sinφᵢ + F*sinφᵢ*sum([(-1)^n * a₂ₙ₋₁[n+1] * cos((2n-1)θ) for n in 0:N])
            
            try
                θᵢ[i] = bisection(f, 0, π/2, 1e-3)
            catch err
                print("Cannot do it here. Sorry! \n")
                return
            end
        end
        #println(θᵢ)
        # Generate the matrices and solve them to find new coefficients
        A = zeros((N+1, N+1))
        B = zeros((N+1, 1))

        for i in 0:N
            for j in 0:N
                A[i+1, j+1] = (-1)^j * sum([cos((2i - 2j) * θᵢ[k]) for k in 1:I])
            end
            B[i+1] = sum([-xᵢ[k]sin((2i-1)θᵢ[k]) + yᵢ[k]cos((2i-1)θᵢ[k]) for k in 1:I])
        end

        #A[end-1, :] = [(-1)^n for n in 0:N]
        #A[end, :] = ones((1, N+1))

        #B[end - 1] = section.D
        #B[end] = section.B / 2
        
        X = zeros((N+1, 1))

        try 
            X = A \ B
        catch err
            print("Im dumb. Sorry! \n")
            return
        end

        F = X[1]
        a₂ₙ₋₁ = X / F

        # Calculate new points from new coefficients
        error = sum([
            (xᵢ[i] + sum([(-1)^n * F * a₂ₙ₋₁[n+1] * sin((2n-1)θᵢ[i]) for n in 0:N]))^2 +
            (yᵢ[i] - sum([(-1)^n * F * a₂ₙ₋₁[n+1] * cos((2n-1)θᵢ[i]) for n in 0:N]))^2
            for i in 1:I])

        # Calculate the error between calculated and actual points
        loopCount += 1
        #print("Loopcount: ", loopCount)
        #print(" | Error :", error)
        #print("\n")
    end

    # Return the coefficients
    print("N: ", N)
    #print("F: ", F, ' ')
    #print("a_2n1: ", a₂ₙ₋₁, ' ')
    print(" error: ", error, ' ')
    print(" Loopcount: ", loopCount, '\n')

    #thetas = LinRange(0, π/2, 100);

    x0s = [-F * sum([(-1)^n*a₂ₙ₋₁[n+1]*sin((2*n - 1)*θᵢ[i]) for n in 0:N]) for i in 1:I]
    y0s = [ F * sum([(-1)^n*a₂ₙ₋₁[n+1]*cos((2*n - 1)*θᵢ[i]) for n in 0:N]) for i in 1:I]

    #x0s = [-F * sum([(-1)^n*a₂ₙ₋₁[n+1]*sin((2*n - 1)*thetas[i]) for n in 0:N]) for i in 1:100]
    #y0s = [ F * sum([(-1)^n*a₂ₙ₋₁[n+1]*cos((2*n - 1)*thetas[i]) for n in 0:N]) for i in 1:100]

    #print("original points: x:", xᵢ, "y:", yᵢ, '\n')
    #print("mapped points:", "x:", x0s, "y:", y0s, '\n')

    p = scatter(xᵢ, -yᵢ, label="Original points")
    plot!(x0s, -y0s, label="Mapped points")
    # display(p)
    filename = string(N, "_image.png")
    savefig(p, filename)
end

function lewisCoefficients(B::Float64, D::Float64, A::Float64)
    H₀ = (B/2) / D
    σₛ = A / (D*B)

    c₁ = 3 + 4*σₛ / π + (1 - 4*σₛ/π)*((H₀ - 1)/(H₀ + 1))^2

    a₃ = (-c₁ + 3 + sqrt(9 - 2*c₁)) / c₁
    a₁ = (H₀ - 1) / (H₀ + 1) * (a₃ + 1)
    
    return [a₁, a₃]
end

function bisection(f, a, b, tol)
    if sign(f(a)) == sign(f(b))
        # print(f(a), " ", f(b), '\n')
        error("function has the same sign at given endpoints")
    end
    mid = (a + b)/2

    while abs(f(mid)) > tol
        sign(f(mid)) == sign(f(a)) ? a=mid : b=mid
        mid = (a + b)/2
    end
    mid
end