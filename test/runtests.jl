using Brooglie
using Test
using LinearAlgebra

@testset "Function introspection" begin
    for N in 1:15
        args = Expr(:tuple, [Symbol("x_$i") for i in 1:N]...)
        f = eval(:($args -> 42))
        @test Brooglie.numberofarguments(f) == N
        # Test also f manually, calling it.
        @test f(ones(N)...) == 42 # Correct num of args
        @test_throws MethodError f(ones(N-1)...)
        @test_throws MethodError f(ones(N+1)...)
    end
end

@testset "Combinations" begin
    @test Brooglie.combinations(1:3, 0) == []
    @test Brooglie.combinations(1:3, 1) == [[1], [2], [3]]
    @test Brooglie.combinations(1:1, 4) == [[1, 1, 1, 1]]
    @test length(Brooglie.combinations([], 9)) == 0
    for n in 1:5
        for k in 1:5
            @test length(Brooglie.combinations(1:n, k)) == n^k
        end
    end
end

@testset "Hamiltonian generation" begin
    # Examples made by hand. Note the recursive nature.
    𝕀(n) = diagm(0 => ones(n))

    expected_1D = [+2 -1 +0
                   -1 +2 -1
                   +0 -1 +2]
    expected_2D = [expected_1D -𝕀(3)       0𝕀(3)
                   -𝕀(3)       expected_1D -𝕀(3)
                   0*𝕀(3)      -𝕀(3)       expected_1D] + 2𝕀(9)
    expected_3D = [expected_2D -𝕀(9) 0𝕀(9)
                   -𝕀(9)  expected_2D -𝕀(9)
                   0*𝕀(9)  -𝕀(9)  expected_2D] + 2𝕀(27)
    V1(x) = 0
    V2(x,y) = 0
    V3(x,y,z) = 0
    @test buildH(V1, N=3, a=-1, b=1, m=1) == expected_1D
    @test buildH(V2, N=3, a=-1, b=1, m=1) == expected_2D
    @test buildH(V3, N=3, a=-1, b=1, m=1) == expected_3D
    W1(x) = 1
    W2(x,y) = 2
    W3(x,y,z) = 3
    ε = 2*1*(2/3)^2 # Carefull with the adimensionalization of V
    @test buildH(W1, N=3, a=-1, b=1, m=1) == expected_1D + ε*𝕀(3)
    @test buildH(W2, N=3, a=-1, b=1, m=1) == expected_2D + ε*2𝕀(9)
    @test buildH(W3, N=3, a=-1, b=1, m=1) == expected_3D + ε*3𝕀(27)
end


@testset "Normalization" begin
    # If we integrate sin(x) from 0 to π the area should be 2. Also,
    # ensure that normalizewf normalizes |φ|² and not φ: it should
    # divide by the square root of the area of |φ|², not by the area.
    # In this case, ∫|sin(x)|²dx ≃ π/2.
    N = 100 # Will grow quickly in 3D
    rr = LinRange(0,π,N)
    f(r) = prod(sin.(r)) # f(x,y,...) = sin(x)*sin(y)*...
    # 1D:
    mesh = [[x] for x in rr]
    φ = f.(mesh)
    A = Brooglie.integrate(φ, π)
    @test isapprox(A, 2, atol=0.1)
    nφ = Brooglie.normalizewf(φ,π)
    @test all(filter(!isnan, φ ./ nφ) .≈ √Brooglie.integrate(φ.^2,π))
    # 2D:
    mesh = [[x,y] for x in rr, y in rr]
    φ = f.(mesh)
    A = Brooglie.integrate(φ, π)
    @test isapprox(A, 2^2, atol=0.1)
    nφ = Brooglie.normalizewf(φ,π)
    @test all(filter(!isnan, φ ./ nφ) .≈ √Brooglie.integrate(φ.^2,π))
    # 3D (precission starts to fall short):
    mesh = [[x,y,z] for x in rr, y in rr, z in rr]
    φ = f.(mesh)
    A = Brooglie.integrate(φ, π)
    @test isapprox(A, 2^3, atol=0.3)
    nφ = Brooglie.normalizewf(φ,π)
    @test all(filter(!isnan, φ ./ nφ) .≈ √Brooglie.integrate(φ.^2,π))
end

# The following tests check mainly that the energy levels and
# wavefunctions are the expected ones. Before them, some quantum
# boilerplate:

"""
     hermite(n,x) ≡ Hₙ(x)

Physicists' Hermite polinomial
"""
function hermite(n,x)
    if iseven(n)
        ζ = n÷2
        return factorial(n) * sum( (-1)^(ζ-ℓ) /
                                   factorial(2ℓ) /
                                   factorial(ζ-ℓ) *
                                   (2x)^(2ℓ) for ℓ in 0:ζ)
    else
        ζ = (n-1)÷2
        return factorial(n) * sum( (-1)^(ζ-ℓ) /
                                   factorial(2ℓ+1) /
                                   factorial(ζ-ℓ) *
                                   (2x)^(2ℓ+1) for ℓ in 0:ζ)
    end
end

"""
Return the `n`∈{0,⋯} wavefunction φₙ(`x`) of an harmonic oscillator of
frequency `ω` and mass `m`. Uses ħ=1.
"""
function QAO(ω,m,n,x)
    return 1/sqrt(2^n * factorial(n)) *
        (m*ω/π)^(1/4) *
        exp(-m*ω*x^2 /2) *
        hermite(n,sqrt(m*ω)*x)
end

"""
Return the `n`∈{1,⋯} wavefunction φₙ(`x`) of a particle of mass `m` in a
box of length `L`."
"""
function box(L,m,n,x)
    A = sqrt(2/L)
    k = n*π/L
    return A*sin(k*(x+L/2))
end

absmaximum(x) = maximum(abs.(x))

@testset "Harmonic oscillator (1D)" begin
    N = 1000
    nev = 10
    m = 1
    ω = 1
    a, b = -10, 10
    V(x) = 1/2 * m * ω^2 * x^2
    EineV, v = solve(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(v) == nev
    for i in 1:nev
        # Test if E ∼ ħω(n+½). Remember, ħ=1.
        n = i-1
        expectedE = ω*(n + 1/2)
        @test isapprox(expectedE, EineV[i],
                       rtol=1e-3)
        # Test if the wavefunction is the expected one. Check only
        # |φ|², neglecting arbitrary phase factors.
        expectedwf = QAO.(ω,m,n,LinRange(a,b,N))
        @test isapprox(abs2.(expectedwf), abs2.(v[i]),
                       norm=absmaximum,
                       rtol=0.01)
    end
end

@testset "Particle in a box (1D)" begin
    N = 1000
    nev = 10
    m = 1
    a, b = -1, 1
    V(x) = 0 # The boundary conditions serve as a box.
    EineV, v = solve(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(v) == nev
    L = b-a
    for n in 1:nev
        # Test if E ≃ n²π²ħ² / 2mL². Remember the atomic units!
        expectedE = n^2*π^2 / (2m*L^2)
        @test isapprox(expectedE, EineV[n]; rtol=1e-2)
        # Test if the wavefunction is the expected one.
        expectedwf = box.(L,m,n,range(a,stop=b,length=N))
        @test isapprox(abs2.(expectedwf), abs2.(v[n]),
                       norm=absmaximum,
                       rtol=0.05)
    end
end

@testset "Harmonic oscillator (2D)" begin
    N = 100
    nev = 10
    m = 1
    ω = 1
    a, b = -10, 10
    V(x,y) = 1/2 * m * ω^2 * (x^2+y^2)
    EineV, vv = solve(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(vv) == nev
    # E ∼ ħω(∑nᵢ+N/2). Remember, ħ=1.
    expectedElist = [ω*(n1+n2+2/2) for n1 in 0:nev-1, n2 in 0:nev-1]
    expectedElist = sort(vec(expectedElist))[1:nev]
    for i in 1:nev
        # Test the energy.
        expectedE = expectedElist[i]
        @test isapprox(expectedE, EineV[i]; rtol=1e-2)
    end
    # Test if the ground wavefunction is the expected one. Only
    # checking ground state because of degeneracy: is very dificult to
    # distinguish between |12⟩ and |21⟩, for example.
    ll = range(a, stop=b, length=N)
    expectedwf = [QAO.(ω,m,0,x)*QAO.(ω,m,0,y) for x in ll, y in ll]
    @test isapprox(abs2.(expectedwf), abs2.(vv[1]),
                   norm=absmaximum,
                   rtol=0.05)
end

@testset "Particle in a box (2D)" begin
    N = 100
    nev = 10
    m = 1
    a, b = -1, 1
    V(x,y) = 0
    EineV, vv = solve(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(vv) == nev
    # E ∼ n²π²ħ² / 2mL², with ħ=1
    L = b-a
    expectedElist = [(n1^2+n2^2)*π^2/(2*m*L^2) for n1 in 1:nev, n2 in 1:nev]
    expectedElist = sort(vec(expectedElist))[1:nev]
    for i in 1:nev
        # Test the energy.
        @test isapprox(expectedElist[i], EineV[i]; rtol=5e-2)
    end
    # Test, again, if the ground wavefunction is the expected one.
    ll = range(a, stop=b, length=N)
    expectedwf = [box.(L,m,1,x)*box.(L,m,1,y) for x in ll, y in ll]
    @test isapprox(abs2.(expectedwf), abs2.(vv[1]),
                   norm=absmaximum,
                   rtol=0.05)
end


@testset "Harmonic oscillator (3D)" begin
    N = 30 # Sloooooowwww
    nev = 10
    m = 1
    ω = 1
    a, b = -10, 10
    V(x,y,z) = 1/2 * m * ω^2 * (x^2+y^2+z^2)
    EineV, vv = solve(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(vv) == nev
    # E ∼ ħω(∑nᵢ+N/2). Remember, ħ=1.
    expectedElist = [ω*(n1+n2+n3+3/2)
                     for n1 in 0:nev-1,
                     n2 in 0:nev-1,
                     n3 in 0:nev-1]
    expectedElist = sort(vec(expectedElist))[1:nev]
    for i in 1:nev
        # Test the energy.
        expectedE = expectedElist[i]
        @test isapprox(expectedE, EineV[i]; rtol=1e-1)
    end
    # Test if the ground wavefunction is the expected one.
    ll = range(a, stop=b, length=N)
    expectedwf = [QAO.(ω,m,0,x)*QAO.(ω,m,0,y)*QAO.(ω,m,0,z)
                  for x in ll, y in ll, z in ll]
    @test isapprox(abs2.(expectedwf), abs2.(vv[1]),
                   norm=absmaximum,
                   atol=0.03)
end

@testset "Particle in a box (3D)" begin
    N = 30
    nev = 10
    m = 1
    a, b = -1, 1
    V(x,y,z) = 0
    EineV, vv = solve(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(vv) == nev
    # E ∼ n²π²ħ² / 2mL², with ħ=1
    L = b-a
    expectedElist = [(n1^2+n2^2+n3^2)*π^2/(2*m*L^2)
                     for n1 in 1:nev, n2 in 1:nev, n3 in 1:nev]
    expectedElist = sort(vec(expectedElist))[1:nev]
    for i in 1:nev
        # Test the energy.
        expectedE = expectedElist[i]
        @test isapprox(expectedE, EineV[i]; rtol=1e-1)
    end
    # Test if the ground wavefunction is the expected one.
    ll = range(a, stop=b, length=N)
    expectedwf = [box.(L,m,1,x)*box.(L,m,1,y)*box.(L,m,1,z)
                  for x in ll, y in ll, z in ll]
    @test isapprox(abs2.(expectedwf), abs2.(vv[1]),
                   norm=absmaximum,
                   atol=0.1)
end
