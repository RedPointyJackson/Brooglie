using Brooglie
using Base.Test

@testset "Normalization" begin
    # If we integrate sin(x) from 0 to π the area should be 2. Also,
    # ensure that normalizewf normalizes |φ|² and not φ: it should
    # divide by the square root of the area of |φ|², not by the area.
    # In this case, ∫|sin(x)|²dx ≃ π/2.
    N = 100 # Will grow quickly in 3D
    rr = linspace(0,π,N)
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

"Physicists' Hermite polinomial `hermite(n,x)` ≡ Hₙ(x)"
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

"Return the `n`∈{0,⋯} wavefunction φₙ(`x`) of an harmonic oscillator of
frequency `ω` and mass `m`. Uses ħ=1."
function QAO(ω,m,n,x)
    return 1/sqrt(2^n * factorial(n)) *
        (m*ω/π)^(1/4) *
        exp(-m*ω*x^2 /2) *
        hermite(n,sqrt(m*ω)*x)
end

"Return the `n`∈{1,⋯} wavefunction φₙ(`x`) of a particle of mass `m` in a
box of length `L`."
function box(L,m,n,x)
    A = sqrt(2/L)
    k = n*π/L
    return A*sin(k*(x+L/2))
end

@testset "Harmonic oscillator (1D)" begin
    N = 1000
    nev = 10
    m = 1
    ω = 1
    a, b = -10, 10
    V(x) = 1/2 * m * ω^2 * x^2
    EineV, v = solve1D(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(v) == nev
    for i in 1:nev
        # Test if E ∼ ħω(n+½). Remember, ħ=1.
        n = i-1
        expectedE = ω*(n + 1/2)
        @test isapprox(expectedE, EineV[i]; rtol=1e-3)
        # Test if the wavefunction is the expected one. Check only
        # |φ|², neglecting arbitrary phase factors.
        expectedwf = QAO.(ω,m,n,linspace(a,b,N))
        @test isapprox(abs2.(expectedwf), abs2.(v[i]), rtol=0.01)
    end
end

@testset "Particle in a box (1D)" begin
    N = 1000
    nev = 10
    m = 1
    a, b = -1, 1
    V(x) = 0 # The boundary conditions serve as a box.
    EineV, v = solve1D(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(v) == nev
    L = b-a
    for n in 1:nev
        # Test if E ≃ n²π²ħ² / 2mL². Remember the atomic units!
        expectedE = n^2*π^2 / (2m*L^2)
        @test isapprox(expectedE, EineV[n]; rtol=1e-2)
        # Test if the wavefunction is the expected one.
        expectedwf = box.(L,m,n,linspace(a,b,N))
        @test isapprox(abs2.(expectedwf), abs2.(v[n]), rtol=0.1)
    end
end

@testset "Harmonic oscillator (2D)" begin
    N = 250
    nev = 10
    m = 1
    ω = 1
    a, b = -10, 10
    V(x,y) = 1/2 * m * ω^2 * (x^2+y^2)
    EineV, v = solve2D(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(v) == nev
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
    ll = linspace(a,b,N)
    expectedwf = [QAO.(ω,m,0,x)*QAO.(ω,m,0,y) for x in ll, y in ll]
    @test isapprox(abs2.(expectedwf), abs2.(v[1]), rtol=0.01)
end

@testset "Particle in a box (2D)" begin
    N = 250
    nev = 10
    m = 1
    a, b = -1, 1
    V(x,y) = 0
    EineV, v = solve2D(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(v) == nev
    # E ∼ n²π²ħ² / 2mL², with ħ=1
    L = b-a
    expectedElist = [(n1^2+n2^2)*π^2/(2*m*L^2) for n1 in 1:nev, n2 in 1:nev]
    expectedElist = sort(vec(expectedElist))[1:nev]
    for i in 1:nev
        # Test the energy.
        expectedE = expectedElist[i]
        @test isapprox(expectedE, EineV[i]; rtol=1e-2)
    end
    # Test, again, if the ground wavefunction is the expected one.
    ll = linspace(a,b,N)
    expectedwf = [box.(L,m,1,x)*box.(L,m,1,y) for x in ll, y in ll]
    @test isapprox(abs2.(expectedwf), abs2.(v[1]), rtol=0.01)
end


@testset "Harmonic oscillator (3D)" begin
    N = 30 # Sloooooowwww
    nev = 10
    m = 1
    ω = 1
    a, b = -10, 10
    V(x,y,z) = 1/2 * m * ω^2 * (x^2+y^2+z^2)
    EineV, v = solve3D(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(v) == nev
    # E ∼ ħω(∑nᵢ+N/2). Remember, ħ=1.
    expectedElist = [ω*(n1+n2+n3+3/2) for n1 in 0:nev-1, n2 in 0:nev-1, n3 in 0:nev-1]
    expectedElist = sort(vec(expectedElist))[1:nev]
    for i in 1:nev
        # Test the energy.
        expectedE = expectedElist[i]
        @test isapprox(expectedE, EineV[i]; rtol=1e-1)
    end
    # Test if the ground wavefunction is the expected one.
    ll = linspace(a,b,N)
    expectedwf = [QAO.(ω,m,0,x)*QAO.(ω,m,0,y)*QAO.(ω,m,0,z) for x in ll, y in ll, z in ll]
    @test_broken isapprox(abs2.(expectedwf), abs2.(v[1]), rtol=0.1)
end

@testset "Particle in a box (3D)" begin
    N = 30
    nev = 10
    m = 1
    a, b = -1, 1
    V(x,y,z) = 0
    EineV, v = solve3D(V, N=N, a=a, b=b, m=m, nev=nev)
    @test EineV == sort(EineV)
    @test length(EineV) == nev
    @test length(v) == nev
    # E ∼ n²π²ħ² / 2mL², with ħ=1
    L = b-a
    expectedElist = [(n1^2+n2^2+n3^2)*π^2/(2*m*L^2) for n1 in 1:nev, n2 in 1:nev, n3 in 1:nev]
    expectedElist = sort(vec(expectedElist))[1:nev]
    for i in 1:nev
        # Test the energy.
        expectedE = expectedElist[i]
        @test isapprox(expectedE, EineV[i]; rtol=1e-1)
    end
    # Test if the ground wavefunction is the expected one.
    ll = linspace(a,b,N)
    expectedwf = [box.(L,m,1,x)*box.(L,m,1,y)*box.(L,m,1,z) for x in ll, y in ll, z in ll]
    @test_broken isapprox(abs2.(expectedwf), abs2.(v[1]), rtol=0.1)
end
