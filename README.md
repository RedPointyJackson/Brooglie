# Brooglie

This software solves the time-independent Schödinger equation for the
1D, 2D and 3D cases employing a finite difference method. For example,
to solve the 1D quantum harmonic oscillator one can use the following
snippet:

```julia
using Brooglie

m = 1
ω = 1
a,b = -10,10 # Limits of the box

# Use atomic units!
V(x) = 0.5*m*ω*x^2

# Find the first 5 eigenstates
E, v = solve1D(V, a=a, b=b, m=m, nev=5)
```

`E` will be in Hartree (atomic units), remember that one Hartree is
equivalent to 27.21 eV. As expected, we obtain E≃ω(n+½) noticing that
ℏ=1 in these units:

|  Analytical     | Numerical  |
| :-------------: | :--------: |
|  0.5            |  0.500952  |
|  1.5            |  1.50275   |
|  2.5            |  2.50436   |
|  3.5            |  3.50576   |
|  4.5            |  4.50696   |



In the array `v` we find 5 vectors, one for each eigenstate:

![alt text](https://github.com/RedPointyJackson/Brooglie/blob/master/harmonic_eigen.png "Eigenstates of the harmonic oscillator")


Read the help of `solve1D`, `solve2D` and `solve3D` to see all the
available options.

### TODO

- [ ] Fix all the tests, some will need to tweak `N` or `maxiters` after the last commits.
- [ ] Add an hydrogen atom to the tests.
