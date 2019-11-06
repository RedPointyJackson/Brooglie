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
V(x) = 0.5*m*ω^2*x^2

# Find the first 5 eigenstates
E, vv = solve(V, a=a, b=b, m=m, nev=5)
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

In the array `vv` we find 5 vectors (1 dimensional), one for each
eigenstate:

![alt text](https://github.com/RedPointyJackson/Brooglie/blob/master/harmonic_eigen.png "Eigenstates of the harmonic oscillator")

Read the help of `solve` to see all the available options. It accepts
potentials with arbitrary number of dimensions, like for example
`V(x,y,z,k,w)`, and creates the Hamiltonian accordingly (read the
documentation for the methods used under the `doc/` directory),
reshaping the `vv` vectors into the correct number of dimensions. For
example, if `V(x,y) = x+y` is used as potential the array `vv` will
contain matrices.

If the Hamiltonian is needed for a more detailed calculation, it can
be returned skipping the diagonalization with the `buildH` function.

### TODO

- [x] Fix all the tests, some will need to tweak `N` or `maxiters` after the last commits.
- [ ] Add an hydrogen atom to the tests.
