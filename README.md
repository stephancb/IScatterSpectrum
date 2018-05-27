# IScatterSpectrum

[![Build Status](https://travis-ci.org/stephancb/IScatterSpectrum.jl.svg?branch=master)](https://travis-ci.org/stephancb/IScatterSpectrum.jl)

[![Coverage Status](https://coveralls.io/repos/stephancb/IScatterSpectrum.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/stephancb/IScatterSpectrum.jl?branch=master)

[![codecov.io](http://codecov.io/github/stephancb/IScatterSpectrum.jl/coverage.svg?branch=master)](http://codecov.io/github/stephancb/IScatterSpectrum.jl?branch=master)

IScatterSpectrum.jl
===================

This package provides tools to calculate the incoherent scatter (IS) spectrum
for a ionospheric plasma.

Facilities like [EISCAT_3D](https://eiscat3d.se//) and the [Advanced
Modular IS Radar](http://amisr.com/amisr/about/amisr-overview/) form multiple
beams, and so greatly increase the amount of data to process. The aim of this
package is to develop and test software for calculating and fitting the IS
spectrum using the high-performance Julia system and its parallel execution.

A reference is the article

Swartz, W. E. (1978), Analytic partial derivatives for least‐squares fitting
incoherent scatter data, Radio Sci., 13(3), 581–589, doi:10.1029/RS013i003p00581.

and FORTRAN IV code written by the same author,
[iscatspe.for](https://github.com/stephancb/IScatterSpectrum.jl/blob/master/src/iscatspe.for).
[ISfortran.jl](https://github.com/stephancb/IScatterSpectrum.jl/blob/master/src/ISfortran.jl)
shows how the reference FORTRAN code can be called from Julia. Presently only
the functions and COMMON BLOCK addresses to call the plasma dispersion function
are provided. In the FORTRAN code the plasm dispersion
function is calculated following the prescription in the book

The Plasma Dispersion Function: The Hilbert Transform of the Gaussian, by Fried
and Conte, Academic Press, 1961.

Over most of the part which is interesting for the IS spectrum the prescription
integrates a differential equation, which does not lend itself easily to
parallel execution.

In
[IScatterSpectrum.jl](https://github.com/stephancb/IScatterSpectrum.jl/blob/master/src/IScatterSpectrum.jl)
`erfcx` from the
[SpecialFunctions](http://juliamath.github.io/SpecialFunctions.jl/stable/index.html)
package is used. Also `faddeeva` using the method of

J.A.C. Weideman, "Computation of the Complex Error Function," SIAM
J. Numerical Analysis, pp. 1497-1518, No. 5, Vol. 31, Oct., 1994 

is provided for verification and experimenting. For physical constants the
unregistered
[PhysicalConstants](https://github.com/DataWookie/PhysicalConstants.jl) is
required.

`struct ScatterVolume` precalculates parameters that change when processing along a
beam or switching beams, but are constant in a fit. `struct Plasma` and
`struct CollisionalPlasma` precalculate parameters that change in the fit
(depending on N<sub>e</sub>, $T_e$, $T_i$, $\nu_i$, composition), but are constant at
all the frequencies.

So far only the power spectrum can be calculated: 
`pwrspec(freq, p::Plasma, s::ScatterVolume)` and `pwrspec(p::Plasma,
s::ScatterVolume)`, the latter for zero frequency which requires special
numerical treatment in the admittance formulation used here. The same functions
are available for `CollisionalPlasma`.

[plottest.jl](./test/plottest.jl) (in the test folder) creates plots of some
example spectra, using the [PlotlyJS](http://spencerlyon.com/PlotlyJS.jl/)
package.

To do
-----

* add the Fourier transform to ACFs. Swartz (1978) prefers a non-fast numerical 
  integration (because only few ACF lags are needed, but many frequencies). But
  this seems not ideal for parallel execution. Several articles on methods for fast
  calculation of Fourier integrals have in the meantime been published.

* add the analytical derivatives with respect to  $N_e$, $T_e$, $T_i$, $\nu_i$,
  composition. The conclusion in 1978 was that numerical derivatives are roughly
  equally fast to calculate, but the analytical formulas seem still preferable,
  after they have been laboriously calculated.

* explore possibilities for parallel execution.

* add fitting to real data, with constraints from models, output of results etc.
