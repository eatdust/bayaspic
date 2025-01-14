# Bayaspic: A fast Bayesian evidence estimator for ASPIC-

### Summary

Bayaspic is a modern fortran code performing parallel estimation of
Bayesian evidences and marginalized posteriors for all inflationary
models encoded in the ASPIC library. It is using a machine-learned
likelihood and parallelisation is made in the space of models
(MPI-based). It relies on either [Multinest]()r [Polychord]() the base
sampler. Although Bayaspic is Free Software and released under the GNU
GPLv3+ public license, both Multinest and Polychord are not. As such,
they cannot be redistributed within the present source code. You are
invited to visit their respective website and abide to their
licensing conditions.

### Compilation

Users of this code are expected to be fluent in using a strongly-typed
language compiler. Please ensure that you have a working installation
of the **gfortran** and **gcc** compilers (or alternatives), the
[aspic](), [lapack](https://heasarc.gsfc.nasa.gov/fitsio/), [blas]()
and [fann]() libraries (the latter being optional). You also need the
source codes of [mlearn](), but only its "fastlike" branch for using
the data (as opposed to learning it).

The source files of [Polychord]() or/and [Multinest]() should be
present in their respective directories ("polychord/" and
"multinest/"). Notice that for building the two nested samplers
simultaneously, two source files have been symlinked.

Within "multinest/"
```
utilsMN.f90 -> utils.f90
priorsMN.f90 -> priors.f90
```

Within "polychord/src/polychord"
```
priorsPC.f90 -> priors.f90
utilsPC.F90 -> utils.F90
```

Disabling dependencies for the fann library can be set by defining by
letting the LDFANN variable empty in the Makefile (or defining the
variable -DNOFANN). Similary, dependencies to Multinest and Polychord
are disabled letting their respective variables NESTDIR and CHORDDIR
undefined.

The flags -DLIKESHEP, -DLIKERBF, -DLIKEFANN allow for setting which
machine-learned likelihoods you are using. Notice that their
respective data must be present (see below).

Editing the provided Makefile might also be needed to specify install
locations of all the needed libraries if they cannot be dynamically
resolved by the linker.

---

### Machine-learned data

The (empty) directories "shepdata/", "rbfdata/" and "fnndata/" should
be filled with specific data files encoding the weights of the machine
learned likelihood (LIKESHEP, LIKERBF and LIKEFANN). At the time if
this writing, they are made available there:
[bayaspic-data](https://curl.irmp.ucl.ac.be/~chris/bayaspic-data)

---

### Output

Bayaspic outputs the nested chains of all models within the directory
"chains/". Ensure that you have enough disk space for that as they can
fill a few hundred GB of data output.

For thinning and analysing these chains, you can either use standard
packages such as [GetDist](), or our dedicated python package
[infdistbayes](https://github.com/cosmicinflation/infdistbayes).


