# Bayaspic: fast bayesian evidences for ASPIC

### Summary

Bayaspic is a modern fortran code performing parallel estimation of
Bayesian evidences and marginalized posteriors for all inflationary
models encoded in the ASPIC library. It is using a machine-learned
likelihood and parallelisation is made in the space of models
(MPI-based). It relies on either
[Multinest](https://github.com/farhanferoz/MultiNest) or
[Polychord](https://github.com/PolyChord/) as the base sampler. Although
Bayaspic is Free Software and released under the GNU GPLv3+ public
license, both Multinest and Polychord are not. As such, they cannot be
redistributed within the present source code. You are invited to visit
their respective website and abide to their licensing restrictions.

### Compilation

Users of this code are expected to be fluent in using a strongly-typed
language compiler. Please ensure that you have a working installation
of the **gfortran** and **gcc** compilers (or alternatives), the
[aspic](https://github.com/cosmicinflation/aspic),
[lapack/blas](https://github.com/Reference-LAPACK/lapack) and
[fann](https://github.com/libfann/fann) libraries (the latter being
optional). You also need the source files of
[mlearn](https://github.com/eatdust/mlearn/tree/fastlikes), but only
its "fastlike" branch for using the data (as opposed to learning it).

The source files of
[Polychord](https://github.com/PolyChord/PolyChordLite.git) or/and
[Multinest](https://github.com/farhanferoz/MultiNest/tree/master/MultiNest_v3.12)
should be present in their respective directories, "polychord/" and
"multinest/". Notice that for building the two nested samplers
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

Disabling dependencies for the fann library can be set by letting the
LDFANN variable empty in the Makefile (or setting the flag -DNOFANN at
build time). Similarly, dependencies to Multinest and Polychord are
disabled by letting their respective variables NESTDIR and CHORDDIR
undefined.

The flags -DLIKESHEP, -DLIKERBF, -DLIKEFANN allow for setting the type
of machine-learned likelihood you are using. Notice that their
respective data must then be present (see below).

Editing the provided Makefile might also be needed for specific
install locations if shared libraries cannot be dynamically resolved.

---

### Machine-learned data

The (empty) directories "shepdata/", "rbfdata/" and "fnndata/" should
be filled with specific data files encoding the weights of the
machine- learned likelihood (LIKESHEP, LIKERBF and LIKEFANN). At the
time if this writing, they are made available there:
[bayaspic-data](https://curl.irmp.ucl.ac.be/~chris/upload/bayaspic-data)

---

### Output

Bayaspic outputs the nested chains of all models within the directory
"chains/". Ensure that you have enough disk space for that as they can
fill a few hundred GB of data output.

For thinning and analysing these chains, you can either use standard
packages such as [GetDist](https://github.com/cmbant/getdist), or our dedicated python package
[infdistbayes](https://github.com/cosmicinflation/infdistbayes).


