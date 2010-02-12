% RNP Documentation -- Compilation and file dependencies
% Victor Liu (vkl@stanford.edu)
% February 12, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# Compilation and file dependencies

Compilation is straightforward: copy the files you need into your project and compile them, that's all there is to it.
The number of files needed is typically minimal, and there are no additional dependencies.

## Dependencies

Usage of any function from one of the modules may entail including additional files into your project.
The dependencies are listed below.

### Header only

* IO
* Random
* TBLAS
* Sparse
* LinearSolve - Also requires TBLAS.h and TLASupport.h

### Sources required

Obviously, the corresponding header files are needed, so those will not be listed.

* Eigensystems - Needs Eigensystems.cpp, TLASupport.h, TBLAS.h
* Generalized Eigensystems - Needs GeneralizedEigensystems.cpp, TLASupport.h, TBLAS.h
* Iterative Linear Solvers - Needs TLASupport and TBLAS.h, as well as the corresponding cpp file from the src/ directory corresponding to the routine name.
* Iterative Generalized Eigensystems - Needs IterativeLinearSolvers.h, Random.h, TLASupport, TBLAS.h, GeneralizedEigensystems.h, GeneralizedEigensystems.cpp, and the corresponding cpp files from the src/ directory corresponding to the routine names used (JDQZ and whatever linear solver).
