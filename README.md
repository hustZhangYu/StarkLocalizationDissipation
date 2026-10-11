# Stark-localization-dissipation
  An external electric field can lead to the localization of the eigensttaes named stark localization. 
  The eigenstates and eigenvalues can be solved analytically. 
  Stark localization show interesting dynamic behaviors like bloch oscillation and dynamic localization. I want to investigate what will happen if we add a dissiapation term to the system.
  Liouvillan gap can be used to determine the timesacle that an open system arrive the steady state. 
  The first part ,we confirm the Liouvillan gap is independent on the particle number for free fermions.(L_N,L_N1)

## Entry points and runtime

Use MATLAB with support for local functions in scripts (R2016b or later).
From the repository root, add the root to the path and run either script:

```matlab
addpath(pwd)
L_N
% Or, in a fresh workspace:
L_N1
```

`L_N.m` builds the full occupation basis; `L_N1.m` builds sectors up to
a chosen particle cutoff. They retain their own parameter defaults and should
not be treated as identical experiments. The default size is L=5. Matrices
and their Liouvillian Kronecker products grow rapidly with L.

`BoundaryLossOperator.m` is the shared endpoint builder used by `L_N1.m`.
Its site argument replaces the duplicated La/Lb routines while preserving
their sector order and occupation-basis signs. The second dissipator remains
commented out in the experiment, as in the original.

```matlab
results = runtests('tests');
assert(all([results.Passed]));
```

The regression test checks exact matrix equality with the original endpoint
builders for L=1 through 6. It does not establish the reported physical gap
independence or validate the dissipative model.

