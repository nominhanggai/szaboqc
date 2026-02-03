# Appendix B Source Code

**Short summary:**  
This repository contains two Fortran implementations of the Appendix B example used for HeH<sup>+</sup> SCF.

- `AppendixBcode.f90` original, Fortran-77 style reference source.  
- `AppendixBcode_revised.f90` Fortran-90 refactor by the author that preserves original variables and algorithms, rewrites control flow (removes GOTOs; uses `DO`/`IF`), and performs diagonalization using LAPACK.

---

## Key changes (changelog)
The file `AppendixBcode_revised.f90` includes the following major modifications relative to `AppendixBcode.f90`:

1. **Removed jumps and rewrote control flow**  
   All unstructured jumps (e.g., `GOTO`) were removed and rewritten using structured `DO` and `IF` blocks.

2. **Retained original variables and integral calculations**  
   The original variable definitions and electronic integral evaluation procedures were preserved to allow direct comparison with the reference code.

3. **LAPACK-based diagonalization**  
   Matrix diagonalization is performed using LAPACK routines (e.g., `DSYEV`).  
   Some variables are redundant and have intentionally not been removed to facilitate validation against the original implementation.

---
## Requirements
- Fortran compiler supporting Fortran-90
- LAPACK library (and BLAS, depending on system configuration)

---

## Build / Compilation

### Linux (example)
```bash
gfortran -o HeH AppendixBcode_revised.f90 -llapack



