# Ice Dynamics II: Force Balance and Models of Glacier Flow

This lecture starts from the general equations of motion (Cauchy/Navier-Stokes) and systematically derives a hierarchy of glacier flow models through scaling analysis. Beginning with the Stokes equations and Glen's flow law as a nonlinear viscosity, it introduces the aspect ratio and shows how different approximations (Full Stokes, Blatter-Pattyn, SSA, SIA) emerge from scaling. The SIA is re-derived formally, the Shallow Shelf Approximation is introduced, and two worked examples -- unconfined ice shelf spreading (the Weertman spreading rate) and laterally confined ice streams -- illustrate the SSA in practice.

## Outline

1. Introduction and motivation (connecting to Ice Dynamics I)
2. Conservation of linear momentum
   - Cauchy's equation of motion
   - The low Reynolds number limit (Stokes approximation)
   - Boundary conditions
3. The Stokes equations for glacier flow
   - Glen's flow law as a nonlinear viscosity
   - The full Stokes system
   - Properties of the Stokes system
4. Scaling analysis and the aspect ratio
   - Characteristic scales
   - Scaling the stress components
   - The scaled momentum equations
   - The model hierarchy from scaling
5. The SIA revisited: formal derivation from the scaled equations
   - When the SIA fails
6. The Shallow Shelf Approximation (SSA)
   - Physical setting and conceptual derivation (depth-averaging)
   - The SSA equations
   - When the SSA applies
   - Diagnosing the stress regime from surface profiles
7. Worked example: unconfined ice shelf spreading (Weertman spreading rate)
8. Worked example: laterally confined ice stream
9. Summary

## Appendices

- Full derivation of the Shallow Shelf Approximation
- Higher-order and hybrid models (Blatter-Pattyn, SIA+SSA, DIVA)
- Thickness profiles: the sliding-modified SIA and the SSA

## Files

- `dynamics_fundamentals_2.tex` — Lecture notes (LaTeX source)
- `dynamics_fundamentals_2.pdf` — Compiled lecture notes
- `dynamics_fundamentals_2_figures.ipynb` — Jupyter notebook generating all lecture figures
- `dynamics_fundamentals_1_figures.ipynb` — Additional figure notebook (carried from Dynamics I)
- `dynamics_fundamentals_2_exercises.ipynb` — Student exercises
- `dynamics_fundamentals_2_exercises_solutions.ipynb` — Exercise solutions (excluded from git)
- `figures/` — Generated figures (PNG and PDF)
- `profiles/` — Ice shelf and East Antarctica profile data (HDF5) and extraction script
