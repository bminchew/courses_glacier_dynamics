# Ice Dynamics II: Force Balance and Models of Glacier Flow

**Ge 193 — Glacier Dynamics, Spring 2026**

Lecture outline for dynamics_fundamentals_2. Follows dynamics_fundamentals_1
(strain rates, Glen's law, simple shear, SIA, Vialov profile).

---

## Prerequisites from earlier lectures

- Stress tensor, deviatoric stress, stress decompositions (intro_stress)
- Mass balance, continuity equation, balance velocity (mass_balance)
- Strain rates, Glen's flow law, parallel-slab solution, SIA flux,
  Vialov profile, perfectly plastic limit (dynamics_fundamentals_1)

---

## Section 1: Introduction and Motivation

- Recap: dynamics_fundamentals_1 derived the velocity from a local force
  balance on a parallel slab, then used it to build the SIA
- The SIA works well in the interior of ice sheets but fails at margins,
  ice streams, divides, and ice shelves (Figure 6 from dynamics_fundamentals_1)
- Goal of this lecture: derive the full equations of motion, then show how
  different simplifications produce the hierarchy of glacier flow models
- Preview the model hierarchy: Full Stokes > Blatter-Pattyn > SSA > SIA
- Key question: what physics does each approximation retain and discard?

## Section 2: Conservation of Linear Momentum

*C&P 8.4*

### 2.1 Cauchy's First Law of Motion

- General form: rho * (Dv/Dt) = div(sigma) + rho*g
- Each term: inertia, stress divergence, body force (gravity)
- Written in component form (index notation and expanded)
- This is the starting point for all glacier flow models

### 2.2 The Low Reynolds Number Limit

- Estimate the Reynolds number for glacier flow:
  Re = rho * U * H / eta ~ 10^{-13}
- Inertial terms are negligible: Stokes flow
- The resulting equation: 0 = div(sigma) + rho*g
- Physical meaning: forces are always in balance; glacier flow is
  quasi-static

### 2.3 Boundary Conditions

- **Free surface (z = s):** stress-free, sigma . n = 0
  (atmospheric pressure neglected)
- **Bed (z = b):** velocity boundary condition (no-slip, sliding law,
  or free-slip); stress continuity
- Sliding laws: brief mention of Weertman-type and Coulomb-type
  (details deferred to the basal sliding lecture)
- **Lateral boundaries:** periodic, symmetry, or prescribed stress
  (e.g., calving front)

## Section 3: The Stokes Equations for Glacier Flow

*C&P 8.5.1–8.5.3*

### 3.1 Constitutive Relation: Glen's Flow Law as a Viscosity

- Recap Glen's law: dot_eps_ij = A * tau_e^{n-1} * tau_ij
- Rewrite as: tau_ij = 2 * eta * dot_eps_ij
  where eta = (1/2) * A^{-1/n} * dot_eps_e^{(1-n)/n}
- eta is the effective viscosity: depends on the strain rate itself
  (shear thinning for n > 1)
- The viscosity is isotropic but nonlinear

### 3.2 The Full Stokes System

- Combine momentum balance + constitutive law + incompressibility:
  - div(sigma) + rho*g = 0
  - sigma = -p*I + 2*eta*dot_eps
  - div(v) = 0
- This is a nonlinear elliptic PDE system for (v, p)
- 4 equations, 4 unknowns (u, v, w, p) in 3D
- Discuss: why is this hard to solve? (nonlinear viscosity, free
  surface, complex geometry, sliding law coupling)

### 3.3 Properties of the Stokes System

- No time derivative of velocity: instantaneous force balance
  (time enters only through the evolving geometry via the continuity
  equation for ice thickness)
- Stress transmission is instantaneous (no elastic waves in this model)
- The pressure is not thermodynamic — it is a Lagrange multiplier
  enforcing incompressibility

## Section 4: Scaling Analysis and the Aspect Ratio

*C&P 8.5.5–8.5.6*

### 4.1 Characteristic Scales

- Define scales: [x] = L, [z] = H, [u] = U, [w] = W
- From incompressibility: W = U * H / L = epsilon * U
  where epsilon = H / L << 1
- Typical values: H ~ 1–3 km, L ~ 100–1000 km, epsilon ~ 10^{-2}–10^{-3}

### 4.2 Scaling the Stress Components

- Scale the strain rates from the velocity scales
- Shear strain rate dot_eps_xz ~ U / H (dominates)
- Longitudinal strain rate dot_eps_xx ~ U / L = epsilon * (U / H)
- Therefore: tau_xx / tau_xz ~ epsilon (longitudinal deviatoric stress
  is O(epsilon) smaller than shear stress)
- Scale the pressure: p ~ rho * g * H (lithostatic to leading order)

### 4.3 The Scaled Momentum Equations

- Write out the x and z momentum equations with all terms scaled
  by their characteristic magnitudes
- Identify O(1), O(epsilon), O(epsilon^2) terms
- The leading-order balance in z: dp/dz = -rho*g (lithostatic)
- The leading-order balance in x: d(tau_xz)/dz = -rho*g*ds/dx
  (recovers the SIA)
- Higher-order terms: longitudinal stress gradients, lateral drag

### 4.4 The Model Hierarchy from Scaling

- Table/figure showing which terms are retained at each order:
  - **O(1) in x, O(1) in z:** SIA (vertical shear only)
  - **O(1) including membrane stresses:** SSA (depth-averaged,
    longitudinal + lateral stresses, no vertical shear)
  - **O(epsilon) corrections:** Blatter-Pattyn / first-order
  - **All terms:** Full Stokes
- Key insight: the choice of model is determined by the ratio of
  basal drag to driving stress. When basal drag is large (frozen bed),
  SIA works. When basal drag is small (sliding, ice shelves), membrane
  stresses matter and SSA or higher is needed.

## Section 5: The SIA Revisited

*C&P 8.5.6, 9.3*

### 5.1 Formal Derivation from the Scaled Equations

- Start from the scaled Stokes equations
- Retain only O(1) terms
- Recover: tau_xz = -rho*g*(s-z)*ds/dx (the slab result, now derived
  generally)
- Velocity profile: same as dynamics_fundamentals_1 but now justified
  for non-uniform geometry
- The SIA flux: q = -2A/(n+2) * (rho*g)^n * |ds/dx|^{n-1} * (ds/dx) * H^{n+2}

### 5.2 When the SIA Fails

- Ice divides: ds/dx = 0, so the SIA predicts zero velocity, but
  longitudinal stresses drive flow
- Ice streams: low basal drag, high sliding velocity; membrane stresses
  are O(1)
- Ice shelves: no basal drag at all; the SIA predicts no deformation
- Margins: steep surface slopes violate the small-slope assumption
- Grounding lines: transition from grounded to floating requires
  stress coupling

## Section 6: The Shallow Shelf Approximation (SSA) — Conceptual Derivation

*C&P 8.5.7, 10.4*

### 6.1 Physical Setting

- Ice shelves: floating slabs with no basal drag
- Ice streams: fast-flowing grounded ice with very low basal drag
- In both cases, the velocity is approximately depth-independent
  (plug flow)
- The dominant stresses are longitudinal (tau_xx) and lateral (tau_xy),
  not vertical shear (tau_xz)

### 6.2 Depth-Averaging the Momentum Equations

- Since u is approximately independent of z, integrate the momentum
  equations over depth
- The depth-integrated force balance (conceptual):
  driving stress = longitudinal stress gradient + lateral drag + basal drag
- In 1D for an ice shelf (no basal drag, no lateral drag):
  d/dx (H * tau_xx) = rho*g*H * ds/dx
  (but tau_xx is related to du/dx through the flow law)

### 6.3 The SSA Equations

- State the SSA equations (full derivation in Appendix A):
  d/dx (2*eta_bar*H * du/dx) + d/dy (eta_bar*H * (du/dy + dv/dx)) + tau_bx
  = rho*g*H * ds/dx
  (and similarly for y)
- eta_bar is the depth-averaged effective viscosity
- Discuss each term: what it represents physically
- Note: this is a 2D elliptic PDE for (u, v) — cheaper than full Stokes
  but still couples all points in the domain

### 6.4 When the SSA Applies

- Ice shelves: tau_b = 0, pure membrane flow
- Ice streams with low basal drag: tau_b is small but nonzero
- NOT appropriate for frozen-bed interior ice sheets (SIA is better there)
- Reference to full derivation: see Appendix A

## Section 7: Worked Example — Unconfined Ice Shelf Spreading

### 7.1 Setup

- 1D ice shelf of uniform thickness H, flowing in x
- No lateral confinement, no basal drag
- Calving front at x = L with ocean pressure boundary condition

### 7.2 The Calving Front Boundary Condition

- At the calving front, the ice stress must balance the ocean
  water pressure:
  sigma_xx(x=L) integrated over depth = integral of -rho_w*g*z over
  the submerged portion
- This gives the "longitudinal stress at the calving front" condition
- Derive the net stress: depends on (rho_w * D^2 - rho_i * H^2)
  where D is the draft

### 7.3 Solution

- For a freely spreading ice shelf with uniform H:
  du/dx = A * (rho_i * g * H * (1 - rho_i/rho_w) / 4)^n
- This is the "Weertman solution" for ice shelf spreading
- The strain rate is uniform and set entirely by the thickness
- Discuss: the ice shelf thins as it spreads (through the
  continuity equation)

### 7.4 Discussion

- The spreading rate depends on H^n — very sensitive to thickness
- No backstress: the ice shelf is unconfined
- In reality, ice shelves are laterally confined, producing
  buttressing (reduced spreading rate)
- This sets up the concept of buttressing for later lectures

## Section 8: Worked Example — Laterally Confined Ice Stream

### 8.1 Setup

- Ice stream of width 2W, thickness H, flowing in x
- Lateral margins provide drag; bed provides weak drag (tau_b << tau_d)
- Assume velocity varies across the stream: u = u(y)
- Neglect longitudinal gradients: d/dx = 0

### 8.2 Force Balance

- The SSA in this geometry reduces to:
  d/dy (eta_bar * du/dy) = -rho*g*ds/dx + tau_b/H
  (lateral shear stress gradient balances driving stress minus basal drag)
- With Glen's law: eta_bar * du/dy = (1/2) * A^{-1/n} * |du/dy|^{(1/n)-1} * du/dy

### 8.3 Solution for n = 1 (Newtonian)

- Linear problem: u(y) = (tau_d - tau_b) * W^2 / (2*eta*H) * (1 - (y/W)^2)
- Parabolic velocity profile across the stream
- Discuss: this is analogous to Poiseuille flow between parallel plates

### 8.4 Solution for General n

- For n = 3: the velocity profile is flatter in the center with
  sharp shear margins
- u(y) = (n/(n+1)) * A * (tau_lat)^n * W * (1 - |y/W|^{(n+1)/n})
  where tau_lat is the lateral shear stress
- The shear is concentrated at the margins — same "plug flow with
  marginal shear" as the vertical profile with large n
- Derive the relationship between driving stress, basal drag,
  lateral drag, and ice stream width

### 8.5 Discussion

- Lateral drag controls ice stream velocity when basal drag is small
- The velocity scales as W^{n+1} — wider streams flow much faster
- Shear heating at the margins can weaken the ice, creating a
  feedback (briefly mention; detail in later lectures)
- Compare with SIA: the SIA would give zero velocity for this geometry
  (no surface slope variation across the stream)

## Section 9: Summary

- The full Stokes equations are the complete description of glacier flow
  (given the constitutive law)
- The aspect ratio epsilon = H/L is the key small parameter
- The SIA retains only vertical shear: valid for slow, frozen-bed ice
  with gentle surface slopes
- The SSA retains membrane stresses and is depth-averaged: valid for
  ice shelves and fast-sliding ice streams
- The choice of model depends on the ratio of basal drag to driving
  stress
- Ice shelf spreading and confined ice streams illustrate the SSA
  and the role of boundary stresses
- Next lecture: basal processes and sliding laws

## Section 10: Reading

Cuffey & Paterson, *Physics of Glaciers*, 4th ed.:
- Section 8.4: Force balance (Cauchy's equation of motion)
- Sections 8.5.1–8.5.3: Glen's flow law as effective viscosity, Stokes equations
- Section 8.5.5–8.5.6: Scaling analysis, shallow ice approximation
- Section 8.5.7: Shallow shelf approximation
- Section 8.5.8: Higher-order and full Stokes models
- Section 10.4: Ice shelf flow and the Weertman spreading solution
- Section 10.7: Ice streams and lateral drag

---

## Appendix A: Full Derivation of the Shallow Shelf Approximation

*This appendix provides the complete SSA derivation referenced in Section 6.*

### A.1 Starting Point: The Stokes Equations

- Write the full 3D Stokes equations with Glen's law viscosity
- Boundary conditions: stress-free surface, basal drag tau_b,
  calving front pressure

### A.2 Depth-Independent Velocity Assumption

- Motivate: when basal drag is small, du/dz ~ 0
- More precisely: the vertical shear stress tau_xz is O(epsilon)
  compared to the membrane stresses
- Therefore u and v are approximately independent of z to leading order

### A.3 Vertical Integration

- Integrate the x-momentum equation from z = b to z = s
- Apply Leibniz rule for the derivative of the integral
- Apply boundary conditions at surface and bed
- Define depth-averaged quantities:
  - H = s - b (thickness)
  - Depth-averaged viscosity eta_bar from the vertically integrated
    stress-strain relation

### A.4 The Depth-Integrated Momentum Equations

- x-equation:
  d/dx(2*eta_bar*H*(2*du/dx + dv/dy)) + d/dy(eta_bar*H*(du/dy + dv/dx))
  - tau_bx = rho_i*g*H * ds/dx

- y-equation:
  d/dy(2*eta_bar*H*(2*dv/dy + du/dx)) + d/dx(eta_bar*H*(du/dy + dv/dx))
  - tau_by = rho_i*g*H * ds/dy

- Where eta_bar = (1/2) * A^{-1/n} * dot_eps_e^{(1/n)-1}
  and dot_eps_e is computed from the depth-averaged strain rates

### A.5 Boundary Condition at the Calving Front

- Derive the depth-integrated stress condition at x = x_cf:
  2*eta_bar*H*(2*du/dx + dv/dy) = (1/2)*rho_i*g*H^2*(1 - rho_i/rho_w)
- This is the "ocean pressure" boundary condition
- Physical meaning: the net horizontal force from the ice overburden
  minus the ocean back-pressure

### A.6 Properties of the SSA

- Elliptic PDE: couples all points in the domain (unlike SIA, which
  is local)
- Nonlinear through the viscosity: requires iterative solution
- 2D problem (x, y only): computationally much cheaper than full Stokes
- Captures: longitudinal stress transmission, lateral drag, buttressing
- Does NOT capture: vertical shear, bridging effects, full 3D stress
  state

### A.7 Comparison with the SIA

- Table comparing SIA and SSA:
  | Property | SIA | SSA |
  |----------|-----|-----|
  | Dominant stress | tau_xz | tau_xx, tau_xy |
  | Velocity | u(x,y,z) | u(x,y) |
  | Equation type | local (algebraic) | elliptic PDE |
  | Basal drag | = driving stress | << driving stress |
  | Applicable to | ice sheet interior | ice shelves, ice streams |

---

## Appendix B: Higher-Order and Hybrid Models

*This appendix provides a detailed derivation and discussion of the
Blatter-Pattyn first-order approximation and hybrid models, referenced
in Section 4.4.*

### B.1 Motivation

- The SIA and SSA are complementary approximations valid in different
  limits
- Many real glaciers have regions where neither is adequate:
  grounding lines, transitions from slow to fast flow, ice divides
- Higher-order models retain more terms from the Stokes equations

### B.2 The Blatter-Pattyn (First-Order) Approximation

- Start from the Stokes equations
- Make one key simplification: assume the vertical velocity is
  determined by incompressibility (no independent vertical momentum
  equation)
- Equivalently: neglect "bridging effects" — horizontal gradients
  of vertical resistive stress (d(R_zz)/dx, d(R_zz)/dy)
- The resulting equations:

  d/dx(2*eta*(2*du/dx + dv/dy)) + d/dy(eta*(du/dy + dv/dx))
  + d/dz(eta*du/dz) = rho*g * ds/dx

  d/dy(2*eta*(2*dv/dy + du/dx)) + d/dx(eta*(du/dy + dv/dx))
  + d/dz(eta*dv/dz) = rho*g * ds/dy

- This is a 3D elliptic PDE for (u, v) only (w from incompressibility,
  p eliminated)
- Boundary conditions: stress-free surface, sliding law at bed

### B.3 Properties of the First-Order Approximation

- Retains both vertical shear AND membrane stresses
- Valid across the full range from SIA to SSA limits
- 3D but only two unknowns (u, v) instead of four (u, v, w, p)
- Captures: grounding line dynamics, ice divide flow, transitions
- Error is O(epsilon^2) compared to full Stokes

### B.4 Hybrid Models: SIA + SSA

- Idea: use the SIA velocity for the depth-varying part and the SSA
  for the depth-averaged part
- Total velocity: u(x,y,z) = u_SSA(x,y) + u_SIA(x,y,z)
  where u_SIA is the SIA deformational velocity profile (zero at bed)
- This is not a formal asymptotic result but a practical approximation
- Works well in many cases because:
  - In the interior (high basal drag): u_SSA ~ 0, u_SIA dominates
    (recovers SIA)
  - On ice shelves (zero basal drag): u_SIA ~ 0, u_SSA dominates
    (recovers SSA)
  - In between: both contribute

### B.5 Implementation Considerations

- SIA is local and cheap; SSA requires solving an elliptic PDE
- Hybrid models solve the SSA globally, then add the SIA correction
  column-by-column
- Much cheaper than Blatter-Pattyn or full Stokes
- Used in many large-scale ice sheet models (e.g., PISM, BISICLES)

### B.6 The Model Hierarchy — Summary Table

| Model | Equations | Unknowns | Stresses retained | Basal drag regime | Cost |
|-------|-----------|----------|-------------------|-------------------|------|
| SIA | Algebraic (local) | u(z) per column | tau_xz only | tau_b ~ tau_d | Cheapest |
| SSA | 2D elliptic PDE | u(x,y) | tau_xx, tau_xy | tau_b << tau_d | Moderate |
| SIA+SSA hybrid | SSA + local SIA | u(x,y) + u(z) | tau_xz + tau_xx, tau_xy | Any | Moderate |
| Blatter-Pattyn | 3D elliptic PDE | u(x,y,z) | All except R_zz gradients | Any | Expensive |
| Full Stokes | 3D elliptic PDE | u,v,w,p(x,y,z) | All | Any | Most expensive |

### B.7 When to Use What

- Ice sheet projections (10^5 yr): SIA or hybrid (computational cost
  matters)
- Century-scale projections: hybrid or Blatter-Pattyn
- Grounding line dynamics: at least Blatter-Pattyn; ideally full Stokes
- Process studies (single glaciers): full Stokes is feasible
- The trend in the field is toward hybrid models with adaptive physics
  (use higher-order only where needed)

---

## Figures needed

1. **Schematic of the full Stokes force balance** — small element with
   stress components, gravity, showing the general setup
2. **Scaling diagram** — showing which stress components dominate at
   different epsilon regimes
3. **Model hierarchy diagram** — visual showing SIA, SSA, hybrid,
   Blatter-Pattyn, full Stokes as a spectrum from simple to complex,
   annotated with applicable regimes
4. **Ice shelf spreading schematic** — showing the calving front
   boundary condition, ocean pressure, spreading flow
5. **Ice shelf spreading solution** — plot of the uniform strain rate
   as a function of thickness
6. **Confined ice stream schematic** — plan view and cross-section
   showing lateral shear margins, plug flow in center
7. **Ice stream velocity profiles** — u(y) across the stream for
   n = 1, 3 (parabolic vs. plug-like), analogous to the vertical
   velocity profiles from dynamics_fundamentals_1
8. **SIA vs SSA normalized thickness profiles** — plot h/h_0 vs x/L
   for the Vialov profile (SIA, n=3) and SSA sliding profiles
   (m=0, 1/3, 1) on the same axes. Shows how the profile shape
   encodes different physics (flow law vs sliding law).
9. **SIA vs SSA comparison** — same glacier geometry solved with both
   approximations, showing where they agree and disagree

## Exercises (for dynamics_fundamentals_2_exercises.ipynb)

1. **Stokes to SIA:** Starting from the scaled momentum equations,
   identify and drop the small terms to arrive at the SIA. Verify
   by computing the velocity profile for a parallel slab.

2. **Ice shelf spreading rate:** Compute the spreading rate for an
   unconfined ice shelf as a function of thickness. Plot for
   H = 200–1000 m with n = 3. How does the spreading rate compare
   to observed values?

3. **Confined ice stream velocity:** Solve for u(y) in a confined
   ice stream for n = 1 and n = 3. Plot the velocity profiles.
   How does the centerline velocity depend on stream width?

4. **Model comparison:** For a simple ice sheet geometry, compute
   the surface velocity using (a) SIA, (b) SSA, (c) SIA+SSA hybrid.
   Plot and compare. Where do they agree? Where do they diverge?

5. **Aspect ratio exploration:** For a real glacier geometry (e.g.,
   a profile from the East Antarctic data), compute the local aspect
   ratio and identify regions where the SIA is and is not valid.
