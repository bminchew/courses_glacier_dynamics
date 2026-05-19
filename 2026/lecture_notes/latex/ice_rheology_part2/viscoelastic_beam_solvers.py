"""
Viscoelastic beam bending solvers for ice-shelf tidal flexure.

Implements a Hermite-cubic finite element solver for the 1D Maxwell
viscoelastic flexure equation on a floating beam with clamped boundary
conditions at the grounding line:

    D_v  w''''_dot  +  rho_w g w  +  t_r rho_w g w_dot  =  q  +  t_r q_dot

where t_r = D_v / D is the flexural relaxation time (proportional to
the Maxwell time), D_v is the viscous flexural rigidity, and D is the
elastic flexural rigidity.

The 1D equation assumes cylindrical bending (eps_yy = 0), appropriate
for ice shelves wider than the flexural wavelength.  This is the
correct 1D reduction of the full plate equation: for Newtonian viscous
flow the Poisson coupling terms cancel in the biharmonic operator
(nu_v = 1/2), so D_v = eta h^3 / 3 is the same for beam and plate.
The elastic rigidity D = E h^3 / [12(1-nu^2)] already includes the
plate plane-stress factor.  For power-law viscosity, the cylindrical
bending constraint (eps_yy = 0) keeps the stress state uniaxial,
avoiding the nonlinear biaxial coupling that makes the general 2D
plate intractable (see Appendix E of the lecture notes).

Supports uniform and graded (non-uniform) meshes.  The per-element
matrix storage enables future extension to variable ice thickness h(x).

Also provides the analytical steady-state solution for the Newtonian
viscoelastic case under monochromatic tidal forcing.
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


# =====================================================================
# Analytical solutions
# =====================================================================

def viscoelastic_complex_amplitude(x, w0, D_v, D_el, omega, rho_w, g):
    """Steady-state complex amplitude W(x) for the Newtonian Maxwell beam.

    The viscoelastic flexure equation under harmonic tidal forcing
    q = rho_w g w0 sin(omega t) becomes, after substituting
    w = Im[W(x) exp(i omega t)]:

        i omega D_v W'''' + rho_w g (1 + i omega t_r) W
            = rho_w g (1 + i omega t_r) w0

    Parameters
    ----------
    x : array_like
        Distances from grounding line (m).
    w0 : float
        Tidal amplitude (m).
    D_v : float
        Viscous flexural rigidity (Pa m^3 s).
    D_el : float
        Elastic flexural rigidity (Pa m^3).
    omega : float
        Tidal angular frequency (rad/s).
    rho_w, g : float
        Seawater density (kg/m^3) and gravitational acceleration (m/s^2).

    Returns
    -------
    W : complex ndarray
        Complex amplitude at each x.
    """
    t_r = D_v / D_el
    alpha = 1 + 1j * omega * t_r

    coeff = -rho_w * g * alpha / (1j * omega * D_v)
    base = coeff ** 0.25
    roots = [base * np.exp(1j * j * np.pi / 2) for j in range(4)]

    decaying = sorted([r for r in roots if r.real < 0], key=lambda r: r.imag)
    lam1, lam2 = decaying

    C1 = w0 * lam2 / (lam1 - lam2)
    C2 = -w0 * lam1 / (lam1 - lam2)

    x = np.asarray(x, dtype=float)
    W = w0 + C1 * np.exp(lam1 * x) + C2 * np.exp(lam2 * x)
    return W


def viscoelastic_analytical(x, t, w0, D_v, D_el, omega, rho_w, g):
    """Analytical steady-state deflection for Newtonian Maxwell beam.

    Parameters
    ----------
    x : array_like
        Distances from grounding line (m).
    t : float or array_like
        Time(s) (s).
    w0 : float
        Tidal amplitude (m).
    D_v, D_el : float
        Viscous and elastic flexural rigidities.
    omega : float
        Tidal angular frequency (rad/s).
    rho_w, g : float
        Seawater density and gravitational acceleration.

    Returns
    -------
    w : ndarray
        Deflection w(x, t).
    """
    W = viscoelastic_complex_amplitude(x, w0, D_v, D_el, omega, rho_w, g)
    t = np.asarray(t)
    if t.ndim == 0:
        return np.imag(W * np.exp(1j * omega * float(t)))
    return np.imag(W[None, :] * np.exp(1j * omega * t[:, None]))


def viscoelastic_analytical_curvature(x, t, w0, D_v, D_el, omega, rho_w, g):
    """Analytical steady-state curvature d^2w/dx^2 for Newtonian Maxwell beam."""
    t_r = D_v / D_el
    alpha = 1 + 1j * omega * t_r
    coeff = -rho_w * g * alpha / (1j * omega * D_v)
    base = coeff ** 0.25
    roots = [base * np.exp(1j * j * np.pi / 2) for j in range(4)]
    decaying = sorted([r for r in roots if r.real < 0], key=lambda r: r.imag)
    lam1, lam2 = decaying
    C1 = w0 * lam2 / (lam1 - lam2)
    C2 = -w0 * lam1 / (lam1 - lam2)
    x = np.asarray(x, dtype=float)
    W_pp = C1 * lam1**2 * np.exp(lam1 * x) + C2 * lam2**2 * np.exp(lam2 * x)
    t = np.asarray(t)
    if t.ndim == 0:
        return np.imag(W_pp * np.exp(1j * omega * float(t)))
    return np.imag(W_pp[None, :] * np.exp(1j * omega * t[:, None]))


# =====================================================================
# Graded mesh generation
# =====================================================================

def generate_graded_mesh(L, ell_f, h_min=25.0, h_max=250.0, Nel_max=160):
    """Generate a graded mesh with denser spacing near the grounding line.

    Element spacing follows the elastic curvature envelope:

        h(x) = h_min + (h_max - h_min) * (1 - exp(-x/ell_f))^p

    where p is chosen by binary search so that the total element count
    equals Nel_max.

    Parameters
    ----------
    L : float
        Domain length (m).
    ell_f : float
        Elastic flexural wavelength (m).
    h_min : float
        Minimum element size near x=0 (m).
    h_max : float
        Maximum element size in the far field (m).
    Nel_max : int
        Target number of elements.

    Returns
    -------
    x_nodes : ndarray
        Node positions, shape (Nel+1,).
    """
    def make_nodes(p):
        nodes = [0.0]
        while nodes[-1] < L:
            x = nodes[-1]
            s = 1.0 - np.exp(-x / ell_f)
            spacing = h_min + (h_max - h_min) * s**p
            nodes.append(min(x + spacing, L))
            if len(nodes) > 10 * Nel_max:
                break
        return np.array(nodes)

    # Binary search for the largest p giving Nel <= Nel_max
    p_lo, p_hi = 0.01, 50.0
    for _ in range(100):
        p_mid = (p_lo + p_hi) / 2
        n = make_nodes(p_mid)
        if len(n) - 1 > Nel_max:
            p_hi = p_mid
        else:
            p_lo = p_mid

    nodes = make_nodes(p_lo)

    # If last element is tiny, merge with previous
    if len(nodes) > 2 and (nodes[-1] - nodes[-2]) < h_min / 2:
        nodes = np.delete(nodes, -2)

    return nodes


# =====================================================================
# FEM infrastructure (Hermite cubics, uniform or non-uniform mesh)
# =====================================================================

class HermiteFEMBeam:
    """Hermite cubic finite element discretization for a 1D beam.

    Degrees of freedom at each node: (w, dw/dx).  Two-point Gauss
    quadrature per element.  Clamped BC at x = 0 (DOFs 0 and 1 fixed).

    Supports uniform and non-uniform meshes.  All reference matrices
    are stored per-element to enable future extension to variable
    ice thickness h(x).

    Parameters
    ----------
    L : float, optional
        Domain length (m).  Used with Nel for a uniform mesh.
    Nel : int, optional
        Number of elements for a uniform mesh.
    x_nodes : array_like, optional
        Node positions for a non-uniform mesh.  Overrides L and Nel.
    """

    def __init__(self, L=None, Nel=None, x_nodes=None):
        if x_nodes is not None:
            self.x_nodes = np.asarray(x_nodes, dtype=float)
            self.Nel = len(self.x_nodes) - 1
            self.L = self.x_nodes[-1] - self.x_nodes[0]
        else:
            self.L = L
            self.Nel = Nel
            self.x_nodes = np.linspace(0, L, Nel + 1)
        self.he = np.diff(self.x_nodes)  # per-element sizes, shape (Nel,)
        self.ndof = 2 * (self.Nel + 1)
        self._build_reference()

    def _build_reference(self):
        Nel = self.Nel
        he = self.he  # (Nel,)

        gp = np.array([0.5 - np.sqrt(3) / 6, 0.5 + np.sqrt(3) / 6])
        gw = np.array([0.5, 0.5])

        # Reference shape functions on [0,1] (element-independent)
        N_ref = np.zeros((2, 4))
        d2N_ref = np.zeros((2, 4))
        for q in range(2):
            xi = gp[q]
            N_ref[q] = [1 - 3*xi**2 + 2*xi**3, xi*(1-xi)**2,
                        3*xi**2 - 2*xi**3, xi**2*(xi-1)]
            d2N_ref[q] = [-6 + 12*xi, -4 + 6*xi, 6 - 12*xi, -2 + 6*xi]

        # Per-element scaling: slope DOFs scale with local he
        # scale[e] = [1, he[e], 1, he[e]]
        scale = np.column_stack([np.ones(Nel), he, np.ones(Nel), he])

        # Physical shape functions and second derivatives per element
        # N_phys[e, q, i], B_phys[e, q, i]
        self.N_phys = N_ref[None, :, :] * scale[:, None, :]  # (Nel, 2, 4)
        self.B_phys = (d2N_ref[None, :, :] / he[:, None, None]**2
                       * scale[:, None, :])  # (Nel, 2, 4)

        # Per-element reference matrices, shape (Nel, 2, 4, 4) or (Nel, 2, 4)
        self.K_ref = np.zeros((Nel, 2, 4, 4))
        self.F_ref = np.zeros((Nel, 2, 4))
        self.M_ref_gp = np.zeros((Nel, 2, 4, 4))

        for q in range(2):
            Bq = self.B_phys[:, q, :]   # (Nel, 4)
            Nq = self.N_phys[:, q, :]   # (Nel, 4)
            jac = gw[q] * he            # (Nel,)
            self.K_ref[:, q] = jac[:, None, None] * (
                Bq[:, :, None] * Bq[:, None, :])
            self.F_ref[:, q] = jac[:, None] * Nq
            self.M_ref_gp[:, q] = jac[:, None, None] * (
                Nq[:, :, None] * Nq[:, None, :])

        self.M_ref = self.M_ref_gp[:, 0] + self.M_ref_gp[:, 1]  # (Nel, 4, 4)

        # Connectivity (unchanged from uniform case)
        rows = np.zeros((Nel, 4, 4), dtype=int)
        cols = np.zeros((Nel, 4, 4), dtype=int)
        rhs_rows = np.zeros((Nel, 4), dtype=int)
        for e in range(Nel):
            d = np.array([2*e, 2*e+1, 2*e+2, 2*e+3])
            rows[e] = d[:, None] * np.ones(4, dtype=int)[None, :]
            cols[e] = np.ones(4, dtype=int)[:, None] * d[None, :]
            rhs_rows[e] = d

        self.rows_flat = rows.ravel()
        self.cols_flat = cols.ravel()
        self.rhs_rows = rhs_rows
        self._gp = gp


# =====================================================================
# Viscoelastic FEM solvers
# =====================================================================

def solve_wdot_viscoelastic(fem, w_curr, w_tide, dw_tide_dt, D_v, dt,
                            rho_w, g, t_relax=0.0):
    """Solve for dw/dt at the current time step (Newtonian viscoelastic).

    Crank-Nicolson discretization of:

        D_v w''''_dot + rho_w g w + t_r rho_w g w_dot = q + t_r q_dot

    Evaluating the buoyancy at the midpoint w^{n+1/2} = w^n + dt/2 * wdot
    and the forcing at t + dt/2 gives the second-order system:

        [D_v K + (dt/2 + t_r) rho_w g M] wdot
            = rho_w g [(w_tide^{n+1/2} - w^n) + t_r dw_tide/dt^{n+1/2}] F

    Parameters
    ----------
    fem : HermiteFEMBeam
    w_curr : ndarray, shape (Nel+1,)
        Current nodal deflections.
    w_tide : float
        Tidal elevation at midpoint (m).
    dw_tide_dt : float
        Time derivative of tidal elevation at midpoint (m/s).
    D_v : float or ndarray
        Viscous flexural rigidity.  Scalar for uniform thickness;
        array of shape (Nel,) for variable thickness (future).
    dt : float
        Time step (s).
    rho_w, g : float
    t_relax : float
        Flexural relaxation time D_v / D_el (s).  0 for pure viscous.

    Returns
    -------
    wdot_vec : ndarray, shape (ndof,)
    """
    Nel = fem.Nel
    ndof = fem.ndof

    # Crank-Nicolson: buoyancy at midpoint w^{n+1/2} = w^n + dt/2 * wdot
    coeff_buoy = (dt / 2 + t_relax) * rho_w * g
    # D_v can be scalar or per-element array (for variable thickness)
    D_v_e = np.broadcast_to(D_v, (Nel,))
    K_elem = (D_v_e[:, None, None] * (fem.K_ref[:, 0] + fem.K_ref[:, 1])
              + coeff_buoy * fem.M_ref)

    K_global = csr_matrix(
        (K_elem.ravel(), (fem.rows_flat, fem.cols_flat)),
        shape=(ndof, ndof)
    ).tolil()

    # RHS load at Gauss points
    gp = fem._gp
    w_gp = (w_curr[:-1, None] * (1 - gp[None, :])
            + w_curr[1:, None] * gp[None, :])
    load = rho_w * g * ((w_tide - w_gp) + t_relax * dw_tide_dt)

    F_elem = (load[:, 0, None] * fem.F_ref[:, 0]
              + load[:, 1, None] * fem.F_ref[:, 1])
    F_global = np.zeros(ndof)
    np.add.at(F_global, fem.rhs_rows.ravel(), F_elem.ravel())

    # Clamped BCs at x = 0
    for bc in [0, 1]:
        K_global[bc, :] = 0
        K_global[:, bc] = 0
        K_global[bc, bc] = 1.0
        F_global[bc] = 0.0

    return spsolve(K_global.tocsr(), F_global)


def solve_wdot_viscoelastic_powerlaw(fem, w_curr, w_tide, dw_tide_dt,
                                     D_n, n_nl, kappa_bg, dt,
                                     rho_w, g, D_el=None, t_relax=0.0,
                                     wdot_prev=None, picard_iters=3):
    """Solve for dw/dt with power-law viscosity and optional elasticity.

    For n > 1 the effective viscous rigidity depends on the curvature rate,
    requiring Picard iteration.  The flexural relaxation time is also
    curvature-rate dependent: t_relax_eff = D_eff / D_el at each Gauss
    point.  For n = 1, D_eff = D_v everywhere and t_relax_eff reduces to
    the constant t_relax = D_v / D_el.

    Parameters
    ----------
    fem : HermiteFEMBeam
    w_curr : ndarray, shape (Nel+1,)
    w_tide : float
        Tidal elevation at midpoint.
    dw_tide_dt : float
        Time derivative of tidal elevation at midpoint.
    D_n : float or ndarray
        Nonlinear flexural rigidity parameter.  Scalar for uniform
        thickness; array of shape (Nel,) for variable thickness (future).
    n_nl : float
        Power-law exponent (1 = Newtonian).
    kappa_bg : float or ndarray
        Background curvature rate for regularization.  Scalar for uniform
        thickness; array of shape (Nel,) for variable thickness (future).
    dt : float
    rho_w, g : float
    D_el : float, ndarray, or None
        Elastic flexural rigidity.  Scalar for uniform thickness;
        array of shape (Nel,) for variable thickness (future).
    t_relax : float
        Flexural relaxation time (Newtonian fallback).  0 for pure viscous.
    wdot_prev : ndarray or None
    picard_iters : int

    Returns
    -------
    wdot_vec : ndarray, shape (ndof,)
    """
    Nel = fem.Nel
    ndof = fem.ndof

    # Broadcast scalar parameters to per-element arrays (enables variable h)
    D_n_e = np.broadcast_to(D_n, (Nel,))
    kbg_e = np.broadcast_to(kappa_bg, (Nel,))
    D_el_e = (np.broadcast_to(D_el, (Nel,))
              if D_el is not None else None)

    wdot_vec = wdot_prev.copy() if wdot_prev is not None else np.zeros(ndof)

    # Interpolate w to Gauss points (constant across Picard iterations)
    gp = fem._gp
    w_gp = (w_curr[:-1, None] * (1 - gp[None, :])
            + w_curr[1:, None] * gp[None, :])

    n_iters = picard_iters if n_nl > 1 else 1
    for iteration in range(n_iters):
        # Compute effective D and local relaxation time at each GP
        D_eff = np.column_stack([D_n_e, D_n_e])  # (Nel, 2)
        t_relax_gp = np.full((Nel, 2), t_relax)

        if n_nl > 1:
            for e in range(Nel):
                dv = wdot_vec[2*e:2*e+4]
                for q in range(2):
                    wpp = fem.B_phys[e, q] @ dv
                    eps_eff = np.sqrt(wpp**2 + kbg_e[e]**2)
                    D_eff[e, q] = D_n_e[e] * eps_eff**((1.0/n_nl) - 1)
                    if D_el_e is not None:
                        t_relax_gp[e, q] = D_eff[e, q] / D_el_e[e]

        # RHS load with local t_relax
        load = rho_w * g * ((w_tide - w_gp) + t_relax_gp * dw_tide_dt)
        F_elem = (load[:, 0, None] * fem.F_ref[:, 0]
                  + load[:, 1, None] * fem.F_ref[:, 1])
        F_global = np.zeros(ndof)
        np.add.at(F_global, fem.rhs_rows.ravel(), F_elem.ravel())

        # Assemble: Crank-Nicolson buoyancy at midpoint w^{n+1/2}
        buoy_gp = rho_w * g * (dt / 2 + t_relax_gp)  # (Nel, 2)
        K_elem = (D_eff[:, 0, None, None] * fem.K_ref[:, 0]
                  + D_eff[:, 1, None, None] * fem.K_ref[:, 1]
                  + buoy_gp[:, 0, None, None] * fem.M_ref_gp[:, 0]
                  + buoy_gp[:, 1, None, None] * fem.M_ref_gp[:, 1])

        K_global = csr_matrix(
            (K_elem.ravel(), (fem.rows_flat, fem.cols_flat)),
            shape=(ndof, ndof)
        ).tolil()

        # Clamped BCs
        F_bc = F_global.copy()
        for bc in [0, 1]:
            K_global[bc, :] = 0
            K_global[:, bc] = 0
            K_global[bc, bc] = 1.0
            F_bc[bc] = 0.0

        wdot_vec = spsolve(K_global.tocsr(), F_bc)

    return wdot_vec


# =====================================================================
# Time-stepping driver
# =====================================================================

def run_viscoelastic_tidal(fem, D_v, D_el, w0, omega, rho_w, g,
                           dt=30.0, n_cycles=6, n_nl=1, A_glen=None,
                           h_plate=None, kappa_bg=None, picard_iters=3,
                           save_last_n_cycles=1):
    """Run the viscoelastic FEM solver under tidal forcing.

    Parameters
    ----------
    fem : HermiteFEMBeam
    D_v : float
        Viscous flexural rigidity (Newtonian).
    D_el : float
        Elastic flexural rigidity.  Set to np.inf for pure viscous.
    w0 : float
        Tidal amplitude (m).
    omega : float
        Tidal angular frequency (rad/s).
    rho_w, g : float
    dt : float
        Time step (s).
    n_cycles : int
        Number of tidal cycles to run.
    n_nl : float
        Power-law exponent (1 = Newtonian).
    A_glen : float or None
        Glen rate factor.  Required if n_nl > 1.
    h_plate : float or None
        Plate thickness.  Required if n_nl > 1.
    kappa_bg : float or None
        Background curvature rate.  Required if n_nl > 1.
    picard_iters : int
    save_last_n_cycles : int

    Returns
    -------
    result : dict with keys:
        'w_hist', 't_hist', 'w_final', 'wdot_prev', 'D_n'
    """
    T_tide = 2 * np.pi / omega
    t_relax = D_v / D_el if np.isfinite(D_el) else 0.0

    # Determine D_n for power-law case
    # (For variable thickness, h_plate and hence D_n would be per-element.)
    if n_nl > 1:
        B_coeff = 2.0 / A_glen**(1.0 / n_nl)
        D_n = (2 * n_nl / (1 + 2 * n_nl)) * B_coeff \
              * (h_plate / 2)**((1.0/n_nl) + 2)
    else:
        D_n = D_v

    n_steps = int(n_cycles * T_tide / dt)
    steps_per_cycle = int(T_tide / dt)
    save_every = max(1, steps_per_cycle // 200)
    save_start = (n_cycles - save_last_n_cycles) * steps_per_cycle

    w = np.zeros(fem.Nel + 1)
    wdot_prev = None
    w_hist = []
    t_hist = []

    for step in range(n_steps):
        t_curr = step * dt
        # Crank-Nicolson: evaluate forcing at midpoint t + dt/2
        t_mid = t_curr + dt / 2
        w_tide = w0 * np.sin(omega * t_mid)
        dw_tide_dt = w0 * omega * np.cos(omega * t_mid)

        if n_nl > 1:
            wdot_vec = solve_wdot_viscoelastic_powerlaw(
                fem, w, w_tide, dw_tide_dt,
                D_n, n_nl, kappa_bg, dt, rho_w, g,
                D_el=D_el if np.isfinite(D_el) else None,
                t_relax=t_relax, wdot_prev=wdot_prev,
                picard_iters=picard_iters
            )
        else:
            wdot_vec = solve_wdot_viscoelastic(
                fem, w, w_tide, dw_tide_dt,
                D_n, dt, rho_w, g, t_relax=t_relax
            )

        wdot_prev = wdot_vec
        w = w + dt * wdot_vec[0::2]
        w[0] = 0.0

        if step >= save_start and step % save_every == 0:
            w_hist.append(w.copy())
            t_hist.append(t_curr)

    return {
        'w_hist': np.array(w_hist),
        't_hist': np.array(t_hist),
        'w_final': w.copy(),
        'wdot_prev': wdot_prev,
        'D_n': D_n,
    }


# =====================================================================
# Convenience: Deborah number sweep
# =====================================================================

def deborah_sweep(fem, eta_values, mu, h_plate, E, nu, w0, omega,
                  rho_w, g, dt=30.0, n_cycles=6):
    """Run Newtonian viscoelastic solver across a range of viscosities.

    Parameters
    ----------
    fem : HermiteFEMBeam
    eta_values : array_like
        Viscosities to sweep (Pa s).
    mu : float
        Shear modulus (Pa).
    h_plate : float
    E, nu : float
    w0, omega : float
    rho_w, g : float
    dt, n_cycles : float, int

    Returns
    -------
    results : list of dicts
        Each dict has keys 'eta', 'De', 't_relax', 'w_hist', 't_hist', etc.
    """
    D_el = E * h_plate**3 / (12 * (1 - nu**2))
    T_tide = 2 * np.pi / omega
    results = []

    for eta in eta_values:
        D_v = eta * h_plate**3 / 3
        t_M = eta / mu
        De = t_M / T_tide

        res = run_viscoelastic_tidal(
            fem, D_v, D_el, w0, omega, rho_w, g,
            dt=dt, n_cycles=n_cycles
        )
        res['eta'] = eta
        res['De'] = De
        res['t_relax'] = D_v / D_el
        res['D_v'] = D_v
        results.append(res)

    return results


# =====================================================================
# Self-tests
# =====================================================================

def _run_tests():
    """Quick self-tests: verify limits and analytical solution."""
    rho_w = 1025.0
    g_val = 9.8
    E_ice = 9.33e9
    nu_ice = 0.325
    mu_ice = 3.52e9
    h = 500.0
    w0 = 1.0
    T_tide = 12 * 3600.0
    omega = 2 * np.pi / T_tide

    D_el = E_ice * h**3 / (12 * (1 - nu_ice**2))
    ell_f = (4 * D_el / (rho_w * g_val))**0.25

    L = 10 * ell_f

    # --- Test with uniform mesh ---
    Nel = 80
    fem = HermiteFEMBeam(L, Nel)
    assert fem.he.shape == (Nel,), f'he shape: {fem.he.shape}'
    assert fem.K_ref.shape == (Nel, 2, 4, 4), f'K_ref shape: {fem.K_ref.shape}'
    print(f'Uniform mesh: {Nel} elements, he = {fem.he[0]:.1f} m')

    # --- Test 1: Pure viscous limit ---
    A_3 = 3.5e-25
    tau_ref = 0.1e6
    eta = 1 / (2 * A_3 * tau_ref**2)
    D_v = eta * h**3 / 3

    print(f'Test 1: pure viscous (eta = {eta:.2e}, D_v = {D_v:.2e})')
    res_visc = run_viscoelastic_tidal(
        fem, D_v, np.inf, w0, omega, rho_w, g_val, dt=30.0, n_cycles=5
    )
    phase = omega * res_visc['t_hist']
    idx_peak = np.argmax(np.sin(phase))
    t_peak = res_visc['t_hist'][idx_peak]

    beta = (rho_w * g_val / (omega * D_v))**0.25
    lam1 = beta * np.exp(1j * 5 * np.pi / 8)
    lam2 = beta * np.exp(1j * 9 * np.pi / 8)
    C1 = w0 * lam2 / (lam1 - lam2)
    C2 = -w0 * lam1 / (lam1 - lam2)
    x = fem.x_nodes
    W_visc = w0 + C1 * np.exp(lam1 * x) + C2 * np.exp(lam2 * x)
    w_ana = np.imag(W_visc * np.exp(1j * omega * t_peak))

    err_visc = np.max(np.abs(res_visc['w_hist'][idx_peak] - w_ana))
    print(f'  max error = {err_visc:.4f} m (should be < 0.05)')
    assert err_visc < 0.05, f'Pure viscous limit failed: {err_visc}'

    # --- Test 2: Viscoelastic analytical vs FEM ---
    eta_ve = 1e14
    D_v_ve = eta_ve * h**3 / 3
    t_M = eta_ve / mu_ice
    De = t_M / T_tide
    print(f'Test 2: viscoelastic (eta = {eta_ve:.1e}, De = {De:.2f})')

    res_ve = run_viscoelastic_tidal(
        fem, D_v_ve, D_el, w0, omega, rho_w, g_val, dt=30.0, n_cycles=6
    )
    phase = omega * res_ve['t_hist']
    idx_peak = np.argmax(np.sin(phase))
    t_peak = res_ve['t_hist'][idx_peak]

    w_ana_ve = viscoelastic_analytical(
        fem.x_nodes, t_peak, w0, D_v_ve, D_el, omega, rho_w, g_val
    )
    err_ve = np.max(np.abs(res_ve['w_hist'][idx_peak] - w_ana_ve))
    print(f'  max error = {err_ve:.4f} m (should be < 0.05)')
    assert err_ve < 0.05, f'VE analytical test failed: {err_ve}'

    # --- Test 3: High-De elastic limit ---
    eta_high = 1e17
    D_v_high = eta_high * h**3 / 3
    De_high = (eta_high / mu_ice) / T_tide
    print(f'Test 3: high De (eta = {eta_high:.1e}, De = {De_high:.1f})')

    res_hi = run_viscoelastic_tidal(
        fem, D_v_high, D_el, w0, omega, rho_w, g_val, dt=30.0, n_cycles=8
    )
    phase = omega * res_hi['t_hist']
    idx_peak = np.argmax(np.sin(phase))
    t_peak = res_hi['t_hist'][idx_peak]

    xi = fem.x_nodes / ell_f
    w_elastic = w0 * np.sin(omega * t_peak) * (
        1 - np.exp(-xi) * (np.cos(xi) + np.sin(xi))
    )
    err_el = np.max(np.abs(res_hi['w_hist'][idx_peak] - w_elastic))
    print(f'  max error = {err_el:.4f} m (should be < 0.15)')
    assert err_el < 0.15, f'Elastic limit test failed: {err_el}'

    # --- Test 4: Graded mesh ---
    x_graded = generate_graded_mesh(L, ell_f, h_min=25.0, h_max=250.0,
                                     Nel_max=100)
    fem_g = HermiteFEMBeam(x_nodes=x_graded)
    print(f'Test 4: graded mesh ({fem_g.Nel} elements, '
          f'he = {fem_g.he[0]:.1f}--{fem_g.he[-1]:.1f} m)')

    res_g = run_viscoelastic_tidal(
        fem_g, D_v_ve, D_el, w0, omega, rho_w, g_val, dt=30.0, n_cycles=6
    )
    phase = omega * res_g['t_hist']
    idx_peak = np.argmax(np.sin(phase))
    t_peak = res_g['t_hist'][idx_peak]

    w_ana_g = viscoelastic_analytical(
        fem_g.x_nodes, t_peak, w0, D_v_ve, D_el, omega, rho_w, g_val
    )
    err_g = np.max(np.abs(res_g['w_hist'][idx_peak] - w_ana_g))
    print(f'  max error = {err_g:.4f} m (should be < 0.05)')
    assert err_g < 0.05, f'Graded mesh test failed: {err_g}'

    print('\nAll tests passed.')


if __name__ == '__main__':
    _run_tests()
