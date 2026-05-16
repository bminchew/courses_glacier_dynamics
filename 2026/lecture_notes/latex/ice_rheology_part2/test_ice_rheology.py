"""Unit tests for ice_rheology_part2 figure-generation code.

Tests physical constants, analytical solutions, dimensional consistency,
boundary conditions, limiting cases, and the FEM solver against known results.
"""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

# =====================================================================
# Physical constants (duplicated from notebook for isolation)
# =====================================================================
E = 9.33e9       # Pa, Young's modulus
mu = 3.52e9      # Pa, shear modulus
K = 8.90e9       # Pa, bulk modulus
nu = 0.325       # Poisson's ratio
lam = 6.55e9     # Pa, first Lame parameter
rho_i = 917.0    # kg/m^3
rho_w = 1025.0   # kg/m^3
g = 9.8          # m/s^2
R_gas = 8.314    # J/(mol K)


# =====================================================================
# Functions under test (extracted from notebook)
# =====================================================================
def A_at_T(T):
    A_ref = 3.5e-25
    T_ref = 263.0
    Q = np.where(np.asarray(T) < 263, 90e3, 139e3)
    return A_ref * np.exp(-Q / R_gas * (1.0 / np.asarray(T) - 1.0 / T_ref))


def maxwell_time(T, tau_e):
    A = A_at_T(T)
    return 1.0 / (2 * mu * A * tau_e**2)


def flexural_wavelength(h):
    D = E * h**3 / (12 * (1 - nu**2))
    return (4 * D / (rho_w * g))**0.25


def compute_D_n(A_Glen, n_nl, h_plate):
    """Plate flexural rigidity parameter D_n."""
    return (4 * n_nl / (1 + 2 * n_nl)) * (h_plate / 2)**((1.0 / n_nl) + 2) / A_Glen**(1.0 / n_nl)


def elastic_flexure(x, t, w0, ell_f, omega):
    xi = x / ell_f
    return w0 * np.sin(omega * t) * (1 - np.exp(-xi) * (np.cos(xi) + np.sin(xi)))


def elastic_curvature(x, t, w0, ell_f, omega):
    xi = x / ell_f
    return w0 * np.sin(omega * t) * 2 * np.exp(-xi) * (np.cos(xi) - np.sin(xi)) / ell_f**2


def newtonian_viscous_solution(x, t, w0, D_v, omega, rho_w, g):
    beta = (rho_w * g / (omega * D_v))**0.25
    lam1 = beta * np.exp(1j * 5 * np.pi / 8)
    lam2 = beta * np.exp(1j * 9 * np.pi / 8)
    C1 = w0 * lam2 / (lam1 - lam2)
    C2 = -w0 * lam1 / (lam1 - lam2)
    eta_hat = C1 * np.exp(lam1 * x) + C2 * np.exp(lam2 * x)
    return np.imag((w0 + eta_hat) * np.exp(1j * omega * t))


def eta_profile(z_norm, T_s, T_b=273.0):
    T_bar = (T_s + T_b) / 2
    T_z = T_bar + (T_s - T_b) / 2 * z_norm
    return A_at_T(T_bar) / A_at_T(T_z)


def compute_z0_newtonian(T_s, T_b=273.0):
    num, _ = quad(lambda z: eta_profile(z, T_s, T_b) * z, -1, 1)
    den, _ = quad(lambda z: eta_profile(z, T_s, T_b), -1, 1)
    return num / den


def compute_Dv_ratio(T_s, T_b=273.0):
    z0 = compute_z0_newtonian(T_s, T_b)
    Dv_tilde, _ = quad(lambda z: eta_profile(z, T_s, T_b) * (z - z0)**2, -1, 1)
    Dv_uniform = 2.0 / 3.0
    return Dv_tilde / Dv_uniform


# =====================================================================
# Tests: Elastic moduli consistency
# =====================================================================
class TestElasticModuli:
    def test_E_from_mu_nu(self):
        """E = 2 mu (1 + nu)."""
        E_check = 2 * mu * (1 + nu)
        assert abs(E - E_check) / E < 0.005

    def test_K_from_E_nu(self):
        """K = E / (3(1 - 2nu))."""
        K_check = E / (3 * (1 - 2 * nu))
        assert abs(K - K_check) / K < 0.005

    def test_lam_from_K_mu(self):
        """lambda = K - 2mu/3."""
        lam_check = K - 2 * mu / 3
        assert abs(lam - lam_check) / lam < 0.01


# =====================================================================
# Tests: Arrhenius rate factor A(T)
# =====================================================================
class TestArrhenius:
    def test_reference_point(self):
        """A(263 K) should equal A_ref = 3.5e-25."""
        assert A_at_T(263) == pytest.approx(3.5e-25, rel=1e-10)

    def test_monotonic_increasing(self):
        """A(T) should increase with temperature."""
        T_vals = np.array([233, 243, 253, 263, 268, 273])
        A_vals = A_at_T(T_vals)
        assert np.all(np.diff(A_vals) > 0)

    def test_vectorized(self):
        """Should handle array inputs."""
        T_arr = np.array([243, 263, 270])
        A_arr = A_at_T(T_arr)
        assert A_arr.shape == (3,)
        for i, T in enumerate(T_arr):
            assert A_arr[i] == pytest.approx(A_at_T(float(T)), rel=1e-10)

    def test_Q_jump_at_263(self):
        """Activation energy changes at 263 K, producing a kink in log(A)."""
        T_below = 262.9
        T_above = 263.1
        # Slope d(ln A)/d(1/T) = -Q/R should be steeper above 263 K
        dlogA_below = np.log(A_at_T(T_below)) - np.log(A_at_T(T_below - 0.1))
        dlogA_above = np.log(A_at_T(T_above + 0.1)) - np.log(A_at_T(T_above))
        # Q_high = 139 kJ/mol > Q_low = 90 kJ/mol, so the slope should be steeper above
        assert abs(dlogA_above) > abs(dlogA_below)


# =====================================================================
# Tests: Maxwell time
# =====================================================================
class TestMaxwellTime:
    def test_dimensions(self):
        """Maxwell time should have units of seconds (positive finite)."""
        t_M = maxwell_time(263, 0.1e6)
        assert t_M > 0
        assert np.isfinite(t_M)

    def test_decreases_with_temperature(self):
        """Higher T -> larger A -> shorter Maxwell time."""
        t_cold = maxwell_time(243, 0.1e6)
        t_warm = maxwell_time(270, 0.1e6)
        assert t_warm < t_cold

    def test_decreases_with_stress(self):
        """Higher stress -> shorter Maxwell time (for n=3)."""
        t_low = maxwell_time(263, 0.05e6)
        t_high = maxwell_time(263, 0.2e6)
        assert t_high < t_low

    def test_known_value(self):
        """Check against Table 2: A(263K)=3.5e-25, tau_e=0.1 MPa -> t_M ~ 11 hr."""
        t_M = maxwell_time(263, 0.1e6)
        t_M_hr = t_M / 3600
        assert t_M_hr == pytest.approx(11, rel=0.15)


# =====================================================================
# Tests: Flexural wavelength
# =====================================================================
class TestFlexuralWavelength:
    def test_dimensions(self):
        """Output should be in meters, positive."""
        lf = flexural_wavelength(500)
        assert lf > 0
        assert 1e3 < lf < 10e3  # reasonable range for 500 m ice

    def test_increases_with_thickness(self):
        """Thicker ice -> longer flexural wavelength."""
        lf_thin = flexural_wavelength(200)
        lf_thick = flexural_wavelength(1000)
        assert lf_thick > lf_thin

    def test_h_scaling(self):
        """ell_f ~ h^(3/4)."""
        h1, h2 = 200, 800
        ratio = flexural_wavelength(h2) / flexural_wavelength(h1)
        expected = (h2 / h1)**0.75
        assert ratio == pytest.approx(expected, rel=1e-10)


# =====================================================================
# Tests: Elastic flexure solution
# =====================================================================
class TestElasticFlexure:
    def setup_method(self):
        self.w0 = 1.0
        self.ell_f = flexural_wavelength(500)
        self.omega = 2 * np.pi / (12 * 3600)
        self.t_peak = np.pi / (2 * self.omega)  # peak tide

    def test_clamped_bc_w(self):
        """w(0) = 0 (clamped grounding line)."""
        w = elastic_flexure(0, self.t_peak, self.w0, self.ell_f, self.omega)
        assert w == pytest.approx(0, abs=1e-15)

    def test_clamped_bc_slope(self):
        """dw/dx(0) = 0 (clamped, checked via finite difference)."""
        dx = 0.01
        w0_val = elastic_flexure(0, self.t_peak, self.w0, self.ell_f, self.omega)
        w1_val = elastic_flexure(dx, self.t_peak, self.w0, self.ell_f, self.omega)
        slope = (w1_val - w0_val) / dx
        assert slope == pytest.approx(0, abs=1e-4)

    def test_far_field(self):
        """w -> w0 * sin(omega * t) as x -> infinity."""
        x_far = 20 * self.ell_f
        w_far = elastic_flexure(x_far, self.t_peak, self.w0, self.ell_f, self.omega)
        w_expected = self.w0 * np.sin(self.omega * self.t_peak)
        assert w_far == pytest.approx(w_expected, rel=1e-6)

    def test_curvature_is_second_derivative(self):
        """Analytical curvature matches numerical d2w/dx2."""
        x = np.linspace(0.1 * self.ell_f, 5 * self.ell_f, 1000)
        w = elastic_flexure(x, self.t_peak, self.w0, self.ell_f, self.omega)
        d2w_num = np.gradient(np.gradient(w, x), x)
        d2w_ana = elastic_curvature(x, self.t_peak, self.w0, self.ell_f, self.omega)
        # Trim edges where np.gradient is less accurate
        mask = (x > 0.5 * self.ell_f) & (x < 4 * self.ell_f)
        np.testing.assert_allclose(d2w_num[mask], d2w_ana[mask], rtol=0.02)

    def test_satisfies_ode(self):
        """D w'''' + rho_w g w = rho_w g w0 sin(omega t) at peak tide."""
        h = 500
        D = E * h**3 / (12 * (1 - nu**2))
        x = np.linspace(2 * self.ell_f, 8 * self.ell_f, 2000)
        w = elastic_flexure(x, self.t_peak, self.w0, self.ell_f, self.omega)
        # Fourth derivative via finite differences
        dx = x[1] - x[0]
        d4w = np.diff(w, n=4) / dx**4
        x4 = x[2:-2]
        w4 = w[2:-2]
        residual = D * d4w + rho_w * g * w4 - rho_w * g * self.w0 * np.sin(self.omega * self.t_peak)
        # Residual should be small relative to the buoyancy term
        scale = rho_w * g * self.w0
        np.testing.assert_allclose(residual / scale, 0, atol=0.01)


# =====================================================================
# Tests: Newtonian viscous solution
# =====================================================================
class TestNewtonianViscous:
    def setup_method(self):
        self.w0 = 1.0
        self.h = 500
        eta = 1.43e14
        self.D_v = eta * self.h**3 / 3
        self.omega = 2 * np.pi / (12 * 3600)
        self.t_peak = np.pi / (2 * self.omega)

    def test_clamped_bc(self):
        """w(0, t) = 0."""
        w = newtonian_viscous_solution(0, self.t_peak, self.w0, self.D_v, self.omega, rho_w, g)
        assert abs(w) < 1e-10

    def test_far_field_amplitude(self):
        """Far from grounding line, amplitude -> w0."""
        x_far = 100e3
        # Check amplitude over one cycle
        t_arr = np.linspace(0, 2 * np.pi / self.omega, 100)
        w_arr = np.array([newtonian_viscous_solution(x_far, t, self.w0, self.D_v, self.omega, rho_w, g) for t in t_arr])
        assert max(abs(w_arr)) == pytest.approx(self.w0, rel=0.01)

    def test_phase_lag_exists(self):
        """Near grounding line, the viscous response should lag the tide."""
        ell_f = flexural_wavelength(self.h)
        x_probe = ell_f
        # Find time of peak elastic response
        t_elastic_peak = np.pi / (2 * self.omega)
        # Find time of peak viscous response at x_probe
        t_arr = np.linspace(0, 2 * np.pi / self.omega, 1000)
        w_arr = np.array([newtonian_viscous_solution(x_probe, t, self.w0, self.D_v, self.omega, rho_w, g) for t in t_arr])
        t_visc_peak = t_arr[np.argmax(w_arr)]
        # Viscous peak should come after elastic peak (phase lag > 0)
        assert t_visc_peak > t_elastic_peak


# =====================================================================
# Tests: Power-law flexural rigidity D_n
# =====================================================================
class TestComputeDn:
    def test_newtonian_limit(self):
        """For n=1, D_n should equal D_v = eta h^3 / 3 (plate)."""
        h = 500
        A_1 = 3.5e-25  # doesn't matter for the ratio, but let's use a value
        # For n=1: B = 3*eta = 3/(2A), and D_n = (2/3)*B*(h/2)^3 = (2/3)*(3/(2A))*(h/2)^3
        # = (1/A)*(h/2)^3 = (h/2)^3 / A
        # But compute_D_n uses Glen convention: D_n = (4*1/3) * (h/2)^3 / A^1
        # = (4/3) * (h/2)^3 / A
        # For plate: D_v = eta*h^3/3 = h^3/(6A) = (h/2)^3 * 8/(6A) = (4/3)*(h/2)^3/A
        D_n = compute_D_n(A_1, 1, h)
        eta = 1 / (2 * A_1)
        D_v = eta * h**3 / 3
        assert D_n == pytest.approx(D_v, rel=1e-10)

    def test_thickness_scaling(self):
        """D_n ~ h^(1/n + 2)."""
        A = 3.5e-25
        n = 3
        h1, h2 = 300, 600
        ratio = compute_D_n(A, n, h2) / compute_D_n(A, n, h1)
        expected = (h2 / h1)**(1.0 / n + 2)
        assert ratio == pytest.approx(expected, rel=1e-10)

    def test_no_spurious_factor(self):
        """Regression: compute_D_n should NOT have an extra factor of 2."""
        h = 500
        A = 3.5e-25
        D_n = compute_D_n(A, 1, h)
        eta = 1 / (2 * A)
        D_v = eta * h**3 / 3
        # Should be exactly equal, not 2x
        assert D_n / D_v == pytest.approx(1.0, rel=1e-10)


# =====================================================================
# Tests: Depth-dependent temperature
# =====================================================================
class TestDepthDependentTemperature:
    def test_uniform_temperature_z0(self):
        """When T_s = T_b, z_0 = 0 (no shift)."""
        z0 = compute_z0_newtonian(273.0, 273.0)
        assert z0 == pytest.approx(0, abs=1e-10)

    def test_uniform_temperature_Dv(self):
        """When T_s = T_b, D_tilde = D_v(T_bar), ratio = 1."""
        ratio = compute_Dv_ratio(273.0, 273.0)
        assert ratio == pytest.approx(1.0, rel=1e-10)

    def test_cold_surface_shifts_neutral_axis_up(self):
        """Cold surface (T_s < T_b) -> z_0 > 0."""
        z0 = compute_z0_newtonian(243, 273)
        assert z0 > 0

    def test_z0_increases_with_contrast(self):
        """Larger T_s - T_b contrast -> larger z_0."""
        z0_small = compute_z0_newtonian(263, 273)
        z0_large = compute_z0_newtonian(243, 273)
        assert z0_large > z0_small

    def test_z0_bounded(self):
        """z_0 should be between 0 and 1 (in units of h/2)."""
        for T_s in [233, 243, 253, 263]:
            z0 = compute_z0_newtonian(T_s, 273)
            assert 0 < z0 < 1

    def test_Dv_ratio_less_than_one(self):
        """Depth-dependent eta reduces D_v relative to uniform at mean T."""
        for T_s in [243, 253, 263]:
            ratio = compute_Dv_ratio(T_s, 273)
            assert ratio < 1.0

    def test_eta_profile_normalized(self):
        """eta(z=0) / eta(T_bar) = 1 at midplane."""
        for T_s in [243, 253, 263]:
            assert eta_profile(0, T_s, 273) == pytest.approx(1.0, rel=1e-10)

    def test_eta_profile_symmetry(self):
        """For T_s = T_b, eta(z) = 1 everywhere."""
        z = np.linspace(-1, 1, 100)
        eta = eta_profile(z, 273, 273)
        np.testing.assert_allclose(eta, 1.0, atol=1e-10)

    def test_parallel_axis_theorem(self):
        """D_tilde = int eta (z-z0)^2 dz = int eta z^2 dz - z0^2 int eta dz."""
        T_s = 253
        T_b = 273
        z0 = compute_z0_newtonian(T_s, T_b)
        lhs, _ = quad(lambda z: eta_profile(z, T_s, T_b) * (z - z0)**2, -1, 1)
        rhs_a, _ = quad(lambda z: eta_profile(z, T_s, T_b) * z**2, -1, 1)
        rhs_b, _ = quad(lambda z: eta_profile(z, T_s, T_b), -1, 1)
        rhs = rhs_a - z0**2 * rhs_b
        assert lhs == pytest.approx(rhs, rel=1e-10)


# =====================================================================
# Tests: Dimensional consistency
# =====================================================================
class TestDimensionalConsistency:
    def test_critical_tidal_amplitude(self):
        """w_0* = 2 eps_bg ell_f^2 / (h omega) should have units of meters."""
        h = 500
        eps_bg = 1e-10       # 1/s
        ell_f = flexural_wavelength(h)  # m
        omega = 2 * np.pi / (12 * 3600)  # 1/s
        w0_star = 2 * eps_bg * ell_f**2 / (h * omega)
        # Should be a length in meters, positive and finite
        assert w0_star > 0
        assert np.isfinite(w0_star)
        # Sanity: should be on the order of centimeters for these params
        assert 0.001 < w0_star < 1.0

    def test_kappa_bg_units(self):
        """kappa_bg = eps_bg / (h/2) should have units 1/(m*s)."""
        h = 500
        eps_bg = 1e-10  # 1/s
        kappa_bg = eps_bg / (h / 2)  # 1/(m*s)
        # Curvature rate from tidal bending: w0*omega/ell_f^2
        ell_f = flexural_wavelength(h)
        omega = 2 * np.pi / (12 * 3600)
        kappa_tidal = 1.0 * omega / ell_f**2  # 1/(m*s), for w0=1m
        # Both should be positive and comparable in magnitude sense
        assert kappa_bg > 0
        assert kappa_tidal > 0
        # kappa_tidal >> kappa_bg for these parameters (nonlinearity matters)
        assert kappa_tidal > kappa_bg

    def test_elastic_flexure_ode_units(self):
        """D [N m] * w'''' [1/m^3] = [Pa], rho_w g w = [Pa]."""
        h = 500
        D = E * h**3 / (12 * (1 - nu**2))  # N m = Pa m^3
        # D * w'''' has units Pa m^3 * m^{-3} = Pa ✓
        # rho_w * g * w has units kg/m^3 * m/s^2 * m = Pa ✓
        # Just verify D is positive and finite
        assert D > 0
        assert np.isfinite(D)


# =====================================================================
# Tests: FEM solver (Hermite cubics, Newtonian beam)
# =====================================================================
class TestFEMSolver:
    """Test the FEM solver against the analytical Newtonian viscous solution."""

    def setup_method(self):
        self.h = 500
        self.w0 = 1.0
        self.T_tide = 12 * 3600
        self.omega = 2 * np.pi / self.T_tide
        self.ell_f = flexural_wavelength(self.h)
        A_3 = 3.5e-25
        tau_ref = 0.1e6
        eta = 1 / (2 * A_3 * tau_ref**2)
        self.D_v = eta * self.h**3 / 3
        self.eps_bg = 1e-10
        self.kappa_bg = self.eps_bg / (self.h / 2)

        # Build a small FEM mesh
        self.Lx = 10 * self.ell_f
        self.Nel = 40
        self.he = self.Lx / self.Nel
        self.x_nodes = np.linspace(0, self.Lx, self.Nel + 1)
        self._setup_fem()

    def _setup_fem(self):
        Nel = self.Nel
        he = self.he
        gp = np.array([0.5 - np.sqrt(3) / 6, 0.5 + np.sqrt(3) / 6])
        gw = np.array([0.5, 0.5])

        N_gp = np.zeros((2, 4))
        d2N_gp = np.zeros((2, 4))
        for q in range(2):
            xi = gp[q]
            N_gp[q] = [1 - 3 * xi**2 + 2 * xi**3, xi * (1 - xi)**2,
                        3 * xi**2 - 2 * xi**3, xi**2 * (xi - 1)]
            d2N_gp[q] = [-6 + 12 * xi, -4 + 6 * xi, 6 - 12 * xi, -2 + 6 * xi]

        scale = np.array([1.0, he, 1.0, he])
        self.N_phys = N_gp * scale[None, :]
        self.B_phys = d2N_gp / he**2 * scale[None, :]

        self.K_ref = np.array([gw[q] * he * np.outer(self.B_phys[q], self.B_phys[q]) for q in range(2)])
        self.F_ref = np.array([gw[q] * he * self.N_phys[q] for q in range(2)])
        self.M_ref = sum(gw[q] * he * np.outer(self.N_phys[q], self.N_phys[q]) for q in range(2))

        self.ndof = 2 * (Nel + 1)
        self._gp = gp

        rows = np.zeros((Nel, 4, 4), dtype=int)
        cols = np.zeros((Nel, 4, 4), dtype=int)
        self.rhs_rows = np.zeros((Nel, 4), dtype=int)
        for e in range(Nel):
            d = np.array([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3])
            rows[e] = d[:, None] * np.ones(4, dtype=int)[None, :]
            cols[e] = np.ones(4, dtype=int)[:, None] * d[None, :]
            self.rhs_rows[e] = d
        self.rows_flat = rows.ravel()
        self.cols_flat = cols.ravel()

    def _solve_wdot(self, w_curr, w_tide, D_n, n_nl, dt):
        Nel = self.Nel
        ndof = self.ndof

        D_eff = np.full((Nel, 2), D_n)
        K_elem = (D_eff[:, 0, None, None] * self.K_ref[0]
                  + D_eff[:, 1, None, None] * self.K_ref[1]
                  + dt * rho_w * g * self.M_ref)
        K_global = csr_matrix((K_elem.ravel(), (self.rows_flat, self.cols_flat)),
                              shape=(ndof, ndof)).tolil()

        w_gp = (w_curr[:-1, None] * (1 - self._gp[None, :])
                + w_curr[1:, None] * self._gp[None, :])
        load = rho_w * g * (w_tide - w_gp)
        F_elem = load[:, 0, None] * self.F_ref[0] + load[:, 1, None] * self.F_ref[1]
        F_global = np.zeros(ndof)
        np.add.at(F_global, self.rhs_rows.ravel(), F_elem.ravel())

        for bc in [0, 1]:
            K_global[bc, :] = 0
            K_global[:, bc] = 0
            K_global[bc, bc] = 1.0
            F_global[bc] = 0.0

        return spsolve(K_global.tocsr(), F_global)

    def test_newtonian_steady_state(self):
        """FEM Newtonian solution should match analytical viscous solution after spin-up."""
        dt = 60.0
        n_cycles = 4
        n_steps = int(n_cycles * self.T_tide / dt)
        w = np.zeros(self.Nel + 1)

        for step in range(n_steps):
            t = step * dt
            w_tide = self.w0 * np.sin(self.omega * t)
            wdot = self._solve_wdot(w, w_tide, self.D_v, 1, dt)
            w = w + dt * wdot[0::2]
            w[0] = 0

        # Compare at peak tide in last cycle
        t_final = n_steps * dt
        t_peak_candidates = np.arange(t_final - self.T_tide, t_final, dt)
        # Run one more cycle and find peak
        w_peak = None
        max_far = -1
        for step in range(int(self.T_tide / dt)):
            t = t_final + step * dt
            w_tide = self.w0 * np.sin(self.omega * t)
            wdot = self._solve_wdot(w, w_tide, self.D_v, 1, dt)
            w = w + dt * wdot[0::2]
            w[0] = 0
            if w[-1] > max_far:
                max_far = w[-1]
                w_peak = w.copy()
                t_at_peak = t

        # Analytical solution at the same time
        w_ana = newtonian_viscous_solution(
            self.x_nodes, t_at_peak, self.w0, self.D_v, self.omega, rho_w, g
        )

        # Agreement should be reasonable (coarse mesh, finite dt)
        np.testing.assert_allclose(w_peak, w_ana, atol=0.05)

    def test_clamped_boundary(self):
        """FEM should enforce w(0) = 0 and w'(0) = 0."""
        w = np.zeros(self.Nel + 1)
        dt = 60.0
        for step in range(100):
            t = step * dt
            w_tide = self.w0 * np.sin(self.omega * t)
            wdot = self._solve_wdot(w, w_tide, self.D_v, 1, dt)
            w = w + dt * wdot[0::2]
            w[0] = 0
        assert w[0] == 0
        # w'(0) is DOF index 1 in the wdot vector
        assert wdot[1] == pytest.approx(0, abs=1e-20)


# =====================================================================
# Tests: Elastic-viscous analogy
# =====================================================================
class TestElasticViscousAnalogy:
    def test_viscous_poisson_ratio(self):
        """Under the analogy E/(1-nu^2) -> 4*eta, D_v = 4*eta*h^3/12 = eta*h^3/3.

        The nu terms cancel from the plate equation (as shown in the notes),
        so D_v does NOT contain (1-nu_v^2) in the denominator.
        """
        h = 500
        eta = 1e14
        D_v = eta * h**3 / 3
        # The analogy replaces the entire prefactor E/(1-nu^2) with 4*eta
        D_from_analogy = 4 * eta * h**3 / 12
        assert D_v == pytest.approx(D_from_analogy, rel=1e-10)

    def test_D_v_formula(self):
        """D_v = eta h^3 / 3."""
        h = 500
        eta = 1.43e14
        D_v = eta * h**3 / 3
        # Same thing from the analogy
        D_analogy = 4 * eta * h**3 / 12
        assert D_v == pytest.approx(D_analogy, rel=1e-10)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
