#!/usr/bin/env python3
"""
Extract normalized ice thickness profiles along flowlines in East Antarctica.

Strategy: seed flowlines in fast-flowing ice (high SNR velocity) and trace
UPSTREAM to the divide.  Keep only flowlines that connect to a ridge
(speed ~ 0, high surface elevation).  Normalize thickness by the divide
value h_0, matching the Vialov profile definition.

Uses:
  - BedMachine v4.1 for ice thickness, surface elevation, and ice mask
  - MEaSUREs Antarctic Boundaries for East Antarctica grounded ice mask
  - MEaSUREs velocity (VX, VY) for flowline tracing
    (following the approach in Quantarctica3/BalanceVelocity)

Outputs:
  - east_antarctica_profiles.h5 with normalized thickness profiles
"""

import numpy as np
import netCDF4
import shapefile
from shapely.geometry import shape as shapely_shape, Point
from shapely.prepared import prep
from scipy.ndimage import uniform_filter, binary_dilation, binary_erosion
from scipy.interpolate import RegularGridInterpolator
from shapely.ops import unary_union
import h5py
import time
import os

# ============================================================
# Paths
# ============================================================
quantarctica = '/Users/minchew/Documents/Research/Antarctica/maps/Quantarctica3'
bedmachine_path = os.path.join(
    quantarctica,
    'TerrainModels/BEDMACHINE/'
    'NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc')
boundaries_path = os.path.join(
    quantarctica,
    'Glaciology/MEaSUREs Antarctic Boundaries/IceBoundaries_Antarctica_v2.shp')
velocity_path = os.path.join(
    quantarctica,
    'BalanceVelocity/antarctica_ice_velocity_450m_v2.nc')
output_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(output_dir, 'east_antarctica_profiles.h5')

# ============================================================
# Parameters
# ============================================================
subsample = 10          # subsample factor (500m -> 5km)
seed_spacing = 10       # seed spacing in grid cells (~50 km at 5km)
seed_speed_min = 50.0   # m/yr, minimum speed for downstream seeds
min_profile_len = 60    # minimum steps (~300 km at 5km)
max_steps = 1200        # maximum flowline steps (~6000 km)
step_frac = 0.8         # fraction of grid cell per step

# Ridge detection: a flowline "connects to the ridge" if it reaches
# a point where speed < ridge_speed_max AND surface > ridge_elev_min
ridge_speed_max = 5.0   # m/yr
ridge_elev_min = 2500   # m

# ============================================================
# 1. Load BedMachine (subsampled)
# ============================================================
print("Loading BedMachine...")
t0 = time.time()
with netCDF4.Dataset(bedmachine_path) as ds:
    x_full = ds.variables['x'][:].astype(np.float64)
    y_full = ds.variables['y'][:].astype(np.float64)
    s = subsample
    thickness = ds.variables['thickness'][::s, ::s].astype(np.float32)
    surface = ds.variables['surface'][::s, ::s].astype(np.float32)
    mask_data = ds.variables['mask'][::s, ::s].astype(np.int8)

x = x_full[::s]
y = y_full[::s]
dx = abs(x[1] - x[0])
ny, nx = thickness.shape
print(f"  Grid: {nx} x {ny}, dx = {dx:.0f} m, loaded in {time.time()-t0:.1f}s")

thickness[thickness < 0] = 0
surface[surface < -9000] = 0
grounded = (mask_data == 2) | (mask_data == 4)
any_ice = (mask_data >= 2) & (mask_data <= 4)  # grounded + floating + Vostok

# ============================================================
# 2. Load velocity (subsampled to same grid)
# ============================================================
print("Loading MEaSUREs velocity...")
t0 = time.time()
with netCDF4.Dataset(velocity_path) as ds:
    x_vel = ds.variables['x'][:].astype(np.float64)
    y_vel = ds.variables['y'][:].astype(np.float64)
    dx_vel = abs(x_vel[1] - x_vel[0])
    vel_sub = max(1, int(round(dx / dx_vel)))
    vx_raw = ds.variables['VX'][::vel_sub, ::vel_sub].astype(np.float32)
    vy_raw = ds.variables['VY'][::vel_sub, ::vel_sub].astype(np.float32)
    x_vel_s = x_vel[::vel_sub]
    y_vel_s = y_vel[::vel_sub]

vx_raw = np.where(np.isfinite(vx_raw), vx_raw, 0.0)
vy_raw = np.where(np.isfinite(vy_raw), vy_raw, 0.0)

# Interpolate to BedMachine grid
y_vel_flip = y_vel_s[::-1]
vx_interp = RegularGridInterpolator(
    (y_vel_flip, x_vel_s), vx_raw[::-1, :],
    method='linear', bounds_error=False, fill_value=0.0)
vy_interp = RegularGridInterpolator(
    (y_vel_flip, x_vel_s), vy_raw[::-1, :],
    method='linear', bounds_error=False, fill_value=0.0)

yy, xx = np.meshgrid(y, x, indexing='ij')
pts = np.column_stack([yy.ravel(), xx.ravel()])
vx_grid = vx_interp(pts).reshape(ny, nx)
vy_grid = vy_interp(pts).reshape(ny, nx)
del vx_raw, vy_raw, pts

# Light smooth
vx_grid = uniform_filter(vx_grid, size=3)
vy_grid = uniform_filter(vy_grid, size=3)
vx_grid[~grounded] = 0
vy_grid[~grounded] = 0
speed = np.sqrt(vx_grid**2 + vy_grid**2)
print(f"  Velocity loaded in {time.time()-t0:.1f}s")
print(f"  Speed range on grounded ice: "
      f"{speed[grounded].min():.1f} – {speed[grounded].max():.0f} m/yr")

# ============================================================
# 3. Create East Antarctica mask from MEaSUREs boundaries
# ============================================================
print("Building East Antarctica mask...")
t0 = time.time()

east_polys = []
sf = shapefile.Reader(boundaries_path)
fields = [f[0] for f in sf.fields[1:]]
for sr in sf.shapeRecords():
    rec = dict(zip(fields, sr.record))
    if rec['Regions'] == 'East' and rec['TYPE'] == 'GR':
        east_polys.append(shapely_shape(sr.shape.__geo_interface__))

east_union = unary_union(east_polys)
east_prepared = prep(east_union)

east_mask = np.zeros((ny, nx), dtype=bool)
coarse = 5
for i in range(0, ny, coarse):
    for j in range(0, nx, coarse):
        if east_prepared.contains(Point(x[j], y[i])):
            east_mask[i:i+coarse, j:j+coarse] = True

print("  Refining mask at boundaries...")
border = binary_dilation(east_mask, iterations=coarse) & ~binary_erosion(
    east_mask, iterations=coarse)
for i, j in np.argwhere(border):
    east_mask[i, j] = east_prepared.contains(Point(x[j], y[i]))

east_grounded = east_mask & grounded & (thickness > 50)
print(f"  East Antarctica grounded cells: {east_grounded.sum():,} "
      f"({time.time()-t0:.1f}s)")

# ============================================================
# 4. Seed flowlines in fast ice, trace UPSTREAM to the divide
# ============================================================
print("Identifying downstream seeds (fast ice)...")

seed_candidates = (
    east_grounded &
    (speed > seed_speed_min) &
    (thickness > 100)
)

# Subsample to seed_spacing
seed_rows, seed_cols = [], []
for i in range(0, ny, seed_spacing):
    for j in range(0, nx, seed_spacing):
        i0, i1 = max(0, i - 2), min(ny, i + 3)
        j0, j1 = max(0, j - 2), min(nx, j + 3)
        window = seed_candidates[i0:i1, j0:j1]
        if window.any():
            # Pick the fastest cell in the window
            sub_spd = speed[i0:i1, j0:j1].copy()
            sub_spd[~window] = -1
            local = np.unravel_index(sub_spd.argmax(), sub_spd.shape)
            seed_rows.append(i0 + local[0])
            seed_cols.append(j0 + local[1])

seed_rows = np.array(seed_rows)
seed_cols = np.array(seed_cols)
print(f"  {len(seed_rows)} downstream seeds identified")

# ============================================================
# 5. Trace flowlines UPSTREAM (against velocity)
# ============================================================
print("Tracing flowlines upstream...")
t0 = time.time()

# Build interpolators
y_flip = y[::-1]
thick_interp = RegularGridInterpolator(
    (y_flip, x), thickness[::-1, :],
    method='linear', bounds_error=False, fill_value=0.0)
surf_interp = RegularGridInterpolator(
    (y_flip, x), surface[::-1, :],
    method='linear', bounds_error=False, fill_value=0.0)
mask_grounded_interp = RegularGridInterpolator(
    (y_flip, x), grounded[::-1, :].astype(np.float32),
    method='nearest', bounds_error=False, fill_value=0.0)
mask_ice_interp = RegularGridInterpolator(
    (y_flip, x), any_ice[::-1, :].astype(np.float32),
    method='nearest', bounds_error=False, fill_value=0.0)
vx_bm_interp = RegularGridInterpolator(
    (y_flip, x), vx_grid[::-1, :],
    method='linear', bounds_error=False, fill_value=0.0)
vy_bm_interp = RegularGridInterpolator(
    (y_flip, x), vy_grid[::-1, :],
    method='linear', bounds_error=False, fill_value=0.0)

step_size = dx * step_frac

profiles = []
n_reached_ridge = 0
margin_h_threshold = 100.0  # m, stop downstream leg when ice thins to this

for k in range(len(seed_rows)):
    seed_px = float(x[seed_cols[k]])
    seed_py = float(y[seed_rows[k]])

    # --- Phase 1: trace UPSTREAM from seed to the divide ---
    px, py = seed_px, seed_py
    up_x = [px]
    up_y = [py]
    up_h = [float(thick_interp((py, px)))]

    reached_ridge = False

    for step in range(max_steps):
        vx_here = float(vx_bm_interp((py, px)))
        vy_here = float(vy_bm_interp((py, px)))
        spd = np.sqrt(vx_here**2 + vy_here**2)

        if spd < 0.5:
            elev = float(surf_interp((py, px)))
            if elev > ridge_elev_min:
                reached_ridge = True
            break

        # Step against the velocity
        ux = -vx_here / spd
        uy = -vy_here / spd
        px += ux * step_size
        py += uy * step_size

        if float(mask_grounded_interp((py, px))) < 0.5:
            break

        h_here = float(thick_interp((py, px)))
        if h_here < 10:
            break

        up_x.append(px)
        up_y.append(py)
        up_h.append(h_here)

        spd_here = np.sqrt(
            float(vx_bm_interp((py, px)))**2 +
            float(vy_bm_interp((py, px)))**2)
        elev_here = float(surf_interp((py, px)))
        if spd_here < ridge_speed_max and elev_here > ridge_elev_min:
            reached_ridge = True
            break

    if not reached_ridge:
        continue

    # --- Phase 2: trace DOWNSTREAM from seed through floating ice ---
    px, py = seed_px, seed_py
    dn_x, dn_y, dn_h, dn_gr = [], [], [], []

    for step in range(max_steps):
        vx_here = float(vx_bm_interp((py, px)))
        vy_here = float(vy_bm_interp((py, px)))
        spd = np.sqrt(vx_here**2 + vy_here**2)

        if spd < 0.5:
            break

        ux = vx_here / spd
        uy = vy_here / spd
        px += ux * step_size
        py += uy * step_size

        # Stop if we leave all ice (ocean or land)
        if float(mask_ice_interp((py, px))) < 0.5:
            break

        h_here = float(thick_interp((py, px)))
        if h_here < margin_h_threshold:
            dn_x.append(px)
            dn_y.append(py)
            dn_h.append(h_here)
            dn_gr.append(float(mask_grounded_interp((py, px))) > 0.5)
            break

        dn_x.append(px)
        dn_y.append(py)
        dn_h.append(h_here)
        dn_gr.append(float(mask_grounded_interp((py, px))) > 0.5)

    # --- Combine: upstream (reversed) + seed + downstream ---
    # Upstream leg is all grounded (traced on grounded ice only)
    n_up = len(up_x)
    up_grounded = [True] * n_up

    xs_arr = np.array(up_x[::-1] + dn_x)
    ys_arr = np.array(up_y[::-1] + dn_y)
    h_arr = np.array(up_h[::-1] + dn_h)
    gr_arr = np.array(up_grounded[::-1] + dn_gr)

    if len(h_arr) < min_profile_len:
        continue

    n_reached_ridge += 1

    ds_arr = np.sqrt(np.diff(xs_arr)**2 + np.diff(ys_arr)**2)
    dist = np.concatenate([[0], np.cumsum(ds_arr)])

    profiles.append({
        'x': xs_arr,
        'y': ys_arr,
        'distance': dist,
        'thickness': h_arr,
        'grounded': gr_arr,
    })

print(f"  {n_reached_ridge} flowlines reached the ridge")
print(f"  {len(profiles)} profiles kept (>= {min_profile_len} steps)"
      f" in {time.time()-t0:.1f}s")

# ============================================================
# 6. Normalize by divide thickness and save to HDF5
# ============================================================
print(f"Saving {len(profiles)} profiles to {output_path}...")

n_pts = 500
xi_norm = np.linspace(0, 1, n_pts)

with h5py.File(output_path, 'w') as f:
    f.attrs['description'] = (
        'Normalized ice thickness profiles along flowlines in East Antarctica. '
        'Flowlines seeded in fast ice and traced upstream to the divide. '
        'Normalized by divide thickness h_0 (Vialov convention).')
    f.attrs['source_thickness'] = 'BedMachine Antarctica v4.1'
    f.attrs['source_velocity'] = 'MEaSUREs ice velocity v2 (450m)'
    f.attrs['source_boundaries'] = 'MEaSUREs Antarctic Boundaries v2'
    f.attrs['n_profiles'] = len(profiles)
    f.attrs['grid_resolution_m'] = float(dx)
    f.attrs['n_normalized_points'] = n_pts

    f.create_dataset('x_norm', data=xi_norm, dtype='f4', compression='gzip')

    h_norm_all = np.zeros((len(profiles), n_pts), dtype=np.float32)
    lengths = np.zeros(len(profiles), dtype=np.float64)
    divide_thicknesses = np.zeros(len(profiles), dtype=np.float64)
    gl_x_norm = np.zeros(len(profiles), dtype=np.float64)

    for i, prof in enumerate(profiles):
        dist = prof['distance']
        h = prof['thickness']
        gr = prof['grounded']
        L = dist[-1]
        h0 = h[0]  # divide thickness (first point after reversal)

        lengths[i] = L
        divide_thicknesses[i] = h0

        # Find grounding line position: last grounded point
        grounded_idx = np.where(gr)[0]
        if len(grounded_idx) > 0:
            gl_x_norm[i] = dist[grounded_idx[-1]] / L
        else:
            gl_x_norm[i] = 1.0

        # Normalize: x/L from divide (0) to margin (1), h/h_0
        x_frac = dist / L
        h_frac = h / h0

        h_norm_all[i, :] = np.interp(xi_norm, x_frac, h_frac)

    f.create_dataset('h_norm', data=h_norm_all,
                     dtype='f4', compression='gzip')
    f.create_dataset('grounding_line_x_norm', data=gl_x_norm,
                     dtype='f4', compression='gzip')
    f.create_dataset('profile_length_m', data=lengths,
                     dtype='f8', compression='gzip')
    f.create_dataset('divide_thickness_m', data=divide_thicknesses,
                     dtype='f8', compression='gzip')

    raw = f.create_group('raw')
    for i, prof in enumerate(profiles):
        g = raw.create_group(f'profile_{i:04d}')
        g.create_dataset('x_epsg3031', data=prof['x'], dtype='f8')
        g.create_dataset('y_epsg3031', data=prof['y'], dtype='f8')
        g.create_dataset('distance_m', data=prof['distance'], dtype='f8')
        g.create_dataset('thickness_m', data=prof['thickness'], dtype='f4')

print(f"\nDone! Wrote {len(profiles)} profiles.")
print(f"  Profile lengths: {lengths.min()/1e3:.0f} – {lengths.max()/1e3:.0f} km "
      f"(median {np.median(lengths)/1e3:.0f} km)")
print(f"  Divide thicknesses: {divide_thicknesses.min():.0f} – "
      f"{divide_thicknesses.max():.0f} m "
      f"(median {np.median(divide_thicknesses):.0f} m)")
