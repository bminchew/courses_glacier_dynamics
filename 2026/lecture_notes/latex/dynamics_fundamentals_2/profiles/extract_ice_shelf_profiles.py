#!/usr/bin/env python3
"""
Extract normalized ice thickness profiles along flowlines through
the Ross, Ronne-Filchner, and Amery ice shelves.

Strategy: seed flowlines at grounding lines of each ice shelf and trace
DOWNSTREAM through the floating ice to the calving front.  Normalize
distance by the total floating-ice flowline length and thickness by the
grounding-line value h_GL.

Uses:
  - BedMachine v4.1 for ice thickness, surface elevation, and ice mask
  - MEaSUREs Antarctic Boundaries for ice shelf polygons
  - MEaSUREs velocity (VX, VY) for flowline tracing

Outputs:
  - ice_shelf_profiles.h5 with normalized thickness profiles per shelf
"""

import numpy as np
import netCDF4
import shapefile
from shapely.geometry import shape as shapely_shape, Point
from shapely.prepared import prep
from shapely.ops import unary_union
from scipy.ndimage import uniform_filter, binary_dilation, binary_erosion
from scipy.interpolate import RegularGridInterpolator
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
output_path = os.path.join(output_dir, 'ice_shelf_profiles.h5')

# ============================================================
# Ice shelves to extract
# ============================================================
shelf_defs = {
    'Ross': ['Ross_East', 'Ross_West'],
    'Ronne-Filchner': ['Ronne', 'Filchner'],
    'Amery': ['Amery'],
}

# ============================================================
# Parameters
# ============================================================
subsample = 10          # subsample factor (500m -> 5km)
seed_spacing = 8        # seed spacing in grid cells (~40 km at 5km)
seed_speed_min = 20.0   # m/yr, minimum speed at grounding line seeds
min_profile_len = 20    # minimum steps (~100 km at 5km)
max_steps = 1200        # maximum flowline steps
step_frac = 0.8         # fraction of grid cell per step
gl_buffer_cells = 3     # cells around GL to search for seeds

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
floating = (mask_data == 3)
any_ice = (mask_data >= 2) & (mask_data <= 4)

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
speed = np.sqrt(vx_grid**2 + vy_grid**2)
print(f"  Velocity loaded in {time.time()-t0:.1f}s")

# ============================================================
# 3. Build interpolators
# ============================================================
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
mask_floating_interp = RegularGridInterpolator(
    (y_flip, x), floating[::-1, :].astype(np.float32),
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

# ============================================================
# 4. Load shelf polygons from MEaSUREs boundaries
# ============================================================
print("Loading ice shelf boundaries...")
sf = shapefile.Reader(boundaries_path)
fields = [f[0] for f in sf.fields[1:]]

shelf_polygons = {}
for shelf_label, names in shelf_defs.items():
    polys = []
    for sr in sf.shapeRecords():
        rec = dict(zip(fields, sr.record))
        if rec['NAME'] in names and rec['TYPE'] == 'FL':
            polys.append(shapely_shape(sr.shape.__geo_interface__))
    shelf_polygons[shelf_label] = unary_union(polys)
    print(f"  {shelf_label}: {len(polys)} polygon(s)")

# ============================================================
# 5. For each shelf: build mask, find GL seeds, trace downstream
# ============================================================
n_pts = 500
xi_norm = np.linspace(0, 1, n_pts)

all_shelf_profiles = {}

for shelf_label, shelf_poly in shelf_polygons.items():
    print(f"\n{'='*60}")
    print(f"Processing {shelf_label} Ice Shelf...")
    print(f"{'='*60}")
    t0 = time.time()

    # Build rasterized shelf mask
    shelf_prepared = prep(shelf_poly)
    shelf_mask = np.zeros((ny, nx), dtype=bool)
    coarse = 5
    for i in range(0, ny, coarse):
        for j in range(0, nx, coarse):
            if shelf_prepared.contains(Point(x[j], y[i])):
                shelf_mask[i:i+coarse, j:j+coarse] = True

    # Refine at boundaries
    border = binary_dilation(shelf_mask, iterations=coarse) & ~binary_erosion(
        shelf_mask, iterations=coarse)
    for i, j in np.argwhere(border):
        shelf_mask[i, j] = shelf_prepared.contains(Point(x[j], y[i]))

    shelf_floating = shelf_mask & floating & (thickness > 50)
    print(f"  Floating cells in shelf: {shelf_floating.sum():,}")

    # Find grounding line: grounded cells adjacent to shelf floating cells
    gl_zone = binary_dilation(shelf_floating, iterations=gl_buffer_cells) & grounded
    gl_zone = gl_zone & (thickness > 100) & (speed > seed_speed_min)

    # Subsample GL seeds
    seed_rows, seed_cols = [], []
    for i in range(0, ny, seed_spacing):
        for j in range(0, nx, seed_spacing):
            i0, i1 = max(0, i - 2), min(ny, i + 3)
            j0, j1 = max(0, j - 2), min(nx, j + 3)
            window = gl_zone[i0:i1, j0:j1]
            if window.any():
                sub_spd = speed[i0:i1, j0:j1].copy()
                sub_spd[~window] = -1
                local = np.unravel_index(sub_spd.argmax(), sub_spd.shape)
                seed_rows.append(i0 + local[0])
                seed_cols.append(j0 + local[1])

    seed_rows = np.array(seed_rows)
    seed_cols = np.array(seed_cols)
    print(f"  {len(seed_rows)} grounding-line seeds")

    # Trace each seed DOWNSTREAM through the shelf
    profiles = []
    for k in range(len(seed_rows)):
        px = float(x[seed_cols[k]])
        py = float(y[seed_rows[k]])

        fl_x = [px]
        fl_y = [py]
        fl_h = [float(thick_interp((py, px)))]

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

            # Stop if we leave all ice
            if float(mask_ice_interp((py, px))) < 0.5:
                break

            h_here = float(thick_interp((py, px)))
            if h_here < 10:
                break

            fl_x.append(px)
            fl_y.append(py)
            fl_h.append(h_here)

        if len(fl_x) < min_profile_len:
            continue

        xs_arr = np.array(fl_x)
        ys_arr = np.array(fl_y)
        h_arr = np.array(fl_h)

        ds_arr = np.sqrt(np.diff(xs_arr)**2 + np.diff(ys_arr)**2)
        dist = np.concatenate([[0], np.cumsum(ds_arr)])

        # Find where the flowline transitions from grounded to floating
        is_gr = np.array([
            float(mask_grounded_interp((ys_arr[m], xs_arr[m]))) > 0.5
            for m in range(len(xs_arr))
        ])
        grounded_idx = np.where(is_gr)[0]
        if len(grounded_idx) > 0:
            gl_idx = grounded_idx[-1]
        else:
            gl_idx = 0

        profiles.append({
            'x': xs_arr,
            'y': ys_arr,
            'distance': dist,
            'thickness': h_arr,
            'grounded': is_gr,
            'gl_idx': gl_idx,
        })

    print(f"  {len(profiles)} profiles kept (>= {min_profile_len} steps)"
          f" in {time.time()-t0:.1f}s")
    all_shelf_profiles[shelf_label] = profiles

# ============================================================
# 6. Normalize and save to HDF5
# ============================================================
print(f"\nSaving to {output_path}...")

with h5py.File(output_path, 'w') as f:
    f.attrs['description'] = (
        'Normalized ice thickness profiles along flowlines through '
        'the Ross, Ronne-Filchner, and Amery ice shelves. '
        'Flowlines seeded at grounding lines and traced downstream. '
        'Distance normalized by total flowline length; thickness '
        'normalized by grounding-line thickness h_GL.')
    f.attrs['source_thickness'] = 'BedMachine Antarctica v4.1'
    f.attrs['source_velocity'] = 'MEaSUREs ice velocity v2 (450m)'
    f.attrs['source_boundaries'] = 'MEaSUREs Antarctic Boundaries v2'
    f.attrs['n_normalized_points'] = n_pts

    f.create_dataset('x_norm', data=xi_norm, dtype='f4', compression='gzip')

    for shelf_label, profiles in all_shelf_profiles.items():
        print(f"  {shelf_label}: {len(profiles)} profiles")
        grp = f.create_group(shelf_label)
        grp.attrs['n_profiles'] = len(profiles)

        if len(profiles) == 0:
            continue

        h_norm_all = np.zeros((len(profiles), n_pts), dtype=np.float32)
        lengths = np.zeros(len(profiles), dtype=np.float64)
        gl_thicknesses = np.zeros(len(profiles), dtype=np.float64)
        gl_x_norm = np.zeros(len(profiles), dtype=np.float64)

        for i, prof in enumerate(profiles):
            dist = prof['distance']
            h = prof['thickness']
            gl_idx = prof['gl_idx']
            L = dist[-1]
            h_gl = h[gl_idx]

            lengths[i] = L
            gl_thicknesses[i] = h_gl
            gl_x_norm[i] = dist[gl_idx] / L if L > 0 else 0.0

            # Normalize: x/L, h/h_GL
            x_frac = dist / L
            h_frac = h / h_gl if h_gl > 0 else h * 0.0

            h_norm_all[i, :] = np.interp(xi_norm, x_frac, h_frac)

        grp.create_dataset('h_norm', data=h_norm_all,
                           dtype='f4', compression='gzip')
        grp.create_dataset('grounding_line_x_norm', data=gl_x_norm,
                           dtype='f4', compression='gzip')
        grp.create_dataset('profile_length_m', data=lengths,
                           dtype='f8', compression='gzip')
        grp.create_dataset('gl_thickness_m', data=gl_thicknesses,
                           dtype='f8', compression='gzip')

        raw = grp.create_group('raw')
        for i, prof in enumerate(profiles):
            g = raw.create_group(f'profile_{i:04d}')
            g.create_dataset('x_epsg3031', data=prof['x'], dtype='f8')
            g.create_dataset('y_epsg3031', data=prof['y'], dtype='f8')
            g.create_dataset('distance_m', data=prof['distance'], dtype='f8')
            g.create_dataset('thickness_m', data=prof['thickness'], dtype='f4')
            g.create_dataset('grounded', data=prof['grounded'], dtype='bool')

print("\nDone!")
for shelf_label, profiles in all_shelf_profiles.items():
    if len(profiles) > 0:
        lengths = np.array([p['distance'][-1] for p in profiles])
        gl_h = np.array([p['thickness'][p['gl_idx']] for p in profiles])
        print(f"  {shelf_label}: {len(profiles)} profiles, "
              f"lengths {lengths.min()/1e3:.0f}–{lengths.max()/1e3:.0f} km, "
              f"GL thickness {gl_h.min():.0f}–{gl_h.max():.0f} m")
