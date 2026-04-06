#!/usr/bin/env python3
"""
Calculate balance velocity for Antarctica using Quantarctica3 data.

Forward flowline integration for Antarctic balance velocity.

Seeds a flowline from every grounded ice cell, traces it downstream
along the observed MEaSUREs velocity field, and deposits that cell's
SMB flux at every grid cell the flowline passes through.  Convergence
is handled naturally — where N flowlines merge into an ice stream,
N× the flux accumulates.  No discrete neighbour routing, no
channelisation, no post-smoothing required.

Algorithm:
    1. Resample MEaSUREs VX, VY (450 m NetCDF) to 1 km grid
    2. For each ice cell, trace a flowline downstream using bilinear
       interpolation of the velocity unit vector
    3. Deposit that cell's ȧ×dx² at every grid cell the flowline
       passes through (visit-stamped to avoid double-counting)
    4. Balance velocity = total accumulated flux / (H × dx)
    5. Fill gaps (slow/stagnant cells) with observed speed

Inputs:
    - BEDMAP2 surface elevation, ice thickness, grounded ice mask (1 km)
    - RACMO surface mass balance (resampled from 35 km)
    - MEaSUREs VX, VY (median of annual velocities, 2017–2022)
    - MEaSUREs ice flow speed (for comparison)

Outputs:
    - BalanceVelocity.tif
"""

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject
from scipy.ndimage import median_filter
import os
import time

# ============================================================
# Paths
# ============================================================
base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
surface_path = os.path.join(base, 'TerrainModels/BEDMAP2/bedmap2_surface.tif')
thickness_path = os.path.join(base, 'TerrainModels/BEDMAP2/bedmap2_thickness.tif')
smb_path = os.path.join(base,
    'Atmosphere/Van Wessem RACMO/RACMO_SurfaceMassBalance_35km.tif')
mask_path = os.path.join(base,
    'TerrainModels/BEDMAP2/bedmap2_icemask_grounded_and_shelves.tif')
velocity_path = os.path.join(base,
    'Glaciology/MEaSUREs Ice Flow Velocity/MEaSUREs_IceFlowSpeed_450m.tif')
vel_nc_path = os.path.join(base,
    'Glaciology/MEaSUREs_AnnualVelocity/MedianVelocity_2017_2022.nc')
out_dir = os.path.join(base, 'BalanceVelocity')
output_tif = os.path.join(out_dir, 'BalanceVelocity.tif')

# ============================================================
# Load data
# ============================================================
print("Loading surface elevation...")
with rasterio.open(surface_path) as src:
    surface = src.read(1).astype(np.float64)
    profile = src.profile.copy()
    bm2_transform = src.transform
    bm2_crs = src.crs
    bm2_shape = src.shape
nrows, ncols = bm2_shape
dx = abs(bm2_transform.a)

print("Loading ice thickness...")
with rasterio.open(thickness_path) as src:
    thickness = src.read(1).astype(np.float64)
    thickness[thickness >= 32767] = 0
    thickness[thickness < 0] = 0

print("Loading ice mask...")
with rasterio.open(mask_path) as src_mask:
    mask_data = np.empty(bm2_shape, dtype=np.float32)
    reproject(rasterio.band(src_mask, 1), mask_data,
              src_transform=src_mask.transform, src_crs=src_mask.crs,
              dst_transform=bm2_transform, dst_crs=bm2_crs,
              resampling=Resampling.nearest)
# 0 = grounded ice, 1 = ice shelf, 127 = ocean/nodata
valid = ((mask_data == 0) | (mask_data == 1)) & (thickness > 10)
grounded = (mask_data == 0) & (thickness > 10)

print("Loading SMB...")
with rasterio.open(smb_path) as src_smb:
    smb = np.empty(bm2_shape, dtype=np.float64)
    reproject(rasterio.band(src_smb, 1), smb,
              src_transform=src_smb.transform, src_crs=src_smb.crs,
              dst_transform=bm2_transform, dst_crs=bm2_crs,
              resampling=Resampling.bilinear)
    smb[smb < -1e30] = 0
smb_myr = smb / 917.0

print(f"  Valid grounded ice cells: {valid.sum():,}")
print(f"  SMB mean: {smb_myr[valid].mean():.3f} m/yr")

# ============================================================
# Step 1: Load median VX, VY (2017–2022, 1 km) → regrid to BEDMAP2
# ============================================================
print("Loading median VX, VY (2017–2022) and regridding to BEDMAP2...")
t0 = time.time()

import netCDF4
from rasterio.transform import from_origin
with netCDF4.Dataset(vel_nc_path) as ds:
    vx_src = np.array(ds.variables['VX'][:], dtype=np.float32)
    vy_src = np.array(ds.variables['VY'][:], dtype=np.float32)
    x_vel = ds.variables['x'][:].astype(np.float64)
    y_vel = ds.variables['y'][:].astype(np.float64)

dx_vel = x_vel[1] - x_vel[0]
dy_vel = y_vel[1] - y_vel[0]
vel_transform = from_origin(x_vel[0] - dx_vel / 2,
                            y_vel[0] - dy_vel / 2,
                            abs(dx_vel), abs(dy_vel))

vx_src = np.where(np.isfinite(vx_src), vx_src, 0.0)
vy_src = np.where(np.isfinite(vy_src), vy_src, 0.0)

vx_grid = np.empty(bm2_shape, dtype=np.float64)
vy_grid = np.empty(bm2_shape, dtype=np.float64)
reproject(vx_src, vx_grid,
          src_transform=vel_transform, src_crs='EPSG:3031',
          dst_transform=bm2_transform, dst_crs=bm2_crs,
          resampling=Resampling.bilinear)
reproject(vy_src, vy_grid,
          src_transform=vel_transform, src_crs='EPSG:3031',
          dst_transform=bm2_transform, dst_crs=bm2_crs,
          resampling=Resampling.bilinear)
del vx_src, vy_src

vx_grid[~valid] = 0; vy_grid[~valid] = 0

print(f"  VX/VY regridded in {time.time()-t0:.1f} s")

# Pre-filter VX, VY with median filter to remove outliers
# that would send flowlines astray
medfilt_win = 11  # 11 km window
print(f"  Median-filtering VX, VY (window = {medfilt_win} km)...")
t0 = time.time()
vx_grid = median_filter(vx_grid, size=medfilt_win).astype(np.float64)
vy_grid = median_filter(vy_grid, size=medfilt_win).astype(np.float64)
vx_grid[~valid] = 0; vy_grid[~valid] = 0
speed = np.sqrt(vx_grid**2 + vy_grid**2)
print(f"  Median filter completed in {time.time()-t0:.1f} s")
print(f"  Speed range on ice: {speed[valid].min():.1f} – {speed[valid].max():.0f} m/yr")

# ============================================================
# Step 3: Forward flowline flux integration (C extension)
# ============================================================
print("Tracing flowlines and accumulating flux...")
t0 = time.time()

import ctypes

lib = ctypes.CDLL(os.path.join(out_dir, 'flux_accumulate.so'))
lib.flowline_integrate.restype = None
lib.flowline_integrate.argtypes = [
    ctypes.c_void_p,   # total_flux (output)
    ctypes.c_void_p,   # smb_flux
    ctypes.c_void_p,   # vx
    ctypes.c_void_p,   # vy
    ctypes.c_void_p,   # speed
    ctypes.c_void_p,   # valid_mask
    ctypes.c_void_p,   # seeds
    ctypes.c_int,       # n_seeds
    ctypes.c_int,       # nrows
    ctypes.c_int,       # ncols
    ctypes.c_int,       # max_steps
]

valid_flat = valid.ravel()
valid_idx = np.where(valid_flat)[0]
n_valid = len(valid_idx)

# Seed from every valid ice cell
seeds_c = np.ascontiguousarray(valid_idx, dtype=np.int32)

local_flux = smb_myr * dx * dx
local_flux[~valid] = 0
smb_flux_c = np.ascontiguousarray(local_flux.ravel(), dtype=np.float64)
total_flux = np.zeros(nrows * ncols, dtype=np.float64)
valid_int = np.ascontiguousarray(valid.ravel().astype(np.int32))
vx_c = np.ascontiguousarray(vx_grid.ravel(), dtype=np.float64)
vy_c = np.ascontiguousarray(vy_grid.ravel(), dtype=np.float64)
sp_c = np.ascontiguousarray(speed.ravel(), dtype=np.float64)

max_steps = 5000  # ~5000 km max flowline length
print(f"  Seeding {n_valid:,} flowlines (max {max_steps} steps each)...")

lib.flowline_integrate(
    total_flux.ctypes.data,
    smb_flux_c.ctypes.data,
    vx_c.ctypes.data, vy_c.ctypes.data, sp_c.ctypes.data,
    valid_int.ctypes.data,
    seeds_c.ctypes.data,
    ctypes.c_int(n_valid), ctypes.c_int(nrows), ctypes.c_int(ncols),
    ctypes.c_int(max_steps),
)

print(f"  Flowline integration completed in {time.time()-t0:.1f} s")
total_smb = (smb_myr[valid] * dx * dx).sum()
print(f"  Total SMB: {total_smb:.3e} m³/yr ({total_smb*917/1e12:.0f} Gt/yr)")

# ============================================================
# Step 4: Balance velocity + post-median filter
# ============================================================
flux_2d = total_flux.reshape(nrows, ncols)
balance_vel = np.abs(flux_2d) / (thickness * dx)
balance_vel[~valid] = 0
balance_vel[thickness <= 10] = 0
balance_vel[~np.isfinite(balance_vel)] = 0

# Light post-filter to remove residual artifacts
post_win = 5
print(f"Post-median filtering balance velocity (window = {post_win} km)...")
t0 = time.time()
balance_vel = median_filter(balance_vel, size=post_win).astype(np.float64)
balance_vel[~valid] = 0
print(f"  Post-filter completed in {time.time()-t0:.1f} s")

# Fill gaps at stagnation points / divides with observed speed
with rasterio.open(velocity_path) as src_vel:
    measures = np.empty(bm2_shape, dtype=np.float64)
    reproject(rasterio.band(src_vel, 1), measures,
              src_transform=src_vel.transform, src_crs=src_vel.crs,
              dst_transform=bm2_transform, dst_crs=bm2_crs,
              resampling=Resampling.bilinear)
    measures[measures < -1e30] = np.nan

gaps = valid & (balance_vel < 0.01) & np.isfinite(measures)
balance_vel[gaps] = measures[gaps]
print(f"  Filled {gaps.sum():,} stagnation-point gaps with observed speed")

print(f"\nBalance velocity statistics (all ice):")
bv = balance_vel[valid & (balance_vel > 0)]
for p in [10, 25, 50, 75, 90, 95, 99]:
    print(f"  {p:>3d}th percentile: {np.percentile(bv, p):>8.1f} m/yr")
print(f"        Max: {bv.max():>8.1f} m/yr")

# ============================================================
# Save GeoTIFF
# ============================================================
print(f"\nSaving {output_tif} ...")
out_profile = profile.copy()
out_profile.update(dtype='float32', count=1, nodata=-9999.0, compress='lzw')
with rasterio.open(output_tif, 'w', **out_profile) as dst:
    out = balance_vel.astype(np.float32)
    out[~valid] = -9999.0
    dst.write(out, 1)

print("\nDone! Output:", output_tif)
