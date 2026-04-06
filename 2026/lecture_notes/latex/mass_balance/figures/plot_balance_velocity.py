#!/usr/bin/env python3
"""
Plot balance velocity, observed velocity, and their difference for Antarctica.

Reads BalanceVelocity.tif (produced by calculate_balance_velocity.py) and
MEaSUREs observed velocity, then generates four maps:

    - BalanceVelocity_map.png        (log-scale balance velocity + elevation contours)
    - ObservedVelocity_map.png       (log-scale MEaSUREs speed + elevation contours + grounding line)
    - BalanceVelocity_difference.png (balance − observed difference)
    - IceThickness_map.png           (BEDMAP2 ice thickness)
"""

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
import cartopy.crs as ccrs
import shapefile as shp
import os

plt.style.use('seaborn-v0_8-poster')

# ============================================================
# Paths
# ============================================================
base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
surface_path = os.path.join(base, 'TerrainModels/BEDMAP2/bedmap2_surface.tif')
thickness_path = os.path.join(base, 'TerrainModels/BEDMAP2/bedmap2_thickness.tif')
mask_path = os.path.join(base,
    'TerrainModels/BEDMAP2/bedmap2_icemask_grounded_and_shelves.tif')
velocity_path = os.path.join(base,
    'Glaciology/MEaSUREs Ice Flow Velocity/MEaSUREs_IceFlowSpeed_450m.tif')
balance_vel_path = os.path.join(base, 'BalanceVelocity/BalanceVelocity.tif')
gl_path = os.path.join(base,
    'Glaciology/MEaSUREs Antarctic Boundaries/GroundingLine_Antarctica_v2')
out_dir = os.path.join(base, 'BalanceVelocity')

# ============================================================
# Projection
# ============================================================
proj = ccrs.SouthPolarStereo(true_scale_latitude=-71)

# ============================================================
# Load data
# ============================================================
print("Loading surface elevation...")
with rasterio.open(surface_path) as src:
    surface = src.read(1).astype(np.float64)
    bm2_transform = src.transform
    bm2_crs = src.crs
    bm2_shape = src.shape
nrows, ncols = bm2_shape

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
valid = ((mask_data == 0) | (mask_data == 1)) & (thickness > 10)

print("Loading balance velocity...")
with rasterio.open(balance_vel_path) as src:
    balance_vel = src.read(1).astype(np.float64)
    balance_vel[balance_vel < -9000] = 0

print("Loading observed velocity (MEaSUREs)...")
with rasterio.open(velocity_path) as src_vel:
    measures = np.empty(bm2_shape, dtype=np.float64)
    reproject(rasterio.band(src_vel, 1), measures,
              src_transform=src_vel.transform, src_crs=src_vel.crs,
              dst_transform=bm2_transform, dst_crs=bm2_crs,
              resampling=Resampling.bilinear)
    measures[measures < -1e30] = np.nan

# ============================================================
# Coordinate grids and projection setup
# ============================================================
x_coords = np.arange(ncols) * bm2_transform.a + bm2_transform.c + bm2_transform.a / 2
y_coords = np.arange(nrows) * bm2_transform.e + bm2_transform.f + bm2_transform.e / 2
step = 5
X_plot, Y_plot = np.meshgrid(x_coords[::step], y_coords[::step])
# Subsampled surface elevation for contours
surf_sub = surface[::step, ::step].copy()
surf_sub[surf_sub <= 0] = np.nan
clevels = np.arange(500, 4500, 500)

# ============================================================
# Grounding line
# ============================================================
print("Loading grounding line...")
gl_sf = shp.Reader(gl_path)

def draw_grounding_line(ax, color='cyan', lw=0.5):
    """Draw grounding line on a cartopy axis, respecting multi-part shapes."""
    for shape in gl_sf.shapes():
        pts = np.array(shape.points)
        parts = list(shape.parts) + [len(pts)]
        for i in range(len(parts) - 1):
            seg = pts[parts[i]:parts[i+1]]
            if len(seg) < 2:
                continue
            ax.plot(seg[:, 0], seg[:, 1], color=color, linewidth=lw, zorder=3,
                    transform=proj)

def setup_ax(ax):
    """Add common map features and gridlines to a cartopy axis."""
    ax.coastlines(resolution='110m', linewidth=0.3, color='0.4')
    gl = ax.gridlines(draw_labels=True, linewidth=1.0, color='gray',
                      alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

# ============================================================
# Plot 1: Balance velocity (log scale, magma)
# ============================================================
print("Plotting balance velocity...")
bv_sub = balance_vel[::step, ::step].copy()
bv_sub[bv_sub <= 0] = np.nan

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': proj})
ax.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
im = ax.pcolormesh(X_plot, Y_plot, bv_sub, norm=LogNorm(vmin=1, vmax=3000),
                   cmap='magma', shading='auto', rasterized=True,
                   transform=proj)
draw_grounding_line(ax, color='cyan', lw=0.4)
cs = ax.contour(X_plot, Y_plot, surf_sub, levels=clevels, colors='white',
                linewidths=0.4, alpha=0.6, transform=proj)
ax.clabel(cs, inline=True, fontsize=6, fmt='%d m')
setup_ax(ax)
cb = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, shrink=0.8)
cb.set_label('Balance Velocity (m/yr)')
ax.set_title('Antarctic Balance Velocity', fontsize=14)
fig.savefig(os.path.join(out_dir, 'BalanceVelocity_map.png'),
            dpi=200, bbox_inches='tight')
plt.close(fig)

# ============================================================
# Plot 2: Observed surface velocity magnitude (MEaSUREs)
# ============================================================
print("Plotting observed velocity...")
measures_sub = measures[::step, ::step].copy()
measures_sub[measures_sub <= 0] = np.nan

fig2, ax2 = plt.subplots(figsize=(10, 10), subplot_kw={'projection': proj})
ax2.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
im2 = ax2.pcolormesh(X_plot, Y_plot, measures_sub, norm=LogNorm(vmin=1, vmax=3000),
                     cmap='magma', shading='auto', rasterized=True,
                     transform=proj)
draw_grounding_line(ax2, color='cyan', lw=0.4)
cs2 = ax2.contour(X_plot, Y_plot, surf_sub, levels=clevels, colors='white',
                  linewidths=0.4, alpha=0.6, transform=proj)
ax2.clabel(cs2, inline=True, fontsize=6, fmt='%d m')
setup_ax(ax2)
cb2 = fig2.colorbar(im2, ax=ax2, orientation='horizontal', pad=0.05, shrink=0.8)
cb2.set_label('Ice Flow Speed (m/yr)')
ax2.set_title('Observed Surface Velocity (MEaSUREs)', fontsize=14)
fig2.savefig(os.path.join(out_dir, 'ObservedVelocity_map.png'),
             dpi=200, bbox_inches='tight')
plt.close(fig2)

# ============================================================
# Plot 3: Difference (balance − observed)
# ============================================================
print("Plotting difference...")
diff = balance_vel - measures
diff[~valid] = np.nan; diff[~np.isfinite(diff)] = np.nan
diff_sub = diff[::step, ::step]
vmax = np.nanpercentile(np.abs(diff_sub[np.isfinite(diff_sub)]), 95)

fig3, ax3 = plt.subplots(figsize=(10, 10), subplot_kw={'projection': proj})
ax3.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
im3 = ax3.pcolormesh(X_plot, Y_plot, diff_sub, cmap='RdBu_r',
                     norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax),
                     shading='auto', rasterized=True,
                     transform=proj)
draw_grounding_line(ax3, color='black', lw=0.4)
setup_ax(ax3)
cb3 = fig3.colorbar(im3, ax=ax3, orientation='horizontal', pad=0.05, shrink=0.8)
cb3.set_label('Balance Velocity − MEaSUREs Velocity (m/yr)')
ax3.set_title('Balance Velocity − Observed Ice Flow Speed', fontsize=14)
fig3.savefig(os.path.join(out_dir, 'BalanceVelocity_difference.png'),
             dpi=200, bbox_inches='tight')
plt.close(fig3)

# ============================================================
# Plot 4: Ice thickness
# ============================================================
print("Plotting ice thickness...")
thick_sub = thickness[::step, ::step].copy()
thick_sub[thick_sub <= 0] = np.nan

fig4, ax4 = plt.subplots(figsize=(10, 10), subplot_kw={'projection': proj})
ax4.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
im4 = ax4.pcolormesh(X_plot, Y_plot, thick_sub, vmin=0, vmax=4500,
                     cmap='viridis', shading='auto', rasterized=True,
                     transform=proj)
draw_grounding_line(ax4, color='cyan', lw=0.4)
cs4 = ax4.contour(X_plot, Y_plot, surf_sub, levels=clevels, colors='white',
                  linewidths=0.4, alpha=0.6, transform=proj)
ax4.clabel(cs4, inline=True, fontsize=6, fmt='%d m')
setup_ax(ax4)
cb4 = fig4.colorbar(im4, ax=ax4, orientation='horizontal', pad=0.05, shrink=0.8)
cb4.set_label('Ice Thickness (m)')
ax4.set_title('Antarctic Ice Thickness (BEDMAP2)', fontsize=14)
fig4.savefig(os.path.join(out_dir, 'IceThickness_map.png'),
             dpi=200, bbox_inches='tight')
plt.close(fig4)

print("\nDone! Plots saved to:", out_dir)
