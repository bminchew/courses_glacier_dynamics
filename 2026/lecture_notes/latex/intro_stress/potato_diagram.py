#!/usr/bin/env python3
"""
Generate a 'potato diagram' for the stress tensor symmetry proof:
an arbitrary body with volume V, surface S, position vector r,
traction T, body force f, and outward normal n-hat.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

fig, ax = plt.subplots(figsize=(5, 4))

# ── Draw an irregular closed body (potato shape) ──
theta = np.linspace(0, 2 * np.pi, 200)
r_base = 1.5
# Add harmonics for irregularity
r = (r_base
     + 0.3 * np.cos(2 * theta)
     + 0.15 * np.sin(3 * theta)
     + 0.1 * np.cos(5 * theta - 0.5))
cx, cy = 2.5, 2.0  # center
xs = cx + r * np.cos(theta)
ys = cy + r * np.sin(theta)

ax.fill(xs, ys, color='lightyellow', edgecolor='black', linewidth=1.5, zorder=2)

# ── Label V (volume) in the interior ──
ax.text(cx - 0.3, cy - 0.3, r'$V$', fontsize=16, ha='center', va='center',
        fontstyle='italic')

# ── Label S (surface) along the boundary ──
ax.annotate(r'$S$', xy=(xs[30], ys[30]), xytext=(xs[30] + 0.5, ys[30] + 0.5),
            fontsize=14, fontstyle='italic',
            arrowprops=dict(arrowstyle='->', color='black', lw=1.0))

# ── Pick a point P inside the body for position vector r ──
px, py = cx + 0.3, cy + 0.2

# Origin for position vector
ox, oy = 0.3, 0.3
ax.plot(ox, oy, 'ko', markersize=4, zorder=5)
ax.annotate('', xy=(px, py), xytext=(ox, oy),
            arrowprops=dict(arrowstyle='->', color='blue', lw=1.8))
ax.text((ox + px) / 2 - 0.25, (oy + py) / 2 + 0.1, r'$\mathbf{r}$',
        fontsize=14, color='blue', fontweight='bold')

# ── Pick a point on the surface for traction and normal ──
idx = 160  # point on the upper-right part of the boundary
sx_pt, sy_pt = xs[idx], ys[idx]

# Outward normal at that point
# Approximate normal from the curve
dx = xs[(idx + 2) % len(xs)] - xs[(idx - 2) % len(xs)]
dy = ys[(idx + 2) % len(ys)] - ys[(idx - 2) % len(ys)]
nx, ny = -dy, dx  # outward normal (rotated 90° clockwise)
norm = np.sqrt(nx**2 + ny**2)
nx, ny = nx / norm, ny / norm

# Draw outward normal n-hat
nlen = 0.7
ax.annotate('', xy=(sx_pt + nlen * nx, sy_pt + nlen * ny),
            xytext=(sx_pt, sy_pt),
            arrowprops=dict(arrowstyle='->', color='green', lw=1.8))
ax.text(sx_pt + nlen * nx + 0.15, sy_pt + nlen * ny + 0.1,
        r'$\hat{\mathbf{n}}$', fontsize=14, color='green', fontweight='bold')

# Draw traction T (at an angle to the surface, not aligned with normal)
tx, ty = 0.6 * nx + 0.5 * (-ny), 0.6 * ny + 0.5 * nx  # mix of normal and tangential
tlen = 0.9
tnorm = np.sqrt(tx**2 + ty**2)
tx, ty = tx / tnorm * tlen, ty / tnorm * tlen
ax.annotate('', xy=(sx_pt + tx, sy_pt + ty),
            xytext=(sx_pt, sy_pt),
            arrowprops=dict(arrowstyle='->', color='red', lw=1.8))
ax.text(sx_pt + tx + 0.15, sy_pt + ty - 0.15,
        r'$\mathbf{T}$', fontsize=14, color='red', fontweight='bold')

# ── Body force f (downward arrow inside the body) ──
fx_pt, fy_pt = cx + 0.6, cy - 0.5
ax.annotate('', xy=(fx_pt, fy_pt - 0.7), xytext=(fx_pt, fy_pt),
            arrowprops=dict(arrowstyle='->', color='purple', lw=1.8))
ax.text(fx_pt + 0.15, fy_pt - 0.35, r'$\mathbf{f}$',
        fontsize=14, color='purple', fontweight='bold')

# ── Small surface element dS near the traction point ──
# Draw a small arc to suggest dS
t_start, t_end = idx - 8, idx + 8
ax.plot(xs[t_start:t_end], ys[t_start:t_end], color='black', linewidth=3.5, zorder=3)
ax.text(xs[idx] - 0.05, ys[idx] - 0.35, r'$dS$', fontsize=11,
        ha='center', color='black')

# ── Formatting ──
ax.set_xlim(-0.5, 5.2)
ax.set_ylim(-0.5, 4.5)
ax.set_aspect('equal')
ax.axis('off')

fig.tight_layout()
fig.savefig('potato_diagram.pdf', bbox_inches='tight', pad_inches=0.1)
print('Saved potato_diagram.pdf')
