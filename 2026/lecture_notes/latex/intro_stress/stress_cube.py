#!/usr/bin/env python3
"""
Generate a Cauchy stress tensor cube figure.
Outputs stress_cube.pdf for inclusion in LaTeX.
"""

import numpy as np
import matplotlib.pyplot as plt

# ── Cube geometry ──────────────────────────────────────────────
L = 1.0  # side length

# 8 vertices of the cube [0, L]^3
# Index: (x_bit, y_bit, z_bit) → vertex number = 4*x + 2*y + z
verts = np.array([[i, j, k] for i in [0, L] for j in [0, L] for k in [0, L]])

# 12 edges as pairs of vertex indices
edges = [
    (0, 1), (0, 2), (0, 4),  # from (0,0,0)
    (1, 3), (1, 5),          # from (0,0,L)
    (2, 3), (2, 6),          # from (0,L,0)
    (3, 7),                  # from (0,L,L)
    (4, 5), (4, 6),          # from (L,0,0)
    (5, 7),                  # from (L,0,L)
    (6, 7),                  # from (L,L,0)
]

# With azim=45, elev=25: viewer is at (+x,+y,+z).
# Hidden corner is (0,0,0) → vertex index 0.
# Hidden edges are the 3 from vertex 0.
hidden = {(0, 1), (0, 2), (0, 4)}

# ── Stress arrows on visible faces ───────────────────────────
# Visible faces (facing the viewer): +x, +y, +z
# Face centers:
#   +x face (x=L): center at (L, L/2, L/2)  — right front face
#   +y face (y=L): center at (L/2, L, L/2)  — left front face
#   +z face (z=L): center at (L/2, L/2, L)  — top face

h = L / 2.0
a = 0.35 * L  # arrow length

# Convention: on a positive face, σ_{ij} acts in the +x_j direction.
# Labels use i,j = x,y,z.
faces = {
    '+x': {
        'center': np.array([L, h, h]),
        'arrows': [
            (np.array([a, 0, 0]), r'$\sigma_{xx}$', 'blue'),    # normal
            (np.array([0, a, 0]), r'$\sigma_{xy}$', 'red'),     # shear
            (np.array([0, 0, a]), r'$\sigma_{xz}$', 'red'),     # shear
        ]
    },
    '+y': {
        'center': np.array([h, L, h]),
        'arrows': [
            (np.array([a, 0, 0]), r'$\sigma_{yx}$', 'red'),     # shear
            (np.array([0, a, 0]), r'$\sigma_{yy}$', 'blue'),    # normal
            (np.array([0, 0, a]), r'$\sigma_{yz}$', 'red'),     # shear
        ]
    },
    '+z': {
        'center': np.array([h, h, L]),
        'arrows': [
            (np.array([a, 0, 0]), r'$\sigma_{zx}$', 'red'),     # shear
            (np.array([0, a, 0]), r'$\sigma_{zy}$', 'red'),     # shear
            (np.array([0, 0, a]), r'$\sigma_{zz}$', 'blue'),    # normal
        ]
    },
}

# ── Plot ──────────────────────────────────────────────────────
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='3d')

# Draw cube edges
for i, j in edges:
    xs = [verts[i, 0], verts[j, 0]]
    ys = [verts[i, 1], verts[j, 1]]
    zs = [verts[i, 2], verts[j, 2]]
    if (i, j) in hidden:
        ax.plot(xs, ys, zs, color='gray', linewidth=0.8, linestyle='--')
    else:
        ax.plot(xs, ys, zs, color='black', linewidth=1.2)

# Draw stress arrows on the three visible faces
for face_name, face_data in faces.items():
    c = face_data['center']
    for dvec, label, color in face_data['arrows']:
        ax.quiver(c[0], c[1], c[2], dvec[0], dvec[1], dvec[2],
                  color=color, arrow_length_ratio=0.15, linewidth=1.8)
        # Label at arrow tip
        tip = c + dvec * 1.2
        ax.text(tip[0], tip[1], tip[2], label, fontsize=11,
                ha='center', va='center', color=color)

# Coordinate axes from the back lower corner (0,0,0)
origin = np.array([0.0, 0.0, 0.0])
axis_len = L * 0.45
axis_labels = [r'$x$', r'$y$', r'$z$']
for i, label in enumerate(axis_labels):
    direction = np.zeros(3)
    direction[i] = axis_len
    ax.quiver(origin[0], origin[1], origin[2],
              direction[0], direction[1], direction[2],
              color='black', arrow_length_ratio=0.18, linewidth=1.5)
    tip = origin + direction * 1.35
    ax.text(tip[0], tip[1], tip[2], label, fontsize=14,
            ha='center', va='center', fontweight='bold')

# ── View and formatting ──────────────────────────────────────
# azim=45: camera in the +x, +y quadrant → sees +x and +y faces as front faces
# elev=25: slightly above → sees +z (top) face
ax.view_init(elev=25, azim=45)

ax.set_box_aspect([1, 1, 1])
ax.set_axis_off()

fig.tight_layout()
fig.savefig('stress_cube.pdf', bbox_inches='tight', pad_inches=0.1)
print('Saved stress_cube.pdf')
