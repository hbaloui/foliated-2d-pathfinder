# Zero-Table O(1) Grid Routing via the Geometry of Integers (`foliated-2d-pathfinder`)
Reference Python 3 implementation for the manuscript:
> **Zero-Table O(1) Grid Routing and Obstacle Pivoting via the Geometry of Integers**
> **Author:** Hassan Baloui (https://www.new-mathematical-concepts.org)
Building upon the discrete framework of [*Geometry of Integers*](https://doi.org/10.5281/zenodo.23456718), this package implements the **Aligned Quadratic Foliation (`A₂`)** and the **Global Vector Transport Operator (`∇_A`)**. Operating directly on 1D non-negative scalar addresses `n, m ∈ ℕ` with **zero precomputed lookup tables**, it computes exact 8-connected Chebyshev distances, stateless next hops, canonical 1-zigzag turning waypoints, and axis-aligned bounding-box (AABB) obstacle bypasses in `O(1)` Word-RAM operations:
    N(g, d) = (|g| + d)² + (|g| + d) + g
## Key Features
- **O(1) Bijective Coordinate Algebra:** Exact conversion between 1D integers `n ∈ ℕ` and 2D half-plane coordinates `(g, d) ∈ ℤ × ℕ` via hardware integer square roots (`math.isqrt`).
- **O(1) Distance, Waypoint & Stateless Next-Hop Oracle:** Evaluates Chebyshev distance `D*(n, m)`, canonical 1-zigzag turning points `X₁(n, m)` and `X₂(n, m)`, and memoryless forwarding `H(x, m)` for scalar addresses exceeding `10¹⁸` in microseconds.
- **O(1) Single-AABB Bypass & Multi-AABB Visibility Routing:** Closed-form 2-zigzag track-interval testing, exterior boundary pivoting, `O(1)` Moore-isolation certification, and exact sparse multi-AABB visibility routing.
- **Publication-Grade Grid Visualization:** Exports 300-DPI `.png` and vector `.pdf` diagrams of the `A₂` address lattice, V-shell contours, obstacles, and geodesic trajectories.
---
## Contents
- `foliated_2d_pathfinder.py` — Self-contained Python 3 routing engine, CLI, and plotter.
- `foliated_2d_pathfinder_manual.pdf` / `.tex` — Technical reference and API manual.
---
## Quickstart
Core routing uses **only the Python 3.10+ standard library** (`math`, `dataclasses`, `heapq`). Install `matplotlib` only if you want to export visual grid plots:
~~~bash
pip install matplotlib
python foliated_2d_pathfinder.py
~~~
### Programmatic API Example
~~~python
import foliated_2d_pathfinder as f2d
# 1. Zero-table O(1) unobstructed routing at 10^18 scale
n = f2d.encode_a2(-4, 0)                            # 16
m = f2d.encode_a2(-4 + pow(10, 9), 6 * pow(10, 8))  # 2,559,999,990,400,000,012
info = f2d.one_zigzag_oracle(n, m, stream_path=False)
print(f"D* = {info['distance']}, Turning Waypoint X1 = {info['waypoint']}")
# 2. O(1) zero-penalty 2-zigzag AABB obstacle bypass (16 -> 45 around {11, 19})
boxes = f2d.parse_obstacles("11, 19")
res = f2d.solve_path_2d(16, 45, boxes)
f2d.print_route_summary({**res, "K": 4, "start": 16, "target": 45})
~~~
---
## Citation
~~~bibtex
@misc{Baloui2026Pathfinder,
  author = {Baloui, Hassan},
  title  = {Zero-Table {O(1)} Grid Routing and Obstacle Pivoting via the Geometry of Integers},
  year   = {2026},
  url    = {https://www.new-mathematical-concepts.org}
}

@misc{Baloui2026GOI,
  author    = {Baloui, Hassan},
  title     = {Geometry of Integers},
  publisher = {Zenodo},
  year      = {2026},
  doi       = {10.5281/zenodo.23456718},
  url       = {https://doi.org/10.5281/zenodo.23456718}
}
~~~
## License
Released under the MIT License.
