Zero-Table O(1) Grid Routing via the Geometry of Integers (foliated-2d-pathfinder)Reference Python 3 implementation for the manuscript:Zero-Table $\mathcal{O}(1)$ Grid Routing and Obstacle Pivoting via the Geometry of IntegersAuthor: Hassan Baloui (new-mathematical-concepts.org)Building upon the discrete framework of Geometry of Integers, this package implements the Aligned Quadratic Foliation ($\mathbb{A}_2$) and the Global Vector Transport Operator ($\nabla_A$). Operating directly on 1D non-negative scalar addresses $n, m \in \mathbb{N}$ with zero precomputed lookup tables, it computes exact 8-connected Chebyshev distances, stateless next hops, canonical 1-zigzag turning waypoints, and axis-aligned bounding-box (AABB) obstacle bypasses in $\mathcal{O}(1)$ Word-RAM operations.Key Features$\mathcal{O}(1)$ Bijective Coordinate Algebra: Exact conversion between 1D integers $n \in \mathbb{N}$ and 2D half-plane coordinates $(g, d) \in \mathbb{Z} \times \mathbb{N}$ via hardware integer square roots (math.isqrt):$$N(g, d) = (\vert{}g\vert{} + d)^2 + (\vert{}g\vert{} + d) + g$$$\mathcal{O}(1)$ Distance, Waypoint & Stateless Next-Hop Oracle: Evaluates $D^*(n, m)$, canonical 1-zigzag turning points $X_1(n, m), X_2(n, m)$, and memoryless forwarding $\mathcal{H}(x, m)$ up to $10^{18}$+ scalar addresses in microseconds.$\mathcal{O}(1)$ Single-AABB Bypass & Multi-AABB Visibility Routing: Closed-form 2-zigzag track-interval testing, exterior boundary pivoting, $\mathcal{O}(1)$ Moore-isolation certification, and exact sparse multi-AABB visibility routing.Publication-Grade Grid Visualization: Exports 300-DPI .png and vector .pdf diagrams of the $\mathbb{A}_2$ address lattice, V-shell contours, obstacles, and geodesic trajectories.Contentsfoliated_2d_pathfinder.py — Self-contained Python 3 routing engine, CLI, and plotter.foliated_2d_pathfinder_manual.pdf / .tex — Technical reference and API manual.QuickstartCore routing uses only the Python 3.10+ standard library (math, dataclasses, heapq). Install matplotlib only if you want to export visual grid plots:pip install matplotlib
python foliated_2d_pathfinder.py
Programmatic API Exampleimport foliated_2d_pathfinder as f2d

# 1. Zero-table O(1) unobstructed routing at 10^18 scale
n = f2d.encode_a2(-4, 0)                      # 16
m = f2d.encode_a2(-4 + 10**9, 6 * 10**8)      # 2,559,999,990,400,000,012
info = f2d.one_zigzag_oracle(n, m, stream_path=False)
print(f"D* = {info['distance']}, Turning Waypoint X1 = {info['waypoint']}")

# 2. O(1) zero-penalty 2-zigzag AABB obstacle bypass (16 -> 45 around {11, 19})
boxes = f2d.parse_obstacles("11, 19")
res = f2d.solve_path_2d(16, 45, boxes)
f2d.print_route_summary({**res, "K": 4, "start": 16, "target": 45})
Citation@misc{Baloui2026Pathfinder,
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
LicenseReleased under the MIT License.
