"""Zero-table O(1) geodesic routing and AABB obstacle pivoting on A_2 (Passes 1-3)."""

from dataclasses import dataclass
import heapq
import math


def sgn(x: int) -> int:
    """Return the three-way sign of an integer x in {-1, 0, +1}."""
    return (x > 0) - (x < 0)


def decode_a2(n: int) -> tuple[int, int, int]:
    """Decode scalar address n in N into aligned coordinates (g, d, k) in O(1) time."""
    if n < 0:
        raise ValueError(f"Scalar address n must be non-negative, got {n}.")
    k = math.isqrt(n)
    g = n - (k * k + k)
    d = k - abs(g)
    return g, d, k


def encode_a2(g: int, d: int) -> int:
    """Encode aligned coordinates (g, d) in Z x N into scalar address N_2(g, d)."""
    if d < 0:
        raise ValueError(f"Aligned depth d must be non-negative, got {d}.")
    k = abs(g) + d
    return k * k + k + g


def chebyshev_dist(n: int, m: int) -> int:
    """Compute the exact Aligned Chebyshev distance D*(n, m) = max(|dg|, |dd|)."""
    gn, dn, _ = decode_a2(n)
    gm, dm, _ = decode_a2(m)
    return max(abs(gm - gn), abs(dm - dn))


def transport_a2(n: int, dg: int, dd: int) -> int:
    """Evaluate the Master Aligned Vector Transport Operator nabla_A(n, (dg, dd))."""
    gn, dn, _ = decode_a2(n)
    if dn + dd < 0:
        raise ValueError("Target aligned depth dn + dd must be non-negative.")
    gm = gn + dg
    return (
        n
        + dg * (2 * gn + dg + 1)
        + dd * (2 * dn + dd + 2 * abs(gm) + 1)
        + (2 * dn + 1) * (abs(gm) - abs(gn))
    )


def next_hop_a2(x: int, m: int) -> int:
    """Stateless O(1) next-hop operator H(x, m) strictly decreasing D*(x, m) by 1."""
    if x == m:
        return m
    gx, dx, _ = decode_a2(x)
    gm, dm, _ = decode_a2(m)
    return transport_a2(x, sgn(gm - gx), sgn(dm - dx))


def moore_neighbors_a2(n: int) -> list[int]:
    """Return the sheared Moore neighborhood M_a(n) (8 interior, 5 at d=0) in O(1) time."""
    gn, dn, _ = decode_a2(n)
    neighbors: list[int] = []
    for eg in (-1, 0, 1):
        for ed in (-1, 0, 1):
            if (eg != 0 or ed != 0) and dn + ed >= 0:
                neighbors.append(encode_a2(gn + eg, dn + ed))
    return neighbors


def one_zigzag_oracle(n: int, m: int, stream_path: bool = True) -> dict:
    """Synthesize D*(n, m), canonical 1-zigzag waypoints, and optional hop stream."""
    gn, dn, _ = decode_a2(n)
    gm, dm, _ = decode_a2(m)
    dg, dd = gm - gn, dm - dn
    eg, ed = sgn(dg), sgn(dd)
    s_diag = min(abs(dg), abs(dd))
    s_card = max(abs(dg), abs(dd)) - s_diag
    waypoint = transport_a2(n, eg * s_diag, ed * s_diag)  # X_1(n, m)
    v_card = (dg - eg * s_diag, 0) if abs(dg) >= abs(dd) else (0, dd - ed * s_diag)
    alt_waypoint = transport_a2(n, v_card[0], v_card[1])  # X_2(n, m)
    path: list[int] = []
    if stream_path:
        cur = n
        path.append(cur)
        while cur != m:
            cur = next_hop_a2(cur, m)
            path.append(cur)
    return {
        "distance": s_diag + s_card,
        "waypoint": waypoint,
        "alt_waypoint": alt_waypoint,
        "s_diag": s_diag,
        "s_card": s_card,
        "path": path,
    }


@dataclass(frozen=True)
class AABB:
    """Axis-aligned bounding box [g_min, g_max] x [d_min, d_max] on A_2."""

    g_min: int
    g_max: int
    d_min: int
    d_max: int

    def contains(self, n: int) -> bool:
        """Return True iff scalar state n lies inside this bounding box."""
        g, d, _ = decode_a2(n)
        return self.g_min <= g <= self.g_max and self.d_min <= d <= self.d_max

    def blocked_track_interval(self, u: int, v: int) -> tuple[int, int, int]:
        """Return (rho_lo, rho_hi, s_card) of blocked 2-zigzag tracks in [0, s_card]."""
        gu, du, _ = decode_a2(u)
        gv, dv, _ = decode_a2(v)
        dg, dd = gv - gu, dv - du
        eg, ed = sgn(dg) or 1, sgn(dd) or 1
        s_diag = min(abs(dg), abs(dd))
        s_card = max(abs(dg), abs(dd)) - s_diag
        g_box = (self.g_min - gu, self.g_max - gu) if eg > 0 else (gu - self.g_max, gu - self.g_min)
        d_box = (self.d_min - du, self.d_max - du) if ed > 0 else (du - self.d_max, du - self.d_min)
        maj, minr = (g_box, d_box) if abs(dg) >= abs(dd) else (d_box, g_box)
        u1, u2 = max(0, maj[0]), min(s_card + s_diag, maj[1])
        v1, v2 = max(0, minr[0]), min(s_diag, minr[1])
        if u1 > u2 or v1 > v2:
            return 1, 0, s_card
        if s_diag == 0:
            return 0, s_card, s_card
        lo = 0 if v2 == s_diag else max(0, u1 - v2)
        hi = s_card if v1 == 0 else min(s_card, u2 - v1)
        return (lo, hi, s_card) if lo <= hi else (1, 0, s_card)

    def intersects_1zigzag(self, u: int, v: int, cardinal_first: bool = False) -> bool:
        """Test in O(1) time whether u -> X_1(u,v) -> v (or X_2) intersects the box."""
        lo, hi, s_card = self.blocked_track_interval(u, v)
        rho = s_card if cardinal_first else 0
        return lo <= rho <= hi

    def admissible_rim_pivots(self) -> list[int]:
        """Generate the <= 12 exterior rim pivots C_adm(B_loc) with d >= 0."""
        g0, g1, d0, d1 = self.g_min, self.g_max, self.d_min, self.d_max
        coords = [
            (g0, d1 + 1), (g1, d1 + 1), (g0, d0 - 1), (g1, d0 - 1),
            (g0 - 1, d0), (g1 + 1, d0), (g0 - 1, d1), (g1 + 1, d1),
            (g0 - 1, d1 + 1), (g1 + 1, d1 + 1), (g0 - 1, d0 - 1), (g1 + 1, d0 - 1),
        ]
        seen: list[int] = []
        for g, d in coords:
            if d >= 0:
                c = encode_a2(g, d)
                if not self.contains(c) and c not in seen:
                    seen.append(c)
        return seen


def parse_obstacles(spec: str) -> list[AABB]:
    """Parse integer IDs (merging adjacent cells into AABBs) or 'g0:g1,d0:d1' boxes."""
    spec = spec.strip()
    if not spec or spec.lower() == "none":
        return []
    boxes: list[AABB] = []
    cells: set[tuple[int, int]] = set()
    for group in spec.split(";"):
        group = group.strip()
        if not group:
            continue
        if ":" in group:
            g_part, d_part = group.split(",")
            g0, g1 = map(int, g_part.split(":"))
            d0, d1 = map(int, d_part.split(":"))
            boxes.append(AABB(min(g0, g1), max(g0, g1), max(0, min(d0, d1)), max(d0, d1)))
        else:
            for tok in group.replace(",", " ").split():
                g, d, _ = decode_a2(int(tok))
                cells.add((g, d))
    while cells:
        g0, d0 = min(cells, key=lambda c: (c[1], c[0]))
        g1 = g0
        while (g1 + 1, d0) in cells:
            g1 += 1
        d1 = d0
        while all((g, d1 + 1) in cells for g in range(g0, g1 + 1)):
            d1 += 1
        for d in range(d0, d1 + 1):
            for g in range(g0, g1 + 1):
                cells.discard((g, d))
        boxes.append(AABB(g0, g1, d0, d1))
    return boxes


def unobstructed_1zigzag_leg(u: int, v: int, boxes: list[AABB]) -> dict | None:
    """Return an unobstructed 1-zigzag leg u -> W -> v (X_1 or X_2), or None if blocked."""
    if not any(b.intersects_1zigzag(u, v, cardinal_first=False) for b in boxes):
        return one_zigzag_oracle(u, v, stream_path=True)
    if not any(b.intersects_1zigzag(u, v, cardinal_first=True) for b in boxes):
        info = one_zigzag_oracle(u, v, stream_path=False)
        w = info["alt_waypoint"]
        leg1 = one_zigzag_oracle(u, w, stream_path=True)["path"]
        leg2 = one_zigzag_oracle(w, v, stream_path=True)["path"]
        return {**info, "waypoint": w, "path": leg1 + leg2[1:]}
    return None


def solve_path_2d(start: int, target: int, boxes: list[AABB]) -> dict:
    """Exact O(1) geodesic / single-AABB pivot solver and multi-AABB visibility router."""
    if any(b.contains(start) or b.contains(target) for b in boxes):
        return {"reachable": False, "reason": "Start or target lies inside an obstacle.", "boxes": boxes}
    if start == target:
        return {
            "reachable": True, "regime": "Trivial 0-Hop (start == target)", "distance": 0,
            "pivots": [], "waypoints": [], "path": [start], "boxes": boxes,
        }
    if all(any(b.contains(u) for b in boxes) for u in moore_neighbors_a2(start)):
        return {
            "reachable": False,
            "reason": f"Start n={start} is isolated (all Moore neighbors in M_a({start}) are blocked).",
            "boxes": boxes,
        }
    if all(any(b.contains(u) for b in boxes) for u in moore_neighbors_a2(target)):
        return {
            "reachable": False,
            "reason": f"Target m={target} is isolated (all Moore neighbors in M_a({target}) are blocked).",
            "boxes": boxes,
        }
    direct = unobstructed_1zigzag_leg(start, target, boxes)
    if direct is not None:
        return {
            "reachable": True, "regime": "Unobstructed 1-Zigzag O(1)", "distance": direct["distance"],
            "pivots": [], "waypoints": [direct["waypoint"]], "path": direct["path"], "boxes": boxes,
        }
    if len(boxes) == 1:
        box = boxes[0]
        lo, hi, s_card = box.blocked_track_interval(start, target)
        if lo > 0 or hi < s_card:
            rho = 0 if lo > 0 else hi + 1
            gn, dn, _ = decode_a2(start)
            gm, dm, _ = decode_a2(target)
            dg, dd = gm - gn, dm - dn
            eg, ed, s_diag = sgn(dg), sgn(dd), min(abs(dg), abs(dd))
            uc = (eg, 0) if abs(dg) >= abs(dd) else (0, ed)
            y1 = transport_a2(start, rho * uc[0], rho * uc[1])
            y2 = transport_a2(y1, eg * s_diag, ed * s_diag)
            p1 = one_zigzag_oracle(start, y1)["path"]
            p2 = one_zigzag_oracle(y1, y2)["path"]
            p3 = one_zigzag_oracle(y2, target)["path"]
            return {
                "reachable": True, "regime": f"O(1) 2-Zigzag Track Bypass (rho*={rho})",
                "distance": chebyshev_dist(start, target), "pivots": [],
                "waypoints": [y1, y2], "path": p1 + p2[1:] + p3[1:], "boxes": boxes,
            }
    nodes = [start, target]
    for b in boxes:
        for c in b.admissible_rim_pivots():
            if c not in nodes and not any(ob.contains(c) for ob in boxes):
                nodes.append(c)
    dist_map: dict[int, int] = {start: 0}
    prev_map: dict[int, tuple[int, dict]] = {}
    pq: list[tuple[int, int]] = [(0, start)]
    while pq:
        d_u, u = heapq.heappop(pq)
        if u == target:
            break
        if d_u > dist_map.get(u, math.inf):
            continue
        for v in nodes:
            if v == u:
                continue
            nd = d_u + chebyshev_dist(u, v)
            if nd < dist_map.get(v, math.inf):
                leg = unobstructed_1zigzag_leg(u, v, boxes)
                if leg is not None:
                    dist_map[v] = nd
                    prev_map[v] = (u, leg)
                    heapq.heappush(pq, (nd, v))
    if target not in dist_map:
        return {"reachable": False, "reason": "No unobstructed path exists.", "boxes": boxes}
    chain: list[tuple[int, int, dict]] = []
    cur = target
    while cur != start:
        p_node, leg = prev_map[cur]
        chain.append((p_node, cur, leg))
        cur = p_node
    chain.reverse()
    pivots = [v for _, v, _ in chain[:-1]]
    waypoints: list[int] = []
    full_path = [start]
    for u, v, leg in chain:
        w = leg["waypoint"]
        if w not in (u, v) and w not in waypoints:
            waypoints.append(w)
        full_path.extend(leg["path"][1:])
    regime = "O(1) Single-AABB Rim Pivot" if len(boxes) == 1 else "Multi-AABB Visibility Graph"
    return {
        "reachable": True, "regime": regime, "distance": dist_map[target],
        "pivots": pivots, "waypoints": pivots + waypoints, "path": full_path, "boxes": boxes,
    }


def prompt_user_inputs() -> dict:
    """Prompt user for K, start, target, and obstacles (with TikZ Example 4.1 defaults)."""
    k_str = input("Grid half-width/height K [default 4]: ").strip()
    start_str = input("Starting integer n [default 16 = N(-4,0)]: ").strip()
    target_str = input("Target integer m [default 45 = N(+3,3)]: ").strip()
    obs_str = input("Obstacles (IDs e.g. '0,2,6,12' or '11,19' or 'g0:g1,d0:d1') [default '0,2,6,12']: ").strip()
    k_disp = int(k_str) if k_str else 4
    start = int(start_str) if start_str else 16
    target = int(target_str) if target_str else 45
    boxes = parse_obstacles(obs_str if obs_str != "" else "0,2,6,12")
    gs, ds, _ = decode_a2(start)
    gt, dt, _ = decode_a2(target)
    k_disp = max(k_disp, abs(gs), ds, abs(gt), dt)
    res = solve_path_2d(start, target, boxes)
    res.update({"K": k_disp, "start": start, "target": target})
    return res


def print_route_summary(res: dict) -> None:
    """Print a concise console report of the foliated 2D routing solution."""
    gs, ds, ks = decode_a2(res["start"])
    gt, dt, kt = decode_a2(res["target"])
    print(f"\n=== A_2 Foliated Pathfinder Summary (Grid K={res['K']}) ===")
    print(f"Start  : n = {res['start']} -> (g={gs:+d}, d={ds}, k={ks})")
    print(f"Target : m = {res['target']} -> (g={gt:+d}, d={dt}, k={kt})")
    if not res["reachable"]:
        print(f"Status : UNREACHABLE ({res['reason']})")
        return
    print(f"Regime : {res['regime']}")
    print(f"Dist D*: {res['distance']} hops | Waypoints/Pivots: {res['waypoints']}")
    print(f"Path   : {' -> '.join(map(str, res['path']))}")


def plot_tikz_style_grid(
    res: dict,
    png_path: str = "foliated_2d_path.png",
    pdf_path: str = "foliated_2d_path.pdf",
) -> None:
    """Render a TikZ-styled figure of the A_2 grid, shells, obstacles, and geodesic path."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyArrowPatch, Patch, Rectangle

    k_max = res["K"]
    boxes: list[AABB] = res["boxes"]
    path: list[int] = res.get("path", [])
    waypoints: set[int] = set(res.get("waypoints", []))
    pivots: set[int] = set(res.get("pivots", []))
    path_set: set[int] = set(path)
    show_labels = k_max <= 18

    fig_w = max(8.5, min(14.5, (2 * k_max + 1) * 0.68 + 2.2))
    fig_h = max(5.8, min(10.5, (k_max + 1) * 0.68 + 2.4))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # 1. Foliated lattice cells (g, d) in [-K, K] x [0, K] (d=0 at top, matching TikZ)
    for d in range(k_max + 1):
        for g in range(-k_max, k_max + 1):
            n_val = encode_a2(g, d)
            in_obs = any(b.contains(n_val) for b in boxes)
            if in_obs:
                fc, ec, lw, hatch = "#334155", "#0F172A", 0.9, "///"
            elif n_val in (res["start"], res["target"]):
                fc, ec, lw, hatch = "#D1FAE5", "#059669", 1.5, None
            elif n_val in waypoints:
                fc, ec, lw, hatch = "#FEF3C7", "#D97706", 1.3, None
            elif n_val in path_set:
                fc, ec, lw, hatch = "#DBEAFE", "#60A5FA", 0.9, None
            elif g == 0:
                fc, ec, lw, hatch = "#F3F6FB", "#C6CDD8", 0.6, None
            else:
                fc, ec, lw, hatch = "#FFFFFF", "#C6CDD8", 0.5, None
            ax.add_patch(
                Rectangle((g - 0.5, -d - 0.5), 1.0, 1.0, facecolor=fc, edgecolor=ec, lw=lw, hatch=hatch, zorder=1)
            )
            if show_labels and n_val not in path_set and not in_obs:
                ax.text(g, -d, str(n_val), ha="center", va="center", fontsize=7.5, color="#64748B", zorder=2)
            elif show_labels and in_obs:
                ax.text(g, -d, str(n_val), ha="center", va="center", fontsize=7.5, weight="bold", color="white", zorder=3)

    # 2. Subtle V-shaped shell guide contours k = |g| + d = const
    for k_sh in range(1, 2 * k_max + 1):
        pts = [(g, -(k_sh - abs(g))) for g in range(-k_max, k_max + 1) if 0 <= k_sh - abs(g) <= k_max]
        if len(pts) >= 2:
            xs, ys = zip(*pts)
            ax.plot(xs, ys, ls=":", lw=0.75, color="#94A3B8", alpha=0.55, zorder=2)

    # 3. Evaluated admissible exterior corners C_adm (subtle hollow amber markers)
    adm_pts = [decode_a2(c)[:2] for b in boxes for c in b.admissible_rim_pivots() if not any(ob.contains(c) for ob in boxes)]
    if adm_pts:
        cg, cd = zip(*[(g, -d) for g, d in adm_pts if abs(g) <= k_max and d <= k_max])
        ax.scatter(cg, cd, s=95, facecolors="none", edgecolors="#D97706", lw=1.1, ls="--", zorder=4)

    # 4. Faint dashed unobstructed reference path if deflected by an obstacle
    ref = one_zigzag_oracle(res["start"], res["target"], stream_path=True)["path"]
    if res["reachable"] and path != ref:
        rx, ry = zip(*[(decode_a2(u)[0], -decode_a2(u)[1]) for u in ref])
        ax.plot(rx, ry, ls="--", lw=1.4, color="#EF4444", alpha=0.55, zorder=3)

    # 5. Optimal geodesic trajectory arrows and node callout badges
    if res["reachable"] and len(path) >= 2:
        for u, v in zip(path[:-1], path[1:]):
            gu, du, _ = decode_a2(u)
            gv, dv, _ = decode_a2(v)
            ax.add_patch(
                FancyArrowPatch(
                    (gu, -du), (gv, -dv), arrowstyle="-|>", mutation_scale=12,
                    lw=2.0, color="#1D4ED8", shrinkA=9, shrinkB=9, zorder=5,
                )
            )
    nodes_to_badge = path if path else [res["start"], res["target"]]
    for u in nodes_to_badge:
        gu, du, _ = decode_a2(u)
        if u in (res["start"], res["target"]):
            bbox = dict(boxstyle="round,pad=0.26", fc="#D1FAE5", ec="#059669", lw=1.5)
            tag = "Start" if u == res["start"] else "Target"
            ax.text(gu, -du + 0.34, tag, ha="center", va="bottom", fontsize=7, weight="bold", color="#065F46", zorder=7)
        elif u in waypoints:
            bbox = dict(boxstyle="circle,pad=0.22", fc="#FEF3C7", ec="#D97706", lw=1.5)
            tag = f"Pivot {u}" if u in pivots else f"Turn {u}"
            ax.text(gu, -du - 0.35, tag, ha="center", va="top", fontsize=7, weight="bold", color="#92400E", zorder=7)
        else:
            bbox = dict(boxstyle="round,pad=0.16", fc="#DBEAFE", ec="#60A5FA", lw=0.9)
        if show_labels or u in (res["start"], res["target"]) or u in waypoints:
            ax.text(gu, -du, str(u), ha="center", va="center", fontsize=8, weight="bold", color="#0F172A", bbox=bbox, zorder=6)

    ax.set_xlim(-k_max - 0.65, k_max + 0.65)
    ax.set_ylim(-k_max - 0.75, 0.75)
    ax.set_xticks(range(-k_max, k_max + 1))
    ax.set_xticklabels([f"{g:+d}" if g != 0 else "0" for g in range(-k_max, k_max + 1)], fontsize=8)
    ax.set_yticks([-d for d in range(k_max + 1)])
    ax.set_yticklabels([str(d) for d in range(k_max + 1)], fontsize=8)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel(r"Signed Generator $g \in \mathbb{Z}$ (Pronic Spine $V_0$ at $g=0$)", fontsize=9.5, weight="bold")
    ax.set_ylabel(r"Aligned Depth $d = k - |g| \in \mathbb{N}$", fontsize=9.5, weight="bold")
    ax.set_aspect("equal")

    hop_str = " -> ".join(map(str, path)) if len(path) <= 14 else f"{' -> '.join(map(str, path[:6]))} -> ... -> {path[-1]}"
    sub = f"D* = {res.get('distance', 'N/A')} hops | Waypoints: {res.get('waypoints', [])} | Path: [{hop_str}]"
    ax.set_title(f"Aligned Quadratic Foliation A_2 — {res.get('regime', 'Unreachable')}\n{sub}", fontsize=10, weight="bold", pad=14)

    legend_items = [
        Patch(fc="#D1FAE5", ec="#059669", lw=1.3, label="Start / Target"),
        Patch(fc="#FEF3C7", ec="#D97706", lw=1.3, label="O(1) Waypoint / Rim Pivot"),
        Patch(fc="#DBEAFE", ec="#60A5FA", lw=0.9, label="Geodesic Hop"),
        Patch(fc="#334155", ec="#0F172A", hatch="///", label="Blocked Obstacle AABB"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#D97706", lw=0, label=r"Admissible Rim $\mathcal{C}_{\mathrm{adm}}$"),
    ]
    ax.legend(handles=legend_items, loc="upper center", bbox_to_anchor=(0.5, -0.04), ncol=3, frameon=True, fontsize=8)
    plt.tight_layout()
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"Saved high-res figures: '{png_path}' (300 DPI) and '{pdf_path}'.")
    plt.show()


if __name__ == "__main__":
    try:
        route_result = prompt_user_inputs()
    except EOFError:
        route_result = {**solve_path_2d(16, 45, parse_obstacles("0,2,6,12")), "K": 4, "start": 16, "target": 45}
    print_route_summary(route_result)
    plot_tikz_style_grid(route_result)
