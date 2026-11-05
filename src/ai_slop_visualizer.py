#!/usr/bin/env python3
# vis_cells.py — закрашивает ячейки как кольцевые сектора
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

MAX_SEGMENTS = 300_000

PAT = re.compile(
    r'^\s*\d+\s+\d+\s*$'
    r'\s*\n\s*uzp\s+([-+0-9.eE]+)\s*$'
    r'\s*\n\s*thetas\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$'
    r'\s*\n\s*radii\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)',
    re.MULTILINE,
)

def parse_fast(text):
    arr = np.asarray(PAT.findall(text), dtype=float)
    if arr.size == 0:
        raise SystemExit("ничего не распарсилось")
    uzp, th0, th1, r0, r1 = arr.T
    return th0, th1, r0, r1, uzp


def visualize_cells(th0, th1, r0, r1, uzp,
                    max_segments=MAX_SEGMENTS,
                    arc_pts=4,
                    r_pad=None):
    """
    Каждая ячейка → закрашенный кольцевой сектор [θ0..θ1] × [r0..r1].

    arc_pts — точек на каждой дуге (4 обычно достаточно,
              для гладкости ставьте 8–16)
    r_pad   — добавка к толщине кольца с каждой стороны. Если None —
              0.5 % от среднего радиуса (для ваших данных r0≈r1,
              без этого ячейки вырождаются в невидимые отрезки).
    """
    N = len(uzp)
    print(f"ячеек: {N:,}")

    if N > max_segments:
        rng = np.random.default_rng(0)
        idx = rng.choice(N, max_segments, replace=False)
        th0, th1, r0, r1, uzp = (a[idx] for a in (th0, th1, r0, r1, uzp))
        N = max_segments
        print(f"даунсэмплировано до {N:,}")

    if r_pad is None:
        r_pad = 0.005 * 0.5 * (np.mean(r0) + np.mean(r1))
    print(f"r_pad = {r_pad:g}")

    # угловые точки: TH имеет форму (N, arc_pts)
    t  = np.linspace(0.0, 1.0, arc_pts)
    TH = th0[:, None] + (th1 - th0)[:, None] * t[None, :]

    r_out = np.maximum(r0, r1) + r_pad                          # (N,)
    r_in  = np.clip(np.minimum(r0, r1) - r_pad, 0.0, None)      # (N,)

    # Внешняя дуга: θ0 → θ1; внутренняя: θ1 → θ0 (замыкаем контур).
    x_out = r_out[:, None] * np.cos(TH)
    y_out = r_out[:, None] * np.sin(TH)
    x_in  = r_in [:, None] * np.cos(TH[:, ::-1])
    y_in  = r_in [:, None] * np.sin(TH[:, ::-1])

    verts = np.concatenate([
        np.stack([x_out, y_out], axis=-1),   # (N, arc_pts, 2)
        np.stack([x_in,  y_in ], axis=-1),   # (N, arc_pts, 2)
    ], axis=1)                                # (N, 2*arc_pts, 2)

    fig, ax = plt.subplots(figsize=(10, 10))
    pc = PolyCollection(verts, cmap="turbo",
                        edgecolors="none", antialiased=False)
    pc.set_array(uzp)          # цвет ячейки = uzp
    ax.add_collection(pc)      # ОДИН artist на все ячейки
    ax.autoscale()
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, ls=":", alpha=0.4)
    fig.colorbar(pc, ax=ax, label="uzp")
    ax.set_title(f"Ячеек: {N:,}  |  uzp ∈ [{uzp.min():.2e}, {uzp.max():.2e}]")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    text = open(sys.argv[1], encoding="utf-8").read()
    # при желании: r_pad=... , max_segments=..., arc_pts=8
    visualize_cells(*parse_fast(text))