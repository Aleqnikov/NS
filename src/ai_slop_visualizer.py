#!/usr/bin/env python3
# vis_cells.py — ячейки (uzp, urp) + векторное поле скорости в декартовой СК
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.path import Path as MplPath

try:
    from scipy.interpolate import griddata
except ImportError:
    griddata = None

MAX_SEGMENTS = 300_000
MAX_ARROWS   =  20_000

PAT = re.compile(
    r'^\s*\d+\s+\d+\s*$'
    r'\s*\n\s*uzp\s+([-+0-9.eE]+)\s*$'
    r'\s*\n\s*urp\s+([-+0-9.eE]+)\s*$'
    r'\s*\n\s*thetas\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$'
    r'\s*\n\s*radii\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)',
    re.MULTILINE,
)


def parse_fast(text):
    arr = np.asarray(PAT.findall(text), dtype=float)
    if arr.size == 0:
        raise SystemExit("ничего не распарсилось")
    uzp, urp, th0, th1, r0, r1 = arr.T
    return th0, th1, r0, r1, uzp, urp


def build_verts(th0, th1, r0, r1, arc_pts, r_pad):
    t = np.linspace(0.0, 1.0, arc_pts)
    TH = th0[:, None] + (th1 - th0)[:, None] * t[None, :]
    r_out = np.maximum(r0, r1) + r_pad
    r_in  = np.clip(np.minimum(r0, r1) - r_pad, 0.0, None)
    x_out = r_out[:, None] * np.cos(TH)
    y_out = r_out[:, None] * np.sin(TH)
    x_in  = r_in [:, None] * np.cos(TH[:, ::-1])
    y_in  = r_in [:, None] * np.sin(TH[:, ::-1])
    return np.concatenate([
        np.stack([x_out, y_out], axis=-1),
        np.stack([x_in,  y_in ], axis=-1),
    ], axis=1)


def _cells_union_path(th0, th1, r0, r1, arc_pts=2, r_pad=0.0):
    """
    Единый compound-path, являющийся объединением всех полигонов ячеек.
    Используется как clip_path для pcolormesh — чтобы заливка была ТОЛЬКО
    внутри ячеек (а «шар» и прочие пустоты оставались пустыми, как в рис. 1–2).
    """
    verts = build_verts(th0, th1, r0, r1, arc_pts=arc_pts, r_pad=r_pad)
    all_verts = []
    all_codes = []
    for v in verts:
        n = len(v)
        if n < 3:
            continue
        v_closed = np.vstack([v, v[:1]])
        codes = np.full(n + 1, MplPath.LINETO, dtype=np.uint8)
        codes[0] = MplPath.MOVETO
        codes[-1] = MplPath.CLOSEPOLY
        all_verts.append(v_closed)
        all_codes.append(codes)
    if not all_verts:
        return None
    return MplPath(np.concatenate(all_verts, axis=0),
                   np.concatenate(all_codes, axis=0))


def plot_field(verts, values, name, n_cells, figsize=(10, 10)):
    fig, ax = plt.subplots(figsize=figsize)
    pc = PolyCollection(verts, cmap="turbo",
                        edgecolors="none", antialiased=False)
    pc.set_array(values)
    ax.add_collection(pc)
    ax.autoscale()
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, ls=":", alpha=0.4)
    fig.colorbar(pc, ax=ax, label=name)
    ax.set_title(
        f"{name}  |  ячеек: {n_cells:,}  |  "
        f"диапазон ∈ [{values.min():.2e}, {values.max():.2e}]"
    )
    plt.tight_layout()
    return fig


def plot_velocity(th0, th1, r0, r1, uzp, urp, n_cells,
                  arrow_frac=1.0,
                  max_arrows=MAX_ARROWS,
                  figsize=(10, 10),
                  uzp_is_x=True,
                  interp=True,
                  interp_res=500,
                  show_arrows=True,
                  show_cells_outline=False,
                  y_min=0.0,
                  clip_to_cells=True):
    """
    Модуль скорости |u| = √(uzp²+urp²) в декартовой СК с интерполяцией.

    - Интерполяция |u| на регулярную декартову сетку (griddata, linear→nearest).
    - Заливка ограничена ОБЪЕДИНЕНИЕМ ячеек через clip_path (по аналогии с рис. 1–2:
      вне ячеек — пусто, «шар» не закрашивается).
    - Дополнительно: область ниже y_min маскируется (по умолчанию y_min=0).
    - Стрелки — только при y ≥ y_min.
    """
    if len(uzp) > max_arrows:
        rng = np.random.default_rng(1)
        idx = rng.choice(len(uzp), max_arrows, replace=False)
        th0, th1, r0, r1, uzp, urp = (
            a[idx] for a in (th0, th1, r0, r1, uzp, urp)
        )

    th_mid = 0.5 * (th0 + th1)
    r_mid  = 0.5 * (r0  + r1 )
    x = r_mid * np.cos(th_mid)
    y = r_mid * np.sin(th_mid)

    if uzp_is_x:
        u_x, u_y = uzp, urp
    else:
        u_x, u_y = urp, uzp

    u_mag = np.hypot(u_x, u_y)

    dr  = np.abs(r1 - r0)
    arc = np.abs(th1 - th0) * r_mid
    cell_size = float(np.median(np.maximum(dr, arc)))
    if not np.isfinite(cell_size) or cell_size <= 0:
        cell_size = 1.0

    max_u = float(u_mag.max()) if len(u_mag) else 0.0
    min_u = float(u_mag.min()) if len(u_mag) else 0.0
    if max_u > 0:
        len_scale = arrow_frac * cell_size / max_u
    else:
        len_scale = 0.0
        print("предупреждение: все |u| == 0 — стрелок не будет")

    fig, ax = plt.subplots(figsize=figsize)

    # ---- интерполированное поле модуля скорости ----
    interp_done = False
    if interp:
        if griddata is None:
            print("предупреждение: scipy не найден — интерполяция недоступна")
        elif len(x) >= 4:
            pad = 0.02 * max(float(np.ptp(x)), float(np.ptp(y)), 1e-12)
            xmin, xmax = float(x.min()) - pad, float(x.max()) + pad
            ymin = max(float(y.min()) - pad, float(y_min))
            ymax = float(y.max()) + pad
            n = int(interp_res)
            xi = np.linspace(xmin, xmax, n)
            yi = np.linspace(ymin, ymax, n)
            XI, YI = np.meshgrid(xi, yi)

            try:
                ZI = griddata((x, y), u_mag, (XI, YI), method="linear")
            except Exception:
                ZI = None
            if ZI is None or np.isnan(ZI).any():
                Zn = griddata((x, y), u_mag, (XI, YI), method="nearest")
                if ZI is None:
                    ZI = Zn
                else:
                    ZI = np.where(np.isnan(ZI), Zn, ZI)

            # маска по y_min (по умолчанию полуплоскость y < 0)
            ZI = np.ma.masked_where(YI < y_min, ZI)

            vmin = min_u if max_u > min_u else 0.0
            pm = ax.pcolormesh(XI, YI, ZI, cmap="turbo",
                               shading="auto",
                               vmin=vmin, vmax=max_u if max_u > 0 else None)

            # ---- КЛЮЧЕВОЕ: обрезка заливки по объединению ячеек ----
            if clip_to_cells:
                r_pad_clip = (0.005 * 0.5 * (float(np.mean(r0)) + float(np.mean(r1)))
                              if len(r0) else 0.0)
                clip_path = _cells_union_path(th0, th1, r0, r1,
                                              arc_pts=2, r_pad=r_pad_clip)
                if clip_path is not None:
                    pm.set_clip_path(clip_path, transform=ax.transData)

            fig.colorbar(pm, ax=ax,
                         label="|u| = √(uzp² + urp²)  (интерп., только внутри ячеек)")
            interp_done = True

    # ---- контуры ячеек (опционально, только для y >= y_min) ----
    if show_cells_outline:
        sel = y >= y_min
        if np.any(sel):
            verts_bg = build_verts(th0[sel], th1[sel], r0[sel], r1[sel],
                                   arc_pts=3, r_pad=0.0)
            pc_bg = PolyCollection(verts_bg,
                                   facecolors="none",
                                   edgecolors="#808080",
                                   linewidths=0.15,
                                   antialiased=False)
            ax.add_collection(pc_bg)

    # ---- стрелки направления (только для y >= y_min) ----
    if max_u > 0 and show_arrows:
        keep = y >= y_min
        if np.any(keep):
            xa, ya = x[keep], y[keep]
            uxa, uya = u_x[keep], u_y[keep]
            uma = u_mag[keep]

            if interp_done:
                ax.quiver(
                    xa, ya,
                    uxa * len_scale, uya * len_scale,
                    angles="xy", scale_units="xy", scale=1,
                    pivot="middle",
                    width=0.0015,
                    color="k", alpha=0.55,
                    )
            else:
                vmin = min_u if max_u > min_u else 0.0
                q = ax.quiver(
                    xa, ya,
                    uxa * len_scale, uya * len_scale,
                    uma,
                    cmap="turbo",
                    clim=(vmin, max_u),
                    pivot="middle",
                    scale_units="xy", scale=1,
                    width=0.0025,
                    )
                fig.colorbar(q, ax=ax, label="|u| = √(uzp² + urp²)")

    ax.set_aspect("equal", adjustable="box")
    ax.autoscale()
    ax.grid(True, ls=":", alpha=0.4)
    mode = "интерполяция |u|" if interp_done else "ячейки/стрелки"
    clip_note = ", заливка ⊆ ячеек (шар пуст)" if (interp_done and clip_to_cells) else ""
    ax.set_title(
        f"модуль скорости |u| — Декарт, {mode}{clip_note}  |  "
        f"точек: {len(u_mag):,} из {n_cells:,}, только y ≥ {y_min:g}\n"
        f"|u| ∈ [{u_mag.min():.2e}, {u_mag.max():.2e}]   "
        f"длина стрелок: max ≈ {arrow_frac:g}×{cell_size:.3g}"
    )
    plt.tight_layout()
    return fig


def visualize_cells(th0, th1, r0, r1, uzp, urp,
                    max_segments=MAX_SEGMENTS,
                    arc_pts=4,
                    r_pad=None,
                    arrow_frac=1.0,
                    uzp_is_x=True,
                    interp_res=500,
                    show_arrows=True,
                    y_min=0.0,
                    clip_to_cells=True):
    N = len(uzp)
    print(f"ячеек: {N:,}")

    if N > max_segments:
        rng = np.random.default_rng(0)
        idx = rng.choice(N, max_segments, replace=False)
        th0, th1, r0, r1, uzp, urp = (
            a[idx] for a in (th0, th1, r0, r1, uzp, urp)
        )
        N = max_segments
        print(f"даунсэмплировано до {N:,}")

    if r_pad is None:
        r_pad = 0.005 * 0.5 * (np.mean(r0) + np.mean(r1))
    print(f"r_pad = {r_pad:g}")

    verts = build_verts(th0, th1, r0, r1, arc_pts, r_pad)

    plot_field(verts, uzp, "uzp", N)
    plot_field(verts, urp, "urp", N)
    plot_velocity(th0, th1, r0, r1, uzp, urp, N,
                  arrow_frac=arrow_frac,
                  uzp_is_x=uzp_is_x,
                  interp=True,
                  interp_res=interp_res,
                  show_arrows=show_arrows,
                  y_min=y_min,
                  clip_to_cells=clip_to_cells)
    plt.show()


if __name__ == "__main__":
    text = open(sys.argv[1], encoding="utf-8").read()
    visualize_cells(*parse_fast(text))