#!/usr/bin/env python3
"""
vis_cells.py — быстрая визуализация осесимметричного расчёта (z, r) вокруг шара.

Что изменилось по сравнению со старой версией
  * Никакого griddata и интерполяции на пиксельную сетку: поля рисуются по
    самим ячейкам (pcolormesh на структурной сетке R×θ). Поэтому границы шара —
    гладкие дуги, а не «пиксели», и всё считается за секунды.
  * Стрелки ставятся на РЕДКОЙ РАВНОМЕРНОЙ декартовой сетке (а не по ячейкам,
    которые у стенки очень густые), значение берётся из ячейки, куда попала точка.
  * Окно по умолчанию — окрестность шара (в радиусах шара), нижняя половина
    зеркалится (осесимметрия), поэтому след виден целиком.
  * Опционально линии тока (--stream) — по ним сразу видно срыв и зону
    обратного течения.

Запуск
  python3 vis_cells.py text.txt                    # |u| + стрелки
  python3 vis_cells.py text.txt --stream           # |u| + линии тока
  python3 vis_cells.py text.txt --field p          # давление (uz, ur, p, flow, mag)
  python3 vis_cells.py text.txt --all              # сразу все поля
  python3 vis_cells.py text.txt --xlim -3 10 --ylim 4 --arrows 30 14
  python3 vis_cells.py text.txt --out pic.png --dpi 200 --no-show
"""
import argparse
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

NUM = r"[-+0-9.eE]+"
RE_HEAD = re.compile(r"^\s*(\d+)\s+(\d+)\s*$")


# --------------------------------------------------------------------------- #
#  Чтение файла
# --------------------------------------------------------------------------- #
def parse(path):
    """Читает вывод Grid::print(). Лишние строки (mass_imb и т.п.) игнорируются,
    если один и тот же (i, j) встречается несколько раз — берётся последний."""
    cells = {}
    cur = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = RE_HEAD.match(line)
            if m:
                cur = [int(m.group(1)), int(m.group(2))] + [None] * 8
                continue
            if cur is None:
                continue
            s = line.strip()
            if s.startswith("uzp"):
                cur[2] = float(s.split()[1])
            elif s.startswith("urp"):
                cur[3] = float(s.split()[1])
            elif s.startswith("sum flows"):
                cur[4] = float(s.split()[2])
            elif s.startswith("p "):
                cur[5] = float(s.split()[1])
            elif s.startswith("thetas"):
                a = s.split()
                cur[6], cur[7] = float(a[1]), float(a[2])      # theta_n, theta_s
            elif s.startswith("radii"):
                a = s.split()
                cur[8], cur[9] = float(a[1]), float(a[2])      # R_e, R_w
                if None not in cur:
                    cells[(cur[0], cur[1])] = cur[2:]
                cur = None
    if not cells:
        raise SystemExit("ничего не распарсилось: проверьте формат файла")

    N = max(k[0] for k in cells)
    M = max(k[1] for k in cells)
    A = np.full((N, M, 8), np.nan)
    for (i, j), v in cells.items():
        A[i - 1, j - 1] = v
    uz, ur, sf, p, th_n, th_s, R_e, R_w = (A[..., k] for k in range(8))
    return dict(uz=uz, ur=ur, flow=sf, p=p, th_n=th_n, th_s=th_s, R_e=R_e, R_w=R_w)


def build_geometry(d):
    """Рёбра ячеек (узлы сетки) и центры."""
    th_edges = np.concatenate([d["th_s"][0, :1], d["th_n"][0, :]])     # M+1
    R_edges = np.concatenate([d["R_w"][:1, 0], d["R_e"][:, 0]])        # N+1
    th_c = 0.5 * (th_edges[:-1] + th_edges[1:])
    R_c = 0.5 * (R_edges[:-1] + R_edges[1:])
    X = R_edges[:, None] * np.cos(th_edges[None, :])                   # (N+1, M+1)
    Y = R_edges[:, None] * np.sin(th_edges[None, :])
    return th_edges, R_edges, th_c, R_c, X, Y


# --------------------------------------------------------------------------- #
#  Рисование
# --------------------------------------------------------------------------- #
def field_values(d, name, U):
    if name == "mag":
        return np.hypot(d["uz"], d["ur"]), r"$|u|$, м/с", "turbo", (0.0, 1.6 * U)
    if name == "uz":
        return d["uz"], r"$u_z$, м/с", "RdBu_r", (-1.6 * U, 1.6 * U)
    if name == "ur":
        return d["ur"], r"$u_r$, м/с", "RdBu_r", (-1.0 * U, 1.0 * U)
    if name == "p":
        q = 0.5 * U * U
        return d["p"], r"$p$ (кинематическое), Па/ρ", "RdBu_r", (-1.5 * q, 1.5 * q)
    if name == "flow":
        v = d["flow"]
        m = np.nanmax(np.abs(v)) or 1.0
        return v, r"$F_e-F_w+F_n-F_s$", "RdBu_r", (-m, m)
    raise SystemExit(f"неизвестное поле {name}")


def sample_cells(d, th_edges, R_edges, x, y):
    """Для точек (x, y) возвращает значения (uz, ur) ячейки, в которую они попали."""
    R = np.hypot(x, y)
    th = np.arctan2(y, x)                      # 0..π для y >= 0
    i = np.searchsorted(R_edges, R) - 1
    j = np.searchsorted(th_edges, th) - 1
    ok = (i >= 0) & (i < d["uz"].shape[0]) & (j >= 0) & (j < d["uz"].shape[1])
    i = np.clip(i, 0, d["uz"].shape[0] - 1)
    j = np.clip(j, 0, d["uz"].shape[1] - 1)
    return d["uz"][i, j], d["ur"][i, j], ok


# --- draw() ---
def draw(d, geo, name, args, U, Rs):
    th_edges, R_edges, th_c, R_c, X, Y = geo
    C, label, cmap, (vmin, vmax) = field_values(d, name, U)

    # --- границы всей расчётной области ---
    # углы сетки — это крайние ячейки; берём по всем узлам (X, Y),
    # плюс небольшой запас, чтобы граница не упиралась в рамку
    X_min, X_max = float(np.nanmin(X)), float(np.nanmax(X))
    Y_max = float(np.nanmax(np.abs(Y)))

    if args.xlim is None:
        span = X_max - X_min
        pad = 0.01 * span if span > 0 else 0.0
        x0, x1 = X_min - pad, X_max + pad
    else:
        x0, x1 = args.xlim[0] * Rs, args.xlim[1] * Rs

    if args.ylim is None:
        yl = Y_max * 1.01
    else:
        yl = args.ylim * Rs
    # ----------------------------------------------------------

    fig, ax = plt.subplots(figsize=(args.fig_w, args.fig_h))
    kw = dict(cmap=cmap, vmin=vmin, vmax=vmax, shading="flat",
              edgecolors="face", linewidth=0, antialiased=False, rasterized=True)
    pm = ax.pcolormesh(X, Y, C, **kw)
    ax.pcolormesh(X, -Y, C, **kw)
    cb = fig.colorbar(pm, ax=ax, label=label, shrink=0.8)

    ax.add_patch(Circle((0, 0), Rs, facecolor="#444444", edgecolor="k",
                        linewidth=0.8, zorder=5))

    ax.set_xlim(x0, x1)
    ax.set_ylim(-yl, yl)
    ax.set_aspect("equal", adjustable="box")

    # дальше — без изменений (стрелки/линии тока, подписи, тики)


    # ---- стрелки / линии тока на равномерной сетке ----
    if name == "mag" or args.stream_on_all:
        nx, ny = args.arrows
        xs = np.linspace(x0, x1, nx)
        ys = np.linspace(-yl, yl, ny)
        XS, YS = np.meshgrid(xs, ys)
        ya = np.abs(YS)
        u, v, ok = sample_cells(d, th_edges, R_edges, XS, ya)
        v = np.where(YS < 0, -v, v)            # зеркальная половина
        ok &= np.hypot(XS, YS) > Rs * 1.002     # не рисуем внутри шара

        if args.stream:
            nxs, nys = 600, int(600 * (2 * yl) / (x1 - x0)) or 300
            xg = np.linspace(x0, x1, nxs)
            yg = np.linspace(-yl, yl, nys)
            XG, YG = np.meshgrid(xg, yg)
            ug, vg, okg = sample_cells(d, th_edges, R_edges, XG, np.abs(YG))
            vg = np.where(YG < 0, -vg, vg)
            okg &= np.hypot(XG, YG) > Rs * 1.002
            ug = np.ma.masked_where(~okg, ug)
            vg = np.ma.masked_where(~okg, vg)
            ax.streamplot(xg, yg, ug, vg, density=args.density, color="k",
                          linewidth=0.6, arrowsize=0.8, zorder=6)
        else:
            dx = (x1 - x0) / max(nx - 1, 1)
            if args.normalize:
                mag = np.hypot(u, v)
                mag[mag == 0] = 1.0
                u, v = u / mag, v / mag
                scale = 1.0 / (0.8 * dx)
            else:
                scale = U / (0.8 * dx)
            ax.quiver(XS[ok], YS[ok], u[ok], v[ok],
                      angles="xy", scale_units="xy", scale=scale,
                      pivot="mid", width=0.0018, color="k", alpha=0.75, zorder=6)

    ax.set_xlabel("z / R_шара" if args.units else "z, м")
    ax.set_ylabel("r / R_шара" if args.units else "r, м")
    if args.units:
        ticks_x = ax.get_xticks()
        ax.set_xticklabels([f"{t / Rs:g}" for t in ticks_x])
        ticks_y = ax.get_yticks()
        ax.set_yticklabels([f"{t / Rs:g}" for t in ticks_y])
    ax.set_title(f"{label}   (U≈{U:.4g} м/с, ячеек {d['uz'].size:,})")
    ax.grid(True, ls=":", alpha=0.3)
    fig.tight_layout()
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--field", default="mag", choices=["mag", "uz", "ur", "p", "flow"])
    ap.add_argument("--all", action="store_true", help="нарисовать все поля")
    ap.add_argument("--xlim", nargs=2, type=float, default=None,
                    help="окно по z в радиусах шара (по умолчанию — вся область)")
    ap.add_argument("--ylim", type=float, default=None,
                    help="окно по r (±) в радиусах шара (по умолчанию — вся область)")
    ap.add_argument("--arrows", nargs=2, type=int, default=[34, 16],
                    help="сколько стрелок по z и по r (редкая равномерная сетка)")
    ap.add_argument("--normalize", action="store_true", help="стрелки одной длины")
    ap.add_argument("--stream", action="store_true", help="линии тока вместо стрелок")
    ap.add_argument("--density", type=float, default=1.6, help="плотность линий тока")
    ap.add_argument("--stream-on-all", action="store_true",
                    help="рисовать стрелки/линии тока и на остальных полях")
    ap.add_argument("--units", action="store_true", help="оси в радиусах шара")
    ap.add_argument("--fig-w", type=float, default=13)
    ap.add_argument("--fig-h", type=float, default=7)
    ap.add_argument("--out", default=None, help="сохранить в файл (для --all: суффикс _<поле>)")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()

    d = parse(args.file)
    geo = build_geometry(d)
    Rs = float(geo[1][0])                                  # радиус шара = внутренняя граница
    U = float(np.nanmedian(d["uz"][-3:, :]))
    if not np.isfinite(U) or U <= 0:
        U = float(np.nanmax(np.hypot(d["uz"], d["ur"])))
    print(f"ячеек: {d['uz'].size:,}  R_шара={Rs:g}  U≈{U:g}")

    names = ["mag", "uz", "ur", "p", "flow"] if args.all else [args.field]
    for n in names:
        fig = draw(d, geo, n, args, U, Rs)
        if args.out:
            path = args.out if not args.all else re.sub(r"(\.\w+)?$", f"_{n}\\1", args.out, count=1)
            fig.savefig(path, dpi=args.dpi)
            print("сохранено:", path)
    if not args.no_show:
        plt.show()

if __name__ == "__main__":
    main()