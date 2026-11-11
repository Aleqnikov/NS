#!/usr/bin/env python3
# cda_sphere.py
#
# CDA для осесимметричного обтекания сферы.
# CDA = F_drag / (0.5 * rho * U^2)
# C_D = CDA / A_ref,  A_ref = pi * R^2
#
# Константы по умолчанию совпадают с C++ ModelConfig:
#   rho = 1.225 кг/м^3, nu = 1.46e-5 м^2/с, nut = 4.38e-5 м^2/с, speed = 100 м/с
#   mu = rho * nu ≈ 1.789e-5 Па*с
#
# Использование:
#   python cda_sphere.py "1 1.txt"
#   python cda_sphere.py "1 1.txt" --U 100 --viscous
#   python cda_sphere.py "1 1.txt" --sa --viscous          # учитывать nut в вязком трении
#   python cda_sphere.py "1 1.txt" --p-sign -1
#
# Важно:
# - p_sign=+1 => Fp = -sum(p*cos(theta)*dA). Если CDA<0, запусти с --p-sign -1.
# - Вязкая часть --viscous чувствительна к первой ячейке.
#   Пока нет SA, используй ламинарный mu (по умолчанию). --sa включает nut.

import re
import math
import argparse

FLOAT = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
BLOCK_RE = re.compile(r'^\s*(\d+)\s+(\d+)\s*$')

KEY_RES = {
    'uzp':       re.compile(r'^\s*uzp\s+(' + FLOAT + r')'),
    'urp':       re.compile(r'^\s*urp\s+(' + FLOAT + r')'),
    'sum flows': re.compile(r'^\s*sum\s+flows\s+(' + FLOAT + r')'),
    'p':         re.compile(r'^\s*p\s+(' + FLOAT + r')'),
    'thetas':    re.compile(r'^\s*thetas\s+(' + FLOAT + r')\s+(' + FLOAT + r')'),
    'radii':     re.compile(r'^\s*radii\s+(' + FLOAT + r')\s+(' + FLOAT + r')'),
}


def parse_file(path):
    cells = {}
    cur = None
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            m = BLOCK_RE.match(line)
            if m:
                i, j = map(int, m.groups())
                cur = {'i': i, 'j': j}
                cells[(i, j)] = cur
                continue
            if cur is None:
                continue
            for key, rx in KEY_RES.items():
                m = rx.match(line)
                if m:
                    if key in ('thetas', 'radii'):
                        cur[key] = [float(m.group(1)), float(m.group(2))]
                    else:
                        cur[key] = float(m.group(1))
                    break
    return cells


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('file')

    # === Параметры, совпадающие с C++ ModelConfig ===
    ap.add_argument('--rho', type=float, default=1.225,
                    help='плотность, кг/м^3  (ModelConfig::rho)')
    ap.add_argument('--nu',  type=float, default=1.46e-5,
                    help='кинематическая вязкость, м^2/с  (ModelConfig::nu)')
    ap.add_argument('--nut', type=float, default=4.38e-5,
                    help='турбулентная вязкость, м^2/с  (ModelConfig::nut)')
    ap.add_argument('--U',   type=float, default=100.0,
                    help='скорость набегания, м/с  (ModelConfig::speed)')
    ap.add_argument('--R',   type=float, default=None,
                    help='радиус сферы, м; по умолчанию min(r) на i=1')
    ap.add_argument('--sa', action='store_true',
                    help='учитывать nut в вязком трении: mu_eff = rho*(nu+nut)')

    # === Управление знаком давления и вязкой частью ===
    ap.add_argument('--p-sign', type=float, default=1.0,
                    help='Fp = -p_sign * sum(p*cos(theta)*dA). Если CDA<0, попробуй -1')
    ap.add_argument('--viscous', action='store_true',
                    help='добавить грубую оценку вязкого трения по первой ячейке')
    args = ap.parse_args()

    cells = parse_file(args.file)
    if not cells:
        raise SystemExit('Не удалось распарсить файл')

    i_min = min(i for (i, j) in cells.keys())
    rows = sorted(
        [c for (i, j), c in cells.items() if i == i_min],
        key=lambda c: c['j']
    )

    if args.R is None:
        R = min(min(c['radii']) for c in rows)
    else:
        R = args.R

    rho = args.rho
    nu  = args.nu
    nut = args.nut
    mu_laminar = rho * nu
    mu_eff     = rho * (nu + nut) if args.sa else mu_laminar

    q = 0.5 * rho * args.U * args.U
    A_ref = math.pi * R * R
    Re_d  = args.U * (2.0 * R) / nu

    # SA: c_w1
    c_b1, c_b2, k, sigma = 0.1355, 0.622, 0.41, 2.0 / 3.0
    c_w1 = c_b1 / (k * k) + (1.0 + c_b2) / sigma

    Fp_plus = 0.0
    Fp_minus = 0.0
    Fv = 0.0
    skipped_visc = 0

    for c in rows:
        th_a, th_b = sorted(c['thetas'])
        th = 0.5 * (th_a + th_b)

        # dA = 2*pi*R^2 * sin(theta) dtheta
        dA = 2.0 * math.pi * R * R * (math.cos(th_a) - math.cos(th_b))

        # Давление: Fz = - integral p * n_z dA, n_z = cos(theta)
        Fp_plus  += -c['p'] * math.cos(th) * dA
        Fp_minus +=  c['p'] * math.cos(th) * dA

        if args.viscous:
            r_cell = 0.5 * (c['radii'][0] + c['radii'][1])
            h = r_cell - R
            if h <= 1e-15:
                skipped_visc += 1
                continue

            # uzp = u_z, urp = u_rho (цилиндрическая радиальная)
            # e_theta = cos(theta)*e_rho - sin(theta)*e_z
            u_theta = c['urp'] * math.cos(th) - c['uzp'] * math.sin(th)

            # tau = mu_eff * du_theta/dn,  du/dn ~ u_theta / h
            tau = mu_eff * u_theta / h

            # e_theta . e_z = -sin(theta)
            Fv += tau * (-math.sin(th)) * dA

    Fp = -args.p_sign * Fp_plus if False else args.p_sign * 0.0
    # аккуратный выбор знака:
    if args.p_sign >= 0:
        Fp = Fp_plus
    else:
        Fp = Fp_minus

    Ftot = Fp + (Fv if args.viscous else 0.0)
    CDAp = Fp / q
    CDp  = CDAp / A_ref

    print(f'R              = {R:.12g} m')
    print(f'U              = {args.U:.12g} m/s')
    print(f'rho            = {rho:.12g} kg/m^3')
    print(f'nu             = {nu:.12g} m^2/s')
    print(f'nut            = {nut:.12g} m^2/s  ({ "used" if args.sa else "ignored" })')
    print(f'mu_laminar     = {mu_laminar:.12g} Pa*s')
    print(f'mu_eff         = {mu_eff:.12g} Pa*s')
    print(f'Re_d           = {Re_d:.12g}')
    print(f'q = 0.5 rho U^2= {q:.12g} Pa')
    print(f'A_ref = pi R^2 = {A_ref:.12g} m^2')
    print(f'c_w1 (SA)      = {c_w1:.12g}')
    print()

    print(f'Fp(sign=+1) = {Fp_plus:.12g} N,  '
          f'CDA = {Fp_plus/q:.12g} m^2,  CD = {Fp_plus/q/A_ref:.12g}')
    print(f'Fp(sign=-1) = {Fp_minus:.12g} N,  '
          f'CDA = {Fp_minus/q:.12g} m^2,  CD = {Fp_minus/q/A_ref:.12g}')
    print(f'Chosen p_sign={args.p_sign:+.0f}:  '
          f'Fp = {Fp:.12g} N,  CDAp = {CDAp:.12g} m^2,  CDp = {CDp:.12g}')

    if args.viscous:
        CDAv = Fv / q
        CDv  = CDAv / A_ref
        CDA_tot = Ftot / q
        CD_tot  = CDA_tot / A_ref
        print()
        print(f'Viscous (first cell):  Fv = {Fv:.12g} N,  '
              f'CDAv = {CDAv:.12g} m^2,  CDv = {CDv:.12g}')
        if skipped_visc:
            print(f'  skipped viscous cells: {skipped_visc}')
        print(f'Total:  F = {Ftot:.12g} N,  '
              f'CDA = {CDA_tot:.12g} m^2,  CD = {CD_tot:.12g}')
    else:
        print()
        print('Viscous part not included. Добавь --viscous для грубой оценки.')

    print()
    print('Если CDA отрицательный — поменяй знак давления: --p-sign -1')


if __name__ == '__main__':
    main()