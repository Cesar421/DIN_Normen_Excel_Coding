# -*- coding: utf-8 -*-
"""
Grafiken der Bemessung (matplotlib).

Alle Abbildungen sind fuer die Einbettung in die Oberflaeche (Tkinter) und in
das Jupyter-Notebook vorgesehen. Titel bzw. Beschriftungen nennen die
zugehoerige Normstelle.

Farbpalette: blau / orange / aqua (kategorial, farbfehlsichtigkeitsgeprueft);
Zustaende OK / NEIN zusaetzlich mit Textbeschriftung, nicht nur ueber die Farbe.
"""

import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Polygon, FancyArrow

# --- Palette ---------------------------------------------------------------
C1 = "#2a78d6"   # blau    - Reihe 1
C2 = "#eb6834"   # orange  - Reihe 2
C3 = "#1baf7a"   # aqua    - Reihe 3
C4 = "#4a3aa7"   # violett - Reihe 4
OK = "#0ca30c"   # Zustand: erfuellt
NOK = "#d03b3b"  # Zustand: nicht erfuellt
SURF = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUT = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BETON = "#e6e4dd"


def _fig(w=9.0, h=5.0, n=1, m=1, **kw):
    fig, ax = plt.subplots(n, m, figsize=(w, h), facecolor=SURF, **kw)
    fig.patch.set_facecolor(SURF)
    return fig, ax


def _stil(ax, grid="y"):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=MUT, labelsize=8.5, length=3, width=0.8)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_color(INK2)
    if grid:
        ax.grid(True, axis=grid, color=GRID, lw=0.7, zorder=0)
        ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.title.set_color(INK)
    return ax


def _unter(ax, text):
    """Untertitel mit der Normstelle; der Titel braucht pad >= 20."""
    ax.annotate(text, xy=(0, 1.0), xycoords="axes fraction",
                xytext=(0, 4), textcoords="offset points",
                fontsize=7.8, color=MUT, ha="left", va="bottom")


# ===========================================================================
# 1. BEWEHRTER QUERSCHNITT
# ===========================================================================
def bild_querschnitt(b, ax=None):
    """Querschnitt mit der gewaehlten Bewehrung (EC2 8.2 / 9.2 / 4.4)."""
    from .baustoffe import betondeckung
    e = b["eingabe"]
    qs = b["qs_feld"]
    bd = betondeckung(e.expositionsklasse, e.phi_laengs, e.phi_buegel, e.d_g,
                      beton=b["beton"])
    cw = bd["c_nom_w"]
    bw, h = e.b, e.h

    if ax is None:
        fig, ax = _fig(6.4, 6.8)
    else:
        fig = ax.figure
    ax.set_facecolor(SURF)

    if e.querschnittstyp == "plattenbalken":
        be, hf = e.b_eff, e.hf
        x0 = (be - bw) / 2.0
        kontur = [(0, 0), (be, 0), (be, hf), (x0 + bw, hf), (x0 + bw, h),
                  (x0, h), (x0, hf), (0, hf)]
        ax.add_patch(Polygon(kontur, closed=True, fc=BETON, ec=INK2, lw=1.4, zorder=1))
        xL, breite_ges = x0, be
    else:
        ax.add_patch(Rectangle((0, 0), bw, h, fc=BETON, ec=INK2, lw=1.4, zorder=1))
        xL, breite_ges = 0.0, bw

    # Buegel
    pb = e.phi_buegel
    ax.add_patch(Rectangle((xL + cw, cw), bw - 2 * cw, h - 2 * cw, fc="none",
                           ec=C2, lw=2.0, zorder=3, joinstyle="round"))
    ax.add_patch(Rectangle((xL + cw + pb, cw + pb), bw - 2 * (cw + pb),
                           h - 2 * (cw + pb), fc="none", ec=C2, lw=2.0,
                           alpha=0.35, zorder=3))

    def _staebe(n, phi, y, farbe=C1):
        if n <= 0:
            return
        xi = xL + cw + pb + phi / 2.0
        xf = xL + bw - cw - pb - phi / 2.0
        xs = [0.5 * (xi + xf)] if n == 1 else np.linspace(xi, xf, n)
        for xb in xs:
            ax.add_patch(Circle((xb, y), phi / 2.0, fc=farbe, ec=SURF, lw=0.8,
                                zorder=5))

    n_u, phi_u = b["n_unten"], e.phi_laengs
    y_u = h - qs.d
    lagen = b.get("lagen_unten", 1)
    if lagen > 1:
        n1 = int(math.ceil(n_u / 2.0))
        sv = max(phi_u, e.d_g + 5.0, 20.0)
        _staebe(n1, phi_u, cw + pb + phi_u / 2.0)
        _staebe(n_u - n1, phi_u, cw + pb + phi_u / 2.0 + phi_u + sv)
    else:
        _staebe(n_u, phi_u, y_u)
    _staebe(b["n_oben"], e.phi_laengs_oben, h - cw - pb - e.phi_laengs_oben / 2.0)

    # Torsionslaengsbewehrung an den Seiten
    tor = b.get("torsion")
    if tor is not None and abs(e.T_Ed) > 1e-9 and tor.erforderlich:
        n_seite = max(1, (tor.n_laengsstaebe - 4) // 2)
        phi_t = 12.0
        for ys in np.linspace(h * 0.30, h * 0.70, n_seite):
            for xs in (xL + cw + pb + phi_t / 2.0,
                       xL + bw - cw - pb - phi_t / 2.0):
                ax.add_patch(Circle((xs, ys), phi_t / 2.0, fc=C4, ec=SURF,
                                    lw=0.8, zorder=5))
        ax.text(xL + bw + 0.08 * breite_ges, h * 0.5,
                "Torsions-\nlaengsstaebe", color=C4, fontsize=8.5,
                ha="left", va="center")

    # Bemassung
    ax.annotate("", xy=(xL, -0.10 * h), xytext=(xL + bw, -0.10 * h),
                arrowprops=dict(arrowstyle="<->", color=MUT, lw=0.9))
    ax.text(xL + bw / 2, -0.135 * h, "b$_w$ = %.0f mm" % bw, ha="center",
            va="top", color=INK2, fontsize=9)
    ax.annotate("", xy=(-0.13 * breite_ges, 0), xytext=(-0.13 * breite_ges, h),
                arrowprops=dict(arrowstyle="<->", color=MUT, lw=0.9))
    ax.text(-0.155 * breite_ges, h / 2, "h = %.0f mm" % h, ha="right",
            va="center", color=INK2, fontsize=9, rotation=90)
    ax.plot([xL + bw + 0.06 * breite_ges] * 2, [h, h - qs.d], color=C1, lw=1.2)
    ax.plot([xL + bw, xL + bw + 0.08 * breite_ges], [h - qs.d, h - qs.d],
            color=C1, lw=0.8, ls=":")
    ax.text(xL + bw + 0.08 * breite_ges, h - qs.d / 2, "d = %.0f mm" % qs.d,
            ha="left", va="center", color=C1, fontsize=9)

    if e.querschnittstyp == "plattenbalken":
        ax.annotate("", xy=(0, h + 0.07 * h), xytext=(e.b_eff, h + 0.07 * h),
                    arrowprops=dict(arrowstyle="<->", color=MUT, lw=0.9))
        ax.text(e.b_eff / 2, h + 0.09 * h, "b$_{eff}$ = %.0f mm" % e.b_eff,
                ha="center", va="bottom", color=INK2, fontsize=9)
        ax.text(e.b_eff * 0.02, e.hf / 2, "h$_f$ = %.0f" % e.hf, ha="left",
                va="center", color=INK2, fontsize=8)

    ax.text(xL + bw / 2, y_u - 0.055 * h,
            "%d $\\varnothing$%.0f  (%.0f mm$^2$)" % (n_u, phi_u, b["As1_vorh"]),
            ha="center", va="top", color=C1, fontsize=9.5, weight="bold")
    ax.text(xL + bw / 2, h - cw - pb + 0.035 * h,
            "%d $\\varnothing$%.0f  (%.0f mm$^2$)"
            % (b["n_oben"], e.phi_laengs_oben, b["As_oben_vorh"]),
            ha="center", va="bottom", color=C1, fontsize=9.5, weight="bold")
    z0 = b["querkraft"].bereiche[0]
    geschl = " (geschlossen)" if abs(e.T_Ed) > 1e-9 else ""
    ax.text(xL + bw + 0.08 * breite_ges, cw + 0.02 * h,
            "Buegel $\\varnothing$%.0f/%.0f mm\n(%d Schenkel)%s"
            % (pb, z0["s"], e.n_schenkel, geschl), ha="left", va="bottom",
            color=C2, fontsize=9)
    ax.annotate("c$_{nom}$ = %.0f mm" % cw, xy=(xL + bw - cw / 2, cw * 1.5),
                xytext=(xL + bw + 0.10 * breite_ges, -0.135 * h),
                arrowprops=dict(arrowstyle="->", color=MUT, lw=0.8, shrinkA=0,
                                shrinkB=2),
                ha="left", va="center", color=MUT, fontsize=8.5)

    ax.set_xlim(-0.30 * breite_ges, 1.48 * breite_ges)
    ax.set_ylim(-0.22 * h, 1.18 * h)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Bewehrter Querschnitt", fontsize=11.5, color=INK, loc="left",
                 pad=22)
    _unter(ax, "Betondeckung: EC2 4.4.1, Gl. (4.1)/(4.2) + NA Tab. 4.4DE | "
               "Stababstaende: EC2 8.2 (2)")
    fig.tight_layout()
    return fig


# ===========================================================================
# 2. BAUSTOFFDIAGRAMME
# ===========================================================================
def bild_baustoffe(b):
    """Bemessungs-Spannungs-Dehnungs-Linien (EC2 3.1.7 und 3.2.7)."""
    C, S = b["beton"], b["stahl"]
    fig, (a1, a2) = _fig(10.0, 4.4, 1, 2)

    _stil(a1)
    e = np.linspace(0, C.eps_cu2, 300)
    a1.plot(e, [C.sigma_c(v) for v in e], color=C1, lw=2.4)
    a1.fill_between(e, 0, [C.sigma_c(v) for v in e], color=C1, alpha=0.13)
    a1.axvline(C.eps_c2, color=MUT, lw=0.9, ls=":")
    a1.axhline(C.fcd, color=MUT, lw=0.9, ls=":")
    a1.annotate("$f_{cd}$ = %.2f N/mm$^2$" % C.fcd, xy=(C.eps_c2 * 0.15, C.fcd),
                xytext=(0, 6), textcoords="offset points", color=INK2, fontsize=9)
    a1.annotate("$\\varepsilon_{c2}$ = %.1f" % C.eps_c2, xy=(C.eps_c2, 0),
                xytext=(3, 8), textcoords="offset points", color=MUT, fontsize=8.5)
    a1.annotate("$\\varepsilon_{cu2}$ = %.1f" % C.eps_cu2, xy=(C.eps_cu2, 0),
                xytext=(-3, 8), textcoords="offset points", color=MUT,
                fontsize=8.5, ha="right")
    a1.set_xlabel("$\\varepsilon_c$ [‰]")
    a1.set_ylabel("$\\sigma_c$ [N/mm$^2$]")
    a1.set_title("Beton %s - Parabel-Rechteck" % C.klasse, fontsize=10,
                 color=INK, loc="left", pad=20)
    _unter(a1, "EC2 3.1.7 (1), Bild 3.3, Gl. (3.17)/(3.18)")

    _stil(a2)
    es = np.linspace(0, S.eps_ud, 300)
    a2.plot(es, [S.sigma_s(v) for v in es], color=C2, lw=2.4)
    a2.axhline(S.fyd, color=MUT, lw=0.9, ls=":")
    a2.annotate("$f_{yd}$ = %.1f N/mm$^2$" % S.fyd, xy=(S.eps_ud * 0.35, S.fyd),
                xytext=(0, 6), textcoords="offset points", color=INK2, fontsize=9)
    a2.annotate("$\\varepsilon_{yd}$ = %.2f" % S.eps_yd, xy=(S.eps_yd, 0),
                xytext=(4, 8), textcoords="offset points", color=MUT, fontsize=8.5)
    a2.annotate("$\\varepsilon_{ud}$ = %.0f" % S.eps_ud, xy=(S.eps_ud, 0),
                xytext=(-3, 8), textcoords="offset points", color=MUT,
                fontsize=8.5, ha="right")
    a2.set_xlabel("$\\varepsilon_s$ [‰]")
    a2.set_ylabel("$\\sigma_s$ [N/mm$^2$]")
    a2.set_title("Betonstahl %s - bilinear" % S.sorte, fontsize=10, color=INK,
                 loc="left", pad=20)
    _unter(a2, "EC2 3.2.7 (2)b), Bild 3.8 | $\\varepsilon_{ud}$: NA NDP zu 3.2.7 (2)")
    fig.tight_layout()
    return fig


# ===========================================================================
# 3. SCHNITTGROESSEN
# ===========================================================================
def bild_schnittgroessen(b):
    """Einhuellende von M und V im GZT + Lastbild (EC0 6.10 / EC2 5.1.3)."""
    e = b["eingabe"]
    einh, traeger = b["einhuellende"], b["traeger"]
    x = einh["x"]

    fig, axs = _fig(10.5, 8.2, 3, 1, sharex=True,
                    gridspec_kw={"height_ratios": [0.6, 1.0, 1.0]})
    a0, a1, a2 = axs

    a0.set_facecolor(SURF)
    a0.plot([0, e.L], [0, 0], color=INK, lw=3, solid_capstyle="butt")
    for ap in traeger.auflager:
        if ap.typ == "gelenkig":
            a0.plot([ap.x], [0], marker="^", ms=12, color=INK2, clip_on=False)
        else:
            a0.plot([ap.x, ap.x], [-0.35, 0.35], color=INK2, lw=4)
    for xa in np.linspace(0.02 * e.L, 0.98 * e.L, 22):
        a0.annotate("", xy=(xa, 0.05), xytext=(xa, 0.75),
                    arrowprops=dict(arrowstyle="->", color=C1, lw=1.0, alpha=0.75))
    a0.plot([0, e.L], [0.75, 0.75], color=C1, lw=1.6)
    q_ed = e.gamma_G * (e.g_k + (25.0 * e.b * e.h / 1e6 if e.eigengewicht else 0.0)) \
        + e.gamma_Q * e.q_k
    a0.text(e.L / 2, 0.95,
            "$q_{Ed}$ = %.2f kN/m   (1,35 $G_k$ + 1,50 $Q_k$)" % q_ed,
            ha="center", color=C1, fontsize=9.5)
    for el in e.einzellasten:
        a0.annotate("", xy=(el[0], 0.05), xytext=(el[0], 1.15),
                    arrowprops=dict(arrowstyle="->", color=C2, lw=2.0))
        a0.text(el[0], 1.2, "%.0f kN" % el[1], ha="center", color=C2, fontsize=9)
    if abs(e.T_Ed) > 1e-9:
        a0.text(e.L * 0.5, -0.45, "$T_{Ed}$ = %.1f kNm" % e.T_Ed, ha="center",
                color=C4, fontsize=9.5, weight="bold")
    a0.set_ylim(-0.7, 1.7)
    a0.axis("off")
    a0.set_title("Statisches System und Bemessungslasten", fontsize=10.5,
                 color=INK, loc="left")

    _stil(a1)
    a1.fill_between(x, 0, einh["Mmax"], color=C1, alpha=0.20, lw=0)
    a1.fill_between(x, 0, einh["Mmin"], color=C1, alpha=0.20, lw=0)
    a1.plot(x, einh["Mmax"], color=C1, lw=2.0)
    a1.plot(x, einh["Mmin"], color=C1, lw=2.0)
    a1.axhline(0, color=AXIS, lw=0.9)
    i = int(np.argmax(einh["Mmax"]))
    a1.plot([x[i]], [einh["Mmax"][i]], "o", color=C1, ms=8)
    a1.annotate("$M_{Ed,max}$ = %.1f kNm" % einh["Mmax"][i],
                xy=(x[i], einh["Mmax"][i]), xytext=(0, 10),
                textcoords="offset points", ha="center", color=C1,
                fontsize=9.5, weight="bold")
    j = int(np.argmin(einh["Mmin"]))
    if einh["Mmin"][j] < -1e-6:
        a1.plot([x[j]], [einh["Mmin"][j]], "o", color=C1, ms=8)
        a1.annotate("$M_{Ed,min}$ = %.1f kNm" % einh["Mmin"][j],
                    xy=(x[j], einh["Mmin"][j]), xytext=(0, -14),
                    textcoords="offset points", ha="center", color=C1,
                    fontsize=9.5, weight="bold")
    a1.invert_yaxis()
    a1.set_ylabel("$M_{Ed}$ [kNm]")
    a1.set_title("Einhuellende der Biegemomente (Zug unten nach unten aufgetragen)",
                 fontsize=10.5, color=INK, loc="left", pad=20)
    _unter(a1, "EC0 6.4.3.2 Gl. (6.10) | feldweise Laststellung EC2 5.1.3")

    _stil(a2)
    a2.fill_between(x, 0, einh["Vmax"], color=C2, alpha=0.20, lw=0)
    a2.fill_between(x, 0, einh["Vmin"], color=C2, alpha=0.20, lw=0)
    a2.plot(x, einh["Vmax"], color=C2, lw=2.0)
    a2.plot(x, einh["Vmin"], color=C2, lw=2.0)
    a2.axhline(0, color=AXIS, lw=0.9)
    k = int(np.argmax(np.abs(einh["Vmax"])))
    a2.annotate("$V_{Ed,max}$ = %.1f kN" % abs(einh["Vmax"][k]),
                xy=(x[k], einh["Vmax"][k]), xytext=(8, 0),
                textcoords="offset points", ha="left", va="center",
                color=C2, fontsize=9.5, weight="bold")
    a2.set_ylabel("$V_{Ed}$ [kN]")
    a2.set_xlabel("x [m]")
    a2.set_title("Einhuellende der Querkraefte", fontsize=10.5, color=INK,
                 loc="left")
    fig.tight_layout()
    return fig


# ===========================================================================
# 4. DEHNUNGEN UND SPANNUNGEN
# ===========================================================================
def bild_dehnungen(b):
    """Dehnungs- und Spannungsverteilung im Bruchzustand (EC2 6.1 / 3.1.7)."""
    e, C, S = b["eingabe"], b["beton"], b["stahl"]
    fx = b.get("biegung_feld")
    qs = b["qs_feld"]
    if fx is None:
        fig, a = _fig(9, 3)
        a.axis("off")
        a.text(0.5, 0.5, "Kein Feldmoment vorhanden", ha="center", va="center",
               color=MUT)
        return fig

    h, d, bw = e.h, qs.d, e.b
    x, ec, es = fx.x, fx.eps_c, fx.eps_s1

    fig, (a0, a1, a2) = _fig(11.0, 5.4, 1, 3,
                             gridspec_kw={"width_ratios": [1.0, 1.1, 1.1]})

    a0.set_facecolor(SURF)
    a0.add_patch(Rectangle((0, 0), bw, h, fc=BETON, ec=INK2, lw=1.2))
    a0.add_patch(Rectangle((0, h - x), bw, x, fc=C1, alpha=0.16, ec="none"))
    a0.plot([0, bw], [h - x, h - x], color=C1, lw=1.6)
    a0.text(bw * 0.5, h - x + 0.02 * h, "Nulllinie  x = %.0f mm" % x,
            ha="center", va="bottom", color=C1, fontsize=8.5)
    a0.plot(np.linspace(bw * 0.15, bw * 0.85, max(b["n_unten"], 1)),
            [h - d] * max(b["n_unten"], 1), "o", color=C1, ms=6)
    a0.set_xlim(-0.3 * bw, 1.3 * bw)
    a0.set_ylim(-0.08 * h, 1.08 * h)
    a0.set_aspect("equal")
    a0.axis("off")
    a0.set_title("Querschnitt", fontsize=10, color=INK, loc="left")

    _stil(a1, grid="x")
    a1.plot([0, 0], [0, h], color=AXIS, lw=0.9)
    a1.plot([-ec, es], [h, h - d], color=C1, lw=2.2, solid_capstyle="round")
    a1.fill_betweenx([h - x, h], [0.0, -ec], [0.0, 0.0], color=C1, alpha=0.12)
    a1.fill_betweenx([h - d, h - x], [es, 0.0], [0.0, 0.0], color=C1, alpha=0.12)
    a1.plot([-ec], [h], "o", color=C1, ms=7)
    a1.plot([es], [h - d], "o", color=C1, ms=7)
    a1.annotate("$\\varepsilon_c$ = %.2f ‰" % ec, xy=(-ec, h),
                xytext=(6, -4), textcoords="offset points", fontsize=9,
                color=C1, ha="left", va="top")
    a1.annotate("$\\varepsilon_{s1}$ = %.2f ‰" % es, xy=(es, h - d),
                xytext=(-6, 8), textcoords="offset points", fontsize=9,
                color=C1, ha="right", va="bottom")
    a1.axhline(h - x, color=MUT, lw=0.8, ls="--")
    a1.set_ylim(-0.05 * h, 1.12 * h)
    a1.set_xlabel("Dehnung $\\varepsilon$ [‰]   (Druck links)")
    a1.set_yticks([0, h - d, h - x, h])
    a1.set_yticklabels(["0", "$A_{s1}$", "x", "h"])
    a1.set_title("Dehnungen (Bemessungspunkt %s)" % fx.punkt, fontsize=10,
                 color=INK, loc="left")

    _stil(a2, grid="x")
    ys = np.linspace(h - x, h, 220)
    sg = np.array([C.sigma_c(ec * (y - (h - x)) / x) for y in ys])
    a2.fill_betweenx(ys, 0, sg, color=C2, alpha=0.22, lw=0)
    a2.plot(sg, ys, color=C2, lw=2.2)
    a2.plot([0, 0], [0, h], color=AXIS, lw=0.9)
    a2.axhline(h - x, color=MUT, lw=0.8, ls="--")
    y_res = h - (qs.d - fx.z)
    a2.plot([0, C.fcd], [y_res, y_res], color=C2, lw=0.8, ls=":")
    a2.annotate("$F_{cd}$ = %.0f kN" % fx.Fc, xy=(C.fcd * 0.5, y_res),
                xytext=(4, 6), textcoords="offset points", fontsize=9, color=C2)
    a2.annotate("", xy=(0, y_res), xytext=(0, h - qs.d),
                arrowprops=dict(arrowstyle="<->", color=C3, lw=1.4))
    a2.text(C.fcd * 0.06, (y_res + h - qs.d) / 2, "z = %.0f mm" % fx.z,
            color=C3, fontsize=9, va="center")
    a2.plot([0], [h - qs.d], "o", color=C1, ms=7)
    a2.annotate("$F_{sd}$ = %.0f kN" % (fx.As1 * fx.sigma_s1 / 1000.0),
                xy=(0, h - qs.d), xytext=(6, -12), textcoords="offset points",
                fontsize=9, color=C1)
    a2.set_xlim(-0.12 * C.fcd, 1.45 * C.fcd)
    a2.set_ylim(-0.05 * h, 1.12 * h)
    a2.set_yticks([0, h - qs.d, h - x, h])
    a2.set_yticklabels(["0", "$A_{s1}$", "x", "h"])
    a2.set_xlabel("Spannung $\\sigma_c$ [N/mm$^2$]   ($f_{cd}$ = %.1f)" % C.fcd)
    a2.set_title("Spannungen und Resultierende", fontsize=10, color=INK, loc="left")

    fig.suptitle("Bruchzustand unter Biegung  |  DIN EN 1992-1-1, 6.1 (2)P "
                 "Bild 6.1  und  3.1.7 (1) Bild 3.3", fontsize=11, color=INK,
                 x=0.012, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


# ===========================================================================
# 5. ALLGEMEINES BEMESSUNGSDIAGRAMM
# ===========================================================================
def bild_bemessungsdiagramm(b, ax=None):
    """Allgemeines Bemessungsdiagramm mu - omega - xi - zeta (EC2 6.1)."""
    from .biegung import bemessungsdiagramm, xi_grenz
    C, S = b["beton"], b["stahl"]
    d = bemessungsdiagramm(C, S, n=400, xi_max=0.62)
    fx = b.get("biegung_feld")

    if ax is None:
        fig, ax = _fig(9.0, 5.6)
    else:
        fig = ax.figure
    _stil(ax)
    m = d["mu"]
    ax.plot(m, d["omega"], color=C1, lw=2.2)
    ax.plot(m, d["xi"], color=C2, lw=2.2)
    ax.plot(m, d["zeta"], color=C3, lw=2.2)

    def _lab(mu_t, reihe, text, farbe, dy, ha="left"):
        y = float(np.interp(mu_t, m, d[reihe]))
        ax.annotate(text, xy=(mu_t, y), xytext=(6 if ha == "left" else -6, dy),
                    textcoords="offset points", color=farbe, fontsize=9.5,
                    ha=ha, weight="bold")

    _lab(m[-1] * 0.95, "omega", "$\\omega$ = $A_{s1}f_{yd}/(b\\,d\\,f_{cd})$",
         C1, -16, "right")
    _lab(m[-1] * 0.95, "xi", "$\\xi$ = x/d", C2, 8, "right")
    _lab(m[-1] * 0.30, "zeta", "$\\zeta$ = z/d", C3, 8, "left")

    lim = xi_grenz(C, b["eingabe"].delta)["xi_lim"]
    mu_lim = float(np.interp(lim, d["xi"], m))
    ax.axvline(mu_lim, color=MUT, lw=1.0, ls="--")
    ax.annotate("$\\mu_{lim}$ = %.3f\n($x_u/d$ = %.2f)" % (mu_lim, lim),
                xy=(mu_lim, 0.05), xytext=(6, 0), textcoords="offset points",
                color=INK2, fontsize=8.8, va="bottom")

    if fx is not None and fx.mu_Eds > 0:
        mu = min(fx.mu_Eds, m[-1])
        for val, col in ((fx.omega, C1), (fx.xi, C2), (fx.zeta, C3)):
            ax.plot([mu], [val], "o", color=col, ms=10, mec=SURF, mew=1.6, zorder=6)
        ax.axvline(mu, color=INK2, lw=1.0, ls=":")
        ax.annotate("$\\mu_{Eds}$ = %.4f\n$\\omega$ = %.4f\n$\\xi$ = %.4f\n"
                    "$\\zeta$ = %.4f" % (fx.mu_Eds, fx.omega, fx.xi, fx.zeta),
                    xy=(mu, 0.60), xytext=(10, 0), textcoords="offset points",
                    fontsize=9, color=INK,
                    bbox=dict(fc=SURF, ec=GRID, lw=0.8, boxstyle="round,pad=0.4"))

    ax.set_xlabel("$\\mu_{Eds}$ = $M_{Eds}/(b\\,d^2 f_{cd})$")
    ax.set_ylabel("$\\omega$ , $\\xi$ , $\\zeta$   [-]")
    ax.set_xlim(0, max(m[-1], 0.30))
    ax.set_ylim(0, 1.05)
    ax.set_title("Allgemeines Bemessungsdiagramm", fontsize=11, color=INK,
                 loc="left", pad=22)
    _unter(ax, "DIN EN 1992-1-1, 6.1 + 3.1.7 (1) | Grenze $x_u/d$: NA NDP zu 5.5 (4) | "
               "Beton %s, Stahl %s" % (C.klasse, S.sorte))
    fig.tight_layout()
    return fig


# ===========================================================================
# 6. QUERKRAFT
# ===========================================================================
def bild_querkraft(b):
    """Einhuellende V gegen die Tragfaehigkeiten und Buegelbewehrung (EC2 6.2)."""
    quer = b["querkraft"]
    e = b["eingabe"]
    x, V = quer.x, quer.V_Ed

    fig, (a1, a2) = _fig(10.5, 7.4, 2, 1, sharex=True)

    _stil(a1)
    a1.fill_between(x, 0, V, color=C2, alpha=0.20, lw=0)
    a1.plot(x, V, color=C2, lw=2.2)
    a1.axhline(quer.V_Rdc, color=C3, lw=1.8, ls="--")
    a1.annotate("$V_{Rd,c}$ = %.0f kN  [EC2 6.2.2, Gl. (6.2a/b)+NA]" % quer.V_Rdc,
                xy=(e.L * 0.5, quer.V_Rdc), xytext=(0, 5),
                textcoords="offset points", color=C3, fontsize=9, ha="center")
    if quer.V_Rdmax < 1.6 * float(np.max(V)):
        a1.axhline(quer.V_Rdmax, color=NOK, lw=1.8, ls="-.")
        a1.annotate("$V_{Rd,max}$ = %.0f kN  [Gl. (6.9)]" % quer.V_Rdmax,
                    xy=(e.L * 0.5, quer.V_Rdmax), xytext=(0, 5),
                    textcoords="offset points", color=NOK, fontsize=9, ha="center")
    a1.plot(x, np.zeros_like(x), color=AXIS, lw=0.9)
    k = int(np.argmax(V))
    a1.annotate("$|V_{Ed}|$ = %.0f kN" % V[k], xy=(x[k], V[k]), xytext=(8, -4),
                textcoords="offset points", color=C2, fontsize=9.5, weight="bold")
    a1.set_ylabel("$|V_{Ed}|$ [kN]")
    a1.set_ylim(0, max(float(np.max(V)) * 1.35,
                       quer.V_Rdmax * 1.1 if quer.V_Rdmax < 1.6 * float(np.max(V))
                       else float(np.max(V)) * 1.35))
    a1.set_title("Einhuellende Querkraft und Tragfaehigkeiten  |  "
                 "cot$\\theta$ = %.2f ($\\theta$ = %.1f Grad) [NA Gl. (6.7aDE)]"
                 % (quer.cot_theta, quer.theta), fontsize=10.5, color=INK,
                 loc="left")

    _stil(a2)
    from .baustoffe import stabflaeche
    A_stab = stabflaeche(e.phi_buegel)
    a2.fill_between(x, 0, quer.asw_schenkel / 100.0, color=C1, alpha=0.20, lw=0)
    a2.plot(x, quer.asw_schenkel / 100.0, color=C1, lw=2.0)
    if quer.asw_torsion > 0:
        a2.axhline(quer.asw_torsion / 100.0, color=C4, lw=1.6, ls="--")
        a2.annotate("Anteil Torsion $a_{sw,T}$ = %.2f cm$^2$/m  [EC2 6.3.2 (3)]"
                    % (quer.asw_torsion / 100.0),
                    xy=(e.L * 0.5, quer.asw_torsion / 100.0), xytext=(0, 5),
                    textcoords="offset points", color=C4, fontsize=8.8,
                    ha="center")
    a2.axhline(quer.asw_min / quer.n_schenkel / 100.0, color=MUT, lw=1.4, ls=":")
    a2.annotate("$a_{sw,min}$/Schenkel = %.2f cm$^2$/m  [NA NDP zu 9.2.2 (5)]"
                % (quer.asw_min / quer.n_schenkel / 100.0),
                xy=(e.L * 0.80, quer.asw_min / quer.n_schenkel / 100.0),
                xytext=(0, -14), textcoords="offset points", color=MUT,
                fontsize=8.8, ha="center")
    for z in quer.bereiche:
        asw_z = A_stab * 1000.0 / z["s"] / 100.0
        a2.plot([z["x1"], z["x2"]], [asw_z, asw_z], color=C2, lw=3.0,
                solid_capstyle="butt")
        a2.plot([z["x1"], z["x1"]], [0, asw_z], color=C2, lw=1.0, alpha=0.5)
        a2.plot([z["x2"], z["x2"]], [0, asw_z], color=C2, lw=1.0, alpha=0.5)
        a2.text(0.5 * (z["x1"] + z["x2"]), asw_z * 1.04,
                "$\\varnothing$%.0f/%.0f" % (e.phi_buegel, z["s"]), ha="center",
                va="bottom", color=C2, fontsize=8.8, weight="bold")
    a2.annotate("erforderlich je Schenkel", xy=(x[int(0.08 * len(x))],
                                                quer.asw_schenkel[int(0.08 * len(x))] / 100.0),
                xytext=(8, 12), textcoords="offset points", color=C1,
                fontsize=9, weight="bold")
    a2.set_ylabel("$a_{sw}$ je Schenkel [cm$^2$/m]")
    a2.set_xlabel("x [m]")
    a2.set_ylim(0, max(float(np.max(quer.asw_schenkel)) / 100.0 * 1.45, 1.0))
    a2.set_title("Erforderliche und gewaehlte Buegelbewehrung", fontsize=10.5,
                 color=INK, loc="left", pad=20)
    _unter(a2, "EC2 6.2.3 (3), Gl. (6.8) + Torsion 6.3.2 (3) | "
               "$s_{max}$ = %.0f mm [NA Tab. NA.9.1 / EC2 9.2.3]" % quer.s_max)
    fig.tight_layout()
    return fig


# ===========================================================================
# 7. TORSION
# ===========================================================================
def bild_torsion(b):
    """Ersatzhohlquerschnitt und Interaktion Torsion/Querkraft (EC2 6.3)."""
    e = b["eingabe"]
    tor = b["torsion"]
    quer = b["querkraft"]
    ehq = tor.detail.get("ehq", {})

    fig, (a1, a2) = _fig(11.0, 5.2, 1, 2,
                         gridspec_kw={"width_ratios": [1.0, 1.15]})

    # --- Ersatzhohlquerschnitt
    a1.set_facecolor(SURF)
    bw, h, t = e.b, e.h, tor.t_ef
    a1.add_patch(Rectangle((0, 0), bw, h, fc=BETON, ec=INK2, lw=1.4))
    if t > 0:
        aussen = [(0, 0), (bw, 0), (bw, h), (0, h)]
        innen = [(t, t), (bw - t, t), (bw - t, h - t), (t, h - t)]
        a1.add_patch(Polygon(aussen, closed=True, fc=C4, alpha=0.22, ec="none"))
        a1.add_patch(Polygon(innen, closed=True, fc=SURF, ec="none"))
        a1.add_patch(Rectangle((t / 2, t / 2), bw - t, h - t, fc="none",
                               ec=C4, lw=2.0, ls="--"))
        # Schubfluss
        for (xa, ya, dx, dy) in [(bw * 0.5, t / 2, bw * 0.18, 0),
                                 (bw - t / 2, h * 0.5, 0, h * 0.14),
                                 (bw * 0.5, h - t / 2, -bw * 0.18, 0),
                                 (t / 2, h * 0.5, 0, -h * 0.14)]:
            a1.add_patch(FancyArrow(xa - dx / 2, ya - dy / 2, dx, dy,
                                    width=bw * 0.012, color=C2,
                                    length_includes_head=True,
                                    head_width=bw * 0.055, zorder=6))
    a1.text(bw * 0.5, h * 0.5, "$A_k$ = %.0f mm$^2$\n$u_k$ = %.0f mm\n"
                               "$t_{ef}$ = %.0f mm" % (tor.A_k, tor.u_k, t),
            ha="center", va="center", color=INK, fontsize=9.5,
            bbox=dict(fc=SURF, ec=GRID, boxstyle="round,pad=0.4"))
    a1.annotate("Schubfluss\n$\\tau_t t_{ef}$ = $T_{Ed}/(2A_k)$",
                xy=(bw * 0.5, -0.08 * h), ha="center", va="top", color=C2,
                fontsize=9)
    a1.set_xlim(-0.25 * bw, 1.25 * bw)
    a1.set_ylim(-0.22 * h, 1.10 * h)
    a1.set_aspect("equal")
    a1.axis("off")
    a1.set_title("Ersatzhohlquerschnitt", fontsize=10.5, color=INK, loc="left",
                 pad=20)
    _unter(a1, "EC2 6.3.1 (3) + 6.3.2 (1), Gl. (6.26)")

    # --- Interaktionsdiagramm
    _stil(a2, grid="both")
    tt = np.linspace(0, 1, 100)
    a2.plot(tt, 1 - tt, color=NOK, lw=2.2)
    a2.fill_between(tt, 0, 1 - tt, color=C3, alpha=0.13)
    a2.annotate("$T_{Ed}/T_{Rd,max}$ + $V_{Ed}/V_{Rd,max}$ = 1,0\n[EC2 6.3.2 (4), "
                "Gl. (6.29)]", xy=(0.52, 0.50), xytext=(0, 14),
                textcoords="offset points", color=NOK, fontsize=9, ha="center")
    a2.text(0.16, 0.16, "zulaessiger\nBereich", color=C3, fontsize=10,
            ha="center", weight="bold")
    if tor.T_Rdmax > 0 and quer.V_Rdmax > 0:
        xt = tor.T_Ed / tor.T_Rdmax
        yv = tor.V_Ed / quer.V_Rdmax
        farbe = OK if tor.interaktion <= 1.0 else NOK
        a2.plot([xt], [yv], "o", color=farbe, ms=13, mec=SURF, mew=2, zorder=6)
        a2.annotate("Bemessungspunkt\n$T_{Ed}/T_{Rd,max}$ = %.3f\n"
                    "$V_{Ed}/V_{Rd,max}$ = %.3f\nSumme = %.3f  (%s)"
                    % (xt, yv, tor.interaktion,
                       "OK" if tor.interaktion <= 1.0 else "NICHT ERFUELLT"),
                    xy=(xt, yv), xytext=(16, 16), textcoords="offset points",
                    fontsize=9, color=INK,
                    bbox=dict(fc=SURF, ec=GRID, lw=0.8, boxstyle="round,pad=0.4"),
                    arrowprops=dict(arrowstyle="->", color=MUT, lw=0.9))
    a2.set_xlim(0, 1.15)
    a2.set_ylim(0, 1.15)
    a2.set_xlabel("$T_{Ed}$ / $T_{Rd,max}$   [-]")
    a2.set_ylabel("$V_{Ed}$ / $V_{Rd,max}$   [-]")
    a2.set_title("Interaktion Torsion - Querkraft (Druckstreben)", fontsize=10.5,
                 color=INK, loc="left", pad=20)
    _unter(a2, "$T_{Rd,max}$ = %.1f kNm [Gl. (6.30), NA $\\nu$ = %.3f] | "
               "$V_{Rd,max}$ = %.0f kN"
           % (tor.T_Rdmax, tor.detail.get("T_Rdmax", {}).get("nu", 0.0),
              quer.V_Rdmax))
    fig.tight_layout()
    return fig


# ===========================================================================
# 8. ZUGKRAFTDECKUNGSLINIE
# ===========================================================================
def bild_zugkraftdeckung(b, ax=None):
    """Zugkraftdeckungslinie mit Versatzmass (EC2 9.2.1.3)."""
    from .konstruktion import zugkraft
    from .baustoffe import stabflaeche
    e, S = b["eingabe"], b["stahl"]
    einh, quer, tor = b["einhuellende"], b["querkraft"], b.get("torsion")
    x = einh["x"]
    Mpos = np.maximum(einh["Mmax"], 0.0)
    F_tor = 0.0
    if tor is not None and abs(e.T_Ed) > 1e-9 and tor.erforderlich:
        F_tor = tor.asl_gesamt * S.fyd / 1000.0 * (tor.detail["ehq"]["b_k"]
                                                   / tor.u_k)
    zk = zugkraft(x, Mpos, quer.z, quer.a_l, 0.0, F_tor)

    if ax is None:
        fig, ax = _fig(10.5, 5.2)
    else:
        fig = ax.figure
    _stil(ax)
    Fs0 = np.abs(Mpos) * 1.0e6 / quer.z / 1.0e3
    ax.plot(x, Fs0, color=MUT, lw=1.4, ls="--")
    ax.annotate("$F_s$ = M/z (ohne Versatz)", xy=(x[len(x) // 2], Fs0[len(x) // 2]),
                xytext=(0, -18), textcoords="offset points", color=MUT,
                fontsize=8.8, ha="center")
    ax.fill_between(x, 0, zk["Fs"], color=C1, alpha=0.18, lw=0)
    ax.plot(x, zk["Fs"], color=C1, lw=2.2)
    txt = "$F_s$ mit Versatzmass $a_l$ = %.0f mm" % quer.a_l
    if F_tor > 0:
        txt += "\n+ Torsionslaengskraft %.0f kN" % F_tor
    ax.annotate(txt, xy=(x[len(x) // 2], zk["Fs"][len(x) // 2]), xytext=(0, 10),
                textcoords="offset points", color=C1, fontsize=9.5, ha="center",
                weight="bold")

    n = b["n_unten"]
    if n > 0:
        F_stab = stabflaeche(e.phi_laengs) * S.fyd / 1000.0
        for k in range(1, n + 1):
            ax.axhline(k * F_stab, color=C2, lw=1.0, alpha=0.45)
        ax.axhline(n * F_stab, color=C2, lw=2.2)
        ax.annotate("Tragfaehigkeit %d $\\varnothing$%.0f = %.0f kN"
                    % (n, e.phi_laengs, n * F_stab), xy=(x[-1], n * F_stab),
                    xytext=(-6, 8), textcoords="offset points", color=C2,
                    fontsize=9.5, ha="right", weight="bold")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("$F_{s}$ [kN]")
    ax.set_ylim(0, None)
    ax.set_title("Zugkraftdeckungslinie", fontsize=11, color=INK, loc="left",
                 pad=22)
    _unter(ax, "Versatzmass $a_l$ = z cot$\\theta$/2 : DIN EN 1992-1-1, "
               "9.2.1.3 (2), Gl. (9.2)")
    fig.tight_layout()
    return fig


# ===========================================================================
# 9. DURCHBIEGUNG
# ===========================================================================
def bild_durchbiegung(b, ax=None):
    """Durchbiegung unter quasi-staendiger Kombination (EC2 7.4.3)."""
    db = b["durchbiegung"]
    e = b["eingabe"]
    if ax is None:
        fig, ax = _fig(10.5, 4.6)
    else:
        fig = ax.figure
    _stil(ax)
    grenz = None
    for f in db:
        if len(f["x"]) == 0:
            continue
        ax.fill_between(f["x"], 0, f["w"], color=C1, alpha=0.18, lw=0)
        ax.plot(f["x"], f["w"], color=C1, lw=2.2)
        j = int(np.argmax(np.abs(f["w"])))
        ax.plot([f["x"][j]], [f["w"][j]], "o", color=C1, ms=8)
        ax.annotate("$w_{max}$ = %.2f mm" % f["w"][j], xy=(f["x"][j], f["w"][j]),
                    xytext=(0, 10), textcoords="offset points", ha="center",
                    color=C1, fontsize=9.5, weight="bold")
        grenz = f["w_grenz"]
        ax.plot([f["feld"][0], f["feld"][1]], [grenz, grenz], color=NOK, lw=1.6,
                ls="--")
    if grenz:
        ax.annotate("Grenzwert L/%.0f = %.2f mm  [EC2 7.4.1 (4)]"
                    % (e.grenze_durchbiegung, grenz), xy=(e.L * 0.5, grenz),
                    xytext=(0, 6), textcoords="offset points", color=NOK,
                    fontsize=9, ha="center")
    ax.axhline(0, color=AXIS, lw=0.9)
    ax.invert_yaxis()
    ax.set_xlabel("x [m]")
    ax.set_ylabel("Durchbiegung w [mm]  (nach unten)")
    ax.set_title("Durchbiegung im Grenzzustand der Gebrauchstauglichkeit",
                 fontsize=11, color=INK, loc="left", pad=22)
    _unter(ax, "Interpolation Zustand I/II: EC2 7.4.3 (3), Gl. (7.18)/(7.19) | "
               "quasi-staendige Kombination EC0 Gl. (6.16b)")
    fig.tight_layout()
    return fig


# ===========================================================================
# 10. AUSNUTZUNGSGRADE
# ===========================================================================
def bild_ausnutzung(b, ax=None):
    """Ausnutzungsgrad eta jedes Nachweises."""
    nw = b["nachweise"]
    if ax is None:
        fig, ax = _fig(9.5, max(3.2, 0.62 * len(nw) + 1.6))
    else:
        fig = ax.figure
    _stil(ax, grid="x")
    namen, etas, oks = [], [], []
    for c in nw:
        if c["vergleich"] == "<=":
            et = c["wert"] / c["grenzwert"] if c["grenzwert"] else 0.0
        else:
            et = c["grenzwert"] / c["wert"] if c["wert"] else 0.0
        namen.append(c["name"])
        etas.append(et)
        oks.append(c["ok"])
    y = np.arange(len(namen))
    for i, (et, ok) in enumerate(zip(etas, oks)):
        ax.barh(y[i], et, height=0.58, color=(OK if ok else NOK),
                hatch=None if ok else "///", edgecolor=SURF, linewidth=2.0,
                zorder=3)
        ax.text(et + 0.02, y[i], "%.2f  %s" % (et, "OK" if ok else "NICHT ERFUELLT"),
                va="center", ha="left", fontsize=9,
                color=(INK2 if ok else NOK), weight="normal" if ok else "bold")
    ax.axvline(1.0, color=INK2, lw=1.4, ls="--")
    ax.annotate("$\\eta$ = 1,00", xy=(1.0, -0.55), xytext=(4, 0),
                textcoords="offset points", color=INK2, fontsize=9,
                va="center", annotation_clip=False)
    ax.set_yticks(y)
    ax.set_yticklabels(namen, fontsize=9)
    ax.set_xlim(0, max(1.25, max(etas) * 1.30 if etas else 1.25))
    ax.set_xlabel("Ausnutzungsgrad $\\eta$ = $E_d$ / $R_d$   [-]")
    ax.invert_yaxis()
    ax.set_title("Zusammenstellung der Nachweise", fontsize=11, color=INK,
                 loc="left", pad=22)
    _unter(ax, "Zustand zusaetzlich durch Textbeschriftung und Schraffur "
               "gekennzeichnet, nicht nur durch die Farbe")
    fig.tight_layout()
    return fig


# ===========================================================================
# Erzeugung aller Bilder
# ===========================================================================
BILDER = [
    ("Querschnitt", bild_querschnitt),
    ("Baustoffe", bild_baustoffe),
    ("Schnittgroessen", bild_schnittgroessen),
    ("Biegung: Dehnungen", bild_dehnungen),
    ("Bemessungsdiagramm", bild_bemessungsdiagramm),
    ("Querkraft", bild_querkraft),
    ("Torsion", bild_torsion),
    ("Zugkraftdeckung", bild_zugkraftdeckung),
    ("Durchbiegung", bild_durchbiegung),
    ("Ausnutzung", bild_ausnutzung),
]


def alle_bilder(b):
    """Liefert [(Name, Figure), ...] mit allen Grafiken des Berichts."""
    out = []
    for name, fn in BILDER:
        try:
            out.append((name, fn(b)))
        except Exception as exc:
            fig, ax = _fig(8, 3)
            ax.axis("off")
            ax.text(0.5, 0.5, "Fehler beim Erzeugen von '%s':\n%s" % (name, exc),
                    ha="center", va="center", color=NOK, fontsize=9, wrap=True)
            out.append((name, fig))
    return out


def pdf_export(b, pfad):
    """Exportiert alle Grafiken in eine PDF-Datei."""
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(pfad) as pdf:
        for _, f in alle_bilder(b):
            pdf.savefig(f, facecolor=SURF)
            plt.close(f)
    return pfad
