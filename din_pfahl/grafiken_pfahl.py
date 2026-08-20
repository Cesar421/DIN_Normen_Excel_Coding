# -*- coding: utf-8 -*-
"""
Grafiken der Pfahlbemessung (matplotlib).

Farbpalette und Stil wie im Modul din_balken.grafiken; Zustaende werden
zusaetzlich durch Text und Schraffur gekennzeichnet, nicht nur durch die Farbe.
"""

import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Polygon, Wedge

from din_balken.grafiken import (C1, C2, C3, C4, OK, NOK, SURF, INK, INK2,
                                 MUT, GRID, AXIS, BETON, _fig, _stil, _unter)

BODEN = ["#d9d2c5", "#cfc4ae", "#c2b79c", "#b6a88a", "#a99a78"]


# ===========================================================================
# 1. PFAHLQUERSCHNITT
# ===========================================================================
def bild_querschnitt(b, ax=None):
    """Kreisquerschnitt mit Laengsbewehrung und Wendel (DIN EN 1536, 7.6)."""
    qs = b["querschnitt"]
    qk = b["querkraft"]
    D, R = qs.D, qs.R
    if ax is None:
        fig, ax = _fig(6.4, 6.6)
    else:
        fig = ax.figure
    ax.set_facecolor(SURF)

    ax.add_patch(Circle((0, 0), R, fc=BETON, ec=INK2, lw=1.6, zorder=1))
    # Wendel
    r_w = R - qs.c_nom - qs.phi_w / 2.0
    ax.add_patch(Circle((0, 0), r_w, fc="none", ec=C2, lw=2.2, zorder=3))
    ax.add_patch(Circle((0, 0), r_w - qs.phi_w, fc="none", ec=C2, lw=2.2,
                        alpha=0.35, zorder=3))
    # Laengsstaebe
    rs = qs.D_s / 2.0
    for i in range(qs.n_l):
        a = 2 * math.pi * i / qs.n_l
        ax.add_patch(Circle((rs * math.sin(a), rs * math.cos(a)), qs.phi_l / 2.0,
                            fc=C1, ec=SURF, lw=0.9, zorder=5))

    # Bemassung
    ax.annotate("", xy=(-R, -1.22 * R), xytext=(R, -1.22 * R),
                arrowprops=dict(arrowstyle="<->", color=MUT, lw=0.9))
    ax.text(0, -1.30 * R, "D = %.0f mm" % D, ha="center", va="top", color=INK2,
            fontsize=10)
    ax.annotate("", xy=(0, 0), xytext=(rs * math.sin(math.radians(35)),
                                       rs * math.cos(math.radians(35))),
                arrowprops=dict(arrowstyle="<->", color=C1, lw=1.0))
    ax.text(rs * 0.35, rs * 0.62, "$D_s$/2 = %.0f" % rs, color=C1, fontsize=8.5,
            rotation=-35)
    ax.annotate("$c_{nom}$ = %.0f mm\n[DIN EN 1536, 7.6.2]" % qs.c_nom,
                xy=(R * math.sin(math.radians(215)) * 0.97,
                    R * math.cos(math.radians(215)) * 0.97),
                xytext=(-1.55 * R, -0.85 * R),
                arrowprops=dict(arrowstyle="->", color=MUT, lw=0.8),
                color=MUT, fontsize=8.5, ha="left")

    ax.text(0, 1.12 * R, "%d $\\varnothing$%.0f  =  %.0f mm$^2$   ($\\rho$ = %.2f %%)"
            % (qs.n_l, qs.phi_l, qs.As_ges, 100 * qs.rho_l),
            ha="center", va="bottom", color=C1, fontsize=10, weight="bold")
    ax.text(1.10 * R, -0.15 * R,
            "Wendel $\\varnothing$%.0f\nGanghoehe %.0f mm" % (qs.phi_w,
                                                              qk["s_gewaehlt"]),
            ha="left", va="center", color=C2, fontsize=9.5)
    ax.text(-1.55 * R, 0.95 * R,
            "lichter Stababstand\n%.0f mm  (min. %.0f mm)"
            % (qs.lichter_stababstand(), 80.0), color=INK2, fontsize=8.5,
            ha="left", va="top")

    ax.set_xlim(-1.75 * R, 1.75 * R)
    ax.set_ylim(-1.45 * R, 1.35 * R)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Pfahlquerschnitt mit Bewehrungskorb", fontsize=11.5,
                 color=INK, loc="left", pad=22)
    _unter(ax, "DIN EN 1536, 7.6.2 / 7.6.3 / 7.6.4 | EC2 9.5.3")
    fig.tight_layout()
    return fig


# ===========================================================================
# 2. LAENGSSCHNITT MIT BAUGRUND
# ===========================================================================
def bild_laengsschnitt(b, ax=None):
    """Pfahllaengsschnitt mit Schichtenfolge und Widerstandsanteilen."""
    e = b["eingabe"]
    trag = b["tragfaehigkeit"]
    D_m = e.D / 1000.0
    if ax is None:
        fig, ax = _fig(7.6, 8.0)
    else:
        fig = ax.figure
    ax.set_facecolor(SURF)

    breite = max(4.0 * D_m, 2.0)
    for i, s in enumerate(e.schichten):
        ax.add_patch(Rectangle((-breite / 2, -s.z_u), breite, s.dicke,
                               fc=BODEN[i % len(BODEN)], ec=MUT, lw=0.6, zorder=1))
        ax.text(-breite / 2 + 0.06, -(s.z_o + s.z_u) / 2,
                "%s\n$E_s$ = %.0f kN/m$^2$\n$q_{s,k}$ = %.0f kN/m$^2$"
                % (s.name[:20], s.E_s, s.q_s_k), fontsize=7.5, color=INK2,
                va="center", ha="left", linespacing=1.5)

    # Pfahl
    ax.add_patch(Rectangle((-D_m / 2, -e.L), D_m, e.L, fc=BETON, ec=INK,
                           lw=1.6, zorder=3))
    lk = e.laenge_bewehrungskorb if e.laenge_bewehrungskorb > 0 else e.L
    ax.add_patch(Rectangle((-D_m / 2 * 0.72, -lk), D_m * 0.72, lk, fc="none",
                           ec=C1, lw=1.4, ls="--", zorder=4))
    ax.annotate("Bewehrungskorb\n%.1f m" % lk, xy=(-D_m / 2 * 0.72, -lk * 0.72),
                xytext=(-breite / 2 * 0.62, -lk * 0.72),
                arrowprops=dict(arrowstyle="->", color=C1, lw=0.8),
                ha="right", va="center", color=C1, fontsize=8.5)

    # Mantelreibung als Pfeile
    for a in trag.anteile:
        if a["dicke"] <= 0:
            continue
        for z in np.linspace(a["z_o"] + 0.25, a["z_u"] - 0.25,
                             max(2, int(a["dicke"]))):
            ax.annotate("", xy=(D_m / 2, -z),
                        xytext=(D_m / 2 + 0.10 + a["q_s_k"] / 400.0, -z),
                        arrowprops=dict(arrowstyle="->", color=C3, lw=1.1))
    ax.text(D_m / 2 + 0.55, -e.L * 0.45,
            "Mantelreibung\n$R_{s,k}$ = %.0f kN" % trag.R_s_k, color=C3,
            fontsize=9.5, ha="left", va="center", weight="bold")

    # Spitzendruck
    for xx in np.linspace(-D_m / 2 * 0.8, D_m / 2 * 0.8, 5):
        ax.annotate("", xy=(xx, -e.L), xytext=(xx, -e.L - 0.55),
                    arrowprops=dict(arrowstyle="->", color=C2, lw=1.4))
    ax.text(0, -e.L - 0.75, "Spitzendruck  $R_{b,k}$ = %.0f kN" % trag.R_b_k,
            ha="center", va="top", color=C2, fontsize=9.5, weight="bold")

    # Kopflasten
    ax.annotate("", xy=(0, 0), xytext=(0, 1.30),
                arrowprops=dict(arrowstyle="->", color=INK, lw=2.4))
    ax.text(0.06, 1.35, "$N_{Ed}$ = %.0f kN" % e.N_Ed, color=INK, fontsize=9.5,
            ha="left")
    if abs(e.H_Ed) > 1e-9:
        ax.annotate("", xy=(0, -0.12), xytext=(-1.05, -0.12),
                    arrowprops=dict(arrowstyle="->", color=C4, lw=2.2))
        ax.text(-1.10, -0.12, "$H_{Ed}$ = %.0f kN" % e.H_Ed, color=C4,
                fontsize=9.5, ha="right", va="center")

    ax.set_xlim(-breite / 2 - 0.2, breite / 2 + 1.5)
    ax.set_ylim(-e.L - 1.6, 1.9)
    ax.set_xlabel("[m]")
    ax.set_ylabel("Tiefe z [m]")
    ax.set_yticks(np.arange(0, -e.L - 1, -2))
    ax.set_yticklabels(["%.0f" % abs(v) for v in np.arange(0, -e.L - 1, -2)])
    _stil(ax, grid=None)
    ax.set_title("Pfahllaengsschnitt und Baugrund", fontsize=11.5, color=INK,
                 loc="left", pad=22)
    _unter(ax, "DIN EN 1997-1, 7.6.2.3 Gl. (7.8) | $q_{b,k}$, $q_{s,k}$ nach "
               "EA-Pfaehle Tab. 5.12-5.15 (Eingabewerte)")
    fig.tight_layout()
    return fig


# ===========================================================================
# 3. M-N-INTERAKTIONSDIAGRAMM
# ===========================================================================
def bild_interaktion(b, ax=None):
    """M-N-Interaktionsdiagramm des Kreisquerschnitts (EC2 6.1)."""
    e = b["eingabe"]
    dg = b["interaktion"]
    qs = b["querschnitt"]
    if ax is None:
        fig, ax = _fig(8.6, 6.2)
    else:
        fig = ax.figure
    _stil(ax, grid="both")

    M = np.concatenate([dg["M"], -dg["M"][::-1]])
    N = np.concatenate([dg["N"], dg["N"][::-1]])
    ax.fill(M, N, color=C1, alpha=0.13, zorder=1)
    ax.plot(dg["M"], dg["N"], color=C1, lw=2.4, zorder=3)
    ax.plot(-dg["M"], dg["N"], color=C1, lw=2.4, alpha=0.5, zorder=3)
    ax.axhline(0, color=AXIS, lw=1.0)
    ax.axvline(0, color=AXIS, lw=1.0)

    i_max = int(np.argmax(dg["M"]))
    ax.plot([dg["M"][i_max]], [dg["N"][i_max]], "o", color=C1, ms=7)
    ax.annotate("$M_{Rd,max}$ = %.0f kNm\nbei $N$ = %.0f kN"
                % (dg["M"][i_max], dg["N"][i_max]),
                xy=(dg["M"][i_max], dg["N"][i_max]), xytext=(10, 0),
                textcoords="offset points", color=C1, fontsize=8.8, va="center")
    ax.annotate("zentrischer Druck\n$N_{Rd}$ = %.0f kN" % dg["N_druck_max"],
                xy=(0, dg["N_druck_max"]), xytext=(12, 8),
                textcoords="offset points", color=INK2, fontsize=8.8)
    ax.annotate("zentrischer Zug\n$N_{Rd}$ = %.0f kN" % dg["N_zug"],
                xy=(0, dg["N_zug"]), xytext=(12, -6), textcoords="offset points",
                color=INK2, fontsize=8.8)

    M_Rd = b["M_Rd"]
    farbe = OK if b["M_Ed_max"] <= M_Rd else NOK
    ax.plot([b["M_Ed_max"]], [e.N_Ed], "o", color=farbe, ms=13, mec=SURF, mew=2,
            zorder=6)
    ax.plot([M_Rd], [e.N_Ed], "s", color=C2, ms=9, mec=SURF, mew=1.5, zorder=6)
    ax.plot([b["M_Ed_max"], M_Rd], [e.N_Ed, e.N_Ed], color=C2, lw=1.2, ls=":")
    ax.annotate("Bemessungspunkt\n$M_{Ed}$ = %.0f kNm\n$N_{Ed}$ = %.0f kN\n"
                "$M_{Rd}$ = %.0f kNm  ($\\eta$ = %.2f)"
                % (b["M_Ed_max"], e.N_Ed, M_Rd, b["M_Ed_max"] / max(M_Rd, 1e-9)),
                xy=(b["M_Ed_max"], e.N_Ed), xytext=(20, -55),
                textcoords="offset points", fontsize=9, color=INK,
                bbox=dict(fc=SURF, ec=GRID, lw=0.8, boxstyle="round,pad=0.4"),
                arrowprops=dict(arrowstyle="->", color=MUT, lw=0.9))

    ax.set_xlabel("$M_{Rd}$ [kNm]")
    ax.set_ylabel("$N_{Rd}$ [kN]   (Druck negativ)")
    ax.set_title("M-N-Interaktionsdiagramm  |  D = %.0f mm, %d $\\varnothing$%.0f"
                 % (qs.D, qs.n_l, qs.phi_l), fontsize=11, color=INK, loc="left",
                 pad=22)
    _unter(ax, "DIN EN 1992-1-1, 6.1 (2)P Bild 6.1 (Bemessungspunkte A/B/C) | "
               "Beton %s, Stahl %s, $\\rho$ = %.2f %%"
           % (b["beton"].klasse, b["stahl"].sorte, 100 * qs.rho_l))
    fig.tight_layout()
    return fig


# ===========================================================================
# 4. HORIZONTAL BELASTETER PFAHL
# ===========================================================================
def bild_horizontal(b):
    """Verlaeufe w(z), M(z) und V(z) aus dem Bettungsmodulverfahren."""
    hor = b["horizontal"]
    e = b["eingabe"]
    fig, (a1, a2, a3) = _fig(11.0, 6.6, 1, 3, sharey=True)

    for ax, wert, farbe, label, einheit in (
            (a1, hor.w, C1, "Verschiebung w", "mm"),
            (a2, hor.M, C2, "Biegemoment M", "kNm"),
            (a3, hor.V, C3, "Querkraft V", "kN")):
        _stil(ax, grid="both")
        ax.plot(wert, -hor.z, color=farbe, lw=2.2)
        ax.fill_betweenx(-hor.z, 0, wert, color=farbe, alpha=0.18)
        ax.axvline(0, color=AXIS, lw=1.0)
        j = int(np.argmax(np.abs(wert)))
        ax.plot([wert[j]], [-hor.z[j]], "o", color=farbe, ms=8)
        ax.annotate("%.1f %s\nz = %.2f m" % (wert[j], einheit, hor.z[j]),
                    xy=(wert[j], -hor.z[j]), xytext=(8, 0),
                    textcoords="offset points", color=farbe, fontsize=9,
                    va="center", weight="bold")
        ax.set_xlabel("%s [%s]" % (label, einheit))
    a1.set_ylabel("Tiefe z [m]")
    a1.set_yticks(np.arange(0, -e.L - 1, -2))
    a1.set_yticklabels(["%.0f" % abs(v) for v in np.arange(0, -e.L - 1, -2)])

    fig.suptitle("Horizontal belasteter Pfahl - Bettungsmodulverfahren  |  "
                 "$H_{Ed}$ = %.0f kN, $M_{Ed}$ = %.0f kNm, Kopf %s  "
                 "[EA-Pfaehle, 6.3]"
                 % (e.H_Ed, e.M_Ed_kopf, e.kopf), fontsize=10.5, color=INK,
                 x=0.012, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


# ===========================================================================
# 5. WIDERSTANDS-SETZUNGS-LINIE
# ===========================================================================
def bild_wsl(b, ax=None):
    """Widerstands-Setzungs-Linie nach EA-Pfaehle, 5.4.5."""
    wsl = b["wsl"]
    trag = b["tragfaehigkeit"]
    e = b["eingabe"]
    if ax is None:
        fig, ax = _fig(8.6, 5.4)
    else:
        fig = ax.figure
    _stil(ax, grid="both")

    ax.plot(wsl["R_s"], wsl["s_mm"], color=C3, lw=2.0)
    ax.plot(wsl["R_b"], wsl["s_mm"], color=C2, lw=2.0)
    ax.plot(wsl["R"], wsl["s_mm"], color=C1, lw=2.6)
    ax.annotate("$R_{s,k}$ (Mantel)", xy=(wsl["R_s"][-1], wsl["s_mm"][-1]),
                xytext=(-8, -10), textcoords="offset points", color=C3,
                fontsize=9, ha="right", weight="bold")
    ax.annotate("$R_{b,k}$ (Spitze)", xy=(wsl["R_b"][-1], wsl["s_mm"][-1]),
                xytext=(-8, -10), textcoords="offset points", color=C2,
                fontsize=9, ha="right", weight="bold")
    ax.annotate("$R_{c,k}$ (gesamt)", xy=(wsl["R"][-1], wsl["s_mm"][-1]),
                xytext=(-8, -10), textcoords="offset points", color=C1,
                fontsize=9.5, ha="right", weight="bold")

    ax.axhline(wsl["s_sg"], color=MUT, lw=1.0, ls=":")
    ax.text(wsl["R"][-1] * 0.02, wsl["s_sg"],
            "$s_{sg}$ = %.1f mm" % wsl["s_sg"], color=MUT, fontsize=8.5,
            va="bottom")
    ax.axhline(wsl["s_g"], color=MUT, lw=1.0, ls=":")
    ax.text(wsl["R"][-1] * 0.02, wsl["s_g"], "$s_g$ = 0,10 D = %.1f mm" % wsl["s_g"],
            color=MUT, fontsize=8.5, va="bottom")

    N_k = abs(min(e.N_k, 0.0))
    ax.plot([N_k], [b["setzung"]], "o", color=C1, ms=12, mec=SURF, mew=2, zorder=6)
    ax.annotate("$N_k$ = %.0f kN\ns = %.1f mm" % (N_k, b["setzung"]),
                xy=(N_k, b["setzung"]), xytext=(14, 14),
                textcoords="offset points", fontsize=9, color=INK,
                bbox=dict(fc=SURF, ec=GRID, lw=0.8, boxstyle="round,pad=0.4"),
                arrowprops=dict(arrowstyle="->", color=MUT, lw=0.9))
    ax.axvline(trag.R_c_d, color=NOK, lw=1.6, ls="--")
    ax.annotate("$R_{c,d}$ = %.0f kN\n[DIN 1054, Tab. A 2.3]" % trag.R_c_d,
                xy=(trag.R_c_d, wsl["s_mm"][-1] * 0.75), xytext=(6, 0),
                textcoords="offset points", color=NOK, fontsize=8.8)

    ax.set_xlabel("Pfahlwiderstand R [kN]")
    ax.set_ylabel("Setzung s [mm]")
    ax.invert_yaxis()
    ax.set_title("Widerstands-Setzungs-Linie", fontsize=11, color=INK,
                 loc="left", pad=22)
    _unter(ax, "EA-Pfaehle, 5.4.5 (vereinfacht linear) | "
               "$s_{sg}$ = 0,5 $R_{s,k}$[MN] + 0,5 cm <= 3 cm ; $s_g$ = 0,10 D")
    fig.tight_layout()
    return fig


# ===========================================================================
# 6. TRAGFAEHIGKEITSANTEILE
# ===========================================================================
def bild_tragfaehigkeit(b, ax=None):
    """Aufteilung des Pfahlwiderstands auf Schichten und Spitze."""
    trag = b["tragfaehigkeit"]
    e = b["eingabe"]
    if ax is None:
        fig, ax = _fig(9.0, max(3.4, 0.6 * (len(trag.anteile) + 3)))
    else:
        fig = ax.figure
    _stil(ax, grid="x")

    namen = ["%s (%.1f-%.1f m)" % (a["name"], a["z_o"], a["z_u"])
             for a in trag.anteile] + ["Spitzendruck $R_{b,k}$"]
    werte = [a["R_s_k"] for a in trag.anteile] + [trag.R_b_k]
    farben = [C3] * len(trag.anteile) + [C2]
    y = np.arange(len(namen))
    for i, (v, f) in enumerate(zip(werte, farben)):
        ax.barh(y[i], v, height=0.6, color=f, edgecolor=SURF, lw=2.0, zorder=3)
        ax.text(v + max(werte) * 0.015, y[i], "%.0f kN  (%.0f %%)"
                % (v, 100 * v / max(trag.R_c_k, 1e-9)), va="center", ha="left",
                fontsize=9, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels(namen, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, max(werte) * 1.35)
    ax.set_xlabel("charakteristischer Widerstand [kN]")
    ax.set_title("Anteile am Pfahlwiderstand  |  $R_{c,k}$ = %.0f kN , "
                 "$R_{c,d}$ = %.0f kN , $F_{c,d}$ = %.0f kN"
                 % (trag.R_c_k, trag.R_c_d, trag.F_c_d), fontsize=10.5,
                 color=INK, loc="left", pad=22)
    _unter(ax, "DIN EN 1997-1, 7.6.2.3 Gl. (7.8) | Teilsicherheitsbeiwerte "
               "DIN 1054 Tab. A 2.3 (%s)" % trag.situation)
    fig.tight_layout()
    return fig


# ===========================================================================
# 7. AUSNUTZUNGSGRADE
# ===========================================================================
def bild_ausnutzung(b, ax=None):
    """Ausnutzungsgrad eta jedes Nachweises."""
    from din_balken.grafiken import bild_ausnutzung as _ba
    return _ba(b, ax)


# ===========================================================================
BILDER = [
    ("Laengsschnitt", bild_laengsschnitt),
    ("Querschnitt", bild_querschnitt),
    ("Tragfaehigkeit", bild_tragfaehigkeit),
    ("Widerstands-Setzungs-Linie", bild_wsl),
    ("Horizontallast", bild_horizontal),
    ("M-N-Interaktion", bild_interaktion),
    ("Ausnutzung", bild_ausnutzung),
]


def alle_bilder(b):
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
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(pfad) as pdf:
        for _, f in alle_bilder(b):
            pdf.savefig(f, facecolor=SURF)
            plt.close(f)
    return pfad
