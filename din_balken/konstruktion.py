# -*- coding: utf-8 -*-
"""
BEWEHRUNGS- UND KONSTRUKTIONSREGELN.

Norm: DIN EN 1992-1-1:2011-01, Abschnitte 8 und 9 + DIN EN 1992-1-1/NA:2013-04

  8.2      Stababstaende
  8.4      Verankerung (fbd, lb,rqd, lbd)
  8.7.3    Uebergreifungsstoesse l0
  9.2.1.1  Mindest- und Hoechstbewehrung fuer Biegung     Gl. (9.1N)
  9.2.1.1  Robustheitsbewehrung (NA)
  9.2.1.3  Versatzmassregel / Zugkraftdeckung
  9.2.1.4  Verankerung am Endauflager                     Gl. (9.3)
  9.2.2    Querkraftbewehrung
  9.2.3    Torsionsbewehrung
"""

import math

import numpy as np

from .normen import ref
from .baustoffe import stabflaeche


# ---------------------------------------------------------------------------
# 9.2.1.1  Mindest- und Hoechstbewehrung
# ---------------------------------------------------------------------------
def mindestbewehrung_biegung(querschnitt, beton, stahl, bt=None):
    """
    DIN EN 1992-1-1, 9.2.1.1 (1), Gl. (9.1N):

        As,min = max( 0,26 (fctm/fyk) bt d ; 0,0013 bt d )

    bt = mittlere Breite der Zugzone (bw beim Plattenbalken mit gezogener Platte).
    """
    bt = querschnitt.bw if bt is None else bt
    d = querschnitt.d
    a1 = 0.26 * beton.fctm / stahl.fyk * bt * d
    a2 = 0.0013 * bt * d
    return dict(As_min=max(a1, a2), term_fctm=a1, term_0013=a2, bt=bt,
                massgebend="0,26 fctm/fyk bt d" if a1 >= a2 else "0,0013 bt d",
                norm=ref("As_min"))


def hoechstbewehrung(querschnitt):
    """DIN EN 1992-1-1/NA, NDP zu 9.2.1.1 (3): As,max = 0,04 Ac."""
    return dict(As_max=0.04 * querschnitt.Ac, Ac=querschnitt.Ac,
                norm=ref("As_max"))


def robustheitsbewehrung(querschnitt, beton, stahl, faktor_z=0.9):
    """
    Mindestbewehrung fuer duktiles Bauteilverhalten (Robustheitsbewehrung).
    DIN EN 1992-1-1/NA, NDP zu 9.2.1.1 (1): die Bewehrung muss das mit fctm
    berechnete Rissmoment aufnehmen koennen:

        As,rob = Mcr / (fyk z)      mit  Mcr = fctm W

    (Es wird fyk und nicht fyd angesetzt, da es sich um einen
    Robustheitsnachweis handelt.)
    """
    W = querschnitt.widerstandsmoment_unten()
    Mcr = beton.fctm * W
    z = faktor_z * querschnitt.d
    return dict(As_rob=Mcr / (stahl.fyk * z), Mcr=Mcr / 1.0e6, W=W, z=z,
                norm=ref("robustheit"))


# ---------------------------------------------------------------------------
# 8.4  Verankerung
# ---------------------------------------------------------------------------
def verbundspannung(beton, guter_verbund=True, phi=16.0):
    """
    DIN EN 1992-1-1, 8.4.2 (2), Gl. (8.2):

        fbd = 2,25 eta1 eta2 fctd

    eta1 = 1,0 (gute Verbundbedingungen) / 0,7 (sonst)
    eta2 = 1,0 fuer phi <= 32 mm ; (132 - phi)/100 fuer phi > 32 mm
    """
    eta1 = 1.0 if guter_verbund else 0.7
    eta2 = 1.0 if phi <= 32.0 else (132.0 - phi) / 100.0
    return dict(fbd=2.25 * eta1 * eta2 * beton.fctd, eta1=eta1, eta2=eta2,
                fctd=beton.fctd, norm=ref("fbd"))


def verankerungslaenge(beton, stahl, phi, sigma_sd=None, guter_verbund=True,
                       alpha1=1.0, alpha2=1.0, alpha3=1.0, alpha4=1.0,
                       alpha5=1.0, zug=True, As_erf=None, As_vorh=None):
    """
    Verankerungslaenge.
    DIN EN 1992-1-1, 8.4.3 (2), Gl. (8.3) und 8.4.4 (1), Gl. (8.4)/(8.6):

        lb,rqd = (phi/4)(sigma_sd/fbd)
        lbd    = alpha1 alpha2 alpha3 alpha4 alpha5 lb,rqd >= lb,min
        lb,min = max(0,3 alpha1 lb,rqd ; 10 phi ; 100 mm)   (Zug)
        lb,min = max(0,6 lb,rqd ; 10 phi ; 100 mm)          (Druck)

    alpha1 = 1,0 gerade / 0,7 mit Haken (Querdeckung >= 3 phi)
    alpha2 = 1 - 0,15 (cd - phi)/phi,  0,7 <= alpha2 <= 1,0
    """
    ss = stahl.fyd if sigma_sd is None else sigma_sd
    if As_erf and As_vorh:
        ss = ss * As_erf / As_vorh
    fb = verbundspannung(beton, guter_verbund, phi)
    lb_rqd = (phi / 4.0) * ss / fb["fbd"]
    prod = alpha1 * alpha2 * alpha3 * alpha4 * alpha5
    if zug:
        lb_min = max(0.3 * alpha1 * lb_rqd, 10.0 * phi, 100.0)
    else:
        lb_min = max(0.6 * lb_rqd, 10.0 * phi, 100.0)
    return dict(lbd=max(prod * lb_rqd, lb_min), lb_rqd=lb_rqd, lb_min=lb_min,
                fbd=fb["fbd"], sigma_sd=ss, alphas=prod, eta1=fb["eta1"],
                eta2=fb["eta2"],
                normen=[ref("fbd"), ref("lb_rqd"), ref("lbd")])


def uebergreifungslaenge(beton, stahl, phi, sigma_sd=None, alpha6=1.4,
                         guter_verbund=True, **kw):
    """
    DIN EN 1992-1-1, 8.7.3 (1), Gl. (8.10):
        l0 = alpha1 alpha2 alpha3 alpha5 alpha6 lb,rqd >= l0,min
        l0,min = max(0,3 alpha6 lb,rqd ; 15 phi ; 200 mm)

    alpha6 = sqrt(rho1/25) mit 1,0 <= alpha6 <= 1,5
    (rho1 = Anteil der im gleichen Schnitt gestossenen Staebe;
     1,4 bei 50 %, 1,5 bei 100 %).
    """
    an = verankerungslaenge(beton, stahl, phi, sigma_sd, guter_verbund, **kw)
    l0 = max(alpha6 * an["lbd"], 0.3 * alpha6 * an["lb_rqd"], 15.0 * phi, 200.0)
    return dict(l0=l0, alpha6=alpha6, lb_rqd=an["lb_rqd"], fbd=an["fbd"],
                norm=ref("l0"))


def verankerung_endauflager(beton, stahl, phi, F_Ed_kN, As_vorh,
                            guter_verbund=True, direkte_lagerung=True):
    """
    Verankerung am Endauflager.
    DIN EN 1992-1-1, 9.2.1.4: zu verankern ist F_Ed (Gl. 9.3) mit
    alpha1 = 1,0 (gerade Staebe); bei direkter Lagerung genuegen 2/3 lbd.
    """
    sigma = min(F_Ed_kN * 1.0e3 / As_vorh if As_vorh > 0 else stahl.fyd,
                stahl.fyd)
    an = verankerungslaenge(beton, stahl, phi, sigma_sd=sigma,
                            guter_verbund=guter_verbund)
    return dict(l_verankerung=(2.0 / 3.0) * an["lbd"] if direkte_lagerung
                else an["lbd"], lbd=an["lbd"], sigma_sd=sigma,
                lb_rqd=an["lb_rqd"], F_Ed=F_Ed_kN,
                normen=[ref("auflagerkraft"), ref("lbd")])


# ---------------------------------------------------------------------------
# 8.2 / 9.2.2  Stababstaende
# ---------------------------------------------------------------------------
def lichter_mindestabstand(phi, d_g=16.0):
    """DIN EN 1992-1-1, 8.2 (2): s >= max(phi ; d_g + 5 mm ; 20 mm)."""
    return dict(s_min=max(phi, d_g + 5.0, 20.0), norm=ref("stababstand"))


def platznachweis(b, n_staebe, phi, c_nom_w, phi_buegel, d_g=16.0, n_lagen=1):
    """
    Prueft, ob die Staebe in die verfuegbare Breite passen.
    DIN EN 1992-1-1, 8.2 (2) + 4.4.1:

        b_erf = 2 c_nom,w + 2 phi_buegel + n phi + (n-1) s_min
    """
    s_min = max(phi, d_g + 5.0, 20.0)
    n_lage = math.ceil(n_staebe / n_lagen)
    b_erf = 2.0 * c_nom_w + 2.0 * phi_buegel + n_lage * phi + (n_lage - 1) * s_min
    s_vorh = ((b - 2.0 * c_nom_w - 2.0 * phi_buegel - n_lage * phi)
              / (n_lage - 1)) if n_lage > 1 else float("inf")
    return dict(ok=b_erf <= b + 1e-6, b_erf=b_erf, b=b, s_min=s_min,
                s_vorh=s_vorh, n_je_lage=n_lage, n_lagen=n_lagen,
                normen=[ref("stababstand"), ref("c_nom")])


def d1_schaetzung(c_nom_w, phi_buegel, phi, n_lagen=1, d_g=16.0):
    """
    Randabstand d1 (Zugrand -> Schwerpunkt der Bewehrung).
    Eine Lage:   d1 = c_nom,w + phi_buegel + phi/2
    Zwei Lagen:  zuzueglich der halben lichten Lagenhoehe.
    """
    d1 = c_nom_w + phi_buegel + phi / 2.0
    if n_lagen > 1:
        s_v = max(phi, d_g + 5.0, 20.0)
        d1 += (n_lagen - 1) * (phi + s_v) / 2.0
    return d1


# ---------------------------------------------------------------------------
# Zugkraftdeckung
# ---------------------------------------------------------------------------
def zugkraft(x_m, M_kNm, z_mm, a_l_mm, N_Ed=0.0, F_torsion_kN=0.0):
    """
    Zugkraft in der Laengsbewehrung mit Versatzmass.
    DIN EN 1992-1-1, 9.2.1.3 (2) (Versatzmassregel):

        Fs(x) = |M(x versetzt)| / z + N_Ed + F_Torsion

    Die Momentenlinie wird um a_l in die unguenstige Richtung verschoben:
    an jeder Stelle wird das groesste |M| im Fenster [x-a_l ; x+a_l] angesetzt.
    Die Torsionslaengskraft ist ueber die Laenge konstant und wird addiert
    (EC2 6.3.2 (3)).
    """
    x = np.asarray(x_m, float)
    M = np.asarray(M_kNm, float)
    a = a_l_mm / 1000.0
    M_links = np.interp(x - a, x, M, left=M[0], right=M[-1])
    M_rechts = np.interp(x + a, x, M, left=M[0], right=M[-1])
    M_vers = np.maximum.reduce([np.abs(M), np.abs(M_links), np.abs(M_rechts)])
    Fs = M_vers * 1.0e6 / z_mm / 1.0e3 + N_Ed + F_torsion_kN    # [kN]
    return dict(x=x, Fs=Fs, M_versetzt=M_vers, a_l=a_l_mm, z=z_mm,
                F_torsion=F_torsion_kN, norm=ref("versatzmass"))
