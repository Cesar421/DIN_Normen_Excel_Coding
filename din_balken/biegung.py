# -*- coding: utf-8 -*-
"""
BIEGEBEMESSUNG im Grenzzustand der Tragfaehigkeit.

Norm: DIN EN 1992-1-1:2011-01, Abschnitt 6.1 (Biegung mit Laengskraft)
      + DIN EN 1992-1-1/NA:2013-04

Bemessungsannahmen (EC2, 6.1 (2)P):
  - Ebenbleiben der Querschnitte (Bernoulli).
  - Vollstaendiger Verbund zwischen Beton und Stahl.
  - Zugfestigkeit des Betons wird nicht angesetzt.
  - Beton: Parabel-Rechteck-Diagramm, EC2 3.1.7 (1), Bild 3.3.
  - Betonstahl: bilinear mit waagerechtem oberem Ast, EC2 3.2.7 (2) b).
  - Grenzdehnungen (EC2, 6.1 (2)P, Bild 6.1 - Bemessungspunkte):
        eps_cu2 = 3,5 Promille   (Beton <= C50/60)
        eps_ud  = 25 Promille    (Stahl, NA NDP zu 3.2.7 (2))

Duktilitaetsgrenze (NA, NDP zu 5.5 (4)):
        xu/d <= 0,45  fuer <= C50/60      (delta = 1,0, ohne Umlagerung)
        xu/d <= 0,35  fuer >= C55/67
Wird die Grenze ueberschritten, ist Druckbewehrung anzuordnen.

Interne Einheiten: N, mm, N*mm.   Schnittstelle: kN, kNm, mm, mm2.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .normen import ref


# ---------------------------------------------------------------------------
# Duktilitaetsgrenze
# ---------------------------------------------------------------------------
def xi_grenz(beton, delta=1.0):
    """
    Zulaessiges xu/d nach DIN EN 1992-1-1/NA, NDP zu 5.5 (4), Gl. (5.10a):

        delta >= k1 + k2 * xu/d      mit  k1 = 0,64 ; k2 = 0,80   (<= C50/60)
        delta >= k3 + k4 * xu/d      mit  k3 = 0,72 ; k4 = 0,80   (>= C55/67)

    delta = Umlagerungsgrad (1,0 = keine Umlagerung).
    """
    if beton.fck <= 50.0:
        k1, k2, obergrenze = 0.64, 0.80, 0.45
    else:
        k1, k2, obergrenze = 0.72, 0.80, 0.35
    xi = (delta - k1) / k2
    return dict(xi_lim=min(max(xi, 0.0), obergrenze), k1=k1, k2=k2, delta=delta,
                obergrenze=obergrenze, normen=[ref("xi_lim"), ref("umlagerung")])


# ---------------------------------------------------------------------------
# Dehnungszustand und Betondruckkraft
# ---------------------------------------------------------------------------
def dehnungszustand(x, d, beton, stahl):
    """
    Dehnungen (eps_c am staerkst gedrueckten Rand, eps_s in As1) [Promille]
    fuer eine Druckzonenhoehe x [mm], nach dem Dehnungsdiagramm von
    DIN EN 1992-1-1, 6.1 (2)P, Bild 6.1.

    x <= x_A  ->  Punkt A: der Stahl ist massgebend, eps_s = eps_ud
    x >  x_A  ->  Punkt B: der Beton ist massgebend, eps_c = eps_cu2
    """
    ecu, eud = beton.eps_cu2, stahl.eps_ud
    x_A = d * ecu / (ecu + eud)
    if x <= 1e-12:
        return 0.0, eud, "A"
    if x <= x_A:
        eps_c = eud * x / max(d - x, 1e-9)
        return min(eps_c, ecu), eud, "A"
    return ecu, ecu * (d - x) / x, "B"


def betondruckkraft(querschnitt, beton, x, eps_c_rand, n=300):
    """
    Resultierende Betondruckkraft.

    Integriert das Parabel-Rechteck-Diagramm (EC2, 3.1.7 (1)) ueber die
    tatsaechliche Druckzone (Rechteck oder Plattenbalken).

    Rueckgabe
    ---------
    Fc : float      Resultierende [N] (positiv = Druck)
    yc : float      Lage der Resultierenden ab oberem Rand [mm]
    """
    x_eff = min(x, querschnitt.h)
    if x_eff <= 1e-9 or eps_c_rand <= 1e-12:
        return 0.0, 0.0
    ys = set(np.linspace(0.0, x_eff, n + 1).tolist())
    if 0.0 < querschnitt.hf < x_eff:            # Sprung der Breite am Plattenrand
        ys.update((querschnitt.hf - 1e-9, querschnitt.hf + 1e-9))
    ys = np.array(sorted(ys))
    eps = eps_c_rand * (1.0 - ys / x)
    sig = np.array([beton.sigma_c(e) for e in eps])
    br = np.array([querschnitt.breite(y) for y in ys])
    f = sig * br                                 # [N/mm]
    Fc = float(np.trapz(f, ys))
    if Fc <= 1e-9:
        return 0.0, 0.0
    return Fc, float(np.trapz(f * ys, ys) / Fc)


def _moment_beton(querschnitt, beton, stahl, x, d):
    """Moment der Betondruckkraft um As1 [N*mm]."""
    eps_c, eps_s, pkt = dehnungszustand(x, d, beton, stahl)
    Fc, yc = betondruckkraft(querschnitt, beton, x, eps_c)
    return Fc * (d - yc), Fc, yc, eps_c, eps_s, pkt


# ---------------------------------------------------------------------------
# Ergebnis
# ---------------------------------------------------------------------------
@dataclass
class ErgebnisBiegung:
    M_Ed: float = 0.0          # [kNm]  Bemessungsmoment
    N_Ed: float = 0.0          # [kN]   Bemessungslaengskraft (Druck negativ)
    M_Eds: float = 0.0         # [kNm]  auf As1 bezogenes Moment
    x: float = 0.0             # [mm]   Druckzonenhoehe
    xi: float = 0.0            # x/d
    z: float = 0.0             # [mm]   innerer Hebelarm
    zeta: float = 0.0          # z/d
    eps_c: float = 0.0         # [Promille] gedrueckter Rand
    eps_s1: float = 0.0        # [Promille] Zugbewehrung
    eps_s2: float = 0.0        # [Promille] Druckbewehrung
    sigma_s1: float = 0.0      # [N/mm2]
    sigma_s2: float = 0.0      # [N/mm2]
    Fc: float = 0.0            # [kN]
    mu_Eds: float = 0.0        # bezogenes Moment
    omega: float = 0.0         # mechanischer Bewehrungsgrad
    As1: float = 0.0           # [mm2] erforderliche Zugbewehrung
    As2: float = 0.0           # [mm2] erforderliche Druckbewehrung
    xi_lim: float = 0.45
    mu_lim: float = 0.0
    punkt: str = "B"
    b_bezug: float = 0.0       # in mu/omega verwendete Bezugsbreite [mm]
    platte_reicht: bool = True  # (Plattenbalken) Druckzone innerhalb der Platte
    ok: bool = True
    hinweise: list = field(default_factory=list)
    normen: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Bemessung
# ---------------------------------------------------------------------------
def bemessung_biegung(querschnitt, beton, stahl, M_Ed, N_Ed=0.0, delta=1.0,
                      druckbewehrung_zulassen=True, tol=1e-8):
    """
    Biegebemessung mit oder ohne Laengskraft.   DIN EN 1992-1-1, 6.1.

    Parameter
    ---------
    querschnitt : Querschnitt  so orientiert, dass As1 auf der Zugseite liegt
    M_Ed : float               Bemessungsmoment [kNm], >= 0
    N_Ed : float               Laengskraft [kN]; DRUCK NEGATIV
    delta : float              Umlagerungsgrad (1,0 = keine Umlagerung)
    """
    d, d2 = querschnitt.d, querschnitt.d2
    fcd, fyd = beton.fcd, stahl.fyd

    r = ErgebnisBiegung(M_Ed=M_Ed, N_Ed=N_Ed)
    r.normen = [ref("biegung"), ref("eps_grenzen"), ref("sigma_eps_c"),
                ref("betonstahl_ec2"), ref("eps_ud"), ref("xi_lim")]

    # --- 1) auf den Schwerpunkt der Zugbewehrung bezogenes Moment ---------
    #        M_Eds = M_Ed - N_Ed * z_s1        (EC2, 6.1)
    M_Eds_Nmm = M_Ed * 1.0e6 - (N_Ed * 1.0e3) * querschnitt.z_s1
    r.M_Eds = M_Eds_Nmm / 1.0e6
    if M_Eds_Nmm <= 0.0:
        r.hinweise.append(
            "M_Eds <= 0: keine Zugbewehrung aus Biegung erforderlich; "
            "als Druckglied nachweisen (EC2, 6.1).")
        r.As1 = 0.0
        return r

    # --- 2) Duktilitaetsgrenze --------------------------------------------
    grenz = xi_grenz(beton, delta)
    r.xi_lim = grenz["xi_lim"]
    x_lim = r.xi_lim * d

    b_bezug = querschnitt.b_eff \
        if (querschnitt.typ == "plattenbalken" and querschnitt.platte_gedrueckt) \
        else querschnitt.b
    r.b_bezug = b_bezug
    r.mu_Eds = M_Eds_Nmm / (b_bezug * d ** 2 * fcd)

    Mc_lim, Fc_lim, yc_lim, ec_lim, es_lim, pkt_lim = \
        _moment_beton(querschnitt, beton, stahl, x_lim, d)
    r.mu_lim = Mc_lim / (b_bezug * d ** 2 * fcd)

    # --- 3) Fall A: einfache Bewehrung (xu <= xu,lim) ---------------------
    if M_Eds_Nmm <= Mc_lim * (1.0 + 1e-9):
        lo, hi = 1e-6 * d, x_lim
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            Mc, *_ = _moment_beton(querschnitt, beton, stahl, mid, d)
            if Mc < M_Eds_Nmm:
                lo = mid
            else:
                hi = mid
            if hi - lo < tol * d:
                break
        x = 0.5 * (lo + hi)
        Mc, Fc, yc, eps_c, eps_s, pkt = \
            _moment_beton(querschnitt, beton, stahl, x, d)
        sigma_s1 = stahl.sigma_s(eps_s)
        As1 = (Fc + N_Ed * 1.0e3) / sigma_s1
        As2 = eps_s2 = sigma_s2 = 0.0
    else:
        # --- 4) Fall B: mit Druckbewehrung (EC2, 6.1) --------------------
        x = x_lim
        Fc, yc, eps_c, eps_s, pkt = Fc_lim, yc_lim, ec_lim, es_lim, pkt_lim
        eps_s2 = eps_c * (x - d2) / x if x > d2 else 0.0
        sigma_s2 = stahl.sigma_s(eps_s2)
        sigma_s1 = stahl.sigma_s(eps_s)
        dM = M_Eds_Nmm - Mc_lim
        if sigma_s2 <= 1.0 or d - d2 <= 0:
            r.ok = False
            r.hinweise.append(
                "d2 zu gross: die Druckbewehrung wird nicht ausgenutzt. "
                "Hoehe vergroessern oder d2 verringern.")
            sigma_s2 = max(sigma_s2, 1.0)
        As2 = dM / ((d - d2) * sigma_s2)
        As1 = (Fc + As2 * sigma_s2 + N_Ed * 1.0e3) / sigma_s1
        if not druckbewehrung_zulassen:
            r.ok = False
        r.hinweise.append(
            "mu_Eds = {:.3f} > mu_lim = {:.3f} (xu/d = {:.2f}): Druckbewehrung "
            "As2 erforderlich [EC2 6.1 + NA NDP zu 5.5 (4)]."
            .format(r.mu_Eds, r.mu_lim, r.xi_lim))
        Mc = Mc_lim

    # --- 5) Ergebnisse ----------------------------------------------------
    r.x, r.xi = x, x / d
    r.z, r.zeta = d - yc, (d - yc) / d
    r.eps_c, r.eps_s1, r.eps_s2 = eps_c, eps_s, eps_s2
    r.sigma_s1, r.sigma_s2 = sigma_s1, sigma_s2
    r.Fc = Fc / 1.0e3
    r.omega = Fc / (b_bezug * d * fcd)
    r.As1, r.As2 = max(As1, 0.0), max(As2, 0.0)
    r.punkt = pkt
    r.platte_reicht = (querschnitt.typ != "plattenbalken") \
        or (not querschnitt.platte_gedrueckt) or (x <= querschnitt.hf + 1e-9)

    if querschnitt.typ == "plattenbalken" and querschnitt.platte_gedrueckt \
            and not r.platte_reicht:
        r.hinweise.append(
            "Druckzone reicht in den Steg (x = {:.0f} mm > hf = {:.0f} mm): "
            "die Berechnung integriert den vollen Plattenbalken [EC2 6.1]."
            .format(x, querschnitt.hf))
    if r.eps_s1 < stahl.eps_yd:
        r.hinweise.append(
            "eps_s1 = {:.2f} < eps_yd = {:.2f} Promille: die Bewehrung fliesst "
            "NICHT (sproedes Versagen). Hoehe vergroessern [EC2 6.1]."
            .format(r.eps_s1, stahl.eps_yd))
    return r


# ---------------------------------------------------------------------------
# Nachweis: Momententragfaehigkeit eines bewehrten Querschnitts
# ---------------------------------------------------------------------------
def momententragfaehigkeit(querschnitt, beton, stahl, As1, As2=0.0, N_Ed=0.0,
                           tol=1e-9):
    """
    Momententragfaehigkeit M_Rd [kNm] bei bekannter Bewehrung.
    Gleichgewicht der Laengskraefte:  Fc + Fs2 + N_Ed - Fs1 = 0   (EC2, 6.1).
    """
    d, d2 = querschnitt.d, querschnitt.d2
    N = N_Ed * 1.0e3

    def rest(x):
        eps_c, eps_s, _ = dehnungszustand(x, d, beton, stahl)
        Fc, _ = betondruckkraft(querschnitt, beton, x, eps_c)
        eps_s2 = eps_c * (x - d2) / max(x, 1e-9)
        return Fc + As2 * stahl.sigma_s(eps_s2) + N - As1 * stahl.sigma_s(eps_s)

    lo, hi = 1e-6 * d, 0.9999 * d
    if rest(lo) * rest(hi) > 0:
        return dict(M_Rd=float("nan"), x=float("nan"), ok=False,
                    hinweis="Kein Gleichgewicht in 0 < x < d "
                            "(Bewehrung/Laengskraft pruefen)",
                    norm=ref("biegung"))
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if rest(lo) * rest(mid) <= 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < tol * d:
            break
    x = 0.5 * (lo + hi)
    eps_c, eps_s, pkt = dehnungszustand(x, d, beton, stahl)
    Fc, yc = betondruckkraft(querschnitt, beton, x, eps_c)
    eps_s2 = eps_c * (x - d2) / x
    Fs2 = As2 * stahl.sigma_s(eps_s2)
    M_s1 = Fc * (d - yc) + Fs2 * (d - d2)               # bezogen auf As1
    return dict(M_Rd=(M_s1 + N * querschnitt.z_s1) / 1.0e6, x=x, xi=x / d,
                z=d - yc, eps_c=eps_c, eps_s1=eps_s, eps_s2=eps_s2,
                sigma_s1=stahl.sigma_s(eps_s), sigma_s2=stahl.sigma_s(eps_s2),
                punkt=pkt, ok=True, hinweis="", norm=ref("biegung"))


# ---------------------------------------------------------------------------
# Allgemeines Bemessungsdiagramm
# ---------------------------------------------------------------------------
def bemessungsdiagramm(beton, stahl, n=400, xi_max=0.80):
    """
    Tabelle mu_Eds - omega - xi - zeta fuer den RECHTECKQUERSCHNITT, Grundlage
    der deutschen Bemessungstafeln (kd-Tafeln / omega-Tafeln).

        mu_Eds = M_Eds/(b d^2 fcd)   omega = As1 fyd/(b d fcd)   xi = x/d
    """
    from .querschnitt import Querschnitt
    q = Querschnitt(b=1000.0, h=1100.0, d1=100.0, d2=50.0)   # d = 1000 mm
    d = q.d
    xs = np.linspace(1e-4, xi_max, n) * d
    mu, om, ze, ec, es, ss = [], [], [], [], [], []
    for x in xs:
        eps_c, eps_s, _ = dehnungszustand(x, d, beton, stahl)
        Fc, yc = betondruckkraft(q, beton, x, eps_c)
        mu.append(Fc * (d - yc) / (q.b * d ** 2 * beton.fcd))
        om.append(Fc / (q.b * d * beton.fcd))
        ze.append((d - yc) / d)
        ec.append(eps_c)
        es.append(eps_s)
        ss.append(stahl.sigma_s(eps_s))
    return dict(xi=xs / d, mu=np.array(mu), omega=np.array(om),
                zeta=np.array(ze), eps_c=np.array(ec), eps_s=np.array(es),
                sigma_s=np.array(ss), normen=[ref("biegung"), ref("sigma_eps_c")])
