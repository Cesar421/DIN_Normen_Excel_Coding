# -*- coding: utf-8 -*-
"""
QUERKRAFTBEMESSUNG.

Norm: DIN EN 1992-1-1:2011-01, Abschnitt 6.2
      + DIN EN 1992-1-1/NA:2013-04 (NDP zu 6.2.2 und 6.2.3)

Ablauf
  1. Innerer Hebelarm z                   NA, NDP zu 6.2.3 (1)
  2. V_Rd,c (ohne Querkraftbewehrung)     Gl. (6.2a)/(6.2b) + NA Gl. (6.3aDE)
  3. V_Ed <= V_Rd,c  ->  nur Mindestbuegel  NA NDP zu 9.2.2 (5)
  4. V_Ed >  V_Rd,c  ->  Fachwerkmodell:
        - Druckstrebenneigung cot(theta)  NA Gl. (6.7aDE)/(6.7bDE)
        - Druckstrebe V_Rd,max            Gl. (6.9)
        - Buegelbewehrung Asw/s           Gl. (6.8)
  5. Groesste Buegelabstaende              NA Tab. NA.9.1
  6. Versatzmass a_l                       Gl. (9.2)

Einheiten: N, mm (V in kN).
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .normen import ref
from .baustoffe import stabflaeche

# uebliche Verlegemasse [mm]
ABSTAENDE = [50, 60, 70, 75, 80, 90, 100, 110, 125, 130, 150, 160, 175,
             180, 200, 220, 225, 250, 275, 300]


def innerer_hebelarm(d, c_nom_l):
    """
    Innerer Hebelarm fuer den Querkraftnachweis.
    DIN EN 1992-1-1/NA, NDP zu 6.2.3 (1):

        z = min(0,9 d ; d - 2 c_v,l ; d - c_v,l - 30 mm)

    c_v,l = Nennmass der Betondeckung der Laengsbewehrung in der Druckzone.
    """
    z = min(0.9 * d, d - 2.0 * c_nom_l, d - c_nom_l - 30.0)
    return dict(z=max(z, 0.5 * d), d=d, c_v_l=c_nom_l, z_09d=0.9 * d,
                z_2c=d - 2.0 * c_nom_l, z_c30=d - c_nom_l - 30.0,
                norm=ref("z_innen"))


def V_Rd_c(beton, bw, d, As_l, N_Ed=0.0, Ac=None):
    """
    Querkrafttragfaehigkeit ohne Querkraftbewehrung.
    DIN EN 1992-1-1, 6.2.2 (1), Gl. (6.2a)/(6.2b) mit den deutschen NDP:

        C_Rd,c = 0,15/gamma_C          (NA, NDP zu 6.2.2 (1))
        k_1    = 0,12
        v_min  = (kappa_1/gamma_C) k^(3/2) fck^(1/2)   (NA, Gl. (6.3aDE))
                 kappa_1 = 0,0525  fuer d <= 600 mm
                 kappa_1 = 0,0375  fuer d >= 800 mm  (dazwischen linear)

    As_l : verankerte Laengszugbewehrung [mm2]
    N_Ed : Laengskraft [kN], DRUCK NEGATIV
    """
    fck, gc = beton.fck, beton.gamma_c
    k = min(1.0 + math.sqrt(200.0 / d), 2.0)
    rho_l = min(As_l / (bw * d), 0.02) if As_l > 0 else 0.0

    # sigma_cp nach EC2: DRUCK POSITIV, <= 0,2 fcd
    sigma_cp = 0.0
    if Ac and N_Ed < 0.0:
        sigma_cp = min(-N_Ed * 1.0e3 / Ac, 0.2 * beton.fcd)
    elif Ac and N_Ed > 0.0:
        sigma_cp = -N_Ed * 1.0e3 / Ac

    C_Rdc, k1 = 0.15 / gc, 0.12
    if d <= 600.0:
        kappa1 = 0.0525
    elif d >= 800.0:
        kappa1 = 0.0375
    else:
        kappa1 = 0.0525 + (0.0375 - 0.0525) * (d - 600.0) / 200.0
    v_min = (kappa1 / gc) * k ** 1.5 * math.sqrt(fck)

    v_a = C_Rdc * k * (100.0 * rho_l * fck) ** (1.0 / 3.0) + k1 * sigma_cp
    v_b = v_min + k1 * sigma_cp
    v = max(v_a, v_b)
    return dict(V_Rdc=v * bw * d / 1.0e3, k=k, rho_l=rho_l, sigma_cp=sigma_cp,
                C_Rdc=C_Rdc, k1=k1, kappa1=kappa1, v_min=v_min, v_a=v_a, v_b=v_b,
                massgebend="Gl.(6.2a)" if v_a >= v_b else "Gl.(6.2b)",
                normen=[ref("VRdc"), ref("vmin")])


def cot_theta_NA(beton, bw, z, V_Ed, N_Ed=0.0, Ac=None, c_beiwert=0.5):
    """
    Druckstrebenneigung nach DIN EN 1992-1-1/NA, NDP zu 6.2.3 (2):

        cot(theta) = (1,2 + 1,4 sigma_cd/fcd) / (1 - V_Rd,cc/V_Ed)
        V_Rd,cc    = c * 0,48 * fck^(1/3) * (1 - 1,2 sigma_cd/fcd) * bw * z
        mit c = 0,5   und   1,0 <= cot(theta) <= 3,0  (Normalbeton)

    sigma_cd = N_Ed/Ac : Betonlaengsspannung im Schwerpunkt, DRUCK NEGATIV.
    """
    fck, fcd = beton.fck, beton.fcd
    sigma_cd = (N_Ed * 1.0e3 / Ac) if Ac else 0.0
    V_Rdcc = c_beiwert * 0.48 * fck ** (1.0 / 3.0) \
        * (1.0 - 1.2 * sigma_cd / fcd) * bw * z / 1.0e3
    zaehler = 1.2 + 1.4 * sigma_cd / fcd
    nenner = 1.0 - V_Rdcc / max(abs(V_Ed), 1e-9)
    if nenner <= 1e-6:
        cot, hinweis = 3.0, "cot(theta) = 3,0 (V_Ed <= V_Rd,cc)"
    else:
        cot, hinweis = zaehler / nenner, ""
    cot_roh = cot
    cot = min(max(cot, 1.0), 3.0)
    return dict(cot_theta=cot, theta_grad=math.degrees(math.atan(1.0 / cot)),
                V_Rdcc=V_Rdcc, sigma_cd=sigma_cd, cot_roh=cot_roh,
                hinweis=hinweis, normen=[ref("cot_theta"), ref("VRdcc")])


def V_Rd_max(beton, bw, z, cot_theta, alpha_cw=1.0):
    """
    Druckstrebentragfaehigkeit.
    DIN EN 1992-1-1, 6.2.3 (3), Gl. (6.9); lotrechte Buegel (alpha = 90 Grad):

        V_Rd,max = alpha_cw bw z nu_1 fcd / (cot(theta) + tan(theta))

    nu_1 = 0,75 nu_2  mit nu_2 = 1,0 fuer <= C50/60  (NA, NDP zu 6.2.3 (3)).
    """
    nu2 = 1.0 if beton.fck <= 50.0 else min(1.1 - beton.fck / 500.0, 1.0)
    nu1 = 0.75 * nu2
    V = alpha_cw * bw * z * nu1 * beton.fcd / (cot_theta + 1.0 / cot_theta) / 1.0e3
    return dict(V_Rdmax=V, nu1=nu1, nu2=nu2, alpha_cw=alpha_cw,
                normen=[ref("VRdmax"), ref("nu1")])


def asw_erforderlich(V_Ed, z, fywd, cot_theta):
    """
    Erforderliche Querkraftbewehrung aus Gl. (6.8):
        Asw/s = V_Ed / (z fywd cot(theta))     [mm2/mm], SUMME aller Schenkel
    """
    return abs(V_Ed) * 1.0e3 / (z * fywd * cot_theta)


def asw_mindest(beton, stahl, bw, alpha=90.0):
    """
    Mindestquerkraftbewehrung.  DIN EN 1992-1-1/NA, NDP zu 9.2.2 (5):
        rho_w,min = 0,16 fctm / fyk        (Balken)
        Asw/s|min = rho_w,min bw sin(alpha)
    """
    rho_min = 0.16 * beton.fctm / stahl.fyk
    return dict(rho_w_min=rho_min,
                asw_min=rho_min * bw * math.sin(math.radians(alpha)),
                norm=ref("rho_w_min"))


def groesster_buegelabstand(V_Ed, V_Rdmax, h, fck=30.0):
    """
    Groesste Buegelabstaende.
    DIN EN 1992-1-1/NA, NDP zu 9.2.2 (6)/(8), Tab. NA.9.1  (<= C50/60):

        V_Ed <= 0,30 V_Rd,max  ->  s <= 0,7 h <= 300 mm ; s_q <= h <= 800 mm
        <= 0,60 V_Rd,max       ->  s <= 0,5 h <= 300 mm ; s_q <= h <= 600 mm
        >  0,60 V_Rd,max       ->  s <= 0,25 h <= 200 mm; s_q <= h <= 600 mm
    """
    r = abs(V_Ed) / V_Rdmax if V_Rdmax > 0 else 1.0
    if r <= 0.30:
        s_l, s_q, zeile = min(0.7 * h, 300.0), min(h, 800.0), "1"
    elif r <= 0.60:
        s_l, s_q, zeile = min(0.5 * h, 300.0), min(h, 600.0), "2"
    else:
        s_l, s_q, zeile = min(0.25 * h, 200.0), min(h, 600.0), "3"
    if fck > 50.0 and r <= 0.30:
        s_l = min(s_l, 200.0)
    return dict(s_max=s_l, s_max_q=s_q, ausnutzung=r, tabellenzeile=zeile,
                norm=ref("s_max"))


def versatzmass(z, cot_theta, cot_alpha=0.0):
    """
    Versatzmass der Zugkraftlinie.
    DIN EN 1992-1-1, 9.2.1.3 (2), Gl. (9.2):  a_l = z (cot(theta) - cot(alpha))/2
    """
    return dict(a_l=z * (cot_theta - cot_alpha) / 2.0, norm=ref("versatzmass"))


# ---------------------------------------------------------------------------
# Bemessung ueber die gesamte Traegerlaenge
# ---------------------------------------------------------------------------
@dataclass
class ErgebnisQuerkraft:
    x: np.ndarray = None            # [m]
    V_Ed: np.ndarray = None         # [kN] Einhuellende (Betrag)
    V_Rdc: float = 0.0              # [kN]
    V_Rdmax: float = 0.0            # [kN]
    V_Rdcc: float = 0.0             # [kN]
    cot_theta: float = 1.0
    theta: float = 45.0
    z: float = 0.0                  # [mm]
    asw_erf: np.ndarray = None      # [mm2/m] Summe aller Schenkel (Querkraft)
    asw_min: float = 0.0            # [mm2/m]
    asw_torsion: float = 0.0        # [mm2/m] je AUSSENSCHENKEL aus Torsion
    asw_schenkel: np.ndarray = None  # [mm2/m] massgebender Aussenschenkel
    s_erf: np.ndarray = None        # [mm] erforderlicher Abstand
    s_max: float = 300.0
    s_max_q: float = 800.0
    bereiche: list = field(default_factory=list)
    phi_buegel: float = 8.0
    n_schenkel: int = 2
    a_l: float = 0.0                # [mm] Versatzmass
    ok: bool = True
    ausnutzung_druckstrebe: float = 0.0
    hinweise: list = field(default_factory=list)
    detail: dict = field(default_factory=dict)
    normen: list = field(default_factory=list)


def bemessung_querkraft(querschnitt, beton, stahl, x_m, V_einh_kN, As_l,
                        c_nom_l, phi_buegel=8.0, n_schenkel=2, N_Ed=0.0,
                        V_bem_kN=None, asw_torsion=0.0, s_max_torsion=None,
                        min_bereichslaenge=0.5):
    """
    Querkraftbemessung des gesamten Traegers.

    Parameter
    ---------
    x_m : array          Abszissen [m]
    V_einh_kN : array    Einhuellende von |V_Ed| [kN]
    As_l : float         vorhandene Laengszugbewehrung [mm2]
    c_nom_l : float      Nennmass der Betondeckung der Laengsbewehrung [mm]
    V_bem_kN : float     Querkraft, die cot(theta) bestimmt (Standard: max |V|)
    asw_torsion : float  [mm2/m] Torsionsbuegel JE AUSSENSCHENKEL (EC2 6.3.2)
    s_max_torsion : float  zusaetzliche Abstandsbegrenzung aus EC2 9.2.3
    """
    bw, d, h = querschnitt.bw, querschnitt.d, querschnitt.h
    fywd = stahl.fyd
    r = ErgebnisQuerkraft(x=np.asarray(x_m, float),
                          V_Ed=np.abs(np.asarray(V_einh_kN, float)),
                          phi_buegel=phi_buegel, n_schenkel=n_schenkel,
                          asw_torsion=asw_torsion)
    r.normen = [ref("querkraft_allg"), ref("VRdc"), ref("vmin"), ref("z_innen"),
                ref("cot_theta"), ref("VRdcc"), ref("VRdmax"), ref("nu1"),
                ref("VRds"), ref("rho_w_min"), ref("s_max"), ref("versatzmass")]

    hb = innerer_hebelarm(d, c_nom_l)
    r.z = hb["z"]

    vc = V_Rd_c(beton, bw, d, As_l, N_Ed, querschnitt.Ac)
    r.V_Rdc = vc["V_Rdc"]

    V_bem = float(np.max(r.V_Ed)) if V_bem_kN is None else abs(V_bem_kN)
    ct = cot_theta_NA(beton, bw, r.z, V_bem, N_Ed, querschnitt.Ac)
    r.cot_theta, r.theta, r.V_Rdcc = ct["cot_theta"], ct["theta_grad"], ct["V_Rdcc"]

    vm = V_Rd_max(beton, bw, r.z, r.cot_theta)
    r.V_Rdmax = vm["V_Rdmax"]
    r.ausnutzung_druckstrebe = V_bem / r.V_Rdmax if r.V_Rdmax > 0 else 9.99
    if r.ausnutzung_druckstrebe > 1.0:
        r.ok = False
        r.hinweise.append(
            "V_Ed = {:.0f} kN > V_Rd,max = {:.0f} kN: Druckstrebenversagen. "
            "bw, h oder Betonfestigkeitsklasse erhoehen "
            "[EC2 6.2.3 (3), Gl. (6.9)].".format(V_bem, r.V_Rdmax))

    am = asw_mindest(beton, stahl, bw)
    r.asw_min = am["asw_min"] * 1000.0                       # mm2/m
    asw = np.array([asw_erforderlich(v, r.z, fywd, r.cot_theta)
                    for v in r.V_Ed]) * 1000.0               # mm2/m
    # Wo V_Ed <= V_Rd,c genuegt die Mindestbewehrung (EC2 6.2.1 (4) + NA 9.2.2)
    asw = np.where(r.V_Ed <= r.V_Rdc, r.asw_min, np.maximum(asw, r.asw_min))
    r.asw_erf = asw

    # --- massgebende Flaeche JE AUSSENSCHENKEL:
    #     Querkraft verteilt sich auf alle Schenkel, Torsion belastet nur die
    #     beiden Aussenschenkel des geschlossenen Buegels (EC2 6.3.2 (3))
    r.asw_schenkel = asw / n_schenkel + asw_torsion

    sm = groesster_buegelabstand(V_bem, r.V_Rdmax, h, beton.fck)
    r.s_max, r.s_max_q = sm["s_max"], sm["s_max_q"]
    if s_max_torsion is not None:
        r.s_max = min(r.s_max, s_max_torsion)

    A_stab = stabflaeche(phi_buegel)                          # ein Schenkel
    r.s_erf = np.minimum(A_stab * 1000.0 / np.maximum(r.asw_schenkel, 1e-9),
                         r.s_max)
    r.bereiche = _bereiche_bilden(r.x, r.s_erf, min_bereichslaenge)
    r.a_l = versatzmass(r.z, r.cot_theta)["a_l"]

    r.detail = dict(V_Rdc=vc, cot=ct, V_Rdmax=vm, asw_min=am, s_max=sm,
                    V_bem=V_bem, A_stab=A_stab, fywd=fywd,
                    A_buegel=n_schenkel * A_stab)
    if float(np.max(r.V_Ed)) <= r.V_Rdc and asw_torsion <= 0:
        r.hinweise.append(
            "V_Ed,max = {:.0f} kN <= V_Rd,c = {:.0f} kN: rechnerisch keine "
            "Querkraftbewehrung erforderlich; es wird die Mindestbewehrung "
            "angeordnet [EC2 6.2.2 (1) + NA NDP zu 9.2.2 (5)]."
            .format(float(np.max(r.V_Ed)), r.V_Rdc))
    return r


def _bereiche_bilden(x, s_erf, l_min=0.5):
    """
    Fasst den erforderlichen Buegelabstand zu Verlegebereichen mit
    genormten Abstaenden zusammen; in jedem Bereich gilt stets s <= s_erf.
    """
    s_verf = np.array(ABSTAENDE, float)
    gewaehlt = np.array([s_verf[s_verf <= sr].max() if (s_verf <= sr).any()
                         else s_verf.min() for sr in s_erf])
    bereiche, anf = [], 0
    for i in range(1, len(x) + 1):
        if i == len(x) or gewaehlt[i] != gewaehlt[anf]:
            bereiche.append([float(x[anf]),
                             float(x[i - 1] if i == len(x) else x[i]),
                             float(gewaehlt[anf])])
            anf = i

    # kurze Bereiche mit dem strengeren Nachbarn verschmelzen
    aenderung = True
    while aenderung and len(bereiche) > 1:
        aenderung = False
        for i, b in enumerate(bereiche):
            if b[1] - b[0] < l_min:
                if i == 0:
                    bereiche[1][0] = b[0]
                    bereiche[1][2] = min(bereiche[1][2], b[2])
                elif i == len(bereiche) - 1:
                    bereiche[-2][1] = b[1]
                    bereiche[-2][2] = min(bereiche[-2][2], b[2])
                else:
                    j = i - 1 if bereiche[i - 1][2] <= bereiche[i + 1][2] else i + 1
                    bereiche[j][2] = min(bereiche[j][2], b[2])
                    if j < i:
                        bereiche[j][1] = b[1]
                    else:
                        bereiche[j][0] = b[0]
                bereiche.pop(i)
                aenderung = True
                break

    ber = _verschmelzen(bereiche)

    # SYMMETRIERUNG: bei symmetrischer Beanspruchung muss auch die Verlegung
    # symmetrisch sein (das Verschmelzen laeuft von links nach rechts).
    x0, x1 = float(x[0]), float(x[-1])
    if np.allclose(s_erf, np.interp(x0 + x1 - x, x, s_erf), rtol=0.02, atol=1.0):
        mitte = 0.5 * (x0 + x1)
        links = [b[:] for b in ber if b[0] < mitte - 1e-9]
        if links:
            links[-1][1] = mitte
            rechts = [[x0 + x1 - b[1], x0 + x1 - b[0], b[2]]
                      for b in reversed(links)]
            ber = _verschmelzen(links + rechts)

    # SICHERHEITSDURCHGANG: jeder Bereich erhaelt den groessten genormten
    # Abstand, der im GESCHLOSSENEN Bereichsintervall s <= s_erf einhaelt.
    for b in ber:
        m = (x >= b[0] - 1e-9) & (x <= b[1] + 1e-9)
        if not np.any(m):
            continue
        grenz = float(np.min(s_erf[m]))
        kand = s_verf[s_verf <= grenz + 1e-9]
        b[2] = float(kand.max()) if kand.size else float(s_verf.min())

    ber = _verschmelzen(ber)
    return [dict(x1=b[0], x2=b[1], s=b[2],
                 n=max(1, int(math.ceil((b[1] - b[0]) * 1000.0 / b[2]))))
            for b in ber]


def _verschmelzen(bereiche):
    """Verbindet benachbarte Bereiche mit gleichem Abstand."""
    out = [list(bereiche[0])]
    for b in bereiche[1:]:
        if abs(b[2] - out[-1][2]) < 1e-9:
            out[-1][1] = b[1]
        else:
            out.append(list(b))
    return out


def verankerungskraft_endauflager(V_Ed, a_l, z, N_Ed=0.0):
    """
    Zu verankernde Zugkraft am Endauflager.
    DIN EN 1992-1-1, 9.2.1.4 (2), Gl. (9.3):  F_Ed = |V_Ed| a_l/z + N_Ed
    """
    return dict(F_Ed=abs(V_Ed) * a_l / z + N_Ed, norm=ref("auflagerkraft"))
