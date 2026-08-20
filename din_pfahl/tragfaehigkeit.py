# -*- coding: utf-8 -*-
"""
AXIALE PFAHLTRAGFAEHIGKEIT (Geotechnik).

Normen: DIN EN 1997-1:2009-09, 7.6  (Nachweis der Tragfaehigkeit von Pfaehlen)
        DIN 1054:2010-12, A 7.6      (ergaenzende deutsche Regelungen)
        EA-Pfaehle (DGGT), Abschnitt 5  (Erfahrungswerte, WSL)

Charakteristischer Widerstand (DIN EN 1997-1, 7.6.2.3, Gl. (7.8)):

    R_c,k = R_b,k + R_s,k = q_b,k A_b + sum( q_s,k,i A_s,i )

Bemessungswert (DIN 1054, A 7.6.2.2):

    R_c,d = R_b,k/gamma_b + R_s,k/gamma_s

Teilsicherheitsbeiwerte (DIN 1054, Tab. A 2.3, Nachweis GEO-2):

    Bemessungssituation   gamma_b   gamma_s   gamma_s,t (Zug)
    BS-P (staendig)        1,10      1,10       1,15
    BS-T (voruebergehend)  1,10      1,10       1,15
    BS-A (aussergewoehnl.) 1,00      1,00       1,10

Nachweis:   F_c,d <= R_c,d        (Druckpfahl)
            F_t,d <= R_t,d        (Zugpfahl, DIN EN 1997-1, 7.6.3)

WICHTIG: q_b,k und q_s,k sind EINGABEWERTE des Anwenders. Sie sind den
Tabellen 5.12 bis 5.15 der EA-Pfaehle (in Abhaengigkeit von q_c bzw. c_u)
oder Probebelastungen zu entnehmen. Dieses Modul enthaelt bewusst KEINE
Erfahrungswerte-Tabellen.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .normen_pfahl import ref

# Teilsicherheitsbeiwerte der Pfahlwiderstaende - DIN 1054, Tab. A 2.3
GAMMA_R = {
    "BS-P": dict(gamma_b=1.10, gamma_s=1.10, gamma_t=1.10, gamma_s_t=1.15),
    "BS-T": dict(gamma_b=1.10, gamma_s=1.10, gamma_t=1.10, gamma_s_t=1.15),
    "BS-A": dict(gamma_b=1.00, gamma_s=1.00, gamma_t=1.00, gamma_s_t=1.10),
}


@dataclass
class ErgebnisTragfaehigkeit:
    R_b_k: float = 0.0        # [kN] Spitzendruckwiderstand
    R_s_k: float = 0.0        # [kN] Mantelreibungswiderstand
    R_c_k: float = 0.0        # [kN] charakteristischer Gesamtwiderstand
    R_c_d: float = 0.0        # [kN] Bemessungswert (Druck)
    R_t_d: float = 0.0        # [kN] Bemessungswert (Zug)
    F_c_d: float = 0.0        # [kN] Bemessungswert der Einwirkung (Druck)
    ausnutzung: float = 0.0
    ok: bool = True
    anteile: list = field(default_factory=list)   # je Schicht
    gamma: dict = field(default_factory=dict)
    situation: str = "BS-P"
    hinweise: list = field(default_factory=list)
    normen: list = field(default_factory=list)


def axiale_tragfaehigkeit(D, schichten, q_b_k, F_c_d=0.0, F_t_d=0.0,
                          situation="BS-P", L_pfahl=None, mantel_ab_tiefe=0.0):
    """
    Axialer Pfahlwiderstand nach DIN EN 1997-1, 7.6.2.3 + DIN 1054, A 7.6.

    Parameter
    ---------
    D : float          Pfahldurchmesser [m]
    schichten : list   Bodenschichten mit q_s_k [kN/m2] und z_o / z_u [m]
    q_b_k : float      charakteristischer Spitzendruck [kN/m2]
    F_c_d : float      Bemessungswert der Druckeinwirkung [kN]
    F_t_d : float      Bemessungswert der Zugeinwirkung [kN]
    situation : str    "BS-P" | "BS-T" | "BS-A"
    mantel_ab_tiefe : float
                       Tiefe [m], ab der Mantelreibung angesetzt wird
                       (z.B. Auffuellung / Absetzbereich ohne Ansatz)
    """
    g = GAMMA_R.get(situation, GAMMA_R["BS-P"])
    U = math.pi * D                       # Mantelumfang [m]
    A_b = math.pi * D ** 2 / 4.0          # Fusspunktflaeche [m2]

    R_s_k, anteile = 0.0, []
    for s in schichten:
        z_o = max(s.z_o, mantel_ab_tiefe)
        z_u = s.z_u
        dl = max(z_u - z_o, 0.0)
        R_i = s.q_s_k * U * dl
        R_s_k += R_i
        anteile.append(dict(name=s.name, z_o=z_o, z_u=z_u, dicke=dl,
                            q_s_k=s.q_s_k, A_s=U * dl, R_s_k=R_i))
    R_b_k = q_b_k * A_b

    r = ErgebnisTragfaehigkeit(R_b_k=R_b_k, R_s_k=R_s_k, R_c_k=R_b_k + R_s_k,
                               F_c_d=F_c_d, anteile=anteile, gamma=g,
                               situation=situation)
    r.R_c_d = R_b_k / g["gamma_b"] + R_s_k / g["gamma_s"]
    r.R_t_d = R_s_k / g["gamma_s_t"]
    r.ausnutzung = F_c_d / r.R_c_d if r.R_c_d > 0 else 9.99
    r.ok = F_c_d <= r.R_c_d + 1e-6
    if F_t_d > 0:
        r.ok = r.ok and (F_t_d <= r.R_t_d + 1e-6)
    r.normen = [ref("geo_nachweis"), ref("geo_widerstand"), ref("geo_gamma"),
                ref("geo_erfahrungswerte")]
    r.hinweise.append(
        "q_b,k und q_s,k sind Eingabewerte. Sie sind den Tabellen 5.12 bis "
        "5.15 der EA-Pfaehle (abhaengig von q_c bzw. c_u) oder einer "
        "Probebelastung zu entnehmen.")
    if L_pfahl is not None and schichten:
        z_ende = max(s.z_u for s in schichten)
        if z_ende < L_pfahl - 1e-6:
            r.hinweise.append(
                "ACHTUNG: die Schichtenfolge reicht nur bis z = {:.2f} m, der "
                "Pfahl ist jedoch {:.2f} m lang. Fuer den Bereich darunter wurde "
                "KEINE Mantelreibung angesetzt (Ergebnis auf der sicheren "
                "Seite). Schichtenfolge bis Pfahlfuss ergaenzen."
                .format(z_ende, L_pfahl))
        elif z_ende > L_pfahl + 1e-6:
            r.hinweise.append(
                "ACHTUNG: die Schichtenfolge reicht bis z = {:.2f} m und damit "
                "unter den Pfahlfuss bei {:.2f} m. Die Mantelreibung wurde ueber "
                "die volle Schichtdicke angesetzt - Schichtgrenzen an die "
                "Pfahllaenge anpassen.".format(z_ende, L_pfahl))
    if F_t_d > 0:
        r.hinweise.append(
            "Bei Zugpfaehlen ist zusaetzlich der Nachweis gegen Aufschwimmen "
            "bzw. das Herausziehen des Bodenblocks zu fuehren "
            "[DIN EN 1997-1, 7.6.3.1].")
    return r


def widerstands_setzungs_linie(D, R_b_k, R_s_k, n=60):
    """
    Vereinfachte Widerstands-Setzungs-Linie (WSL) nach EA-Pfaehle, 5.4.5.

    Charakteristische Setzungen:
        s_sg [cm] = 0,50 * R_s,k [MN] + 0,50 cm   ,  s_sg <= 3,0 cm
                    (Grenzsetzung fuer die volle Mantelreibung)
        s_g       = 0,10 * D
                    (Grenzsetzung fuer den vollen Spitzendruck)

    Zwischen 0 und den Grenzsetzungen wird linear interpoliert (Vereinfachung;
    die EA-Pfaehle geben gekruemmte Verlaeufe an).
    """
    s_sg = min(0.5 * (R_s_k / 1000.0) + 0.5, 3.0) / 100.0     # [m]
    s_g = 0.10 * D                                            # [m]
    s = np.linspace(0.0, max(s_g * 1.3, s_sg * 1.3), n)
    R_s = R_s_k * np.minimum(s / max(s_sg, 1e-9), 1.0)
    R_b = R_b_k * np.minimum(s / max(s_g, 1e-9), 1.0)
    return dict(s=s * 1000.0, s_mm=s * 1000.0, R_s=R_s, R_b=R_b, R=R_s + R_b,
                s_sg=s_sg * 1000.0, s_g=s_g * 1000.0,
                R_c_k=R_b_k + R_s_k, norm=ref("geo_wsl"))


def pfahlkopfsetzung(wsl, F_k):
    """Setzung [mm] bei charakteristischer Einwirkung F_k [kN] aus der WSL."""
    R, s = wsl["R"], wsl["s_mm"]
    if F_k <= 0:
        return 0.0
    if F_k >= R[-1]:
        return float(s[-1])
    return float(np.interp(F_k, R, s))
