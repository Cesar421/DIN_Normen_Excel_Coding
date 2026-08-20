# -*- coding: utf-8 -*-
"""
Baustoffe nach DIN EN 1992-1-1:2011-01 + DIN EN 1992-1-1/NA:2013-04.

Beton ............. EC2, 3.1  (Tab. 3.1, Gl. 3.15/3.16, Bild 3.3)
Betonstahl ........ EC2, 3.2  + DIN 488-1 (B500A / B500B)
Dauerhaftigkeit ... EC2, 4.2 / 4.4 + NA Tab. 4.3DE / 4.4DE
"""

import math
from dataclasses import dataclass

from .normen import ref

# ---------------------------------------------------------------------------
# Festigkeitsklassen - DIN EN 1992-1-1, 3.1.2, Tab. 3.1
# ---------------------------------------------------------------------------
BETONKLASSEN = [
    "C12/15", "C16/20", "C20/25", "C25/30", "C30/37", "C35/45", "C40/50",
    "C45/55", "C50/60", "C55/67", "C60/75", "C70/85", "C80/95", "C90/105",
    "C100/115",
]

# Lieferbare Stabdurchmesser - DIN 488-2 / DIN 488-4
DURCHMESSER = [6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 25.0, 28.0, 32.0, 40.0]


def stabflaeche(phi):
    """Querschnittsflaeche eines Stabes mit Durchmesser phi [mm] -> [mm2]."""
    return math.pi * phi ** 2 / 4.0


def stabflaeche_n(n, phi):
    """Querschnittsflaeche von n Staeben mit Durchmesser phi [mm] -> [mm2]."""
    return n * stabflaeche(phi)


# ---------------------------------------------------------------------------
# BETON
# ---------------------------------------------------------------------------
@dataclass
class Beton:
    """
    Normalbeton nach DIN EN 1992-1-1, 3.1.

    Alle Spannungen in N/mm2 (= MPa), Dehnungen in Promille.
    """
    klasse: str = "C30/37"
    gamma_c: float = 1.50      # NA, 2.4.2.4 (1), Tab. 2.1DE (staendige Situation)
    alpha_cc: float = 0.85     # NA, NDP zu 3.1.6 (1)P
    alpha_ct: float = 0.85     # NA, NDP zu 3.1.6 (2)P
    d_g: float = 16.0          # Groesstkorndurchmesser [mm] (EC2, 8.2)

    # --- charakteristische Festigkeiten (Tab. 3.1) ------------------------
    @property
    def fck(self):
        """Charakteristische Zylinderdruckfestigkeit [N/mm2]."""
        return float(self.klasse.lstrip("C").split("/")[0])

    @property
    def fck_wuerfel(self):
        """Charakteristische Wuerfeldruckfestigkeit [N/mm2]."""
        return float(self.klasse.split("/")[1])

    @property
    def fcm(self):
        """fcm = fck + 8   (Tab. 3.1)."""
        return self.fck + 8.0

    @property
    def fctm(self):
        """Mittlere Zugfestigkeit (Tab. 3.1)."""
        if self.fck <= 50.0:
            return 0.30 * self.fck ** (2.0 / 3.0)
        return 2.12 * math.log(1.0 + self.fcm / 10.0)

    @property
    def fctk005(self):
        """5 %-Fraktile der Zugfestigkeit (Tab. 3.1)."""
        return 0.7 * self.fctm

    @property
    def fctk095(self):
        """95 %-Fraktile der Zugfestigkeit (Tab. 3.1)."""
        return 1.3 * self.fctm

    @property
    def Ecm(self):
        """Sekantenmodul Ecm = 22000 (fcm/10)^0,3 [N/mm2] (Tab. 3.1)."""
        return 22000.0 * (self.fcm / 10.0) ** 0.3

    # --- Bemessungswerte ---------------------------------------------------
    @property
    def fcd(self):
        """fcd = alpha_cc fck / gamma_C   (EC2, 3.1.6 (1)P, Gl. 3.15)."""
        return self.alpha_cc * self.fck / self.gamma_c

    @property
    def fctd(self):
        """fctd = alpha_ct fctk;0,05 / gamma_C   (EC2, 3.1.6 (2)P, Gl. 3.16)."""
        return self.alpha_ct * self.fctk005 / self.gamma_c

    # --- Parabel-Rechteck-Diagramm (EC2, 3.1.7 (1), Bild 3.3) --------------
    @property
    def eps_c2(self):
        """Dehnung am Beginn des waagerechten Astes [Promille]."""
        if self.fck <= 50.0:
            return 2.0
        return 2.0 + 0.085 * (self.fck - 50.0) ** 0.53

    @property
    def eps_cu2(self):
        """Bruchdehnung des Betons unter Druck [Promille]."""
        if self.fck <= 50.0:
            return 3.5
        return 2.6 + 35.0 * ((90.0 - self.fck) / 100.0) ** 4

    @property
    def n_exp(self):
        """Exponent der Parabel."""
        if self.fck <= 50.0:
            return 2.0
        return 1.4 + 23.4 * ((90.0 - self.fck) / 100.0) ** 4

    def sigma_c(self, eps):
        """
        Betonspannung [N/mm2] bei einer (positiven) Druckdehnung eps
        [Promille].   EC2, 3.1.7 (1), Gl. (3.17)/(3.18).
        """
        if eps <= 0.0:
            return 0.0
        if eps < self.eps_c2:
            return self.fcd * (1.0 - (1.0 - eps / self.eps_c2) ** self.n_exp)
        if eps <= self.eps_cu2 + 1e-9:
            return self.fcd
        return 0.0

    # --- Kennwerte der Betondruckzone --------------------------------------
    def alpha_R(self, eps_c):
        """
        Voelligkeitsbeiwert der Betondruckzone:
            alpha_R = Fc / (b x fcd)
        Analytische Integration des Parabel-Rechteck-Diagramms.
        EC2, 3.1.7 (1) + 6.1.
        """
        e2, n = self.eps_c2, self.n_exp
        if eps_c <= 0.0:
            return 0.0
        if eps_c <= e2:
            r = eps_c / e2
            if r < 0.02:                       # Reihenentwicklung (Ausloeschung)
                return n * r / 2.0 * (1.0 - (n - 1.0) * r / 3.0)
            return 1.0 - (1.0 - (1.0 - r) ** (n + 1.0)) / (r * (n + 1.0))
        a_par = e2 * (1.0 - 1.0 / (n + 1.0))   # Flaeche der vollen Parabel
        return (a_par + (eps_c - e2)) / eps_c

    def k_a(self, eps_c):
        """
        Lage der Betondruckkraft, gemessen vom staerkst gedrueckten Rand:
            a = k_a * x.   EC2, 3.1.7 (1) + 6.1.
        """
        e2, n = self.eps_c2, self.n_exp
        aR = self.alpha_R(eps_c)
        if aR <= 0.0:
            return 1.0 / 3.0
        if eps_c <= e2:
            r = eps_c / e2
            if r < 0.02:
                return 1.0 / 3.0 + n * r / 36.0
            m = 0.5 - (1.0 - (1.0 - r) ** (n + 1.0) * (1.0 + r * (n + 1.0))) \
                / (r ** 2 * (n + 1.0) * (n + 2.0))
            return 1.0 - m / aR
        m_par = (0.5 - 1.0 / ((n + 1.0) * (n + 2.0))) * e2 ** 2
        m_rec = 0.5 * (eps_c ** 2 - e2 ** 2)
        m = (m_par + m_rec) / eps_c ** 2
        return 1.0 - m / aR

    # --- Rissbildung / Kriechen -------------------------------------------
    def fct_eff(self, frueh=False):
        """
        fct,eff fuer den Rissnachweis (EC2, 7.3.2 (2)).
        `frueh=True` -> Zwang in den ersten 28 Tagen: fct,eff = 0,50 fctm.
        """
        return 0.5 * self.fctm if frueh else self.fctm

    def Ec_eff(self, phi_kriech):
        """Wirksamer Modul mit Kriechen: Ec,eff = Ecm/(1+phi)  (EC2, 7.4.3 (5))."""
        return self.Ecm / (1.0 + phi_kriech)

    def kennwerte(self):
        return {
            "Klasse": self.klasse,
            "fck [N/mm2]": self.fck,
            "fcm [N/mm2]": self.fcm,
            "fctm [N/mm2]": self.fctm,
            "fctk;0,05 [N/mm2]": self.fctk005,
            "Ecm [N/mm2]": self.Ecm,
            "gamma_C": self.gamma_c,
            "alpha_cc": self.alpha_cc,
            "fcd [N/mm2]": self.fcd,
            "fctd [N/mm2]": self.fctd,
            "eps_c2 [permil]": self.eps_c2,
            "eps_cu2 [permil]": self.eps_cu2,
            "n": self.n_exp,
            "alpha_R (eps_cu2)": self.alpha_R(self.eps_cu2),
            "k_a (eps_cu2)": self.k_a(self.eps_cu2),
        }

    def normstellen(self):
        return [ref("beton_tab31"), ref("fcd"), ref("alpha_cc"),
                ref("fctd"), ref("alpha_ct"), ref("sigma_eps_c")]


# ---------------------------------------------------------------------------
# BETONSTAHL
# ---------------------------------------------------------------------------
@dataclass
class Betonstahl:
    """
    Betonstahl nach DIN 488-1 und DIN EN 1992-1-1, 3.2.

    Verwendet wird die bilineare Linie mit WAAGERECHTEM oberem Ast
    (EC2, 3.2.7 (2) b), Bild 3.8). Die Dehnung wird dennoch auf
    eps_ud = 25 Promille begrenzt (NA, NDP zu 3.2.7 (2)) - uebliche deutsche
    Praxis der Bemessungstafeln.
    """
    sorte: str = "B500B"
    fyk: float = 500.0         # DIN 488-1, Tab. 4
    gamma_s: float = 1.15      # NA, 2.4.2.4 (1), Tab. 2.1DE
    Es: float = 200000.0       # EC2, 3.2.7 (4)
    eps_ud: float = 25.0       # NA, NDP zu 3.2.7 (2) [Promille]

    @property
    def fyd(self):
        """fyd = fyk / gamma_S   (EC2, 3.2.7 (2))."""
        return self.fyk / self.gamma_s

    @property
    def eps_yd(self):
        """Bemessungswert der Streckgrenzendehnung [Promille]."""
        return 1000.0 * self.fyd / self.Es

    @property
    def k(self):
        """Charakteristisches Verhaeltnis ft/fy (DIN 488-1: A -> 1,05 ; B -> 1,08)."""
        return 1.05 if self.sorte.upper().endswith("A") else 1.08

    def sigma_s(self, eps):
        """
        Stahlspannung [N/mm2] bei der Dehnung eps [Promille]
        (waagerechter oberer Ast, EC2, 3.2.7 (2) b)).
        """
        s = math.copysign(1.0, eps)
        e = abs(eps)
        if e <= self.eps_yd:
            return s * self.Es * e / 1000.0
        return s * self.fyd

    def kennwerte(self):
        return {
            "Sorte": self.sorte,
            "fyk [N/mm2]": self.fyk,
            "gamma_S": self.gamma_s,
            "fyd [N/mm2]": self.fyd,
            "Es [N/mm2]": self.Es,
            "eps_yd [permil]": self.eps_yd,
            "eps_ud [permil]": self.eps_ud,
            "ft/fy (k)": self.k,
        }

    def normstellen(self):
        return [ref("betonstahl_din488"), ref("betonstahl_ec2"), ref("eps_ud")]


# ---------------------------------------------------------------------------
# DAUERHAFTIGKEIT - Expositionsklassen
# EC2, 4.2 Tab. 4.1 ; NA Tab. 4.3DE (Bauteilklasse) und Tab. 4.4DE (c_min,dur)
# ---------------------------------------------------------------------------
# (c_min,dur [mm] fuer Bauteilklasse S4 ; Delta_c_dev [mm] ;
#  Mindestbetonfestigkeitsklasse nach DIN 1045-2, Tab. F.2.1 ;
#  w_max [mm] fuer Stahlbeton nach NA Tab. 7.1DE)
EXPOSITION = {
    "X0":  dict(c_min_dur=10.0, dc_dev=10.0, klasse_min="C12/15", w_max=0.40,
                text="Kein Korrosions- oder Angriffsrisiko"),
    "XC1": dict(c_min_dur=10.0, dc_dev=10.0, klasse_min="C16/20", w_max=0.40,
                text="Trocken oder staendig nass (Karbonatisierung)"),
    "XC2": dict(c_min_dur=20.0, dc_dev=15.0, klasse_min="C16/20", w_max=0.30,
                text="Nass, selten trocken (Karbonatisierung)"),
    "XC3": dict(c_min_dur=20.0, dc_dev=15.0, klasse_min="C20/25", w_max=0.30,
                text="Maessige Feuchte (Karbonatisierung)"),
    "XC4": dict(c_min_dur=25.0, dc_dev=15.0, klasse_min="C25/30", w_max=0.30,
                text="Wechselnd nass und trocken (Karbonatisierung)"),
    "XD1": dict(c_min_dur=40.0, dc_dev=15.0, klasse_min="C30/37", w_max=0.30,
                text="Maessige Feuchte, Chloride ausser Meerwasser"),
    "XD2": dict(c_min_dur=40.0, dc_dev=15.0, klasse_min="C35/45", w_max=0.30,
                text="Nass, selten trocken, Chloride ausser Meerwasser"),
    "XD3": dict(c_min_dur=40.0, dc_dev=15.0, klasse_min="C35/45", w_max=0.30,
                text="Wechselnd nass und trocken, Chloride ausser Meerwasser"),
    "XS1": dict(c_min_dur=40.0, dc_dev=15.0, klasse_min="C30/37", w_max=0.30,
                text="Salzhaltige Luft, kein unmittelbarer Meerwasserkontakt"),
    "XS2": dict(c_min_dur=40.0, dc_dev=15.0, klasse_min="C35/45", w_max=0.30,
                text="Unter Wasser (Meerwasser)"),
    "XS3": dict(c_min_dur=40.0, dc_dev=15.0, klasse_min="C35/45", w_max=0.30,
                text="Tidebereich, Spritzwasser- und Sprühnebelbereich"),
}


def betondeckung(expo_klasse, phi_laengs, phi_buegel=0.0, d_g=16.0,
                 delta_c_dev=None, beton=None):
    """
    Nennmass der Betondeckung nach DIN EN 1992-1-1, 4.4.1 und NA Tab. 4.4DE.

        c_min   = max(c_min,b ; c_min,dur ; 10 mm)      Gl. (4.2)
        c_nom   = c_min + Delta_c_dev                   Gl. (4.1)

    Massgebend fuer die Ausfuehrung ist die Deckung des BUEGELS (Verlegemass);
    daraus folgt die Deckung der Laengsbewehrung.
    """
    if expo_klasse not in EXPOSITION:
        raise ValueError("Unbekannte Expositionsklasse: " + str(expo_klasse))
    e = EXPOSITION[expo_klasse]
    dc = e["dc_dev"] if delta_c_dev is None else float(delta_c_dev)

    c_min_b_w = phi_buegel if phi_buegel > 0 else phi_laengs
    c_min_w = max(c_min_b_w, e["c_min_dur"], 10.0)
    c_nom_w = c_min_w + dc

    c_min_l = max(phi_laengs, e["c_min_dur"], 10.0)
    c_nom_l = max(c_min_l + dc, c_nom_w + phi_buegel)

    hinweis = ""
    if beton is not None:
        idx_req = BETONKLASSEN.index(e["klasse_min"])
        idx_use = BETONKLASSEN.index(beton.klasse)
        if idx_use < idx_req:
            hinweis = ("Betonfestigkeitsklasse " + beton.klasse
                       + " < Mindestklasse " + e["klasse_min"] + " fuer "
                       + expo_klasse + " (DIN 1045-2, Tab. F.2.1)")
    return {
        "expo_klasse": expo_klasse,
        "beschreibung": e["text"],
        "c_min_dur": e["c_min_dur"],
        "Delta_c_dev": dc,
        "c_min_w": c_min_w,
        "c_nom_w": c_nom_w,     # Nennmass am Buegel
        "c_nom_l": c_nom_l,     # Nennmass an der Laengsbewehrung
        "w_max": e["w_max"],
        "klasse_min_beton": e["klasse_min"],
        "hinweis": hinweis,
        "normen": [ref("expos"), ref("c_nom"), ref("c_min"), ref("c_min_dur"),
                   ref("dc_dev")],
    }


def stabwahl(As_erf, b, c_nom_w, phi_buegel, d_g=16.0,
             n_max=10, n_min=2, max_ergebnis=6):
    """
    Schlaegt Kombinationen n x phi vor, die As_erf [mm2] abdecken und in einer
    Lage in die Breite b [mm] passen.

    Lichter Stababstand (EC2, 8.2 (2)): s_l >= max(phi ; d_g + 5 mm ; 20 mm)
    """
    ergebnis = []
    for phi in DURCHMESSER:
        for n in range(n_min, n_max + 1):
            As = stabflaeche_n(n, phi)
            if As < As_erf:
                continue
            s_min = max(phi, d_g + 5.0, 20.0)
            b_erf = 2.0 * c_nom_w + 2.0 * phi_buegel + n * phi + (n - 1) * s_min
            passt = b_erf <= b + 1e-6
            s_licht = ((b - 2.0 * c_nom_w - 2.0 * phi_buegel - n * phi) / (n - 1)
                       if n > 1 else float("inf"))
            ergebnis.append(dict(n=n, phi=phi, As=As, ueberschuss=As / As_erf - 1.0,
                                 passt_1lage=passt, s_licht=s_licht, b_erf=b_erf))
            break
    ergebnis.sort(key=lambda r: (not r["passt_1lage"], r["ueberschuss"]))
    return ergebnis[:max_ergebnis]
