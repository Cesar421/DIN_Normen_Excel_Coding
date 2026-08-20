# -*- coding: utf-8 -*-
"""
Gesamtbemessung eines Stahlbetonbalkens nach DIN EN 1992-1-1 + Nationalem Anhang.

Ablauf (mit der jeweils angewendeten Normstelle):

  1. Baustoffe                        EC2 3.1 / 3.2 + DIN 488-1
  2. Dauerhaftigkeit, Betondeckung    EC2 4.4 + NA Tab. 4.4DE
  3. Querschnitt                      EC2 5.3.2.1
  4. Einwirkungen, Schnittgroessen    EC0 6.4.3.2 Gl. (6.10) / 6.5.3 Gl. (6.16b)
  5. Biegung GZT                      EC2 6.1 + NA NDP zu 5.5 (4)
  6. Mindest-/Hoechstbewehrung        EC2 9.2.1.1 Gl. (9.1N) + NA
  7. Querkraft GZT                    EC2 6.2 + NA NDP zu 6.2.2 / 6.2.3
  8. Torsion GZT                      EC2 6.3 + NA NDP zu 6.3.2 + EC2 9.2.3
  9. Rissbreite GZG                   EC2 7.3 + NA Gl. (7.11DE)
 10. Durchbiegung GZG                 EC2 7.4
 11. Verankerung, Versatzmass         EC2 8.4 / 9.2.1.3 / 9.2.1.4
 12. Zusammenstellung der Nachweise
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .normen import ref, normentabelle
from .baustoffe import (Beton, Betonstahl, betondeckung, stabflaeche,
                        stabflaeche_n)
from .querschnitt import Querschnitt
from .schnittgroessen import (Durchlauftraeger, Auflager, GELENKIG, EINGESPANNT,
                              bemessungsquerkraft)
from .biegung import bemessung_biegung, momententragfaehigkeit, xi_grenz
from .querkraft import bemessung_querkraft, verankerungskraft_endauflager
from .torsion import bemessung_torsion, gleichgewichtstorsion_hinweis
from .gebrauchstauglichkeit import (rissbreite, mindestbewehrung_riss,
                                    nachweis_durchbiegung)
from .konstruktion import (mindestbewehrung_biegung, hoechstbewehrung,
                           robustheitsbewehrung, verankerungslaenge,
                           verankerung_endauflager, platznachweis,
                           d1_schaetzung, uebergreifungslaenge)


# ---------------------------------------------------------------------------
# Eingabedaten
# ---------------------------------------------------------------------------
@dataclass
class EingabeBalken:
    # --- Geometrie -------------------------------------------------------
    L: float = 6.00                 # Gesamtlaenge [m]
    auflager: list = field(default_factory=lambda: [(0.0, GELENKIG, 0.30),
                                                    (6.0, GELENKIG, 0.30)])
    querschnittstyp: str = "rechteck"   # "rechteck" | "plattenbalken"
    b: float = 300.0                # Stegbreite bw [mm]
    h: float = 600.0                # Gesamthoehe [mm]
    b_eff: float = 1200.0           # mitwirkende Plattenbreite [mm]
    hf: float = 150.0               # Plattendicke [mm]

    # --- Baustoffe -------------------------------------------------------
    betonklasse: str = "C30/37"
    stahlsorte: str = "B500B"
    expositionsklasse: str = "XC1"
    d_g: float = 16.0               # Groesstkorn [mm]

    # --- Einwirkungen ----------------------------------------------------
    g_k: float = 15.0               # staendige Streckenlast [kN/m] (ohne EG)
    q_k: float = 20.0               # veraenderliche Streckenlast [kN/m]
    eigengewicht: bool = True       # Eigengewicht ansetzen (25 kN/m3, EC1 Tab. A.1)
    einzellasten: list = field(default_factory=list)  # [(x[m], P[kN], "G"/"Q")]
    psi_2: float = 0.30             # DIN EN 1990/NA, Tab. NA.A.1.1
    gamma_G: float = 1.35
    gamma_Q: float = 1.50
    laststellungen: bool = True     # feldweise Laststellung (EC2, 5.1.3)
    N_Ed: float = 0.0               # Laengskraft [kN], Druck negativ
    T_Ed: float = 0.0               # Torsionsmoment [kNm] (Bemessungswert)
    gleichgewichtstorsion: bool = True   # False -> Vertraeglichkeitstorsion
    kastenquerschnitt: bool = False      # NA: nu = 0,75 nu_2 statt 0,525 nu_2

    # --- Bewehrung -------------------------------------------------------
    phi_laengs: float = 20.0        # Durchmesser untere Laengsbewehrung [mm]
    phi_laengs_oben: float = 16.0   # Durchmesser obere Laengsbewehrung [mm]
    phi_buegel: float = 8.0         # Buegeldurchmesser [mm]
    n_schenkel: int = 2             # Schenkel je Buegel
    delta: float = 1.0              # Umlagerungsgrad

    # --- Gebrauchstauglichkeit -------------------------------------------
    phi_kriech: float = 2.0         # Kriechzahl phi(inf,t0)
    eps_cs: float = 0.0             # Schwinddehnung [-] (z.B. 5e-4)
    grenze_durchbiegung: float = 250.0   # l/250
    verformungsempfindlich: bool = False

    def beton(self):
        return Beton(self.betonklasse, d_g=self.d_g)

    def stahl(self):
        return Betonstahl(self.stahlsorte)


# ---------------------------------------------------------------------------
# Hilfsfunktionen fuer den Bericht
# ---------------------------------------------------------------------------
def _block(titel, normen=None):
    return dict(titel=titel, normen=list(normen or []), zeilen=[], nachweise=[])


def _z(bl, text):
    bl["zeilen"].append(text)


def _nw(bl, name, wert, grenzwert, ok, einheit="", vergleich="<="):
    bl["nachweise"].append(dict(name=name, wert=wert, grenzwert=grenzwert,
                                ok=bool(ok), einheit=einheit,
                                vergleich=vergleich))


def _f(v, n=2):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "-"
    return ("{:." + str(n) + "f}").format(v)


# ---------------------------------------------------------------------------
# Gesamtbemessung
# ---------------------------------------------------------------------------
def bemessung_balken(e):
    """Fuehrt die vollstaendige Bemessung durch und liefert einen Bericht."""
    C = e.beton()
    S = e.stahl()
    b = dict(eingabe=e, beton=C, stahl=S, bloecke=[], ok_gesamt=True)
    B = b["bloecke"]

    # =====================================================================
    # 1. BAUSTOFFE
    # =====================================================================
    bl = _block("1. BAUSTOFFE",
                C.normstellen() + S.normstellen() + [ref("gamma_M")])
    _z(bl, "Beton {}  [EC2 3.1.2, Tab. 3.1]".format(C.klasse))
    _z(bl, "   fck = {} N/mm2   fcm = {} N/mm2   fctm = {} N/mm2   fctk;0,05 = {} N/mm2"
       .format(_f(C.fck, 1), _f(C.fcm, 1), _f(C.fctm, 2), _f(C.fctk005, 2)))
    _z(bl, "   Ecm = {} N/mm2".format(_f(C.Ecm, 0)))
    _z(bl, "   gamma_C = {}   alpha_cc = {}   [NA 2.4.2.4 Tab. 2.1DE / NDP zu 3.1.6(1)P]"
       .format(_f(C.gamma_c, 2), _f(C.alpha_cc, 2)))
    _z(bl, "   fcd = alpha_cc fck/gamma_C = {} * {} / {} = {} N/mm2   [Gl. (3.15)]"
       .format(_f(C.alpha_cc, 2), _f(C.fck, 0), _f(C.gamma_c, 2), _f(C.fcd, 2)))
    _z(bl, "   fctd = alpha_ct fctk;0,05/gamma_C = {} N/mm2            [Gl. (3.16)]"
       .format(_f(C.fctd, 3)))
    _z(bl, "   Parabel-Rechteck-Diagramm [EC2 3.1.7(1), Bild 3.3]:")
    _z(bl, "      eps_c2 = {} permil   eps_cu2 = {} permil   n = {}"
       .format(_f(C.eps_c2, 2), _f(C.eps_cu2, 2), _f(C.n_exp, 2)))
    _z(bl, "      alpha_R = {}   k_a = {}   (Druckzone im Bruchzustand)"
       .format(_f(C.alpha_R(C.eps_cu2), 4), _f(C.k_a(C.eps_cu2), 4)))
    _z(bl, "")
    _z(bl, "Betonstahl {}  [DIN 488-1, Tab. 4 ; EC2 3.2]".format(S.sorte))
    _z(bl, "   fyk = {} N/mm2   gamma_S = {}   fyd = {} N/mm2"
       .format(_f(S.fyk, 0), _f(S.gamma_s, 2), _f(S.fyd, 2)))
    _z(bl, "   Es = {} N/mm2   eps_yd = {} permil   eps_ud = {} permil  [NA NDP zu 3.2.7(2)]"
       .format(_f(S.Es, 0), _f(S.eps_yd, 2), _f(S.eps_ud, 1)))
    B.append(bl)

    # =====================================================================
    # 2. DAUERHAFTIGKEIT UND BETONDECKUNG
    # =====================================================================
    bd = betondeckung(e.expositionsklasse, e.phi_laengs, e.phi_buegel, e.d_g,
                      beton=C)
    bl = _block("2. DAUERHAFTIGKEIT UND BETONDECKUNG", bd["normen"])
    _z(bl, "Expositionsklasse {}: {}   [EC2 4.2, Tab. 4.1]"
       .format(bd["expo_klasse"], bd["beschreibung"]))
    _z(bl, "   c_min,dur = {} mm            [NA Tab. 4.4DE]".format(_f(bd["c_min_dur"], 0)))
    _z(bl, "   c_min,b = phi_Buegel = {} mm  [EC2 4.4.1.2 (3)]".format(_f(e.phi_buegel, 0)))
    _z(bl, "   c_min = max(c_min,b ; c_min,dur ; 10) = {} mm   [Gl. (4.2)]"
       .format(_f(bd["c_min_w"], 0)))
    _z(bl, "   Delta_c_dev = {} mm          [NA NDP zu 4.4.1.3 (1)P]"
       .format(_f(bd["Delta_c_dev"], 0)))
    _z(bl, "   c_nom (Buegel)        = c_min + Delta_c_dev = {} mm   [Gl. (4.1)]"
       .format(_f(bd["c_nom_w"], 0)))
    _z(bl, "   c_nom (Laengsstab)    = {} mm".format(_f(bd["c_nom_l"], 0)))
    if bd["hinweis"]:
        _z(bl, "   !! " + bd["hinweis"])
    _z(bl, "   w_max = {} mm (zulaessige Rissbreite)  [NA Tab. 7.1DE]"
       .format(_f(bd["w_max"], 2)))
    B.append(bl)

    c_w, c_l = bd["c_nom_w"], bd["c_nom_l"]

    # =====================================================================
    # 3. QUERSCHNITT
    # =====================================================================
    d1 = d1_schaetzung(c_w, e.phi_buegel, e.phi_laengs, 1, e.d_g)
    d1_oben = d1_schaetzung(c_w, e.phi_buegel, e.phi_laengs_oben, 1, e.d_g)

    def _qs(feld=True, lagen=1):
        """Querschnitt orientiert: feld=True -> Feldmoment (Zug unten)."""
        dd1 = d1_schaetzung(c_w, e.phi_buegel, e.phi_laengs, lagen, e.d_g) \
            if feld else d1_oben
        dd2 = d1_oben if feld else d1
        if e.querschnittstyp == "plattenbalken" and feld:
            return Querschnitt(b=e.b, h=e.h, d1=dd1, d2=dd2, typ="plattenbalken",
                               b_eff=e.b_eff, hf=e.hf, platte_gedrueckt=True)
        return Querschnitt(b=e.b, h=e.h, d1=dd1, d2=dd2, typ="rechteck")

    qs_feld = _qs(True, 1)
    qs_stuetz = _qs(False)

    bl = _block("3. QUERSCHNITT",
                [ref("b_eff")] if e.querschnittstyp == "plattenbalken" else [])
    _z(bl, qs_feld.beschreibung())
    _z(bl, "   d1 (unten) = c_nom + phi_w + phi_l/2 = {} + {} + {}/2 = {} mm"
       .format(_f(c_w, 0), _f(e.phi_buegel, 0), _f(e.phi_laengs, 0), _f(d1, 1)))
    _z(bl, "   d (Feldbereich)   = h - d1 = {} - {} = {} mm"
       .format(_f(e.h, 0), _f(d1, 1), _f(qs_feld.d, 1)))
    _z(bl, "   d (Stuetzbereich) = h - d1,oben = {} mm".format(_f(qs_stuetz.d, 1)))
    _z(bl, "   Ac = {} mm2 = {} m2".format(_f(qs_feld.Ac, 0), _f(qs_feld.Ac / 1e6, 4)))
    if e.querschnittstyp == "plattenbalken":
        _z(bl, "   b_eff = {} mm  [EC2 5.3.2.1 (3), Gl. (5.7)] - Eingabewert"
           .format(_f(e.b_eff, 0)))
    B.append(bl)

    # =====================================================================
    # 4. EINWIRKUNGEN UND SCHNITTGROESSEN
    # =====================================================================
    g_eg = 25.0 * (e.b * e.h) / 1.0e6 if e.eigengewicht else 0.0
    g_ges = e.g_k + g_eg

    auflager = [Auflager(x=a[0], typ=a[1], breite=a[2]) for a in e.auflager]
    EI = C.Ecm * qs_feld.traegheitsmoment() / 1.0e9   # N/mm2 * mm4 -> kNm2
    traeger = Durchlauftraeger(e.L, auflager, EI=EI)
    traeger.strecke(0.0, e.L, g_ges, art="G")
    traeger.strecke(0.0, e.L, e.q_k, art="Q")
    for el in e.einzellasten:
        traeger.einzel(el[0], el[1], el[2] if len(el) > 2 else "G")

    einh = traeger.einhuellende(e.gamma_G, e.gamma_Q, e.laststellungen)
    qs_last = traeger.quasi_staendig(e.psi_2)
    b["traeger"], b["einhuellende"], b["quasi_staendig"] = traeger, einh, qs_last

    x = einh["x"]
    Mmax = float(np.max(einh["Mmax"]))
    i_max = int(np.argmax(einh["Mmax"]))
    Mmin = float(np.min(einh["Mmin"]))
    i_min = int(np.argmin(einh["Mmin"]))
    Veinh = np.maximum(np.abs(einh["Vmax"]), np.abs(einh["Vmin"]))

    bl = _block("4. EINWIRKUNGEN UND SCHNITTGROESSEN",
                [ref("eigengewicht"), ref("komb_GZT"), ref("gamma_F"),
                 ref("komb_QS"), ref("psi"), ref("lastfaelle")])
    _z(bl, "Charakteristische Einwirkungen:")
    if e.eigengewicht:
        _z(bl, "   Eigengewicht g_EG = 25,0 kN/m3 * {} m2 = {} kN/m   [EC1 Tab. A.1]"
           .format(_f(e.b * e.h / 1e6, 4), _f(g_eg, 2)))
    _z(bl, "   g_k = {} kN/m (+ EG) -> G_k = {} kN/m".format(_f(e.g_k, 2), _f(g_ges, 2)))
    _z(bl, "   q_k = {} kN/m".format(_f(e.q_k, 2)))
    for el in e.einzellasten:
        _z(bl, "   Einzellast {} kN bei x = {} m ({})"
           .format(_f(el[1], 2), _f(el[0], 2), el[2] if len(el) > 2 else "G"))
    if abs(e.T_Ed) > 1e-9:
        _z(bl, "   Torsionsmoment T_Ed = {} kNm ({})"
           .format(_f(e.T_Ed, 2),
                   "Gleichgewichtstorsion" if e.gleichgewichtstorsion
                   else "Vertraeglichkeitstorsion"))
    _z(bl, "")
    _z(bl, "Grundkombination GZT  [EC0 6.4.3.2, Gl. (6.10)]:")
    _z(bl, "   Ed = {} * G_k + {} * Q_k      [EC0/NA Tab. NA.A.1.2(B)]"
       .format(_f(e.gamma_G, 2), _f(e.gamma_Q, 2)))
    _z(bl, "   q_Ed = {} kN/m".format(_f(e.gamma_G * g_ges + e.gamma_Q * e.q_k, 2)))
    if e.laststellungen:
        _z(bl, "   Feldweise Laststellung der veraenderlichen Last  [EC2 5.1.3]")
    _z(bl, "")
    _z(bl, "Quasi-staendige Kombination GZG  [EC0 6.5.3, Gl. (6.16b)]:")
    _z(bl, "   q_qs = G_k + psi_2 Q_k = {} + {} * {} = {} kN/m"
       .format(_f(g_ges, 2), _f(e.psi_2, 2), _f(e.q_k, 2),
               _f(g_ges + e.psi_2 * e.q_k, 2)))
    _z(bl, "")
    _z(bl, "Bemessungsschnittgroessen:")
    _z(bl, "   M_Ed,max (Feld)   = {} kNm  bei x = {} m".format(_f(Mmax, 2), _f(x[i_max], 2)))
    _z(bl, "   M_Ed,min (Stuetze) = {} kNm  bei x = {} m".format(_f(Mmin, 2), _f(x[i_min], 2)))
    _z(bl, "   V_Ed,max = {} kN".format(_f(float(np.max(Veinh)), 2)))
    if abs(e.T_Ed) > 1e-9:
        _z(bl, "   T_Ed     = {} kNm".format(_f(e.T_Ed, 2)))
    for i, (xr, _R, _Mr) in enumerate(einh["faelle"][0].auflagerkraefte):
        _z(bl, "   Auflagerkraft bei x = {} m : R_Ed = {} kN (Einhuellende)"
           .format(_f(xr, 2), _f(float(einh["R"][i]), 2)))
    B.append(bl)

    # =====================================================================
    # 5. BIEGEBEMESSUNG
    # =====================================================================
    grenz = xi_grenz(C, e.delta)
    bl = _block("5. BIEGEBEMESSUNG (GZT)",
                [ref("biegung"), ref("eps_grenzen"), ref("sigma_eps_c"),
                 ref("betonstahl_ec2"), ref("xi_lim")])
    _z(bl, "Duktilitaetsgrenze  [NA NDP zu 5.5 (4), Gl. (5.10a)]:")
    _z(bl, "   delta = {}  ->  xu/d <= (delta - k1)/k2 = ({} - {})/{} = {}"
       .format(_f(e.delta, 2), _f(e.delta, 2), _f(grenz["k1"], 2),
               _f(grenz["k2"], 2), _f(grenz["xi_lim"], 3)))
    _z(bl, "")

    # --- untere Bewehrung (Feldmoment) -----------------------------------
    biege_feld, n_unten, lagen = None, 0, 1
    if Mmax > 1e-6:
        for lagen in (1, 2):
            qs_feld = _qs(True, lagen)
            biege_feld = bemessung_biegung(qs_feld, C, S, Mmax, e.N_Ed, e.delta)
            n_unten = max(2, int(math.ceil(biege_feld.As1 / stabflaeche(e.phi_laengs))))
            platz = platznachweis(e.b, n_unten, e.phi_laengs, c_w, e.phi_buegel,
                                  e.d_g, lagen)
            if platz["ok"]:
                break
        As1_vorh = stabflaeche_n(n_unten, e.phi_laengs)
        _z(bl, "--- UNTERE BEWEHRUNG (M_Ed = {} kNm) ---".format(_f(Mmax, 2)))
        _z(bl, "   M_Eds = M_Ed - N_Ed z_s1 = {} kNm".format(_f(biege_feld.M_Eds, 2)))
        _z(bl, "   mu_Eds = M_Eds/(b d^2 fcd) = {:.4g}*10^6/({} * {}^2 * {}) = {}"
           .format(biege_feld.M_Eds, _f(biege_feld.b_bezug, 0), _f(qs_feld.d, 1),
                   _f(C.fcd, 2), _f(biege_feld.mu_Eds, 4)))
        _z(bl, "   mu_lim (xu/d = {}) = {}".format(_f(biege_feld.xi_lim, 2),
                                                   _f(biege_feld.mu_lim, 4)))
        _z(bl, "   -> xi = x/d = {}   x = {} mm   (Bemessungspunkt {})"
           .format(_f(biege_feld.xi, 4), _f(biege_feld.x, 1), biege_feld.punkt))
        _z(bl, "   -> zeta = z/d = {}   z = {} mm".format(_f(biege_feld.zeta, 4),
                                                          _f(biege_feld.z, 1)))
        _z(bl, "   -> eps_c = {} permil   eps_s1 = {} permil   sigma_s1 = {} N/mm2"
           .format(_f(biege_feld.eps_c, 2), _f(biege_feld.eps_s1, 2),
                   _f(biege_feld.sigma_s1, 1)))
        _z(bl, "   -> omega = {}   Fcd = {} kN".format(_f(biege_feld.omega, 4),
                                                       _f(biege_feld.Fc, 1)))
        _z(bl, "   As1,erf (Biegung) = (Fcd + N_Ed)/sigma_s1 = {} mm2"
           .format(_f(biege_feld.As1, 0)))
        if biege_feld.As2 > 0:
            _z(bl, "   As2,erf = {} mm2 (Druckbewehrung) [EC2 6.1]"
               .format(_f(biege_feld.As2, 0)))
        _z(bl, "   GEWAEHLT: {} phi {} = {} mm2   (As,vorh/As,erf = {})"
           .format(n_unten, _f(e.phi_laengs, 0), _f(As1_vorh, 0),
                   _f(As1_vorh / max(biege_feld.As1, 1e-9), 3)))
        platz = platznachweis(e.b, n_unten, e.phi_laengs, c_w, e.phi_buegel,
                              e.d_g, lagen)
        _z(bl, "   Platznachweis [EC2 8.2 (2)]: b_erf = {} mm {} b = {} mm  ({} Lage/n)"
           .format(_f(platz["b_erf"], 0), "<=" if platz["ok"] else ">",
                   _f(e.b, 0), platz["n_lagen"]))
        _nw(bl, "Platz untere Bewehrung", platz["b_erf"], e.b, platz["ok"], "mm")
        mr = momententragfaehigkeit(qs_feld, C, S, As1_vorh, biege_feld.As2, e.N_Ed)
        _z(bl, "   NACHWEIS: M_Rd = {} kNm >= M_Ed = {} kNm  (eta = {})"
           .format(_f(mr["M_Rd"], 2), _f(Mmax, 2), _f(Mmax / max(mr["M_Rd"], 1e-9), 3)))
        _nw(bl, "Biegung Feld  M_Ed/M_Rd", Mmax, mr["M_Rd"], mr["M_Rd"] >= Mmax, "kNm")
        for hh in biege_feld.hinweise:
            _z(bl, "   !! " + hh)
        b["M_Rd_feld"] = mr["M_Rd"]
    else:
        n_unten = 2
        As1_vorh = stabflaeche_n(n_unten, e.phi_laengs)
        _z(bl, "--- UNTERE BEWEHRUNG ---")
        _z(bl, "   Kein Feldmoment (M_Ed,max = 0): die Unterseite ist gedrueckt.")
        _z(bl, "   Es wird eine konstruktive Montagebewehrung angeordnet.")
        _z(bl, "   GEWAEHLT: 2 phi {} = {} mm2".format(_f(e.phi_laengs, 0),
                                                       _f(As1_vorh, 0)))
        b["M_Rd_feld"] = 0.0

    b.update(biegung_feld=biege_feld, As1_vorh=As1_vorh, n_unten=n_unten,
             qs_feld=qs_feld, lagen_unten=lagen)

    # --- obere Bewehrung (Stuetzmoment) ----------------------------------
    biege_stuetz, n_oben = None, 2
    As_oben_vorh = stabflaeche_n(2, e.phi_laengs_oben)
    if Mmin < -1e-6:
        biege_stuetz = bemessung_biegung(qs_stuetz, C, S, abs(Mmin), e.N_Ed, e.delta)
        n_oben = max(2, int(math.ceil(biege_stuetz.As1 / stabflaeche(e.phi_laengs_oben))))
        As_oben_vorh = stabflaeche_n(n_oben, e.phi_laengs_oben)
        _z(bl, "")
        _z(bl, "--- OBERE BEWEHRUNG (M_Ed = {} kNm) ---".format(_f(Mmin, 2)))
        _z(bl, "   Rechteckquerschnitt bw/h (die Platte liegt in der Zugzone)")
        _z(bl, "   mu_Eds = {}   xi = {}   zeta = {}"
           .format(_f(biege_stuetz.mu_Eds, 4), _f(biege_stuetz.xi, 4),
                   _f(biege_stuetz.zeta, 4)))
        _z(bl, "   eps_s1 = {} permil   As1,erf = {} mm2"
           .format(_f(biege_stuetz.eps_s1, 2), _f(biege_stuetz.As1, 0)))
        if biege_stuetz.As2 > 0:
            _z(bl, "   As2,erf = {} mm2 (Druckbewehrung)".format(_f(biege_stuetz.As2, 0)))
        _z(bl, "   GEWAEHLT: {} phi {} = {} mm2".format(n_oben, _f(e.phi_laengs_oben, 0),
                                                        _f(As_oben_vorh, 0)))
        pl_o = platznachweis(e.b, n_oben, e.phi_laengs_oben, c_w, e.phi_buegel,
                             e.d_g, 1)
        _nw(bl, "Platz obere Bewehrung", pl_o["b_erf"], e.b, pl_o["ok"], "mm")
        mro = momententragfaehigkeit(qs_stuetz, C, S, As_oben_vorh,
                                     biege_stuetz.As2, e.N_Ed)
        _z(bl, "   NACHWEIS: M_Rd = {} kNm >= |M_Ed| = {} kNm"
           .format(_f(mro["M_Rd"], 2), _f(abs(Mmin), 2)))
        _nw(bl, "Biegung Stuetze |M_Ed|/M_Rd", abs(Mmin), mro["M_Rd"],
            mro["M_Rd"] >= abs(Mmin), "kNm")
        for hh in biege_stuetz.hinweise:
            _z(bl, "   !! " + hh)
        b["M_Rd_stuetze"] = mro["M_Rd"]
    else:
        _z(bl, "")
        _z(bl, "--- OBERE BEWEHRUNG ---")
        _z(bl, "   Kein Stuetzmoment: konstruktive Montagebewehrung")
        _z(bl, "   GEWAEHLT: 2 phi {} = {} mm2".format(_f(e.phi_laengs_oben, 0),
                                                       _f(As_oben_vorh, 0)))
        b["M_Rd_stuetze"] = 0.0
    b.update(biegung_stuetze=biege_stuetz, As_oben_vorh=As_oben_vorh,
             n_oben=n_oben, qs_stuetz=qs_stuetz)
    B.append(bl)

    # =====================================================================
    # 6. MINDEST- UND HOECHSTBEWEHRUNG
    # =====================================================================
    amin = mindestbewehrung_biegung(b["qs_feld"], C, S)
    amax = hoechstbewehrung(b["qs_feld"])
    arob = robustheitsbewehrung(b["qs_feld"], C, S)
    ariss = mindestbewehrung_riss(b["qs_feld"], C, S)
    As_noetig = max(amin["As_min"], arob["As_rob"])

    bl = _block("6. MINDEST- UND HOECHSTBEWEHRUNG",
                [ref("As_min"), ref("As_max"), ref("robustheit"), ref("As_min_riss")])
    _z(bl, "As,min (Biegung) = max(0,26 fctm/fyk bt d ; 0,0013 bt d)   "
           "[EC2 9.2.1.1, Gl. (9.1N)]")
    _z(bl, "   = max(0,26 * {}/{} * {} * {} ; 0,0013 * {} * {})"
       .format(_f(C.fctm, 2), _f(S.fyk, 0), _f(amin["bt"], 0), _f(b["qs_feld"].d, 1),
               _f(amin["bt"], 0), _f(b["qs_feld"].d, 1)))
    _z(bl, "   = max({} ; {}) = {} mm2   (massgebend: {})"
       .format(_f(amin["term_fctm"], 0), _f(amin["term_0013"], 0),
               _f(amin["As_min"], 0), amin["massgebend"]))
    _z(bl, "As,rob (Robustheit, Mcr = fctm W = {} kNm) = {} mm2   [NA NDP zu 9.2.1.1 (1)]"
       .format(_f(arob["Mcr"], 2), _f(arob["As_rob"], 0)))
    _z(bl, "As,min (Rissbewehrung, kc={} k={}) = {} mm2   [EC2 7.3.2, Gl. (7.1)]"
       .format(_f(ariss["kc"], 2), _f(ariss["k"], 2), _f(ariss["As_min"], 0)))
    _z(bl, "As,max = 0,04 Ac = {} mm2   [NA NDP zu 9.2.1.1 (3)]".format(_f(amax["As_max"], 0)))
    _z(bl, "")
    _z(bl, "As,vorh (unten) = {} mm2   |   As,vorh (oben) = {} mm2"
       .format(_f(As1_vorh, 0), _f(As_oben_vorh, 0)))
    _z(bl, "Die Mindestbewehrung ist auf der jeweils GEZOGENEN Seite nachzuweisen")
    _z(bl, "[EC2 9.2.1.1 (1)].")
    if Mmax > 1e-6:
        _nw(bl, "As,min Unterseite (Feldzug)", As1_vorh, As_noetig,
            As1_vorh >= As_noetig, "mm2", ">=")
    if Mmin < -1e-6:
        _nw(bl, "As,min Oberseite (Stuetzzug)", As_oben_vorh, As_noetig,
            As_oben_vorh >= As_noetig, "mm2", ">=")
    _nw(bl, "As,vorh <= As,max", As1_vorh + As_oben_vorh, amax["As_max"],
        As1_vorh + As_oben_vorh <= amax["As_max"], "mm2")
    b["mindestbewehrung"] = dict(amin=amin, amax=amax, arob=arob, ariss=ariss)
    B.append(bl)

    # =====================================================================
    # 7. QUERKRAFTBEMESSUNG
    # =====================================================================
    As_l_quer = As1_vorh if Mmin >= -1e-6 else min(As1_vorh, As_oben_vorh)
    vbem = bemessungsquerkraft(einh, auflager, b["qs_feld"].d / 1000.0)
    V_bem = max([v["V_d"] for v in vbem]) if vbem else float(np.max(Veinh))

    # 1. Durchgang: nur Querkraft -> liefert cot(theta), V_Rd,max, V_Rd,c
    quer = bemessung_querkraft(b["qs_feld"], C, S, x, Veinh, As_l_quer, c_l,
                               e.phi_buegel, e.n_schenkel, e.N_Ed,
                               V_bem_kN=V_bem)

    # =====================================================================
    # 8. TORSIONSBEMESSUNG  (vor der endgueltigen Buegelwahl)
    # =====================================================================
    tor = bemessung_torsion(b["qs_feld"], C, S, e.T_Ed, V_bem, quer.cot_theta,
                            quer.V_Rdmax, quer.V_Rdc, e.kastenquerschnitt)
    asw_tor = tor.asw_je_schenkel if (tor.erforderlich and abs(e.T_Ed) > 1e-9) else 0.0
    smax_tor = tor.s_max if abs(e.T_Ed) > 1e-9 else None

    # 2. Durchgang: Buegel fuer Querkraft UND Torsion
    quer = bemessung_querkraft(b["qs_feld"], C, S, x, Veinh, As_l_quer, c_l,
                               e.phi_buegel, e.n_schenkel, e.N_Ed,
                               V_bem_kN=V_bem, asw_torsion=asw_tor,
                               s_max_torsion=smax_tor)
    b["querkraft"], b["torsion"] = quer, tor
    det = quer.detail

    bl = _block("7. QUERKRAFTBEMESSUNG (GZT)", quer.normen)
    _z(bl, "Innerer Hebelarm  [NA NDP zu 6.2.3 (1)]:")
    _z(bl, "   z = min(0,9d ; d-2c_v,l ; d-c_v,l-30) = min({} ; {} ; {}) = {} mm"
       .format(_f(0.9 * b["qs_feld"].d, 1), _f(b["qs_feld"].d - 2 * c_l, 1),
               _f(b["qs_feld"].d - c_l - 30, 1), _f(quer.z, 1)))
    _z(bl, "")
    _z(bl, "Bemessungsquerkraft im Abstand d vom Auflagerrand  [EC2 6.2.1 (8)]:")
    for v in vbem:
        _z(bl, "   Auflager x = {} m ({}): V(Rand) = {} kN -> V_Ed(d) = {} kN"
           .format(_f(v["auflager_x"], 2), v["seite"], _f(v["V_rand"], 1), _f(v["V_d"], 1)))
    _z(bl, "")
    vc = det["V_Rdc"]
    _z(bl, "Tragfaehigkeit ohne Querkraftbewehrung  [EC2 6.2.2 (1), Gl. (6.2a/b) + NA]:")
    _z(bl, "   k = 1 + sqrt(200/d) = {} (<= 2,0)   rho_l = As_l/(bw d) = {} (<= 0,02)"
       .format(_f(vc["k"], 3), _f(vc["rho_l"], 5)))
    _z(bl, "   C_Rd,c = 0,15/gamma_C = {}   k1 = {}   [NA NDP zu 6.2.2 (1)]"
       .format(_f(vc["C_Rdc"], 3), _f(vc["k1"], 2)))
    _z(bl, "   v_min = (kappa1/gamma_C) k^1,5 sqrt(fck) = {} N/mm2 (kappa1 = {})  "
           "[Gl. (6.3aDE)]".format(_f(vc["v_min"], 3), _f(vc["kappa1"], 4)))
    _z(bl, "   V_Rd,c = {} kN   (massgebend {})".format(_f(quer.V_Rdc, 1), vc["massgebend"]))
    _z(bl, "   HINWEIS: rho_l wurde mit As_l = {} mm2 gebildet. Nach EC2 6.2.2 (1)"
       .format(_f(As_l_quer, 0)))
    _z(bl, "   muss diese Bewehrung >= (lbd + d) ueber den betrachteten Schnitt")
    _z(bl, "   hinausgefuehrt werden; in Bereichen mit gestaffelter Bewehrung pruefen.")
    _z(bl, "")
    ct = det["cot"]
    _z(bl, "Fachwerkmodell  [NA NDP zu 6.2.3 (2), Gl. (6.7aDE)/(6.7bDE)]:")
    _z(bl, "   V_Rd,cc = 0,5 * 0,48 * fck^(1/3) (1 - 1,2 sigma_cd/fcd) bw z = {} kN"
       .format(_f(quer.V_Rdcc, 1)))
    _z(bl, "   cot(theta) = (1,2 + 1,4 sigma_cd/fcd)/(1 - V_Rd,cc/V_Ed) = {} -> {} (1,0..3,0)"
       .format(_f(ct["cot_roh"], 3), _f(quer.cot_theta, 3)))
    _z(bl, "   theta = {} Grad".format(_f(quer.theta, 1)))
    _z(bl, "   V_Rd,max = alpha_cw bw z nu1 fcd/(cot+tan) = {} kN   (nu1 = {})  [Gl. (6.9)]"
       .format(_f(quer.V_Rdmax, 1), _f(det["V_Rdmax"]["nu1"], 3)))
    _nw(bl, "Druckstrebe V_Ed/V_Rd,max", det["V_bem"], quer.V_Rdmax,
        det["V_bem"] <= quer.V_Rdmax, "kN")
    _z(bl, "")
    _z(bl, "Querkraftbewehrung  [EC2 6.2.3 (3), Gl. (6.8)]:")
    asw_v = det["V_bem"] * 1e3 / (quer.z * S.fyd * quer.cot_theta) * 1000.0
    _z(bl, "   asw = V_Ed/(z fywd cot(theta)) ;  V_Ed = {} kN -> asw = {} cm2/m (alle Schenkel)"
       .format(_f(det["V_bem"], 1), _f(asw_v / 100.0, 2)))
    _z(bl, "   asw,min = rho_w,min bw = 0,16 fctm/fyk * bw = {} cm2/m   [NA NDP zu 9.2.2 (5)]"
       .format(_f(quer.asw_min / 100.0, 2)))
    _z(bl, "   s_max = {} mm ; s_max,quer = {} mm   [NA Tab. NA.9.1, Zeile {}]"
       .format(_f(quer.s_max, 0), _f(quer.s_max_q, 0), det["s_max"]["tabellenzeile"]))
    _z(bl, "   Versatzmass a_l = z cot(theta)/2 = {} mm   [EC2 9.2.1.3 (2), Gl. (9.2)]"
       .format(_f(quer.a_l, 1)))
    for hh in quer.hinweise:
        _z(bl, "   !! " + hh)
    B.append(bl)

    # --- Torsionsblock ----------------------------------------------------
    ehq = tor.detail.get("ehq", {})
    bl = _block("8. TORSIONSBEMESSUNG (GZT)", tor.normen)
    if abs(e.T_Ed) <= 1e-9:
        _z(bl, "Kein Torsionsmoment angesetzt (T_Ed = 0).")
        _z(bl, "")
        _z(bl, gleichgewichtstorsion_hinweis())
    else:
        art = "Gleichgewichtstorsion" if e.gleichgewichtstorsion \
            else "Vertraeglichkeitstorsion"
        _z(bl, "Art der Torsion: {}   [EC2 6.3.1 (2)]".format(art))
        if not e.gleichgewichtstorsion:
            _z(bl, "   Bei Vertraeglichkeitstorsion darf im GZT auf den Nachweis")
            _z(bl, "   verzichtet werden; es genuegt die Mindestbewehrung nach")
            _z(bl, "   9.2.2 und 9.2.3. Der Nachweis wird hier dennoch gefuehrt.")
        _z(bl, "")
        _z(bl, "Ersatzhohlquerschnitt  [EC2 6.3.1 (3)]:")
        _z(bl, "   A = {} mm2   u = {} mm".format(_f(ehq.get("A", 0), 0),
                                                  _f(ehq.get("u", 0), 0)))
        _z(bl, "   t_ef,i = A/u = {} mm  ->  gewaehlt {} mm  "
               "(>= 2 d1 = {} mm, <= {} mm)"
           .format(_f(ehq.get("t_ef_roh", 0), 1), _f(tor.t_ef, 1),
                   _f(ehq.get("t_min", 0), 1), _f(ehq.get("t_max", 0), 1)))
        _z(bl, "   A_k = (b - t_ef)(h - t_ef) = {} * {} = {} mm2"
           .format(_f(ehq.get("b_k", 0), 1), _f(ehq.get("h_k", 0), 1), _f(tor.A_k, 0)))
        _z(bl, "   u_k = 2[(b-t_ef) + (h-t_ef)] = {} mm".format(_f(tor.u_k, 0)))
        _z(bl, "")
        _z(bl, "Risstorsionsmoment  [EC2 6.3.2 (5), Gl. (6.31)]:")
        _z(bl, "   T_Rd,c = fctd t_ef 2 A_k = {} * {} * 2 * {} = {} kNm"
           .format(_f(C.fctd, 3), _f(tor.t_ef, 1), _f(tor.A_k, 0), _f(tor.T_Rdc, 2)))
        _z(bl, "   T_Ed = {} kNm {} T_Rd,c = {} kNm"
           .format(_f(tor.T_Ed, 2), ">" if tor.T_Ed > tor.T_Rdc else "<=",
                   _f(tor.T_Rdc, 2)))
        _z(bl, "")
        entb = tor.detail.get("entbehrlich", {})
        _z(bl, "Verzicht auf rechnerische Torsionsbewehrung  [NA NDP zu 6.3.2 (5)]:")
        _z(bl, "   Gl. (6.31aDE): T_Ed = {} kNm {} V_Ed bw/4,5 = {} kNm"
           .format(_f(tor.T_Ed, 2), "<=" if entb.get("bedingung_a") else ">",
                   _f(entb.get("grenze_a", 0), 2)))
        _z(bl, "   Gl. (6.31bDE): V_Ed [1+4,5 T_Ed/(V_Ed bw)] = {} kN {} V_Rd,c = {} kN"
           .format(_f(entb.get("V_wirksam", 0), 1),
                   "<=" if entb.get("bedingung_b") else ">", _f(quer.V_Rdc, 1)))
        _z(bl, "   -> Torsionsbewehrung {}"
           .format("ENTBEHRLICH" if entb.get("erfuellt") else "ERFORDERLICH"))
        _z(bl, "")
        trm = tor.detail.get("T_Rdmax", {})
        _z(bl, "Druckstrebe unter Torsion  [EC2 6.3.2 (4), Gl. (6.30)]:")
        _z(bl, "   nu = {} * nu_2 = {}   ({})   [NA NDP zu 6.3.2 (4)]"
           .format("0,75" if e.kastenquerschnitt else "0,525",
                   _f(trm.get("nu", 0), 4),
                   "Kastenquerschnitt" if e.kastenquerschnitt else "Vollquerschnitt"))
        _z(bl, "   T_Rd,max = 2 nu alpha_cw fcd A_k t_ef sin(theta) cos(theta) = {} kNm"
           .format(_f(tor.T_Rdmax, 2)))
        _z(bl, "")
        _z(bl, "Interaktion Torsion / Querkraft  [EC2 6.3.2 (4), Gl. (6.29)]:")
        _z(bl, "   T_Ed/T_Rd,max + V_Ed/V_Rd,max = {}/{} + {}/{} = {} <= 1,0"
           .format(_f(tor.T_Ed, 1), _f(tor.T_Rdmax, 1), _f(tor.V_Ed, 1),
                   _f(quer.V_Rdmax, 1), _f(tor.interaktion, 3)))
        _nw(bl, "Interaktion T/V (Druckstreben)", tor.interaktion, 1.0,
            tor.interaktion <= 1.0, "-")
        _z(bl, "")
        _z(bl, "Erforderliche Torsionsbewehrung:")
        _z(bl, "   Buegel je AUSSENSCHENKEL:")
        _z(bl, "      asw,T = T_Ed/(2 A_k fywd cot(theta)) = {} cm2/m"
           .format(_f(tor.asw_je_schenkel / 100.0, 2)))
        _z(bl, "   Laengsbewehrung  [EC2 6.3.2 (3), Gl. (6.28)]:")
        _z(bl, "      sum(Asl) = T_Ed cot(theta) u_k /(2 A_k fyd) = {} mm2"
           .format(_f(tor.asl_gesamt, 0)))
        _z(bl, "      gleichmaessig ueber u_k verteilen, mindestens ein Stab je Ecke;")
        _z(bl, "      gewaehlt {} Staebe -> {} mm2 je Stab   [EC2 9.2.3]"
           .format(tor.n_laengsstaebe, _f(tor.asl_je_ecke, 0)))
        _z(bl, "   Groesster Buegelabstand [EC2 9.2.3 (3)]:")
        smt = tor.detail.get("s_max", {})
        _z(bl, "      s <= min(u_k/8 ; kleinste Abmessung) = min({} ; {}) = {} mm"
           .format(_f(smt.get("u_k_8", 0), 0), _f(smt.get("min_abmessung", 0), 0),
                   _f(tor.s_max, 0)))
        _z(bl, "      zusammen mit Tab. NA.9.1 massgebend: s <= {} mm".format(_f(quer.s_max, 0)))
        for hh in tor.hinweise:
            _z(bl, "   !! " + hh)
    B.append(bl)

    # --- Buegelbereiche (Querkraft + Torsion) ----------------------------
    bl = _block("9. BUEGELBEWEHRUNG (Querkraft + Torsion)",
                [ref("VRds"), ref("rho_w_min"), ref("s_max"), ref("torsion_konstr")])
    A_stab = det["A_stab"]
    _z(bl, "Massgebend ist der AUSSENSCHENKEL des geschlossenen Buegels:")
    _z(bl, "   asw,Schenkel = asw,V/{} + asw,T".format(e.n_schenkel))
    _z(bl, "   max. erforderlich = {} cm2/m je Schenkel"
       .format(_f(float(np.max(quer.asw_schenkel)) / 100.0, 2)))
    _z(bl, "   Buegel phi {} mit {} Schenkeln: A_Schenkel = {} mm2"
       .format(_f(e.phi_buegel, 0), e.n_schenkel, _f(A_stab, 1)))
    if abs(e.T_Ed) > 1e-9:
        _z(bl, "   Torsion erfordert GESCHLOSSENE Buegel mit Uebergreifung oder "
               "Haken [EC2 9.2.3 (2)]")
    _z(bl, "")
    _z(bl, "   VERLEGEBEREICHE:")
    for zz in quer.bereiche:
        asw_z = A_stab * 1000.0 / zz["s"] / 100.0
        _z(bl, "      x = {} .. {} m :  Buegel phi{}/{:.0f} mm  ({} Stk, "
               "asw = {} cm2/m je Schenkel)"
           .format(_f(zz["x1"], 2), _f(zz["x2"], 2), _f(e.phi_buegel, 0), zz["s"],
                   zz["n"], _f(asw_z, 2)))
    B.append(bl)

    # =====================================================================
    # 10. RISSBREITE (GZG)
    # =====================================================================
    M_qs_pos = float(np.max(qs_last.M))
    M_qs_neg = float(np.min(qs_last.M))
    if abs(M_qs_neg) > M_qs_pos:
        M_riss, As_riss, qs_riss = abs(M_qs_neg), As_oben_vorh, qs_stuetz
        phi_riss, seite = e.phi_laengs_oben, "Oberseite (Stuetzzug)"
    else:
        M_riss, As_riss, qs_riss = M_qs_pos, As1_vorh, b["qs_feld"]
        phi_riss, seite = e.phi_laengs, "Unterseite (Feldzug)"
    riss = rissbreite(qs_riss, C, S, M_riss, As_riss, phi_riss, 0.0,
                      e.phi_kriech, bd["w_max"])
    b["rissbreite"] = riss

    bl = _block("10. RISSBREITENBEGRENZUNG (GZG)",
                [ref("wk"), ref("eps_sm"), ref("sr_max"), ref("hc_eff"),
                 ref("w_max"), ref("komb_QS")])
    _z(bl, "Quasi-staendiges Moment M_qs = {} kNm  [EC0 Gl. (6.16b)]".format(_f(M_riss, 2)))
    _z(bl, "Nachgewiesene Seite: {} ; As = {} mm2 aus phi {}"
       .format(seite, _f(As_riss, 0), _f(phi_riss, 0)))
    _z(bl, "Zustand II, alpha_e = Es/Ec,eff mit phi = {}:".format(_f(e.phi_kriech, 2)))
    _z(bl, "   Ec,eff = Ecm/(1+phi) = {} N/mm2   alpha_e = {}"
       .format(_f(C.Ecm / (1 + e.phi_kriech), 0), _f(riss.get("alpha_e", 0), 2)))
    _z(bl, "   x_II = {} mm   I_II = {:.4g} mm4".format(_f(riss.get("x_II", 0), 1),
                                                        riss.get("I_II", 0)))
    _z(bl, "   sigma_s = alpha_e M (d-x)/I_II = {} N/mm2".format(_f(riss.get("sigma_s", 0), 1)))
    _z(bl, "   hc,ef = min(2,5(h-d) ; (h-x)/3 ; h/2) = {} mm   [EC2 7.3.2 (3)]"
       .format(_f(riss.get("hc_ef", 0), 1)))
    _z(bl, "   Ac,eff = {} mm2   rho_p,eff = {}".format(_f(riss.get("Ac_eff", 0), 0),
                                                        _f(riss.get("rho_p_eff", 0), 5)))
    _z(bl, "   (eps_sm - eps_cm) = {:.5f}   (kt = {})   [Gl. (7.9)]"
       .format(riss.get("eps_sm_cm", 0), _f(riss.get("kt", 0.4), 1)))
    if riss.get("min_massgebend"):
        _z(bl, "      -> der Mindestwert 0,6 sigma_s/Es ist massgebend")
    _z(bl, "   sr,max = min(phi/(3,6 rho_p,eff) ; sigma_s phi/(3,6 fct,eff))")
    _z(bl, "          = min({} ; {}) = {} mm   [NA Gl. (7.11DE)]"
       .format(_f(riss.get("sr_a", 0), 1), _f(riss.get("sr_b", 0), 1),
               _f(riss.get("sr_max", 0), 1)))
    _z(bl, "   wk = sr,max (eps_sm - eps_cm) = {} mm <= w_max = {} mm   [Gl. (7.8)]"
       .format(_f(riss.get("wk", 0), 3), _f(bd["w_max"], 2)))
    _nw(bl, "Rissbreite wk", riss.get("wk", 0), bd["w_max"], riss.get("ok", True), "mm")
    B.append(bl)

    # =====================================================================
    # 11. DURCHBIEGUNG (GZG)
    # =====================================================================
    felder, kragarme = list(traeger.felder), []
    xa = [a.x for a in auflager]
    if xa[-1] < e.L - 1e-6:
        kragarme.append(len(felder))
        felder.append((xa[-1], e.L))
    if not felder:
        kragarme.append(0)
        felder.append((xa[0], e.L))
    K = 0.4 if (kragarme and len(felder) == 1) else (1.0 if len(felder) == 1 else 1.3)
    As_zug = As_oben_vorh if (kragarme and len(felder) == 1) else As1_vorh
    As_druck = As1_vorh if As_zug != As1_vorh else As_oben_vorh
    durchb = nachweis_durchbiegung(qs_last.x, qs_last.M, b["qs_feld"], C, S,
                                   As_zug, As_druck, felder, e.phi_kriech,
                                   e.eps_cs, e.grenze_durchbiegung, K=K,
                                   sigma_s=max(riss.get("sigma_s", 310.0), 1.0),
                                   kragarme=tuple(kragarme))
    b["durchbiegung"] = durchb

    bl = _block("11. VERFORMUNGSNACHWEIS (GZG)",
                [ref("durchbiegung_rech"), ref("durchbiegung_ld"),
                 ref("durchbiegung_NA"), ref("kriechen"), ref("komb_QS")])
    _z(bl, "Berechnung durch Interpolation Zustand I / II  [EC2 7.4.3, Gl. (7.18)/(7.19)]")
    _z(bl, "   phi(inf,t0) = {}   eps_cs = {:.1e}   beta = 0,5 (Langzeiteinwirkung)"
       .format(_f(e.phi_kriech, 2), e.eps_cs))
    for f in durchb:
        _z(bl, "   {} {:.2f} - {:.2f} m (L = {} mm):"
           .format("Kragarm" if f["kragarm"] else "Feld", f["feld"][0], f["feld"][1],
                   _f(f["L"], 0)))
        _z(bl, "      w_max = {} mm  bei x = {} m".format(_f(f["w_max"], 2),
                                                          _f(f["x_max"], 2)))
        _z(bl, "      Grenzwert L/{:.0f} = {} mm  ->  eta = {}"
           .format(e.grenze_durchbiegung, _f(f["w_grenz"], 2), _f(f["ausnutzung"], 3)))
        _z(bl, "      Biegeschlankheit l/d = {} ; zul. Gl. (7.16) = {} ; NA = {}"
           .format(_f(f["ld_vorh"], 1), _f(f["ld_zul"], 1), _f(f["ld_NA"], 1)))
        _nw(bl, "Durchbiegung {:.1f}-{:.1f} m".format(*f["feld"]), f["w_max"],
            f["w_grenz"], f["ok"], "mm")
    B.append(bl)

    # =====================================================================
    # 12. VERANKERUNG
    # =====================================================================
    bl = _block("12. VERANKERUNG UND BEWEHRUNGSFUEHRUNG",
                [ref("fbd"), ref("lb_rqd"), ref("lbd"), ref("l0"),
                 ref("versatzmass"), ref("auflagerkraft")])
    ver = verankerungslaenge(C, S, e.phi_laengs, guter_verbund=True,
                             As_erf=(biege_feld.As1 if biege_feld else None),
                             As_vorh=(As1_vorh if biege_feld else None))
    _z(bl, "Verbund  [EC2 8.4.2 (2), Gl. (8.2)]:")
    _z(bl, "   fbd = 2,25 eta1 eta2 fctd = 2,25 * {} * {} * {} = {} N/mm2"
       .format(_f(ver["eta1"], 1), _f(ver["eta2"], 2), _f(C.fctd, 3), _f(ver["fbd"], 2)))
    _z(bl, "   lb,rqd = (phi/4)(sigma_sd/fbd) = {} mm   (sigma_sd = {} N/mm2)  [Gl. (8.3)]"
       .format(_f(ver["lb_rqd"], 0), _f(ver["sigma_sd"], 1)))
    _z(bl, "   lbd = {} mm  (lb,min = {} mm)   [Gl. (8.4)/(8.6)]"
       .format(_f(ver["lbd"], 0), _f(ver["lb_min"], 0)))
    ue = uebergreifungslaenge(C, S, e.phi_laengs, alpha6=1.4)
    _z(bl, "   l0 (Stoss, 50 % im Schnitt, alpha6 = 1,4) = {} mm   [EC2 8.7.3]"
       .format(_f(ue["l0"], 0)))
    _z(bl, "")
    V_auf = max([v["V_rand"] for v in vbem]) if vbem else float(np.max(Veinh))
    fa = verankerungskraft_endauflager(V_auf, quer.a_l, quer.z, e.N_Ed)
    _z(bl, "Endauflager  [EC2 9.2.1.4 (2), Gl. (9.3)]:")
    _z(bl, "   F_Ed = |V_Ed| a_l/z + N_Ed = {} * {}/{} = {} kN"
       .format(_f(V_auf, 1), _f(quer.a_l, 1), _f(quer.z, 1), _f(fa["F_Ed"], 1)))
    n_auf = max(2, int(math.ceil(0.25 * n_unten)))
    va = verankerung_endauflager(C, S, e.phi_laengs, fa["F_Ed"],
                                 stabflaeche_n(n_auf, e.phi_laengs))
    _z(bl, "   Zum Auflager gefuehrt: {} phi {} (>= 25 % der Feldbewehrung, "
           "EC2 9.2.1.4 (1))".format(n_auf, _f(e.phi_laengs, 0)))
    _z(bl, "   sigma_sd = F_Ed/As = {} N/mm2 -> Verankerungslaenge = {} mm "
           "(2/3 lbd, direkte Lagerung)"
       .format(_f(va["sigma_sd"], 1), _f(va["l_verankerung"], 0)))
    _z(bl, "   Versatzmass a_l = {} mm: die Staebe sind um a_l ueber den "
           "rechnerischen".format(_f(quer.a_l, 1)))
    _z(bl, "   Endpunkt hinaus zu fuehren  [EC2 9.2.1.3 (2), Gl. (9.2)]")
    if abs(e.T_Ed) > 1e-9 and tor.erforderlich:
        _z(bl, "")
        _z(bl, "   TORSION: die Torsionslaengsbewehrung ({} mm2 gesamt) ist ueber die"
           .format(_f(tor.asl_gesamt, 0)))
        _z(bl, "   volle Bauteillaenge zu fuehren und an den Enden voll zu verankern")
        _z(bl, "   [EC2 9.2.3]. Der auf die Zugzone entfallende Anteil ist zur")
        _z(bl, "   Biegebewehrung zu ADDIEREN.")
    b["verankerung"], b["verankerung_auflager"] = ver, va
    B.append(bl)

    # =====================================================================
    # 13. ZUSAMMENFASSUNG
    # =====================================================================
    bl = _block("13. ZUSAMMENSTELLUNG DER NACHWEISE", [])
    alle_nw = []
    for bb in B:
        alle_nw.extend(bb["nachweise"])
    for c in alle_nw:
        eta = (c["wert"] / c["grenzwert"]) if (c["vergleich"] == "<=" and c["grenzwert"]) \
            else ((c["grenzwert"] / c["wert"]) if c["wert"] else 0.0)
        _z(bl, "   [{}]  {:<42s} {:>10s} {} {:>10s} {}   eta = {}"
           .format("OK" if c["ok"] else "NEIN", c["name"], _f(c["wert"], 2),
                   c["vergleich"], _f(c["grenzwert"], 2), c["einheit"], _f(eta, 3)))
        if not c["ok"]:
            b["ok_gesamt"] = False
    _z(bl, "")
    _z(bl, "   GESAMTERGEBNIS: " + ("NACHWEISE ERFUELLT" if b["ok_gesamt"]
                                    else "NACHWEISE NICHT ERFUELLT"))
    _z(bl, "")
    _z(bl, "   GEWAEHLTE BEWEHRUNG")
    _z(bl, "      unten  : {} phi {} = {} mm2 ({} Lage/n)"
       .format(n_unten, _f(e.phi_laengs, 0), _f(As1_vorh, 0), b.get("lagen_unten", 1)))
    _z(bl, "      oben   : {} phi {} = {} mm2"
       .format(n_oben, _f(e.phi_laengs_oben, 0), _f(As_oben_vorh, 0)))
    if abs(e.T_Ed) > 1e-9 and tor.erforderlich:
        _z(bl, "      Torsion: zusaetzlich {} mm2 Laengsbewehrung, gleichmaessig auf "
               "{} Staebe".format(_f(tor.asl_gesamt, 0), tor.n_laengsstaebe))
        _z(bl, "               entlang u_k = {} mm verteilt".format(_f(tor.u_k, 0)))
    for zz in quer.bereiche:
        _z(bl, "      Buegel : x = {} .. {} m  phi{}/{:.0f} mm ({} Schenkel{})"
           .format(_f(zz["x1"], 2), _f(zz["x2"], 2), _f(e.phi_buegel, 0), zz["s"],
                   e.n_schenkel, ", geschlossen" if abs(e.T_Ed) > 1e-9 else ""))
    B.append(bl)

    b["nachweise"] = alle_nw
    return b


# ---------------------------------------------------------------------------
# Textbericht
# ---------------------------------------------------------------------------
def bericht_text(b, mit_normen=True):
    """Erzeugt den vollstaendigen Bericht als Text."""
    e = b["eingabe"]
    out = []
    out.append("=" * 100)
    out.append("BEMESSUNG EINES STAHLBETONBALKENS")
    out.append("Norm: DIN EN 1992-1-1:2011-01 + DIN EN 1992-1-1/NA:2013-04 "
               "(Eurocode 2 mit deutschem NA)")
    out.append("=" * 100)
    out.append("Traeger {:.2f} m | Querschnitt {} | {} | {}"
               .format(e.L, b["qs_feld"].beschreibung(), e.betonklasse, e.stahlsorte))
    out.append("")
    for bl in b["bloecke"]:
        out.append("")
        out.append("-" * 100)
        out.append(bl["titel"])
        if bl["normen"]:
            out.append("   Normgrundlage:")
            gesehen = set()
            for n in bl["normen"]:
                if n.id in gesehen:
                    continue
                gesehen.add(n.id)
                out.append("      - " + n.lang())
        out.append("-" * 100)
        out.extend(bl["zeilen"])
    if mit_normen:
        out.append("")
        out.append("=" * 100)
        out.append(normentabelle())
    out.append("")
    out.append("=" * 100)
    out.append("HINWEIS: automatisch erzeugte Ergebnisse. Sie sind von einem")
    out.append("verantwortlichen Ingenieur gegen die gueltige Ausgabe von")
    out.append("DIN EN 1992-1-1 und ihres Nationalen Anhangs zu pruefen. Nicht")
    out.append("erfasst: Durchstanzen (6.4), Ermuedung (6.8), Vorspannung (5.10),")
    out.append("aussergewoehnliche und seismische Bemessungssituationen, Brandfall.")
    out.append("=" * 100)
    return "\n".join(out)
