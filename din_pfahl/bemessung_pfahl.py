# -*- coding: utf-8 -*-
"""
Gesamtbemessung eines BOHRPFAHLS.

Ablauf (mit der jeweils angewendeten Normstelle):

  1. Baustoffe                      EC2 3.1 / 3.2 + DIN 488-1
  2. Geometrie und Betondeckung     DIN EN 1536, 7.6.2
  3. Einwirkungen                   EC0 6.4.3.2 Gl. (6.10) + DIN 1054 Tab. A 2.1
  4. Axiale Tragfaehigkeit          DIN EN 1997-1 7.6 + DIN 1054 A 7.6
  5. Horizontalbelastung            EA-Pfaehle 6.3 (Bettungsmodulverfahren)
  6. Laengsbewehrung (M-N)          EC2 6.1 (Interaktionsdiagramm)
  7. Mindestbewehrung               DIN EN 1536 7.6.3 Tab. 4 + NA zu EC2 9.5.2
  8. Querbewehrung (Wendel)         EC2 6.2 + 9.5.3 + DIN EN 1536 7.6.4
  9. Knicknachweis                  EA-Pfaehle 4.7 / EC2 5.8
 10. Konstruktive Durchbildung      DIN EN 1536 7.6
"""

import math
from dataclasses import dataclass, field

import numpy as np

from din_balken.baustoffe import Beton, Betonstahl, stabflaeche

from .normen_pfahl import ref, normentabelle
from .kreisquerschnitt import (Kreisquerschnitt, interaktionsdiagramm,
                               M_Rd_bei_N, erforderliche_bewehrung,
                               mindestbewehrung_pfahl,
                               mindestbewehrung_druckglied, betondeckung_pfahl,
                               konstruktive_pruefung, schnittgroessen_bei_x)
from .bettung import (Bodenschicht, pfahl_horizontal, bettungsmodul,
                      knicklast_gebettet)
from .tragfaehigkeit import (axiale_tragfaehigkeit, widerstands_setzungs_linie,
                             pfahlkopfsetzung, GAMMA_R)


# ---------------------------------------------------------------------------
# Eingabedaten
# ---------------------------------------------------------------------------
@dataclass
class EingabePfahl:
    # --- Geometrie -------------------------------------------------------
    D: float = 900.0                # Pfahldurchmesser [mm]
    L: float = 15.0                 # Pfahllaenge [m]
    unter_stuetzfluessigkeit: bool = False

    # --- Baustoffe -------------------------------------------------------
    betonklasse: str = "C25/30"     # DIN EN 1536, 6.3: i.d.R. >= C25/30
    stahlsorte: str = "B500B"
    d_g: float = 16.0               # Groesstkorn [mm]

    # --- Baugrund --------------------------------------------------------
    schichten: list = field(default_factory=list)   # list[Bodenschicht]
    q_b_k: float = 1500.0           # charakt. Spitzendruck [kN/m2] (EA-Pfaehle)
    mantel_ab_tiefe: float = 0.0    # ab dieser Tiefe wird Mantelreibung angesetzt
    situation: str = "BS-P"         # DIN 1054 Bemessungssituation

    # --- Einwirkungen (Bemessungswerte am Pfahlkopf) ---------------------
    N_Ed: float = -2500.0           # [kN] Laengskraft, DRUCK NEGATIV
    N_k: float = -1850.0            # [kN] charakteristisch (fuer die Setzung)
    H_Ed: float = 150.0             # [kN] Horizontalkraft
    M_Ed_kopf: float = 0.0          # [kNm] Kopfmoment
    kopf: str = "frei"              # "frei" | "eingespannt"
    N_Ed_zug: float = 0.0           # [kN] Zugkraft (positiv), falls vorhanden

    # --- Bewehrung -------------------------------------------------------
    phi_l: float = 20.0             # Laengsstabdurchmesser [mm]
    n_l: int = 0                    # 0 -> automatisch ermitteln
    phi_w: float = 10.0             # Wendeldurchmesser [mm]
    laenge_bewehrungskorb: float = 0.0   # 0 -> volle Pfahllaenge

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
                                ok=bool(ok), einheit=einheit, vergleich=vergleich))


def _f(v, n=2):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "-"
    return ("{:." + str(n) + "f}").format(v)


# ---------------------------------------------------------------------------
# Querkraftnachweis des Kreisquerschnitts
# ---------------------------------------------------------------------------
def querkraft_kreis(qs, beton, stahl, V_Ed, N_Ed):
    """
    Querkraftnachweis fuer den Kreisquerschnitt.
    DIN EN 1992-1-1, 6.2 + NA; EC2 enthaelt keine expliziten Regeln fuer
    Kreisquerschnitte, daher werden angesetzt:

        b_w = D                    (Ersatzbreite)
        d   = D/2 + D_s/pi         (Schwerpunkt der gezogenen Stabhaelfte)
        z   = 0,9 d
        A_sl = As/2                (gezogene Stabhaelfte)

    Die Wendel wird wie ein zweischnittiger Buegel behandelt (2 Schenkel).
    """
    D, d = qs.D, qs.d_eff
    bw = D
    z = 0.9 * d
    A_sl = qs.As_ges / 2.0

    k = min(1.0 + math.sqrt(200.0 / d), 2.0)
    rho_l = min(A_sl / (bw * d), 0.02)
    # Laengsdruckspannung (Druck POSITIV nach EC2), <= 0,2 fcd
    sigma_cp = min(max(-N_Ed, 0.0) * 1.0e3 / qs.Ac, 0.2 * beton.fcd)
    C_Rdc, k1 = 0.15 / beton.gamma_c, 0.12
    kappa1 = 0.0525 if d <= 600.0 else (0.0375 if d >= 800.0
                                        else 0.0525 + (0.0375 - 0.0525) * (d - 600.0) / 200.0)
    v_min = (kappa1 / beton.gamma_c) * k ** 1.5 * math.sqrt(beton.fck)
    v_a = C_Rdc * k * (100.0 * rho_l * beton.fck) ** (1.0 / 3.0) + k1 * sigma_cp
    v_b = v_min + k1 * sigma_cp
    V_Rdc = max(v_a, v_b) * bw * d / 1.0e3

    # Druckstrebenneigung nach NA Gl. (6.7aDE)
    sigma_cd = N_Ed * 1.0e3 / qs.Ac                 # Druck negativ
    V_Rdcc = 0.5 * 0.48 * beton.fck ** (1 / 3) * (1.0 - 1.2 * sigma_cd / beton.fcd) \
        * bw * z / 1.0e3
    nenner = 1.0 - V_Rdcc / max(abs(V_Ed), 1e-9)
    cot = 3.0 if nenner <= 1e-6 else (1.2 + 1.4 * sigma_cd / beton.fcd) / nenner
    cot = min(max(cot, 1.0), 3.0)

    nu2 = 1.0 if beton.fck <= 50.0 else min(1.1 - beton.fck / 500.0, 1.0)
    nu1 = 0.75 * nu2
    V_Rdmax = bw * z * nu1 * beton.fcd / (cot + 1.0 / cot) / 1.0e3

    asw_erf = abs(V_Ed) * 1.0e3 / (z * stahl.fyd * cot) * 1000.0    # mm2/m
    rho_w_min = 0.16 * beton.fctm / stahl.fyk
    asw_min = rho_w_min * bw * 1000.0                                # mm2/m
    asw = max(asw_erf, asw_min) if abs(V_Ed) > V_Rdc else asw_min

    A_wendel = stabflaeche(qs.phi_w)
    s_erf = 2.0 * A_wendel * 1000.0 / max(asw, 1e-9)                 # mm

    # Groesste Abstaende: EC2 9.5.3 (Druckglieder) + DIN EN 1536, 7.6.4
    s_max_ec2 = min(20.0 * qs.phi_l, qs.D, 400.0)
    s_max_1536 = 400.0
    s_min_1536 = 100.0
    s_max = min(s_max_ec2, s_max_1536)
    s_gew = max(min(math.floor(min(s_erf, s_max) / 25.0) * 25.0, s_max), s_min_1536)

    return dict(bw=bw, d=d, z=z, A_sl=A_sl, k=k, rho_l=rho_l, sigma_cp=sigma_cp,
                v_min=v_min, V_Rdc=V_Rdc, V_Rdcc=V_Rdcc, cot=cot,
                theta=math.degrees(math.atan(1.0 / cot)), nu1=nu1,
                V_Rdmax=V_Rdmax, asw_erf=asw_erf, asw_min=asw_min, asw=asw,
                A_wendel=A_wendel, s_erf=s_erf, s_max=s_max,
                s_max_ec2=s_max_ec2, s_min=s_min_1536, s_gewaehlt=s_gew,
                ok_druckstrebe=abs(V_Ed) <= V_Rdmax,
                bewehrung_erforderlich=abs(V_Ed) > V_Rdc,
                normen=[ref("qs_querkraft"), ref("qs_VRdc"), ref("qs_cot"),
                        ref("qs_wendel"), ref("pfahl_quer")])


# ---------------------------------------------------------------------------
# Gesamtbemessung
# ---------------------------------------------------------------------------
def bemessung_pfahl(e):
    """Fuehrt die vollstaendige Pfahlbemessung durch und liefert einen Bericht."""
    C, S = e.beton(), e.stahl()
    b = dict(eingabe=e, beton=C, stahl=S, bloecke=[], ok_gesamt=True)
    B = b["bloecke"]
    D_m = e.D / 1000.0

    # =====================================================================
    # 1. BAUSTOFFE
    # =====================================================================
    bl = _block("1. BAUSTOFFE", [ref("qs_fcd"), ref("pfahl_beton")])
    _z(bl, "Beton {}  [EC2 3.1.2, Tab. 3.1 ; DIN EN 1536, 6.3]".format(C.klasse))
    _z(bl, "   fck = {} N/mm2   fctm = {} N/mm2   Ecm = {} N/mm2"
       .format(_f(C.fck, 1), _f(C.fctm, 2), _f(C.Ecm, 0)))
    _z(bl, "   fcd = alpha_cc fck/gamma_C = 0,85 * {} / 1,50 = {} N/mm2   [Gl. (3.15)]"
       .format(_f(C.fck, 0), _f(C.fcd, 2)))
    _z(bl, "   eps_c2 = {} permil   eps_cu2 = {} permil".format(_f(C.eps_c2, 2),
                                                                _f(C.eps_cu2, 2)))
    _z(bl, "Betonstahl {}: fyk = {} N/mm2   fyd = {} N/mm2   eps_ud = {} permil"
       .format(S.sorte, _f(S.fyk, 0), _f(S.fyd, 2), _f(S.eps_ud, 0)))
    _z(bl, "")
    _z(bl, "Hinweis DIN EN 1536, 6.3: Pfahlbeton mit weicher bis fliessfaehiger")
    _z(bl, "Konsistenz, Groesstkorn <= 32 mm, Mindestzementgehalt beachten.")
    B.append(bl)

    # =====================================================================
    # 2. GEOMETRIE UND BETONDECKUNG
    # =====================================================================
    bdk = betondeckung_pfahl(e.D, e.unter_stuetzfluessigkeit)
    c_nom = bdk["c_nom"]
    bl = _block("2. GEOMETRIE UND BETONDECKUNG", [ref("pfahl_deckung")])
    _z(bl, "Bohrpfahl D = {} mm , L = {} m".format(_f(e.D, 0), _f(e.L, 2)))
    _z(bl, "   Ac = pi D^2/4 = {} mm2 = {} m2".format(_f(math.pi * e.D ** 2 / 4, 0),
                                                      _f(math.pi * D_m ** 2 / 4, 4)))
    _z(bl, "   U  = pi D = {} m   A_b = {} m2".format(_f(math.pi * D_m, 3),
                                                      _f(math.pi * D_m ** 2 / 4, 4)))
    _z(bl, "   c_nom = {} mm   [DIN EN 1536, 7.6.2: >= 60 mm bei D >= 0,6 m ; "
           ">= 50 mm bei D < 0,6 m]".format(_f(c_nom, 0)))
    if e.unter_stuetzfluessigkeit:
        _z(bl, "   Betonage unter Stuetzfluessigkeit -> c_nom >= 75 mm")
    B.append(bl)

    # =====================================================================
    # 3. EINWIRKUNGEN
    # =====================================================================
    bl = _block("3. EINWIRKUNGEN AM PFAHLKOPF",
                [ref("komb_GZT"), ref("gamma_F"), ref("geo_einwirkung")])
    _z(bl, "Bemessungswerte (Grundkombination GZT, EC0 Gl. (6.10);")
    _z(bl, "Teilsicherheitsbeiwerte DIN 1054, Tab. A 2.1: gamma_G = 1,35 / "
           "gamma_Q = 1,50):")
    _z(bl, "   N_Ed = {} kN   ({})".format(_f(e.N_Ed, 1),
                                           "Druck" if e.N_Ed < 0 else "Zug"))
    _z(bl, "   H_Ed = {} kN".format(_f(e.H_Ed, 1)))
    _z(bl, "   M_Ed (Kopf) = {} kNm".format(_f(e.M_Ed_kopf, 1)))
    _z(bl, "   Kopfausbildung: {}".format(
        "frei drehbar" if e.kopf == "frei" else "drehstarr eingespannt"))
    if e.N_Ed_zug > 0:
        _z(bl, "   Zugkraft N_Ed,Zug = {} kN".format(_f(e.N_Ed_zug, 1)))
    _z(bl, "   N_k (charakteristisch, fuer die Setzung) = {} kN".format(_f(e.N_k, 1)))
    B.append(bl)

    # =====================================================================
    # 4. AXIALE TRAGFAEHIGKEIT
    # =====================================================================
    trag = axiale_tragfaehigkeit(D_m, e.schichten, e.q_b_k, abs(min(e.N_Ed, 0.0)),
                                 e.N_Ed_zug, e.situation, L_pfahl=e.L,
                                 mantel_ab_tiefe=e.mantel_ab_tiefe)
    wsl = widerstands_setzungs_linie(D_m, trag.R_b_k, trag.R_s_k)
    setzung = pfahlkopfsetzung(wsl, abs(min(e.N_k, 0.0)))
    b["tragfaehigkeit"], b["wsl"], b["setzung"] = trag, wsl, setzung

    bl = _block("4. AXIALE PFAHLTRAGFAEHIGKEIT (GEO-2)", trag.normen)
    _z(bl, "Charakteristischer Widerstand  [DIN EN 1997-1, 7.6.2.3, Gl. (7.8)]:")
    _z(bl, "   R_c,k = R_b,k + R_s,k = q_b,k A_b + sum(q_s,k,i A_s,i)")
    _z(bl, "")
    _z(bl, "   {:<18s}{:>8s}{:>8s}{:>10s}{:>12s}{:>12s}"
       .format("Schicht", "z_o [m]", "z_u [m]", "q_s,k", "A_s [m2]", "R_s,k [kN]"))
    for a in trag.anteile:
        _z(bl, "   {:<18s}{:>8s}{:>8s}{:>10s}{:>12s}{:>12s}"
           .format(a["name"][:18], _f(a["z_o"], 2), _f(a["z_u"], 2),
                   _f(a["q_s_k"], 0), _f(a["A_s"], 3), _f(a["R_s_k"], 0)))
    _z(bl, "   {:<18s}{:>50s}".format("Summe Mantel", _f(trag.R_s_k, 0)))
    _z(bl, "")
    _z(bl, "   R_b,k = q_b,k A_b = {} * {} = {} kN"
       .format(_f(e.q_b_k, 0), _f(math.pi * D_m ** 2 / 4, 4), _f(trag.R_b_k, 0)))
    _z(bl, "   R_c,k = {} + {} = {} kN".format(_f(trag.R_b_k, 0), _f(trag.R_s_k, 0),
                                               _f(trag.R_c_k, 0)))
    _z(bl, "")
    g = trag.gamma
    _z(bl, "Bemessungswert  [DIN 1054, A 7.6.2.2, Tab. A 2.3, {}]:".format(e.situation))
    _z(bl, "   gamma_b = {} ; gamma_s = {} ; gamma_s,t = {}"
       .format(_f(g["gamma_b"], 2), _f(g["gamma_s"], 2), _f(g["gamma_s_t"], 2)))
    _z(bl, "   R_c,d = R_b,k/gamma_b + R_s,k/gamma_s = {}/{} + {}/{} = {} kN"
       .format(_f(trag.R_b_k, 0), _f(g["gamma_b"], 2), _f(trag.R_s_k, 0),
               _f(g["gamma_s"], 2), _f(trag.R_c_d, 0)))
    _z(bl, "")
    _z(bl, "   NACHWEIS: F_c,d = {} kN <= R_c,d = {} kN   (eta = {})"
       .format(_f(trag.F_c_d, 0), _f(trag.R_c_d, 0), _f(trag.ausnutzung, 3)))
    _nw(bl, "Axiale Tragfaehigkeit F_c,d/R_c,d", trag.F_c_d, trag.R_c_d,
        trag.F_c_d <= trag.R_c_d, "kN")
    if e.N_Ed_zug > 0:
        _z(bl, "   ZUGPFAHL: R_t,d = R_s,k/gamma_s,t = {} kN >= F_t,d = {} kN"
           .format(_f(trag.R_t_d, 0), _f(e.N_Ed_zug, 0)))
        _nw(bl, "Zugtragfaehigkeit F_t,d/R_t,d", e.N_Ed_zug, trag.R_t_d,
            e.N_Ed_zug <= trag.R_t_d, "kN")
    _z(bl, "")
    _z(bl, "Setzung aus der Widerstands-Setzungs-Linie  [EA-Pfaehle, 5.4.5]:")
    _z(bl, "   s_sg = 0,5 R_s,k[MN] + 0,5 cm <= 3 cm  ->  {} mm".format(_f(wsl["s_sg"], 1)))
    _z(bl, "   s_g  = 0,10 D = {} mm".format(_f(wsl["s_g"], 1)))
    _z(bl, "   s(N_k = {} kN) = {} mm".format(_f(abs(min(e.N_k, 0.0)), 0), _f(setzung, 1)))
    for hh in trag.hinweise:
        _z(bl, "   !! " + hh)
    B.append(bl)

    # =====================================================================
    # 5. HORIZONTALBELASTUNG
    # =====================================================================
    EI_brutto = C.Ecm * (math.pi * e.D ** 4 / 64.0) / 1.0e9      # kNm2
    EI_riss = 0.5 * EI_brutto        # gerissener Zustand, Naeherung
    hor = pfahl_horizontal(e.L, D_m, EI_riss, e.schichten, e.H_Ed,
                           e.M_Ed_kopf, e.kopf)
    b["horizontal"] = hor
    M_Ed_max = max(hor.M_max, abs(e.M_Ed_kopf))
    V_Ed_max = max(hor.V_max, abs(e.H_Ed))

    bl = _block("5. HORIZONTAL BELASTETER PFAHL (Bettungsmodulverfahren)",
                hor.normen)
    _z(bl, "Bettung  [EA-Pfaehle, 6.3]:  k_s = E_s/D (D <= 1 m) bzw. E_s/1 m")
    for s in e.schichten:
        _z(bl, "   {:<18s} z = {} .. {} m : E_s = {} kN/m2 -> k_s = {} kN/m3"
           .format(s.name[:18], _f(s.z_o, 2), _f(s.z_u, 2), _f(s.E_s, 0),
                   _f(bettungsmodul(s.E_s, D_m), 0)))
    _z(bl, "")
    _z(bl, "   Biegesteifigkeit EI = {} MNm2 (Zustand II angenaehert mit 0,5 EI_b)"
       .format(_f(EI_riss / 1000.0, 1)))
    _z(bl, "   elastische Laenge 1/lambda = {} m ; L/(1/lambda) = {}"
       .format(_f(hor.elastische_laenge, 2),
               _f(e.L / max(hor.elastische_laenge, 1e-9), 2)))
    _z(bl, "")
    _z(bl, "   Kopfverschiebung w = {} mm".format(_f(hor.w_kopf, 2)))
    _z(bl, "   M_max = {} kNm bei z = {} m".format(_f(hor.M_max, 1), _f(hor.z_Mmax, 2)))
    _z(bl, "   V_max = {} kN".format(_f(hor.V_max, 1)))
    for hh in hor.hinweise:
        _z(bl, "   !! " + hh)
    B.append(bl)

    # =====================================================================
    # 6. LAENGSBEWEHRUNG (M-N-INTERAKTION)
    # =====================================================================
    mb_1536_probe = mindestbewehrung_pfahl(
        Kreisquerschnitt(D=e.D, c_nom=c_nom, phi_l=e.phi_l, n_l=6, phi_w=e.phi_w))
    mb_stuetze = mindestbewehrung_druckglied(e.N_Ed, S)
    As_min = max(mb_1536_probe["As_min"], mb_stuetze["As_min"])

    if e.n_l and e.n_l >= 6:
        qs = Kreisquerschnitt(D=e.D, c_nom=c_nom, phi_l=e.phi_l, n_l=e.n_l,
                              phi_w=e.phi_w)
        diag = interaktionsdiagramm(qs, C, S)
        M_Rd = M_Rd_bei_N(diag, e.N_Ed)
        ok_mn = M_Rd >= M_Ed_max - 1e-9
    else:
        qs, diag, M_Rd, ok_mn = erforderliche_bewehrung(
            e.D, c_nom, e.phi_l, e.phi_w, C, S, e.N_Ed, M_Ed_max, As_min=As_min)
    b["querschnitt"], b["interaktion"], b["M_Rd"] = qs, diag, M_Rd

    bl = _block("6. LAENGSBEWEHRUNG - M-N-INTERAKTION",
                [ref("qs_biegung"), ref("qs_punkt_c")])
    _z(bl, "Kreisquerschnitt mit gleichmaessig verteilten Laengsstaeben.")
    _z(bl, "Das Interaktionsdiagramm wird ueber das volle Dehnungsdiagramm nach")
    _z(bl, "EC2 6.1 (2)P, Bild 6.1 (Bemessungspunkte A, B, C) ermittelt.")
    _z(bl, "")
    _z(bl, "   D_s (Bewehrungskreis) = D - 2(c_nom + phi_w) - phi_l = {} mm"
       .format(_f(qs.D_s, 0)))
    _z(bl, "   GEWAEHLT: {} phi {} = {} mm2   (rho = {} %)"
       .format(qs.n_l, _f(qs.phi_l, 0), _f(qs.As_ges, 0), _f(100 * qs.rho_l, 2)))
    _z(bl, "   lichter Stababstand = {} mm".format(_f(qs.lichter_stababstand(), 1)))
    _z(bl, "")
    _z(bl, "   Eckwerte des Interaktionsdiagramms:")
    _z(bl, "      N_Rd (zentrischer Zug)   = {} kN".format(_f(diag["N_zug"], 0)))
    _z(bl, "      N_Rd (zentrischer Druck) = {} kN".format(_f(diag["N_druck_max"], 0)))
    _z(bl, "      M_Rd,max                 = {} kNm".format(_f(diag["M_max"], 1)))
    _z(bl, "")
    _z(bl, "   NACHWEIS bei N_Ed = {} kN:".format(_f(e.N_Ed, 0)))
    _z(bl, "      M_Rd = {} kNm  {}  M_Ed = {} kNm   (eta = {})"
       .format(_f(M_Rd, 1), ">=" if ok_mn else "<", _f(M_Ed_max, 1),
               _f(M_Ed_max / max(M_Rd, 1e-9), 3)))
    _nw(bl, "M-N-Interaktion M_Ed/M_Rd", M_Ed_max, M_Rd, ok_mn, "kNm")
    B.append(bl)

    # =====================================================================
    # 7. MINDESTBEWEHRUNG
    # =====================================================================
    mb = mindestbewehrung_pfahl(qs)
    kon = konstruktive_pruefung(qs, e.d_g)
    bl = _block("7. MINDESTBEWEHRUNG UND KONSTRUKTIVE REGELN",
                [ref("pfahl_As_min"), ref("pfahl_laengs"), ref("qs_As_stuetze")])
    _z(bl, "Mindestlaengsbewehrung  [DIN EN 1536, 7.6.3, Tab. 4]:")
    _z(bl, "   Ac = {} m2  ->  Regel: {}".format(_f(mb["Ac_m2"], 3), mb["regel"]))
    _z(bl, "   As,min (EN 1536) = {} mm2".format(_f(mb["As_min"], 0)))
    _z(bl, "Mindestbewehrung von Druckgliedern  [NA zu EC2, NDP zu 9.5.2 (2)]:")
    _z(bl, "   As,min = 0,15 |N_Ed|/fyd = 0,15 * {} * 10^3 / {} = {} mm2"
       .format(_f(abs(e.N_Ed), 0), _f(S.fyd, 1), _f(mb_stuetze["As_min"], 0)))
    _z(bl, "   massgebend: As,min = {} mm2".format(_f(As_min, 0)))
    _z(bl, "   As,vorh = {} mm2".format(_f(qs.As_ges, 0)))
    _nw(bl, "As,vorh >= As,min", qs.As_ges, As_min, qs.As_ges >= As_min - 1e-6,
        "mm2", ">=")
    _z(bl, "")
    _z(bl, "Konstruktive Anforderungen  [DIN EN 1536, 7.6.3 / 7.6.4]:")
    _z(bl, "   Stabanzahl n = {} {} 6                     -> {}"
       .format(kon["n_l"], ">=" if kon["n_ok"] else "<", "OK" if kon["n_ok"] else "NEIN"))
    _z(bl, "   Stabdurchmesser phi = {} mm {} 16 mm       -> {}"
       .format(_f(kon["phi_l"], 0), ">=" if kon["phi_ok"] else "<",
               "OK" if kon["phi_ok"] else "NEIN"))
    _z(bl, "   lichter Stababstand = {} mm {} {} mm      -> {}"
       .format(_f(kon["s_licht"], 1), ">=" if kon["abstand_ok"] else "<",
               _f(kon["s_min"], 0), "OK" if kon["abstand_ok"] else "NEIN"))
    _z(bl, "   Wendeldurchmesser phi_w = {} mm {} {} mm  -> {}"
       .format(_f(qs.phi_w, 0), ">=" if kon["wendel_ok"] else "<",
               _f(kon["phi_w_min"], 1), "OK" if kon["wendel_ok"] else "NEIN"))
    _nw(bl, "Stabanzahl n >= 6", kon["n_l"], 6, kon["n_ok"], "-", ">=")
    _nw(bl, "Stabdurchmesser >= 16 mm", kon["phi_l"], 16.0, kon["phi_ok"], "mm", ">=")
    _nw(bl, "lichter Stababstand", kon["s_licht"], kon["s_min"],
        kon["abstand_ok"], "mm", ">=")
    b["mindestbewehrung"] = dict(en1536=mb, stuetze=mb_stuetze, As_min=As_min,
                                 konstruktiv=kon)
    B.append(bl)

    # =====================================================================
    # 8. QUERBEWEHRUNG (WENDEL)
    # =====================================================================
    qk = querkraft_kreis(qs, C, S, V_Ed_max, e.N_Ed)
    b["querkraft"] = qk
    bl = _block("8. QUERBEWEHRUNG - WENDEL", qk["normen"])
    _z(bl, "EC2 enthaelt fuer Kreisquerschnitte keine expliziten Querkraftregeln.")
    _z(bl, "Angesetzt werden (uebliche Praxis):")
    _z(bl, "   b_w = D = {} mm     d = D/2 + D_s/pi = {} mm     z = 0,9 d = {} mm"
       .format(_f(qk["bw"], 0), _f(qk["d"], 1), _f(qk["z"], 1)))
    _z(bl, "   A_sl = As/2 = {} mm2 (gezogene Stabhaelfte) -> rho_l = {}"
       .format(_f(qk["A_sl"], 0), _f(qk["rho_l"], 5)))
    _z(bl, "")
    _z(bl, "Tragfaehigkeit ohne Querkraftbewehrung  [EC2 6.2.2 + NA]:")
    _z(bl, "   k = {}   sigma_cp = {} N/mm2 (<= 0,2 fcd)   v_min = {} N/mm2"
       .format(_f(qk["k"], 3), _f(qk["sigma_cp"], 2), _f(qk["v_min"], 3)))
    _z(bl, "   V_Rd,c = {} kN   |   V_Ed = {} kN  ->  {}"
       .format(_f(qk["V_Rdc"], 1), _f(V_Ed_max, 1),
               "Querkraftbewehrung rechnerisch erforderlich"
               if qk["bewehrung_erforderlich"] else "Mindestbewehrung genuegt"))
    _z(bl, "")
    _z(bl, "Fachwerkmodell  [NA Gl. (6.7aDE) ; EC2 Gl. (6.8)/(6.9)]:")
    _z(bl, "   V_Rd,cc = {} kN  ->  cot(theta) = {} (theta = {} Grad)"
       .format(_f(qk["V_Rdcc"], 1), _f(qk["cot"], 3), _f(qk["theta"], 1)))
    _z(bl, "   V_Rd,max = {} kN   (nu1 = {})".format(_f(qk["V_Rdmax"], 1),
                                                     _f(qk["nu1"], 3)))
    _nw(bl, "Druckstrebe V_Ed/V_Rd,max", V_Ed_max, qk["V_Rdmax"],
        qk["ok_druckstrebe"], "kN")
    _z(bl, "")
    _z(bl, "   asw,erf = V_Ed/(z fywd cot) = {} cm2/m".format(_f(qk["asw_erf"] / 100, 2)))
    _z(bl, "   asw,min = 0,16 fctm/fyk * bw = {} cm2/m".format(_f(qk["asw_min"] / 100, 2)))
    _z(bl, "   massgebend asw = {} cm2/m".format(_f(qk["asw"] / 100, 2)))
    _z(bl, "   Wendel phi {} (2 Schenkel, A = {} mm2) -> s_erf = {} mm"
       .format(_f(qs.phi_w, 0), _f(2 * qk["A_wendel"], 1), _f(qk["s_erf"], 0)))
    _z(bl, "   Grenzen: s <= min(20 phi_l ; D ; 400) = {} mm  [EC2 9.5.3]"
       .format(_f(qk["s_max_ec2"], 0)))
    _z(bl, "            100 mm <= s <= 400 mm                 [DIN EN 1536, 7.6.4]")
    _z(bl, "   GEWAEHLT: Wendel phi {} , Ganghoehe s = {} mm"
       .format(_f(qs.phi_w, 0), _f(qk["s_gewaehlt"], 0)))
    _nw(bl, "Wendelabstand s <= s_max", qk["s_gewaehlt"], qk["s_max"],
        qk["s_gewaehlt"] <= qk["s_max"] + 1e-6, "mm")
    B.append(bl)

    # =====================================================================
    # 9. KNICKNACHWEIS
    # =====================================================================
    Es_min = min([s.E_s for s in e.schichten]) if e.schichten else 0.0
    cu_min = min([s.c_u_k for s in e.schichten if s.c_u_k > 0], default=0.0)
    ks_min = bettungsmodul(Es_min, D_m) if Es_min > 0 else 0.0
    kn = knicklast_gebettet(EI_riss, ks_min, D_m)
    bl = _block("9. KNICKNACHWEIS", [ref("geo_knicken"), ref("qs_theta2")])
    _z(bl, "EA-Pfaehle, 4.7: ein Knicknachweis ist nur bei sehr weichen Boeden")
    _z(bl, "erforderlich (Richtwert c_u < 10 kN/m2).")
    if cu_min > 0:
        _z(bl, "   kleinste undraenierte Kohaesion c_u = {} kN/m2".format(_f(cu_min, 1)))
    _z(bl, "   kleinster Steifemodul E_s = {} kN/m2 -> k_s = {} kN/m3"
       .format(_f(Es_min, 0), _f(ks_min, 0)))
    _z(bl, "   ideale Knicklast (Engesser) N_ki = 2 sqrt(EI k) = {} kN"
       .format(_f(kn["N_ki"], 0)))
    _z(bl, "   |N_Ed| = {} kN  ->  N_ki/|N_Ed| = {}"
       .format(_f(abs(e.N_Ed), 0), _f(kn["N_ki"] / max(abs(e.N_Ed), 1e-9), 1)))
    if cu_min > 0 and cu_min < 10.0:
        _z(bl, "   !! c_u < 10 kN/m2: Knicknachweis nach EC2 5.8 zwingend fuehren.")
        _nw(bl, "Knicksicherheit N_ki/|N_Ed| >= 3", kn["N_ki"] / max(abs(e.N_Ed), 1e-9),
            3.0, kn["N_ki"] / max(abs(e.N_Ed), 1e-9) >= 3.0, "-", ">=")
    else:
        _z(bl, "   -> kein Knicknachweis erforderlich (Boden ausreichend fest).")
    B.append(bl)

    # =====================================================================
    # 10. ZUSAMMENFASSUNG
    # =====================================================================
    bl = _block("10. ZUSAMMENSTELLUNG DER NACHWEISE", [])
    alle_nw = []
    for bb in B:
        alle_nw.extend(bb["nachweise"])
    for c in alle_nw:
        eta = (c["wert"] / c["grenzwert"]) if (c["vergleich"] == "<=" and c["grenzwert"]) \
            else ((c["grenzwert"] / c["wert"]) if c["wert"] else 0.0)
        _z(bl, "   [{}]  {:<40s} {:>11s} {} {:>11s} {}   eta = {}"
           .format("OK" if c["ok"] else "NEIN", c["name"], _f(c["wert"], 2),
                   c["vergleich"], _f(c["grenzwert"], 2), c["einheit"], _f(eta, 3)))
        if not c["ok"]:
            b["ok_gesamt"] = False
    _z(bl, "")
    _z(bl, "   GESAMTERGEBNIS: " + ("NACHWEISE ERFUELLT" if b["ok_gesamt"]
                                    else "NACHWEISE NICHT ERFUELLT"))
    _z(bl, "")
    _z(bl, "   GEWAEHLTE BEWEHRUNG")
    _z(bl, "      Laengsbewehrung : {} phi {} = {} mm2 (rho = {} %)"
       .format(qs.n_l, _f(qs.phi_l, 0), _f(qs.As_ges, 0), _f(100 * qs.rho_l, 2)))
    _z(bl, "      Wendel          : phi {} , Ganghoehe {} mm"
       .format(_f(qs.phi_w, 0), _f(qk["s_gewaehlt"], 0)))
    lk = e.laenge_bewehrungskorb if e.laenge_bewehrungskorb > 0 else e.L
    _z(bl, "      Bewehrungskorb  : Laenge {} m von {} m Pfahllaenge"
       .format(_f(lk, 2), _f(e.L, 2)))
    _z(bl, "      Betondeckung    : c_nom = {} mm".format(_f(c_nom, 0)))
    B.append(bl)

    b["nachweise"] = alle_nw
    b["c_nom"] = c_nom
    b["M_Ed_max"], b["V_Ed_max"] = M_Ed_max, V_Ed_max
    return b


# ---------------------------------------------------------------------------
# Textbericht
# ---------------------------------------------------------------------------
def bericht_text(b, mit_normen=True):
    """Erzeugt den vollstaendigen Pfahlbericht als Text."""
    e = b["eingabe"]
    out = []
    out.append("=" * 100)
    out.append("BEMESSUNG EINES BOHRPFAHLS")
    out.append("Normen: DIN EN 1536:2015-10 | DIN EN 1997-1 + DIN 1054:2010-12 | "
               "EA-Pfaehle | DIN EN 1992-1-1 + NA")
    out.append("=" * 100)
    out.append("Bohrpfahl D = {:.0f} mm , L = {:.2f} m | {} | {}"
               .format(e.D, e.L, e.betonklasse, e.stahlsorte))
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
    out.append("HINWEIS: automatisch erzeugte Ergebnisse. q_b,k und q_s,k sind")
    out.append("Eingabewerte des Anwenders (EA-Pfaehle, Tab. 5.12 bis 5.15 oder")
    out.append("Probebelastung). Nicht erfasst: Pfahlgruppenwirkung, negative")
    out.append("Mantelreibung, zyklische und dynamische Einwirkungen, Erdbeben,")
    out.append("Setzungsdifferenzen des Gesamtbauwerks. Alle Ergebnisse sind von")
    out.append("einem verantwortlichen Ingenieur zu pruefen.")
    out.append("=" * 100)
    return "\n".join(out)
