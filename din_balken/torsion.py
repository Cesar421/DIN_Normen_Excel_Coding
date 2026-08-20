# -*- coding: utf-8 -*-
"""
TORSIONSBEMESSUNG.

Norm: DIN EN 1992-1-1:2011-01, Abschnitt 6.3  (Torsion)
      + DIN EN 1992-1-1/NA:2013-04 (NDP zu 6.3.2)
      Konstruktive Durchbildung: EC2, 9.2.3

Modell (EC2, 6.3.1): ERSATZHOHLQUERSCHNITT
--------------------------------------------------------------------------
Der Vollquerschnitt wird durch einen duennwandigen geschlossenen Querschnitt
ersetzt. Der Torsionsschubfluss laeuft in einem geschlossenen Ring um:

    t_ef,i = A / u        (>= 2 d1 ; <= tatsaechliche Wanddicke)
    A_k    = von den Mittellinien der Waende eingeschlossene Flaeche
    u_k    = Umfang von A_k

    tau_t,i * t_ef,i = T_Ed / (2 A_k)                    Gl. (6.26)
    V_Ed,i           = tau_t,i * t_ef,i * z_i            Gl. (6.27)

Bewehrung
--------------------------------------------------------------------------
    Buegel (je AUSSENSCHENKEL):
        Asw/s = T_Ed / (2 A_k fywd cot(theta))
    Laengsbewehrung, gleichmaessig ueber u_k verteilt:
        sum(Asl) fyd / u_k = T_Ed / (2 A_k) cot(theta)   Gl. (6.28)

Druckstrebe
--------------------------------------------------------------------------
    T_Rd,max = 2 nu alpha_cw fcd A_k t_ef,i sin(theta) cos(theta)   Gl. (6.30)
    NA, NDP zu 6.3.2 (4):
        nu = 0,525 nu_2  bei Vollquerschnitten
        nu = 0,75  nu_2  bei Kastenquerschnitten mit Bewehrung an beiden
                         Wandseiten

Interaktion Torsion + Querkraft
--------------------------------------------------------------------------
    T_Ed/T_Rd,max + V_Ed/V_Rd,max <= 1,0                 Gl. (6.29)

Verzicht auf rechnerische Torsionsbewehrung (NA, NDP zu 6.3.2 (5)):
    T_Ed <= V_Ed bw / 4,5                                Gl. (6.31aDE)
    V_Ed [1 + 4,5 T_Ed/(V_Ed bw)] <= V_Rd,c              Gl. (6.31bDE)

HINWEIS zur Anwendung: die Werte nu = 0,525/0,75 und die Gleichungen
(6.31aDE)/(6.31bDE) sind national festgelegte Parameter. Gegen die
gueltige Ausgabe des Nationalen Anhangs pruefen.
"""

import math
from dataclasses import dataclass, field

from .normen import ref


# ---------------------------------------------------------------------------
# Ersatzhohlquerschnitt
# ---------------------------------------------------------------------------
def ersatzhohlquerschnitt(querschnitt, d1=None):
    """
    Geometrie des Ersatzhohlquerschnitts nach DIN EN 1992-1-1, 6.3.1 (3).

        t_ef,i = A/u ,   jedoch >= 2 d1 (doppelte Randabstand der Laengsstaebe)
                         und <= tatsaechliche Wanddicke

    Fuer den Plattenbalken wird auf der sicheren Seite nur der STEG als
    torsionswirksamer Ring angesetzt (die Platte wird nicht mitgerechnet).
    """
    b = querschnitt.bw
    h = querschnitt.h
    d1 = querschnitt.d1 if d1 is None else d1

    A = b * h
    u = 2.0 * (b + h)
    t_min = 2.0 * d1
    t_max = min(b, h) / 2.0
    t_ef = min(max(A / u, t_min), t_max)

    b_k = b - t_ef
    h_k = h - t_ef
    A_k = b_k * h_k
    u_k = 2.0 * (b_k + h_k)
    return dict(A=A, u=u, t_ef=t_ef, t_ef_roh=A / u, t_min=t_min, t_max=t_max,
                b_k=b_k, h_k=h_k, A_k=A_k, u_k=u_k, b=b, h=h,
                normen=[ref("torsion_allg")])


def T_Rd_c(beton, ehq):
    """
    Risstorsionsmoment.  DIN EN 1992-1-1, 6.3.2 (5), Gl. (6.31):

        T_Rd,c = fctd * t_ef,i * 2 * A_k     [kNm]
    """
    return dict(T_Rdc=beton.fctd * ehq["t_ef"] * 2.0 * ehq["A_k"] / 1.0e6,
                norm=ref("TRdc"))


def T_Rd_max(beton, ehq, cot_theta, alpha_cw=1.0, kasten=False):
    """
    Torsionsdruckstrebe.  DIN EN 1992-1-1, 6.3.2 (4), Gl. (6.30):

        T_Rd,max = 2 nu alpha_cw fcd A_k t_ef,i sin(theta) cos(theta)   [kNm]

    NA, NDP zu 6.3.2 (4):  nu = 0,525 nu_2 (Vollquerschnitt)
                           nu = 0,75  nu_2 (Kastenquerschnitt)
    """
    nu2 = 1.0 if beton.fck <= 50.0 else min(1.1 - beton.fck / 500.0, 1.0)
    nu = (0.75 if kasten else 0.525) * nu2
    theta = math.atan(1.0 / cot_theta)
    T = (2.0 * nu * alpha_cw * beton.fcd * ehq["A_k"] * ehq["t_ef"]
         * math.sin(theta) * math.cos(theta)) / 1.0e6
    return dict(T_Rdmax=T, nu=nu, nu2=nu2, kasten=kasten, theta_grad=math.degrees(theta),
                normen=[ref("TRdmax"), ref("nu_torsion")])


def asw_torsion(T_Ed, ehq, fywd, cot_theta):
    """
    Erforderliche Torsionsbuegel JE AUSSENSCHENKEL [mm2/mm]:

        Asw/s = T_Ed / (2 A_k fywd cot(theta))       (aus Gl. (6.26)/(6.28))
    """
    return abs(T_Ed) * 1.0e6 / (2.0 * ehq["A_k"] * fywd * cot_theta)


def asl_torsion(T_Ed, ehq, fyd, cot_theta):
    """
    Erforderliche Torsionslaengsbewehrung (Summe, ueber u_k verteilt) [mm2].
    DIN EN 1992-1-1, 6.3.2 (3), Gl. (6.28):

        sum(Asl) = T_Ed cot(theta) u_k / (2 A_k fyd)
    """
    return abs(T_Ed) * 1.0e6 * cot_theta * ehq["u_k"] / (2.0 * ehq["A_k"] * fyd)


def bewehrung_entbehrlich(T_Ed, V_Ed, bw, V_Rdc):
    """
    Verzicht auf rechnerische Torsionsbewehrung bei Vollquerschnitten.
    DIN EN 1992-1-1/NA, NDP zu 6.3.2 (5), Gl. (6.31aDE)/(6.31bDE):

        T_Ed <= V_Ed bw / 4,5                                (6.31aDE)
        V_Ed [1 + 4,5 T_Ed/(V_Ed bw)] <= V_Rd,c              (6.31bDE)

    T_Ed [kNm], V_Ed [kN], bw [mm], V_Rd,c [kN].
    """
    T = abs(T_Ed) * 1.0e6          # Nmm
    V = abs(V_Ed) * 1.0e3          # N
    grenze_a = V * bw / 4.5        # Nmm
    bed_a = T <= grenze_a + 1e-9
    if V > 1e-9:
        V_wirk = V * (1.0 + 4.5 * T / (V * bw))
    else:
        V_wirk = float("inf") if T > 0 else 0.0
    bed_b = V_wirk <= V_Rdc * 1.0e3 + 1e-6
    return dict(erfuellt=bool(bed_a and bed_b), bedingung_a=bool(bed_a),
                bedingung_b=bool(bed_b), grenze_a=grenze_a / 1.0e6,
                V_wirksam=V_wirk / 1.0e3, V_Rdc=V_Rdc, norm=ref("torsion_ohne"))


def s_max_torsion(ehq, querschnitt):
    """
    Groesster Buegelabstand fuer Torsion.  DIN EN 1992-1-1, 9.2.3 (3):

        s <= u_k/8  und  <= kleinste Querschnittsabmessung
        (zusaetzlich gilt s_max aus dem Querkraftnachweis)
    """
    s = min(ehq["u_k"] / 8.0, min(querschnitt.bw, querschnitt.h))
    return dict(s_max=s, u_k_8=ehq["u_k"] / 8.0,
                min_abmessung=min(querschnitt.bw, querschnitt.h),
                norm=ref("torsion_konstr"))


# ---------------------------------------------------------------------------
# Vollstaendige Torsionsbemessung
# ---------------------------------------------------------------------------
@dataclass
class ErgebnisTorsion:
    T_Ed: float = 0.0            # [kNm]
    V_Ed: float = 0.0            # [kN]
    A_k: float = 0.0             # [mm2]
    u_k: float = 0.0             # [mm]
    t_ef: float = 0.0            # [mm]
    cot_theta: float = 1.0
    T_Rdc: float = 0.0           # [kNm] Risstorsionsmoment
    T_Rdmax: float = 0.0         # [kNm] Druckstrebe
    V_Rdmax: float = 0.0         # [kN]
    interaktion: float = 0.0     # T_Ed/T_Rd,max + V_Ed/V_Rd,max
    asw_je_schenkel: float = 0.0  # [mm2/m] Torsionsbuegel je Aussenschenkel
    asl_gesamt: float = 0.0      # [mm2] Torsionslaengsbewehrung (Summe)
    asl_je_ecke: float = 0.0     # [mm2] rechnerisch je Eckstab (4 Ecken)
    n_laengsstaebe: int = 4
    s_max: float = 300.0         # [mm]
    erforderlich: bool = True
    ok: bool = True
    hinweise: list = field(default_factory=list)
    detail: dict = field(default_factory=dict)
    normen: list = field(default_factory=list)


def bemessung_torsion(querschnitt, beton, stahl, T_Ed, V_Ed, cot_theta,
                      V_Rdmax, V_Rdc, kasten=False, n_laengsstaebe=None,
                      alpha_cw=1.0):
    """
    Torsionsbemessung und Interaktion mit der Querkraft.

    Parameter
    ---------
    T_Ed : float        Bemessungstorsionsmoment [kNm]
    V_Ed : float        gleichzeitig wirkende Bemessungsquerkraft [kN]
    cot_theta : float   Druckstrebenneigung (dieselbe wie beim Querkraft-
                        nachweis, EC2 6.3.2 (2))
    V_Rdmax, V_Rdc : float  aus dem Querkraftnachweis [kN]
    kasten : bool       True bei Kastenquerschnitt (NA: nu = 0,75 nu_2)
    """
    r = ErgebnisTorsion(T_Ed=abs(T_Ed), V_Ed=abs(V_Ed), cot_theta=cot_theta,
                        V_Rdmax=V_Rdmax)
    r.normen = [ref("torsion_allg"), ref("torsion_schubfluss"),
                ref("torsion_laengs"), ref("torsion_interaktion"),
                ref("TRdmax"), ref("nu_torsion"), ref("TRdc"),
                ref("torsion_ohne"), ref("torsion_konstr")]

    ehq = ersatzhohlquerschnitt(querschnitt)
    r.A_k, r.u_k, r.t_ef = ehq["A_k"], ehq["u_k"], ehq["t_ef"]

    trc = T_Rd_c(beton, ehq)
    r.T_Rdc = trc["T_Rdc"]
    trm = T_Rd_max(beton, ehq, cot_theta, alpha_cw, kasten)
    r.T_Rdmax = trm["T_Rdmax"]
    smt = s_max_torsion(ehq, querschnitt)
    r.s_max = smt["s_max"]

    if r.T_Ed <= 1e-9:
        r.erforderlich = False
        r.hinweise.append("Kein Torsionsmoment angesetzt (T_Ed = 0).")
        r.detail = dict(ehq=ehq, T_Rdc=trc, T_Rdmax=trm, s_max=smt)
        return r

    # --- Interaktion der Druckstreben, Gl. (6.29) ------------------------
    r.interaktion = r.T_Ed / r.T_Rdmax + (r.V_Ed / V_Rdmax if V_Rdmax > 0 else 0.0)
    if r.interaktion > 1.0:
        r.ok = False
        r.hinweise.append(
            "T_Ed/T_Rd,max + V_Ed/V_Rd,max = {:.3f} > 1,0: Druckstrebenversagen "
            "unter Torsion und Querkraft. Querschnitt vergroessern "
            "[EC2 6.3.2 (4), Gl. (6.29)].".format(r.interaktion))

    # --- Verzicht auf rechnerische Torsionsbewehrung (NA) ----------------
    entb = bewehrung_entbehrlich(r.T_Ed, r.V_Ed, querschnitt.bw, V_Rdc)
    if entb["erfuellt"]:
        r.erforderlich = False
        r.hinweise.append(
            "Gl. (6.31aDE) und (6.31bDE) erfuellt: rechnerische Torsions"
            "bewehrung entbehrlich; Mindestbewehrung nach 9.2.2/9.2.3 anordnen "
            "[NA NDP zu 6.3.2 (5)].")

    # --- erforderliche Bewehrung -----------------------------------------
    r.asw_je_schenkel = asw_torsion(r.T_Ed, ehq, stahl.fyd, cot_theta) * 1000.0
    r.asl_gesamt = asl_torsion(r.T_Ed, ehq, stahl.fyd, cot_theta)

    # Laengsstaebe: mindestens einer je Ecke, Abstand entlang u_k <= 350 mm
    n_min = max(4, int(math.ceil(r.u_k / 350.0)))
    r.n_laengsstaebe = n_min if n_laengsstaebe is None else max(4, n_laengsstaebe)
    r.asl_je_ecke = r.asl_gesamt / r.n_laengsstaebe

    if r.T_Ed > r.T_Rdc:
        r.hinweise.append(
            "T_Ed = {:.1f} kNm > T_Rd,c = {:.1f} kNm: der Querschnitt reisst "
            "unter Torsion; Torsionsbewehrung ist erforderlich "
            "[EC2 6.3.2 (5), Gl. (6.31)].".format(r.T_Ed, r.T_Rdc))

    r.detail = dict(ehq=ehq, T_Rdc=trc, T_Rdmax=trm, s_max=smt, entbehrlich=entb)
    return r


def gleichgewichtstorsion_hinweis():
    """
    Hinweistext zur Unterscheidung Gleichgewichts-/Vertraeglichkeitstorsion
    nach DIN EN 1992-1-1, 6.3.1 (2).
    """
    return ("EC2 6.3.1 (2): Ist das Torsionsmoment fuer das Gleichgewicht "
            "NICHT erforderlich (Vertraeglichkeitstorsion, z.B. Einspannung "
            "einer Platte in einen Randbalken), darf im GZT auf den "
            "Torsionsnachweis verzichtet werden; es genuegt eine "
            "Mindestbewehrung nach 9.2.2 und 9.2.3 zur Begrenzung der "
            "Rissbildung. Bei GLEICHGEWICHTSTORSION (statisch bestimmte "
            "Abtragung, z.B. auskragende Konsole am Randbalken) ist der "
            "Nachweis zwingend zu fuehren.")
