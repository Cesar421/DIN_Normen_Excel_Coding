# -*- coding: utf-8 -*-
"""
HORIZONTAL BELASTETER PFAHL - Bettungsmodulverfahren.

Grundlage: EA-Pfaehle (DGGT), Abschnitt 6.3, in Verbindung mit
           DIN EN 1997-1, 7.7 (quer belastete Pfaehle)

Modell (Winkler-Bettung)
------------------------
Der Pfahl wird als elastisch gebetteter Balken abgebildet:

    EI w''''(z) + k_s(z) D w(z) = 0

    k_s = E_s / D     fuer D <= 1,0 m
    k_s = E_s / 1,0 m fuer D >  1,0 m      [EA-Pfaehle, 6.3]

Geloest wird mit Balkenelementen (2 Freiheitsgrade je Knoten) und in den
Knoten konzentrierten Bettungsfedern.

Kopfrandbedingungen
    "frei"        - Pfahlkopf frei drehbar (H und M werden aufgebracht)
    "eingespannt" - Pfahlkopf drehstarr in die Kopfplatte eingebunden

GRENZEN DES MODELLS
    Die Bettung ist LINEAR angesetzt. Die tatsaechlich mobilisierbare
    Bettungsspannung ist durch den Erdwiderstand begrenzt; bei grossen
    Kopfverschiebungen ist ein nichtlineares p-y-Verfahren zu verwenden
    (EA-Pfaehle, 6.4). Die Ausnutzung wird als Hinweis ausgegeben.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .normen_pfahl import ref


@dataclass
class Bodenschicht:
    """
    Bodenschicht laengs des Pfahls.

    z_o, z_u : float   Ober- und Unterkante ab Pfahlkopf [m]
    E_s : float        Steifemodul [kN/m2]  (fuer k_s = E_s/D)
    q_s_k : float      charakteristischer Mantelreibungswert [kN/m2]
    gamma : float      Wichte [kN/m3]        (fuer den Erdwiderstand)
    phi_k : float      Reibungswinkel [Grad] (fuer den Erdwiderstand)
    c_u_k : float      undraenierte Kohaesion [kN/m2] (bindige Boeden)
    name : str
    """
    z_o: float
    z_u: float
    E_s: float = 20000.0
    q_s_k: float = 60.0
    gamma: float = 19.0
    phi_k: float = 30.0
    c_u_k: float = 0.0
    name: str = "Schicht"

    @property
    def dicke(self):
        return self.z_u - self.z_o


def bettungsmodul(E_s, D_m):
    """
    Bettungsmodul k_s [kN/m3] nach EA-Pfaehle, 6.3:
        k_s = E_s / D      fuer D <= 1,0 m
        k_s = E_s / 1,0 m  fuer D >  1,0 m
    """
    return E_s / D_m if D_m <= 1.0 else E_s / 1.0


@dataclass
class ErgebnisBettung:
    z: np.ndarray = None        # [m] ab Pfahlkopf
    w: np.ndarray = None        # [mm] Horizontalverschiebung
    M: np.ndarray = None        # [kNm]
    V: np.ndarray = None        # [kN]
    p: np.ndarray = None        # [kN/m] Bettungsreaktion
    w_kopf: float = 0.0         # [mm]
    M_max: float = 0.0          # [kNm]
    z_Mmax: float = 0.0         # [m]
    V_max: float = 0.0          # [kN]
    k_s_mittel: float = 0.0     # [kN/m3]
    elastische_laenge: float = 0.0   # [m]  1/lambda
    hinweise: list = field(default_factory=list)
    normen: list = field(default_factory=list)


def pfahl_horizontal(L, D, EI, schichten, H=0.0, M_kopf=0.0,
                     kopf="frei", n_elem=200):
    """
    Berechnet den horizontal belasteten Pfahl nach dem Bettungsmodulverfahren.

    Parameter
    ---------
    L : float          Pfahllaenge [m]
    D : float          Pfahldurchmesser [m]
    EI : float         Biegesteifigkeit [kNm2]
    schichten : list[Bodenschicht]
    H : float          Horizontalkraft am Kopf [kN]
    M_kopf : float     Kopfmoment [kNm]
    kopf : str         "frei" | "eingespannt"
    """
    n = int(n_elem)
    le = L / n
    z_kn = np.linspace(0.0, L, n + 1)
    ndof = 2 * (n + 1)
    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)

    def _k_s_bei(z):
        for s in schichten:
            if s.z_o - 1e-9 <= z <= s.z_u + 1e-9:
                return bettungsmodul(s.E_s, D)
        return bettungsmodul(schichten[-1].E_s, D) if schichten else 0.0

    # --- Balkensteifigkeit
    ke = EI / le ** 3 * np.array([
        [12.0,     6 * le,      -12.0,   6 * le],
        [6 * le,   4 * le ** 2, -6 * le, 2 * le ** 2],
        [-12.0,    -6 * le,     12.0,    -6 * le],
        [6 * le,   2 * le ** 2, -6 * le, 4 * le ** 2],
    ])
    for e in range(n):
        d = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        K[np.ix_(d, d)] += ke

    # --- Bettungsfedern (auf die Knoten verteilt)
    k_feder = np.zeros(n + 1)
    for i, z in enumerate(z_kn):
        l_trib = le if 0 < i < n else le / 2.0
        k_feder[i] = _k_s_bei(z) * D * l_trib      # [kN/m]
        K[2 * i, 2 * i] += k_feder[i]

    # --- Belastung am Kopf
    # Das Vorzeichen des Kopfmomentes wird so gewaehlt, dass ein positives
    # M_kopf denselben Biegesinn erzeugt wie eine positive Horizontalkraft H
    # (innerer Biegemoment M(0) = M_kopf).
    F[0] += H
    F[1] += -M_kopf

    fest = []
    if kopf == "eingespannt":
        fest.append(1)                              # Verdrehung = 0
    frei = [i for i in range(ndof) if i not in fest]
    u = np.zeros(ndof)
    u[frei] = np.linalg.solve(K[np.ix_(frei, frei)], F[frei])

    # --- Schnittgroessen aus den Elementendkraeften
    Mz = np.zeros(n + 1)
    Vz = np.zeros(n + 1)
    for e in range(n):
        ue = u[[2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]]
        fe = ke @ ue                                # [V1, M1, V2, M2]
        Mz[e] += -fe[1]
        Mz[e + 1] += fe[3]
        Vz[e] += fe[0]
        Vz[e + 1] += -fe[2]
    # an den Innenknoten liegen zwei Elementbeitraege vor -> mitteln
    Mz[1:n] *= 0.5
    Vz[1:n] *= 0.5

    w_mm = u[0::2] * 1000.0
    p = k_feder * u[0::2] / np.where(np.arange(n + 1) == 0, le / 2.0,
                                     np.where(np.arange(n + 1) == n, le / 2.0, le))

    r = ErgebnisBettung(z=z_kn, w=w_mm, M=Mz, V=Vz, p=p,
                        w_kopf=float(w_mm[0]),
                        M_max=float(np.max(np.abs(Mz))),
                        z_Mmax=float(z_kn[int(np.argmax(np.abs(Mz)))]),
                        V_max=float(np.max(np.abs(Vz))))
    ks_m = float(np.mean([_k_s_bei(z) for z in z_kn]))
    r.k_s_mittel = ks_m
    lam = (ks_m * D / (4.0 * EI)) ** 0.25 if EI > 0 and ks_m > 0 else 0.0
    r.elastische_laenge = 1.0 / lam if lam > 0 else float("inf")
    r.normen = [ref("geo_bettung")]

    if lam > 0 and L * lam < 2.5:
        r.hinweise.append(
            "L*lambda = {:.2f} < 2,5: der Pfahl verhaelt sich KURZ (starr). "
            "Das Verfahren des elastisch gebetteten Balkens ist hier nur "
            "eingeschraenkt gueltig [EA-Pfaehle, 6.3].".format(L * lam))
    p_max = float(np.max(np.abs(p)))
    if p_max > 0:
        r.hinweise.append(
            "Groesste Bettungsreaktion p = {:.1f} kN/m entspricht einer "
            "Sohlspannung von {:.1f} kN/m2. Sie ist gegen den moeglichen "
            "Erdwiderstand zu begrenzen [EA-Pfaehle, 6.4]."
            .format(p_max, p_max / D))
    return r


def knicklast_gebettet(EI, k_s, D):
    """
    Ideale Knicklast eines elastisch gebetteten Pfahls (Engesser):

        N_ki = 2 sqrt(EI k)     mit k = k_s D  [kN/m2]

    Nach EA-Pfaehle, 4.7 ist ein Knicknachweis nur bei sehr weichen Boeden
    erforderlich (Richtwert c_u < 10 kN/m2).
    """
    k = k_s * D
    return dict(N_ki=2.0 * math.sqrt(max(EI * k, 0.0)), k=k,
                norm=ref("geo_knicken"))
