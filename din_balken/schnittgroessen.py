# -*- coding: utf-8 -*-
"""
Schnittgroessenermittlung fuer Durchlauftraeger.

Weggroessenverfahren (Euler-Bernoulli-Balken, 2 Freiheitsgrade je Knoten)
+ Integration der Gleichgewichtsbedingungen fuer die Verlaeufe V(x) und M(x).

Einwirkungskombinationen:
    GZT   DIN EN 1990, 6.4.3.2, Gl. (6.10)   ->  1,35 G + 1,50 Q
          Beiwerte nach DIN EN 1990/NA, Tab. NA.A.1.2(B)
    GZG   DIN EN 1990, 6.5.3, Gl. (6.16b)    ->  G + psi_2 Q  (quasi-staendig)
    Feldweise Laststellung nach DIN EN 1992-1-1, 5.1.3.

Vorzeichenregelung
    - Lasten nach UNTEN positiv.
    - Auflagerkraefte nach OBEN positiv.
    - Biegemoment POSITIV = Zug am UNTEREN Rand (Feldmoment).
    - Querkraft V(x) = Summe der nach oben gerichteten Kraefte links von x.
"""

import itertools
from dataclasses import dataclass, field

import numpy as np

from .normen import ref

GELENKIG = "gelenkig"
EINGESPANNT = "eingespannt"


@dataclass
class Auflager:
    x: float                    # Lage [m]
    typ: str = GELENKIG         # "gelenkig" | "eingespannt"
    breite: float = 0.30        # Auflagerbreite [m] (fuer V im Abstand d)


@dataclass
class Streckenlast:
    x1: float
    x2: float
    q1: float                   # [kN/m] nach unten
    q2: float = None            # [kN/m] (None -> konstant)
    art: str = "G"              # "G" staendig | "Q" veraenderlich

    def __post_init__(self):
        if self.q2 is None:
            self.q2 = self.q1


@dataclass
class Einzellast:
    x: float
    P: float                    # [kN] nach unten
    art: str = "G"


@dataclass
class Ergebnis:
    x: np.ndarray               # [m]
    V: np.ndarray               # [kN]
    M: np.ndarray               # [kNm]
    auflagerkraefte: list       # [(x, R [kN], M_Auflager [kNm])]
    w: np.ndarray = None        # elastische Durchbiegung [mm] (informativ)


class Durchlauftraeger:
    """
    Traeger mit einem oder mehreren Feldern und punktfoermigen Auflagern.

    Parameter
    ---------
    L : float             Gesamtlaenge [m]
    auflager : list[Auflager]
    EI : float            Biegesteifigkeit [kNm2] (nur bei statisch unbestimmten
                          Systemen von Einfluss)
    """

    def __init__(self, L, auflager, EI=1.0e5, npkt=1201):
        self.L = float(L)
        self.auflager = sorted(auflager, key=lambda a: a.x)
        self.EI = float(EI)
        self.npkt = int(npkt)
        self.streckenlasten = []
        self.einzellasten = []
        if len(self.auflager) < 1:
            raise ValueError("Mindestens ein Auflager erforderlich")

    # -- Lastdefinition ---------------------------------------------------
    def strecke(self, x1, x2, q1, q2=None, art="G"):
        self.streckenlasten.append(Streckenlast(x1, x2, q1, q2, art))
        return self

    def einzel(self, x, P, art="G"):
        self.einzellasten.append(Einzellast(x, P, art))
        return self

    # -- Feldgeometrie ----------------------------------------------------
    @property
    def felder(self):
        """Liste von (x_anf, x_end) zwischen benachbarten Auflagern."""
        xs = [a.x for a in self.auflager]
        return [(xs[i], xs[i + 1]) for i in range(len(xs) - 1)]

    # -- Berechnung -------------------------------------------------------
    def _knoten(self):
        p = {0.0, self.L}
        p.update(a.x for a in self.auflager)
        p.update(c.x for c in self.einzellasten)
        for c in self.streckenlasten:
            p.update((c.x1, c.x2))
        p = sorted(v for v in p if -1e-9 <= v <= self.L + 1e-9)
        verfeinert = set(p)
        for i in range(len(p) - 1):
            a, b = p[i], p[i + 1]
            for k in range(1, 4):
                verfeinert.add(a + (b - a) * k / 4.0)
        return sorted(verfeinert)

    def _q_bei(self, x, faktor=None):
        """
        Gewichtete Streckenlast [kN/m] an der Stelle x.
        `faktor(last, x) -> float` beruecksichtigt Teilsicherheitsbeiwert und
        Laststellung (die Wichtung haengt von der LAGE ab, weil eine Last
        mehrere Felder ueberspannen kann).
        """
        q = 0.0
        for c in self.streckenlasten:
            if c.x1 - 1e-12 <= x <= c.x2 + 1e-12 and c.x2 > c.x1:
                t = (x - c.x1) / (c.x2 - c.x1)
                qi = c.q1 + (c.q2 - c.q1) * min(max(t, 0.0), 1.0)
                q += qi * (faktor(c, x) if faktor is not None else 1.0)
        return q

    def berechnen(self, faktor_G=1.0, faktor_Q=1.0, laststellung=None):
        """
        Berechnet den Traeger fuer  faktor_G * (G-Lasten) + faktor_Q * (Q-Lasten).

        `laststellung(x) -> bool` gibt an, ob die veraenderliche Last an der
        Stelle x wirkt (feldweise Laststellung, DIN EN 1992-1-1, 5.1.3).
        """
        def fak(c, xpos):
            if c.art.upper().startswith("Q"):
                if laststellung is not None and not laststellung(xpos):
                    return 0.0
                return faktor_Q
            return faktor_G

        knoten = np.array(self._knoten())
        nk = len(knoten)
        ndof = 2 * nk
        K = np.zeros((ndof, ndof))
        F = np.zeros(ndof)

        # --- Assemblierung (w positiv nach OBEN, Kraefte positiv nach oben)
        for e in range(nk - 1):
            le = knoten[e + 1] - knoten[e]
            if le <= 0:
                continue
            k = self.EI / le ** 3 * np.array([
                [12.0,     6 * le,      -12.0,   6 * le],
                [6 * le,   4 * le ** 2, -6 * le, 2 * le ** 2],
                [-12.0,    -6 * le,     12.0,    -6 * le],
                [6 * le,   2 * le ** 2, -6 * le, 4 * le ** 2],
            ])
            dofs = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
            K[np.ix_(dofs, dofs)] += k

            # Jedes Element liegt vollstaendig in einem Feld (an jedem Auflager
            # sitzt ein Knoten). Die Laststellung wird daher EINMAL in der
            # Elementmitte ausgewertet - so erbt ein Element, das genau an
            # einem Auflager beginnt, nicht die Stellung des Vorfeldes.
            xm = 0.5 * (knoten[e] + knoten[e + 1])

            def fak_e(c, _x, xm=xm):
                return fak(c, xm)

            w1 = -self._q_bei(knoten[e] + 1e-9, fak_e)
            w2 = -self._q_bei(knoten[e + 1] - 1e-9, fak_e)
            fe = np.array([
                le * (7 * w1 + 3 * w2) / 20.0,
                le ** 2 * (3 * w1 + 2 * w2) / 60.0,
                le * (3 * w1 + 7 * w2) / 20.0,
                -le ** 2 * (2 * w1 + 3 * w2) / 60.0,
            ])
            F[dofs] += fe

        for c in self.einzellasten:
            i = int(np.argmin(np.abs(knoten - c.x)))
            F[2 * i] += -fak(c, c.x) * c.P

        # --- Randbedingungen
        fest = []
        for a in self.auflager:
            i = int(np.argmin(np.abs(knoten - a.x)))
            fest.append(2 * i)
            if a.typ == EINGESPANNT:
                fest.append(2 * i + 1)
        fest = sorted(set(fest))
        frei = [i for i in range(ndof) if i not in fest]

        u = np.zeros(ndof)
        if frei:
            u[frei] = np.linalg.solve(K[np.ix_(frei, frei)], F[frei])
        Rvec = K @ u - F

        auflagerkraefte = []
        for a in self.auflager:
            i = int(np.argmin(np.abs(knoten - a.x)))
            R = Rvec[2 * i]
            Mr = Rvec[2 * i + 1] if a.typ == EINGESPANNT else 0.0
            auflagerkraefte.append((float(knoten[i]), float(R), float(Mr)))

        # --- Schnittgroessen durch Integration des Gleichgewichts
        eps = 1e-7
        krit = [a.x for a in self.auflager] + [c.x for c in self.einzellasten]
        for c in self.streckenlasten:
            krit += [c.x1, c.x2]
        pkt = set(np.linspace(0.0, self.L, self.npkt).tolist())
        for xc in krit:
            for s in (-eps, 0.0, eps):
                pkt.add(min(max(xc + s, 0.0), self.L))
        x = np.array(sorted(pkt))

        q = np.array([self._q_bei(xi, fak) for xi in x])
        Iq = np.concatenate(([0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(x))))

        V = -Iq
        for xr, R, _ in auflagerkraefte:
            V += np.where(x >= xr - 1e-9, R, 0.0)
        for c in self.einzellasten:
            V -= np.where(x >= c.x - 1e-9, fak(c, c.x) * c.P, 0.0)

        # Der letzte Punkt liegt knapp HINTER dem rechten Endauflager (V wuerde
        # dessen Auflagerkraft enthalten und 0 werden); im Traeger gilt der
        # linksseitige Grenzwert.
        if len(x) > 1 and any(abs(xr - self.L) < 1e-7
                              for xr, _, _ in auflagerkraefte):
            V[-1] = V[-2]

        M = np.concatenate(([0.0], np.cumsum(0.5 * (V[1:] + V[:-1]) * np.diff(x))))
        for xr, _, Mr in auflagerkraefte:
            if abs(Mr) > 1e-12:
                M += np.where(x >= xr - 1e-9, -Mr, 0.0)

        wk = -u[0::2] * 1000.0        # [mm], positiv nach unten
        return Ergebnis(x=x, V=V, M=M, auflagerkraefte=auflagerkraefte,
                        w=np.interp(x, knoten, wk))

    # -- Einhuellende ------------------------------------------------------
    def einhuellende(self, gamma_G=1.35, gamma_Q=1.50, laststellungen=True,
                     gamma_G_inf=1.00):
        """
        Einhuellende von M und V im GZT.

        DIN EN 1990, Gl. (6.10) mit gamma_G = 1,35 / gamma_Q = 1,50
        (DIN EN 1990/NA, Tab. NA.A.1.2(B)).
        Mit `laststellungen=True` werden die 2^n feldweisen Laststellungen der
        veraenderlichen Last durchlaufen (DIN EN 1992-1-1, 5.1.3). Zusaetzlich
        wird jede Stellung mit gamma_G,inf = 1,00 auf der GESAMTEN staendigen
        Last untersucht (die staendige Einwirkung ist nach DIN EN 1990, 6.4.3.2
        eine einzige Einwirkung und wird nicht feldweise variiert).
        """
        faelle = []
        felder = self.felder if self.felder else [(0.0, self.L)]
        n = len(felder)
        if not laststellungen or n > 10:
            faelle.append(self.berechnen(gamma_G, gamma_Q, None))
        else:
            for stellung in itertools.product([False, True], repeat=n):
                def stellung_fn(xpos, stellung=stellung):
                    for i, (a, b) in enumerate(felder):
                        if a - 1e-9 <= xpos <= b + 1e-9:
                            return stellung[i]
                    return True   # Kragarme: immer belastet (unguenstig)
                faelle.append(self.berechnen(gamma_G, gamma_Q, stellung_fn))
                faelle.append(self.berechnen(gamma_G_inf, gamma_Q, stellung_fn))

        x = faelle[0].x
        return dict(
            x=x,
            Mmax=np.max([c.M for c in faelle], axis=0),
            Mmin=np.min([c.M for c in faelle], axis=0),
            Vmax=np.max([c.V for c in faelle], axis=0),
            Vmin=np.min([c.V for c in faelle], axis=0),
            R=np.max([[abs(r[1]) for r in c.auflagerkraefte] for c in faelle],
                     axis=0),
            faelle=faelle,
            normen=[ref("komb_GZT"), ref("gamma_F"), ref("lastfaelle")])

    def quasi_staendig(self, psi2=0.3):
        """Quasi-staendige Kombination G + psi_2 Q (DIN EN 1990, Gl. 6.16b)."""
        return self.berechnen(1.0, psi2, None)

    def charakteristisch(self):
        """Charakteristische Kombination G + Q (DIN EN 1990, Gl. 6.14b)."""
        return self.berechnen(1.0, 1.0, None)


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------
def wert_bei(x_arr, y_arr, x0):
    """Interpoliert y(x0)."""
    return float(np.interp(x0, x_arr, y_arr))


def bemessungsquerkraft(einh, auflager, d_m):
    """
    Bemessungsquerkraft an jedem Auflager im Abstand d vom AUFLAGERRAND
    (DIN EN 1992-1-1, 6.2.1 (8)); gilt bei Lasteintrag am oberen Rand und
    direkter Lagerung.
    """
    x, Vmax, Vmin = einh["x"], einh["Vmax"], einh["Vmin"]
    out = []
    for a in auflager:
        for seite in (-1, +1):
            x_rand = a.x + seite * a.breite / 2.0
            x_d = x_rand + seite * d_m
            if x_rand < -1e-9 or x_rand > x[-1] + 1e-9:
                continue
            x_rand = min(max(x_rand, 0.0), x[-1])
            x_d = min(max(x_d, 0.0), x[-1])
            Vr = max(abs(wert_bei(x, Vmax, x_rand)), abs(wert_bei(x, Vmin, x_rand)))
            Vd = max(abs(wert_bei(x, Vmax, x_d)), abs(wert_bei(x, Vmin, x_d)))
            out.append(dict(auflager_x=a.x, seite="links" if seite < 0 else "rechts",
                            x_rand=x_rand, V_rand=Vr, x_d=x_d, V_d=Vd,
                            norm=ref("V_rand")))
    return out
