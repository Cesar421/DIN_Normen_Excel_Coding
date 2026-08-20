# -*- coding: utf-8 -*-
"""
KREISQUERSCHNITT eines Bohrpfahls - Querschnittsbemessung.

Norm: DIN EN 1992-1-1:2011-01, Abschnitt 6.1 (Biegung mit Laengskraft)
      + DIN EN 1992-1-1/NA:2013-04
      DIN EN 1536:2015-10, 7.6 (Bewehrung von Bohrpfaehlen)

Modell
------
Der Kreisquerschnitt wird mit dem vollen Dehnungsdiagramm nach EC2 6.1 (2)P,
Bild 6.1 abgebildet (Bemessungspunkte A, B und C):

    Punkt A : x <= x_A     -> Stahl massgebend, eps_s = eps_ud (Zugrand)
    Punkt B : x_A < x <= D -> Beton massgebend, eps_c = eps_cu2 (Druckrand)
    Punkt C : x > D        -> ueberdrueckt; im Abstand
                              y_C = (1 - eps_c2/eps_cu2) D vom Druckrand
                              gilt eps = eps_c2   [EC2 6.1 (5)]

Bezugssystem: y wird vom staerkst gedrueckten Rand gemessen (0 .. D).
Die Breite des Kreises betraegt   b(y) = 2 sqrt(y (D - y)).

Vorzeichen: DRUCK NEGATIV fuer N_Ed / N_Rd (Konvention des EC2).
Intern wird mit Druck POSITIV gerechnet und am Ende umgesetzt.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .normen_pfahl import ref


@dataclass
class Kreisquerschnitt:
    """
    Parameter
    ---------
    D : float          Pfahldurchmesser [mm]
    c_nom : float      Nennmass der Betondeckung [mm] (am Wendel gemessen)
    phi_l : float      Durchmesser der Laengsstaebe [mm]
    n_l : int          Anzahl der Laengsstaebe
    phi_w : float      Durchmesser der Wendel/Ringe [mm]
    """
    D: float = 900.0
    c_nom: float = 60.0
    phi_l: float = 20.0
    n_l: int = 10
    phi_w: float = 10.0

    @property
    def R(self):
        return self.D / 2.0

    @property
    def Ac(self):
        """Bruttobetonflaeche [mm2]."""
        return math.pi * self.D ** 2 / 4.0

    @property
    def D_s(self):
        """Durchmesser des Bewehrungskreises (Stabachsen) [mm]."""
        return self.D - 2.0 * (self.c_nom + self.phi_w) - self.phi_l

    @property
    def A_stab(self):
        """Flaeche eines Laengsstabes [mm2]."""
        return math.pi * self.phi_l ** 2 / 4.0

    @property
    def As_ges(self):
        """Gesamte Laengsbewehrung [mm2]."""
        return self.n_l * self.A_stab

    @property
    def rho_l(self):
        """Laengsbewehrungsgrad As/Ac [-]."""
        return self.As_ges / self.Ac

    @property
    def d_eff(self):
        """
        Statische Nutzhoehe des Kreisquerschnitts fuer den Querkraftnachweis:
            d = D/2 + D_s/pi
        (Schwerpunkt der gezogenen Stabhaelfte; uebliche Naeherung, da EC2
        fuer Kreisquerschnitte keine explizite Regel enthaelt.)
        """
        return self.D / 2.0 + self.D_s / math.pi

    @property
    def I_brutto(self):
        """Bruttotraegheitsmoment [mm4]."""
        return math.pi * self.D ** 4 / 64.0

    def breite(self, y):
        """Breite des Kreises in der Tiefe y ab dem gedrueckten Rand [mm]."""
        if y <= 0.0 or y >= self.D:
            return 0.0
        return 2.0 * math.sqrt(y * (self.D - y))

    def stab_lagen(self, offset_grad=0.0):
        """
        Tiefen y_i der Laengsstaebe ab dem gedrueckten Rand [mm].
        Die Staebe liegen gleichmaessig auf dem Bewehrungskreis; der erste Stab
        liegt bei `offset_grad` (0 Grad = staerkst gedrueckter Rand).
        """
        rs = self.D_s / 2.0
        winkel = np.radians(offset_grad) + 2.0 * np.pi * np.arange(self.n_l) / self.n_l
        return self.R - rs * np.cos(winkel)

    def lichter_stababstand(self):
        """Lichter Abstand der Laengsstaebe auf dem Bewehrungskreis [mm]."""
        if self.n_l < 2:
            return float("inf")
        bogen = math.pi * self.D_s / self.n_l
        return bogen - self.phi_l

    def beschreibung(self):
        return ("Bohrpfahl D = {:.0f} mm, {} phi {:.0f} "
                "(As = {:.0f} mm2, rho = {:.2f} %), Wendel phi {:.0f}"
                .format(self.D, self.n_l, self.phi_l, self.As_ges,
                        100 * self.rho_l, self.phi_w))


# ---------------------------------------------------------------------------
# Dehnungszustand nach EC2 6.1 (2)P, Bild 6.1
# ---------------------------------------------------------------------------
def dehnungsebene(qs, x, beton, stahl):
    """
    Liefert (eps_rand, funktion eps(y)) fuer eine Nulllinienlage x [mm].

    x > 0 ab dem gedrueckten Rand gemessen. x -> unendlich entspricht
    zentrischem Druck.
    """
    D = qs.D
    ecu, ec2, eud = beton.eps_cu2, beton.eps_c2, stahl.eps_ud
    y_zug = float(np.max(qs.stab_lagen()))       # tiefster (gezogener) Stab

    if x >= D:
        # Punkt C: eps_c2 im Abstand y_C vom gedrueckten Rand  [EC2 6.1 (5)]
        y_C = (1.0 - ec2 / ecu) * D
        eps_rand = ec2 * x / max(x - y_C, 1e-9)
        eps_rand = min(eps_rand, ecu)
        punkt = "C"
    else:
        x_A = y_zug * ecu / (ecu + eud)          # Grenze Punkt A / Punkt B
        if x <= x_A:
            eps_rand = eud * x / max(y_zug - x, 1e-9)
            eps_rand = min(eps_rand, ecu)
            punkt = "A"
        else:
            eps_rand = ecu
            punkt = "B"

    def eps(y):
        return eps_rand * (1.0 - y / x) if x > 1e-12 else 0.0

    return eps_rand, eps, punkt


def schnittgroessen_bei_x(qs, beton, stahl, x, n_int=200, offset_grad=0.0,
                          beton_verdraengung=True):
    """
    Innere Schnittgroessen (N_Rd, M_Rd) fuer eine Nulllinienlage x.

    Rueckgabe (Druck NEGATIV fuer N_Rd, wie in EC2):
        N_Rd [kN], M_Rd [kNm], Zusatzinformationen
    """
    D, R = qs.D, qs.R
    eps_rand, eps_fn, punkt = dehnungsebene(qs, x, beton, stahl)

    # --- Betondruckkraft ueber die gedrueckte Kreisflaeche ---------------
    y_max = min(x, D)
    Fc, Mc = 0.0, 0.0
    yc = 0.0
    if y_max > 1e-9 and eps_rand > 1e-12:
        ys = np.linspace(0.0, y_max, n_int + 1)
        bs = np.array([qs.breite(y) for y in ys])
        sg = np.array([beton.sigma_c(eps_fn(y)) for y in ys])
        f = sg * bs                                    # [N/mm]
        Fc = float(np.trapz(f, ys))
        if Fc > 1e-9:
            yc = float(np.trapz(f * ys, ys) / Fc)
        Mc = Fc * (R - yc)

    # --- Stahlkraefte -----------------------------------------------------
    Fs, Ms = 0.0, 0.0
    y_staebe = qs.stab_lagen(offset_grad)
    eps_staebe, sig_staebe = [], []
    for y_i in y_staebe:
        e_i = eps_fn(y_i)
        e_i = max(min(e_i, beton.eps_cu2), -stahl.eps_ud)
        s_i = stahl.sigma_s(e_i)
        # verdraengter Beton bei gedrueckten Staeben abziehen
        if beton_verdraengung and e_i > 0:
            s_i -= beton.sigma_c(e_i)
        F_i = qs.A_stab * s_i
        Fs += F_i
        Ms += F_i * (R - y_i)
        eps_staebe.append(e_i)
        sig_staebe.append(s_i)

    N_druck = Fc + Fs                                  # Druck positiv
    return dict(N_Rd=-N_druck / 1.0e3, M_Rd=(Mc + Ms) / 1.0e6,
                x=x, eps_rand=eps_rand, punkt=punkt, Fc=Fc / 1.0e3,
                Fs=Fs / 1.0e3, y_c=yc,
                eps_staebe=np.array(eps_staebe), sigma_staebe=np.array(sig_staebe),
                y_staebe=y_staebe,
                eps_zug=float(np.min(eps_staebe)) if eps_staebe else 0.0)


# ---------------------------------------------------------------------------
# M-N-Interaktionsdiagramm
# ---------------------------------------------------------------------------
def interaktionsdiagramm(qs, beton, stahl, n_punkte=160, offset_grad=0.0):
    """
    Vollstaendiges M-N-Interaktionsdiagramm des Kreisquerschnitts.
    DIN EN 1992-1-1, 6.1.

    Rueckgabe: dict mit N [kN] (Druck negativ) und M [kNm], jeweils fuer die
    positive Momentenrichtung; das Diagramm ist zu M = 0 symmetrisch.
    """
    D = qs.D
    # Nulllinienlagen: von sehr klein (reine Biegung) bis sehr gross (Druck)
    xs = np.concatenate([
        np.linspace(0.02 * D, 1.00 * D, n_punkte // 2),
        np.geomspace(1.01 * D, 60.0 * D, n_punkte - n_punkte // 2),
    ])
    N, M, pkt = [], [], []
    for x in xs:
        r = schnittgroessen_bei_x(qs, beton, stahl, x, offset_grad=offset_grad)
        N.append(r["N_Rd"])
        M.append(r["M_Rd"])
        pkt.append(r["punkt"])

    # reiner Zug: alle Staebe auf fyd
    N_zug = qs.As_ges * stahl.fyd / 1.0e3
    N = np.array([N_zug] + N)
    M = np.array([0.0] + M)
    pkt = ["Zug"] + pkt
    return dict(N=N, M=M, punkt=pkt, x=np.concatenate([[0.0], xs]),
                N_zug=N_zug, N_druck_max=float(np.min(N)),
                M_max=float(np.max(M)),
                normen=[ref("qs_biegung"), ref("qs_punkt_c")])


def M_Rd_bei_N(diagramm, N_Ed):
    """
    Momententragfaehigkeit M_Rd [kNm] bei gegebener Laengskraft N_Ed [kN]
    (Druck negativ) durch Interpolation im Interaktionsdiagramm.
    """
    N, M = diagramm["N"], diagramm["M"]
    if N_Ed > diagramm["N_zug"] or N_Ed < diagramm["N_druck_max"]:
        return 0.0
    # N ist nicht monoton -> ueber die Kurve laufen und den Schnitt suchen
    beste = 0.0
    for i in range(len(N) - 1):
        n1, n2 = N[i], N[i + 1]
        if (n1 - N_Ed) * (n2 - N_Ed) <= 0 and abs(n2 - n1) > 1e-12:
            t = (N_Ed - n1) / (n2 - n1)
            beste = max(beste, M[i] + t * (M[i + 1] - M[i]))
        elif abs(n1 - N_Ed) < 1e-9:
            beste = max(beste, M[i])
    return beste


def erforderliche_bewehrung(D, c_nom, phi_l, phi_w, beton, stahl, N_Ed, M_Ed,
                            n_min=6, n_max=40, As_min=0.0):
    """
    Sucht die kleinste Stabanzahl n (bei festem Stabdurchmesser), fuer die
    M_Rd(N_Ed) >= M_Ed gilt und die Mindestbewehrung eingehalten wird.

    Rueckgabe: (Kreisquerschnitt, Interaktionsdiagramm, M_Rd, ok)
    """
    for n in range(max(6, n_min), n_max + 1):
        qs = Kreisquerschnitt(D=D, c_nom=c_nom, phi_l=phi_l, n_l=n, phi_w=phi_w)
        if qs.As_ges < As_min:
            continue
        dg = interaktionsdiagramm(qs, beton, stahl, n_punkte=120)
        M_Rd = M_Rd_bei_N(dg, N_Ed)
        if M_Rd >= abs(M_Ed) - 1e-9:
            return qs, dg, M_Rd, True
    qs = Kreisquerschnitt(D=D, c_nom=c_nom, phi_l=phi_l, n_l=n_max, phi_w=phi_w)
    dg = interaktionsdiagramm(qs, beton, stahl, n_punkte=120)
    return qs, dg, M_Rd_bei_N(dg, N_Ed), False


# ---------------------------------------------------------------------------
# Mindestbewehrung nach DIN EN 1536, 7.6.3, Tab. 4
# ---------------------------------------------------------------------------
def mindestbewehrung_pfahl(qs):
    """
    Mindestlaengsbewehrung von Bohrpfaehlen.
    DIN EN 1536:2015-10, 7.6.3, Tabelle 4:

        Ac <= 0,5 m2          ->  As >= 0,5  % Ac
        0,5 m2 < Ac <= 1,0 m2 ->  As >= 25 cm2
        Ac > 1,0 m2           ->  As >= 0,25 % Ac

    Zusaetzlich: mindestens 6 Staebe, phi >= 16 mm, lichter Abstand >= 100 mm.
    """
    Ac_m2 = qs.Ac / 1.0e6
    if Ac_m2 <= 0.5:
        As_min, regel = 0.005 * qs.Ac, "0,5 % Ac (Ac <= 0,5 m2)"
    elif Ac_m2 <= 1.0:
        As_min, regel = 2500.0, "25 cm2 (0,5 < Ac <= 1,0 m2)"
    else:
        As_min, regel = 0.0025 * qs.Ac, "0,25 % Ac (Ac > 1,0 m2)"
    return dict(As_min=As_min, regel=regel, Ac_m2=Ac_m2,
                n_min=6, phi_min=16.0, s_licht_min=100.0,
                norm=ref("pfahl_As_min"))


def mindestbewehrung_druckglied(N_Ed, stahl):
    """
    Mindestlaengsbewehrung von Druckgliedern.
    DIN EN 1992-1-1/NA, NDP zu 9.5.2 (2):  As,min = 0,15 |N_Ed| / fyd
    """
    return dict(As_min=0.15 * abs(N_Ed) * 1.0e3 / stahl.fyd,
                norm=ref("qs_As_stuetze"))


def betondeckung_pfahl(D, unter_stuetzfluessigkeit=False):
    """
    Betondeckung von Bohrpfaehlen.
    DIN EN 1536:2015-10, 7.6.2:

        c_nom >= 60 mm  bei D >= 0,6 m
        c_nom >= 50 mm  bei D <  0,6 m
        c_nom >= 75 mm  beim Betonieren unter Stuetzfluessigkeit
    """
    c = 60.0 if D >= 600.0 else 50.0
    if unter_stuetzfluessigkeit:
        c = max(c, 75.0)
    return dict(c_nom=c, D=D, unter_stuetzfluessigkeit=unter_stuetzfluessigkeit,
                norm=ref("pfahl_deckung"))


def konstruktive_pruefung(qs, d_g=16.0):
    """
    Konstruktive Nachweise nach DIN EN 1536, 7.6.3 / 7.6.4 und EC2 9.5.3.
    """
    s_min = 100.0 if d_g > 20.0 else 80.0
    s_licht = qs.lichter_stababstand()
    phi_w_min = max(6.0, qs.phi_l / 4.0)
    return dict(
        n_ok=qs.n_l >= 6, n_l=qs.n_l,
        phi_ok=qs.phi_l >= 16.0, phi_l=qs.phi_l,
        abstand_ok=s_licht >= s_min - 1e-6, s_licht=s_licht, s_min=s_min,
        wendel_ok=qs.phi_w >= phi_w_min - 1e-9, phi_w_min=phi_w_min,
        normen=[ref("pfahl_laengs"), ref("pfahl_quer"), ref("qs_wendel")])
