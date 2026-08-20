# -*- coding: utf-8 -*-
"""
Geometrie des Querschnitts.

Rechteckquerschnitt und Plattenbalken nach DIN EN 1992-1-1, 5.3.2.1
(mitwirkende Plattenbreite b_eff).

Bezugssystem: y wird vom STAERKST GEDRUECKTEN Rand nach unten gemessen.
"""

import math
from dataclasses import dataclass

from .normen import ref


@dataclass
class Querschnitt:
    """
    Parameter
    ---------
    b : float       Stegbreite bw [mm]  (Gesamtbreite bei Rechteckquerschnitt)
    h : float       Gesamthoehe [mm]
    d1 : float      Randabstand der Zugbewehrung [mm]
                    (Zugrand -> Schwerpunkt As1)
    d2 : float      Randabstand der Druckbewehrung [mm]
    typ : str       "rechteck" | "plattenbalken"
    b_eff : float   mitwirkende Plattenbreite [mm]  (nur Plattenbalken)
    hf : float      Plattendicke [mm]               (nur Plattenbalken)
    platte_gedrueckt : bool
                    True  -> Platte liegt in der Druckzone (Feldmoment)
                    False -> Platte gezogen; Bemessung als Rechteck mit bw
                             (Stuetzmoment)
    """
    b: float = 300.0
    h: float = 600.0
    d1: float = 50.0
    d2: float = 50.0
    typ: str = "rechteck"
    b_eff: float = None
    hf: float = 0.0
    platte_gedrueckt: bool = True

    def __post_init__(self):
        if self.typ == "plattenbalken":
            if not self.b_eff or self.b_eff < self.b:
                raise ValueError("b_eff muss >= b (Stegbreite) sein")
            if self.hf <= 0:
                raise ValueError("hf muss > 0 sein beim Plattenbalken")
        else:
            self.b_eff = self.b
            self.hf = 0.0

    # -- Grundgroessen -----------------------------------------------------
    @property
    def d(self):
        """Statische Nutzhoehe d = h - d1 [mm]."""
        return self.h - self.d1

    @property
    def bw(self):
        """Stegbreite [mm]."""
        return self.b

    @property
    def Ac(self):
        """Bruttobetonflaeche [mm2]."""
        if self.typ == "plattenbalken":
            return self.b * self.h + (self.b_eff - self.b) * self.hf
        return self.b * self.h

    @property
    def u_aussen(self):
        """Aeusserer Umfang [mm] (fuer die Torsion, EC2 6.3.1)."""
        if self.typ == "plattenbalken":
            return 2.0 * (self.b_eff + self.h)
        return 2.0 * (self.b + self.h)

    @property
    def z_s1(self):
        """
        Abstand vom Schwerpunkt der Betonflaeche zum Schwerpunkt der
        Zugbewehrung [mm], positiv nach unten.
        Verwendet fuer die Umrechnung von N_Ed: M_Eds = M_Ed - N_Ed z_s1.
        """
        return self.y_schwerpunkt_zug()

    def y_schwerpunkt(self):
        """Schwerpunktlage vom oberen Rand [mm]."""
        if self.typ == "plattenbalken" and self.platte_gedrueckt:
            A1 = self.b_eff * self.hf
            A2 = self.b * (self.h - self.hf)
            return (A1 * self.hf / 2.0
                    + A2 * (self.hf + (self.h - self.hf) / 2.0)) / (A1 + A2)
        return self.h / 2.0

    def y_schwerpunkt_zug(self):
        """Abstand Schwerpunkt -> Schwerpunkt As1 [mm]."""
        return self.d - self.y_schwerpunkt()

    def breite(self, y):
        """Breite [mm] in der Tiefe y ab dem oberen Rand."""
        if self.typ == "plattenbalken" and self.platte_gedrueckt:
            return self.b_eff if y <= self.hf else self.b
        return self.b

    # -- Traegheitsmomente (Zustand I, ungerissen) -------------------------
    def traegheitsmoment(self):
        """Bruttotraegheitsmoment [mm4] um den Schwerpunkt."""
        if self.typ == "plattenbalken" and self.platte_gedrueckt:
            yc = self.y_schwerpunkt()
            I = (self.b_eff * self.hf ** 3 / 12.0
                 + self.b_eff * self.hf * (yc - self.hf / 2.0) ** 2)
            hw = self.h - self.hf
            I += (self.b * hw ** 3 / 12.0
                  + self.b * hw * (self.hf + hw / 2.0 - yc) ** 2)
            return I
        return self.b * self.h ** 3 / 12.0

    def widerstandsmoment_unten(self):
        """W = I / (h - y_s) bezogen auf den Zugrand [mm3]."""
        return self.traegheitsmoment() / (self.h - self.y_schwerpunkt())

    def beschreibung(self):
        if self.typ == "plattenbalken":
            return ("Plattenbalken bw/h = {:.0f}/{:.0f} mm, "
                    "b_eff = {:.0f} mm, hf = {:.0f} mm"
                    .format(self.b, self.h, self.b_eff, self.hf))
        return "Rechteckquerschnitt b/h = {:.0f}/{:.0f} mm".format(self.b, self.h)


def mitwirkende_plattenbreite(bw, l0, b1, b2=None, innenbalken=True):
    """
    Mitwirkende Plattenbreite nach DIN EN 1992-1-1, 5.3.2.1 (3), Gl. (5.7):

        b_eff   = sum(b_eff,i) + bw  <= b
        b_eff,i = 0,2 bi + 0,1 l0 <= 0,2 l0   und  <= bi

    Parameter
    ---------
    bw : float   Stegbreite [mm]
    l0 : float   Abstand der Momentennullpunkte [mm]
                 (Einfeldtraeger l0 = 0,85 l ; Innenfeld 0,70 l ;
                  Kragarm 1,50 l ... DIN EN 1992-1-1, 5.3.2.1 (2), Bild 5.2)
    b1, b2 : float  halbe lichte Plattenbreiten beidseits [mm]
    """
    def bi(b):
        if b is None:
            return 0.0
        return min(0.2 * b + 0.1 * l0, 0.2 * l0, b)

    b_eff = bi(b1) + (bi(b2) if innenbalken else 0.0) + bw
    return dict(b_eff=b_eff, l0=l0, norm=ref("b_eff"))
