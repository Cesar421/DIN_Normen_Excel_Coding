# -*- coding: utf-8 -*-
"""
Verzeichnis der Normstellen fuer die PFAHLBEMESSUNG.
================================================================================
Zugrunde gelegte Normen
--------------------------------------------------------------------------------
[EN1536]  DIN EN 1536:2015-10
          "Ausfuehrung von Arbeiten im Spezialtiefbau - Bohrpfaehle"
          -> Ausfuehrung, Betondeckung, Mindestbewehrung, Bewehrungskorb
[EC7]     DIN EN 1997-1:2009-09 + DIN EN 1997-1/NA:2010-12
          "Eurocode 7: Entwurf, Berechnung und Bemessung in der Geotechnik"
[DIN1054] DIN 1054:2010-12 (+ A1:2012-08)
          "Baugrund - Sicherheitsnachweise im Erd- und Grundbau -
           Ergaenzende Regelungen zu DIN EN 1997-1"
[EA]      EA-Pfaehle (Empfehlungen des Arbeitskreises "Pfaehle" der DGGT),
          2. Auflage - Erfahrungswerte q_b,k / q_s,k und Widerstands-Setzungs-Linien
[EC2]     DIN EN 1992-1-1:2011-01 + /NA:2013-04  (Querschnittsbemessung)
[EC0]     DIN EN 1990:2010-12 + /NA               (Einwirkungskombinationen)
[DIN488]  DIN 488-1:2009-08                       (Betonstahl)
================================================================================
Weitere Ausfuehrungsnormen (hier nicht im Detail umgesetzt):
  DIN EN 12699  Verdraengungspfaehle
  DIN EN 14199  Mikropfaehle
  DIN 4093      Duesenstrahlverfahren
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Normstelle:
    id: str
    norm: str
    abschnitt: str
    titel: str
    gleichung: str = ""

    def __str__(self):
        gl = ", " + self.gleichung if self.gleichung else ""
        return "[" + self.norm.split(":")[0] + ", " + self.abschnitt + gl + "]"

    def lang(self):
        gl = " " + self.gleichung if self.gleichung else ""
        return self.norm + ", Abschnitt " + self.abschnitt + gl + " - " + self.titel


_EN1536 = "DIN EN 1536:2015-10"
_EC7 = "DIN EN 1997-1:2009-09"
_D1054 = "DIN 1054:2010-12"
_EA = "EA-Pfaehle (DGGT), 2. Auflage"
_EC2 = "DIN EN 1992-1-1:2011-01"
_NA2 = "DIN EN 1992-1-1/NA:2013-04"
_EC0 = "DIN EN 1990:2010-12"
_EC0NA = "DIN EN 1990/NA:2010-12"

_N = {}


def _add(id_, norm, abschnitt, titel, gleichung=""):
    _N[id_] = Normstelle(id_, norm, abschnitt, titel, gleichung)


# ---------------------------------------------------------------------------
# AUSFUEHRUNG UND KONSTRUKTION (DIN EN 1536)
# ---------------------------------------------------------------------------
_add("pfahl_allg", _EN1536, "1 / 3",
     "Anwendungsbereich und Begriffe fuer Bohrpfaehle")
_add("pfahl_beton", _EN1536, "6.3",
     "Anforderungen an den Pfahlbeton (Konsistenz, Groesstkorn, Mindestgehalt)")
_add("pfahl_deckung", _EN1536, "7.6.2",
     "Betondeckung: c_nom >= 60 mm bei D >= 0,6 m; >= 50 mm bei D < 0,6 m; "
     ">= 75 mm beim Betonieren unter Stuetzfluessigkeit")
_add("pfahl_As_min", _EN1536, "7.6.3, Tab. 4",
     "Mindestlaengsbewehrung: 0,5 % Ac (Ac <= 0,5 m2) / 25 cm2 "
     "(0,5 < Ac <= 1,0 m2) / 0,25 % Ac (Ac > 1,0 m2)")
_add("pfahl_laengs", _EN1536, "7.6.3",
     "Laengsbewehrung: mindestens 6 Staebe, phi >= 16 mm, lichter Abstand "
     ">= 100 mm (>= 80 mm bei Groesstkorn <= 20 mm)")
_add("pfahl_quer", _EN1536, "7.6.4",
     "Querbewehrung (Wendel oder Ringe): phi >= 6 mm, lichter Abstand "
     ">= 100 mm und <= 400 mm")
_add("pfahl_korb", _EN1536, "7.6.5",
     "Bewehrungskorb: Abstandhalter, Steifigkeit, Transport und Einbau")

# ---------------------------------------------------------------------------
# GEOTECHNIK (DIN EN 1997-1 / DIN 1054 / EA-Pfaehle)
# ---------------------------------------------------------------------------
_add("geo_nachweis", _EC7, "7.6.2.1",
     "Nachweis der axialen Pfahltragfaehigkeit: F_c,d <= R_c,d")
_add("geo_widerstand", _EC7, "7.6.2.3, Gl. (7.8)",
     "R_c,k = R_b,k + R_s,k = q_b,k A_b + sum(q_s,k,i A_s,i)")
_add("geo_gamma", _D1054, "A 7.6.2.2, Tab. A 2.3",
     "Teilsicherheitsbeiwerte der Pfahlwiderstaende: gamma_b = gamma_s = "
     "gamma_t = 1,10 (BS-P) ; Zugpfaehle gamma_s,t = 1,15 (BS-P)")
_add("geo_einwirkung", _D1054, "A 2.4.7.6.1, Tab. A 2.1",
     "Teilsicherheitsbeiwerte der Einwirkungen im Nachweis GEO-2: "
     "gamma_G = 1,35 ; gamma_Q = 1,50 (BS-P)")
_add("geo_erfahrungswerte", _EA, "Tab. 5.12 bis 5.15",
     "Erfahrungswerte q_b,k und q_s,k fuer Bohrpfaehle in nichtbindigen und "
     "bindigen Boeden (Eingabewerte des Anwenders)")
_add("geo_wsl", _EA, "5.4.5",
     "Widerstands-Setzungs-Linie (WSL) aus Erfahrungswerten")
_add("geo_zugpfahl", _EC7, "7.6.3",
     "Nachweis von Zugpfaehlen (Herausziehwiderstand)")
_add("geo_knicken", _EA, "4.7",
     "Knicknachweis nur bei sehr weichen Boeden erforderlich "
     "(Richtwert c_u < 10 kN/m2)")
_add("geo_bettung", _EA, "6.3",
     "Horizontal belastete Pfaehle: Bettungsmodulverfahren, k_s = E_s/D "
     "(D <= 1 m) bzw. k_s = E_s/1 m (D > 1 m)")
_add("geo_gruppen", _EA, "8",
     "Pfahlgruppen: Gruppenwirkung und Verteilung der Pfahlkraefte "
     "(hier NICHT erfasst)")

# ---------------------------------------------------------------------------
# QUERSCHNITTSBEMESSUNG (DIN EN 1992-1-1)
# ---------------------------------------------------------------------------
_add("qs_biegung", _EC2, "6.1",
     "Bemessung fuer Biegung mit Laengskraft - Kreisquerschnitt, "
     "M-N-Interaktionsdiagramm", "Bild 6.1")
_add("qs_punkt_c", _EC2, "6.1 (5), Bild 6.1",
     "Bemessungspunkt C bei ueberdrueckten Querschnitten: eps_c2 im Abstand "
     "(1 - eps_c2/eps_cu2) h vom gedrueckten Rand")
_add("qs_querkraft", _EC2, "6.2",
     "Querkraftnachweis; fuer Kreisquerschnitte enthaelt EC2 keine expliziten "
     "Regeln - es werden bw = D und d = D/2 + D_s/pi angesetzt")
_add("qs_VRdc", _NA2, "NDP zu 6.2.2 (1)",
     "C_Rd,c = 0,15/gamma_C ; k1 = 0,12 ; v_min nach Gl. (6.3aDE)")
_add("qs_cot", _NA2, "NDP zu 6.2.3 (2), Gl. (6.7aDE)",
     "Druckstrebenneigung 1,0 <= cot(theta) <= 3,0")
_add("qs_theta2", _EC2, "5.8", "Nachweis nach Theorie II. Ordnung (Knicken)")
_add("qs_wendel", _EC2, "9.5.3",
     "Querbewehrung von Druckgliedern: phi_w >= max(6 mm ; phi_l/4); "
     "s_cl,t <= min(20 phi_l ; kleinste Abmessung ; 400 mm)")
_add("qs_As_stuetze", _NA2, "NDP zu 9.5.2 (2)",
     "Mindestlaengsbewehrung von Druckgliedern As,min = 0,15 |N_Ed|/fyd")
_add("qs_fcd", _EC2, "3.1.6 (1)P", "fcd = alpha_cc fck/gamma_C", "Gl. (3.15)")

# ---------------------------------------------------------------------------
# EINWIRKUNGEN
# ---------------------------------------------------------------------------
_add("komb_GZT", _EC0, "6.4.3.2, Gl. (6.10)",
     "Grundkombination GZT: sum(gamma_G G) + gamma_Q Q_1 + sum(gamma_Q psi_0 Q_i)")
_add("komb_QS", _EC0, "6.5.3, Gl. (6.16b)",
     "Quasi-staendige Kombination GZG: sum(G) + sum(psi_2 Q)")
_add("gamma_F", _EC0NA, "Tab. NA.A.1.2(B)",
     "gamma_G = 1,35 / gamma_Q = 1,50")


def ref(id_):
    """Gibt die unter `id_` hinterlegte Normstelle zurueck."""
    if id_ not in _N:
        raise KeyError("Normstelle nicht hinterlegt: " + repr(id_))
    return _N[id_]


def alle():
    return list(_N.values())


def normentabelle():
    zeilen = ["VERWENDETE NORMSTELLEN - PFAHLBEMESSUNG", "=" * 96]
    aktuell = None
    for n in sorted(alle(), key=lambda x: (x.norm, x.abschnitt)):
        if n.norm != aktuell:
            aktuell = n.norm
            zeilen.append("")
            zeilen.append("--- " + aktuell + " ---")
        gl = " " + n.gleichung if n.gleichung else ""
        zeilen.append("  {:<30s}{:<18s} {}".format(n.abschnitt, gl, n.titel))
    return "\n".join(zeilen)
