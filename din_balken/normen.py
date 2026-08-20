# -*- coding: utf-8 -*-
"""
Verzeichnis der angewendeten deutschen Normstellen (DIN).
================================================================================
Zugrunde gelegte Normen
--------------------------------------------------------------------------------
[EC2]       DIN EN 1992-1-1:2011-01
            "Eurocode 2: Bemessung und Konstruktion von Stahlbeton- und
             Spannbetontragwerken - Teil 1-1: Allgemeine Bemessungsregeln und
             Regeln fuer den Hochbau"
[NA]        DIN EN 1992-1-1/NA:2013-04 (+ A1:2015-12)
            Nationaler Anhang - National festgelegte Parameter (NDP)
[EC0]       DIN EN 1990:2010-12 + DIN EN 1990/NA:2010-12  (Tragwerksplanung)
[EC1]       DIN EN 1991-1-1:2010-12 + /NA:2010-12         (Einwirkungen, Wichten)
[DIN488]    DIN 488-1:2009-08                             (Betonstahl B500A/B500B)
[DIN1045-2] DIN EN 206-1 / DIN 1045-2                     (Beton, Expositionsklassen)

Historie: DIN 1045-1:2008-08 wurde 2010 zurueckgezogen und durch [EC2]+[NA]
ersetzt. Die entsprechende Stelle in DIN 1045-1 ist informativ angegeben
(Feld `din1045`).
================================================================================
Jede Berechnung des Pakets gibt die angewendete Normstelle zurueck.
Aufruf:  ref("biegung")  ->  Objekt Normstelle mit Norm, Abschnitt, Titel,
                            Gleichung.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Normstelle:
    """Eine einzelne Normstelle."""
    id: str
    norm: str           # z.B. "DIN EN 1992-1-1:2011-01"
    abschnitt: str      # z.B. "6.2.3 (3)"
    titel: str          # Kurzbeschreibung
    gleichung: str = ""  # z.B. "Gl. (6.9)"
    din1045: str = ""   # entsprechender Abschnitt in DIN 1045-1:2008 (informativ)

    def __str__(self):
        gl = ", " + self.gleichung if self.gleichung else ""
        return "[" + self.norm.split(":")[0] + ", " + self.abschnitt + gl + "]"

    def lang(self):
        gl = " " + self.gleichung if self.gleichung else ""
        alt = "  (~ DIN 1045-1, " + self.din1045 + ")" if self.din1045 else ""
        return self.norm + ", Abschnitt " + self.abschnitt + gl + " - " + self.titel + alt


_EC2 = "DIN EN 1992-1-1:2011-01"
_NA = "DIN EN 1992-1-1/NA:2013-04"
_EC0 = "DIN EN 1990:2010-12"
_EC0NA = "DIN EN 1990/NA:2010-12"
_EC1 = "DIN EN 1991-1-1:2010-12"
_D488 = "DIN 488-1:2009-08"

_N = {}


def _add(id_, norm, abschnitt, titel, gleichung="", din1045=""):
    _N[id_] = Normstelle(id_, norm, abschnitt, titel, gleichung, din1045)


# ---------------------------------------------------------------------------
# 1. GRUNDLAGEN DER BEMESSUNG / SICHERHEITSKONZEPT
# ---------------------------------------------------------------------------
_add("gamma_M", _NA, "2.4.2.4 (1), Tab. 2.1DE",
     "Teilsicherheitsbeiwerte der Baustoffe gamma_C = 1,50 / gamma_S = 1,15",
     din1045="5.3.3")
_add("komb_GZT", _EC0, "6.4.3.2, Gl. (6.10)",
     "Grundkombination GZT: sum(gamma_G G) + gamma_Q Q_1 + sum(gamma_Q psi_0 Q_i)")
_add("gamma_F", _EC0NA, "Tab. NA.A.1.2(B)",
     "Teilsicherheitsbeiwerte der Einwirkungen gamma_G = 1,35 / gamma_Q = 1,50")
_add("komb_QS", _EC0, "6.5.3, Gl. (6.16b)",
     "Quasi-staendige Kombination GZG: sum(G) + sum(psi_2 Q)")
_add("psi", _EC0NA, "Tab. NA.A.1.1", "Kombinationsbeiwerte psi_0 / psi_1 / psi_2")
_add("eigengewicht", _EC1, "Tab. A.1",
     "Wichte von Stahlbeton gamma = 25,0 kN/m3")

# ---------------------------------------------------------------------------
# 2. BAUSTOFFE
# ---------------------------------------------------------------------------
_add("beton_tab31", _EC2, "3.1.2, Tab. 3.1",
     "Festigkeitsklassen und Eigenschaften des Betons (fck, fcm, fctm, Ecm)",
     din1045="9.1.2, Tab. 9")
_add("fcd", _EC2, "3.1.6 (1)P",
     "Bemessungswert der Betondruckfestigkeit fcd = alpha_cc fck / gamma_C",
     "Gl. (3.15)", din1045="9.1.6")
_add("alpha_cc", _NA, "NDP zu 3.1.6 (1)P", "alpha_cc = 0,85 (Normalbeton)")
_add("fctd", _EC2, "3.1.6 (2)P",
     "Bemessungswert der Betonzugfestigkeit fctd = alpha_ct fctk;0,05 / gamma_C",
     "Gl. (3.16)")
_add("alpha_ct", _NA, "NDP zu 3.1.6 (2)P", "alpha_ct = 0,85")
_add("sigma_eps_c", _EC2, "3.1.7 (1), Bild 3.3",
     "Parabel-Rechteck-Diagramm fuer die Querschnittsbemessung",
     "Gl. (3.17)/(3.18)", din1045="9.1.5")
_add("betonstahl_din488", _D488, "Tab. 4 / 5",
     "Betonstahl B500A / B500B: fyk = 500 N/mm2")
_add("betonstahl_ec2", _EC2, "3.2.2 / 3.2.7, Bild 3.8",
     "Bemessungs-Spannungs-Dehnungs-Linie des Betonstahls (waagerechter oberer Ast)",
     din1045="9.2.3")
_add("eps_ud", _NA, "NDP zu 3.2.7 (2)", "Grenzdehnung des Betonstahls eps_ud = 25 permil")

# ---------------------------------------------------------------------------
# 3. DAUERHAFTIGKEIT UND BETONDECKUNG
# ---------------------------------------------------------------------------
_add("expos", _EC2, "4.2, Tab. 4.1", "Expositionsklassen (X0, XC, XD, XS, XF, XA)")
_add("c_nom", _EC2, "4.4.1.1 (2)P", "c_nom = c_min + Delta_c_dev", "Gl. (4.1)",
     din1045="6.3")
_add("c_min", _EC2, "4.4.1.2 (2)P", "c_min = max(c_min,b ; c_min,dur ; 10 mm)",
     "Gl. (4.2)")
_add("c_min_dur", _NA, "NDP zu 4.4.1.2 (5), Tab. 4.4DE",
     "Mindestbetondeckung aus Dauerhaftigkeit je Expositionsklasse")
_add("dc_dev", _NA, "NDP zu 4.4.1.3 (1)P", "Delta_c_dev = 15 mm (10 mm bei XC1)")

# ---------------------------------------------------------------------------
# 5. SCHNITTGROESSENERMITTLUNG
# ---------------------------------------------------------------------------
_add("b_eff", _EC2, "5.3.2.1 (3)", "Mitwirkende Plattenbreite von Plattenbalken",
     "Gl. (5.7)/(5.7a)/(5.7b)", din1045="7.3.1")
_add("l_eff", _EC2, "5.3.2.2", "Effektive Stuetzweite l_eff", "Gl. (5.8)")
_add("lastfaelle", _EC2, "5.1.3", "Lastfaelle und Lastkombinationen (Feldweise Anordnung)")
_add("umlagerung", _NA, "NDP zu 5.5 (4)",
     "Momentenumlagerung: delta >= k1 + k2 xu/d mit k1 = 0,64 ; k2 = 0,80 (<= C50/60)",
     "Gl. (5.10a)")
_add("xi_lim", _NA, "NDP zu 5.5 (4)",
     "Duktilitaetsgrenze xu/d <= 0,45 (<= C50/60) / <= 0,35 (>= C55/67) bei delta = 1,0")

# ---------------------------------------------------------------------------
# 6. GRENZZUSTAENDE DER TRAGFAEHIGKEIT (GZT)
# ---------------------------------------------------------------------------
_add("biegung", _EC2, "6.1", "Biegung mit und ohne Laengskraft - Bemessungsannahmen",
     "Bild 6.1", din1045="10.2")
_add("eps_grenzen", _EC2, "6.1 (2)P, Bild 6.1",
     "Zulaessige Dehnungsverteilungen (Dehnungsdiagramm mit Bemessungspunkten A/B/C)")
_add("druckbewehrung", _EC2, "6.1", "Druckbewehrung wenn xu > xu,lim")
_add("querkraft_allg", _EC2, "6.2.1", "Querkraft - allgemeines Nachweisverfahren",
     din1045="10.3")
_add("VRdc", _NA, "NDP zu 6.2.2 (1)",
     "Querkrafttragfaehigkeit ohne Querkraftbewehrung: C_Rd,c = 0,15/gamma_C ; k1 = 0,12",
     "Gl. (6.2a)/(6.2b)")
_add("vmin", _NA, "NDP zu 6.2.2 (1), Gl. (6.3aDE)",
     "v_min = (kappa_1/gamma_C) k^1,5 fck^0,5 ; kappa_1 = 0,0525 (d<=600) / 0,0375 (d>=800)")
_add("z_innen", _NA, "NDP zu 6.2.3 (1)",
     "Innerer Hebelarm z = min(0,9 d ; d - 2 c_v,l ; d - c_v,l - 30 mm)")
_add("cot_theta", _NA, "NDP zu 6.2.3 (2), Gl. (6.7aDE)",
     "Druckstrebenneigung: 1,0 <= cot(theta) <= 3,0 (Normalbeton)")
_add("VRdcc", _NA, "NDP zu 6.2.3 (2), Gl. (6.7bDE)",
     "V_Rd,cc = c 0,48 fck^(1/3) (1 - 1,2 sigma_cd/fcd) bw z mit c = 0,5")
_add("VRdmax", _EC2, "6.2.3 (3)", "Druckstrebentragfaehigkeit V_Rd,max", "Gl. (6.9)",
     din1045="10.3.4")
_add("nu1", _NA, "NDP zu 6.2.3 (3)", "nu_1 = 0,75 nu_2 ; nu_2 = 1,0 fuer <= C50/60")
_add("VRds", _EC2, "6.2.3 (3)",
     "Querkraftbewehrung V_Rd,s = (Asw/s) z fywd cot(theta)", "Gl. (6.8)")
_add("V_rand", _EC2, "6.2.1 (8)",
     "Bemessungsquerkraft im Abstand d vom Auflagerrand (direkte Lagerung)")
_add("versatzmass", _EC2, "9.2.1.3 (2)",
     "Versatzmass der Zugkraftlinie a_l = z cot(theta)/2", "Gl. (9.2)",
     din1045="13.2.2")
_add("auflagerkraft", _EC2, "9.2.1.4 (2)",
     "Zu verankernde Zugkraft am Endauflager F_Ed = |V_Ed| a_l/z + N_Ed", "Gl. (9.3)")

# --- Torsion (EC2 6.3) -----------------------------------------------------
_add("torsion_allg", _EC2, "6.3.1",
     "Torsion - Ersatzhohlquerschnitt duennwandiger geschlossener Querschnitt",
     din1045="10.4")
_add("torsion_schubfluss", _EC2, "6.3.2 (1)",
     "Schubfluss tau_t,i t_ef,i = T_Ed/(2 A_k)", "Gl. (6.26)")
_add("torsion_VEdi", _EC2, "6.3.2 (1)",
     "Querkraft je Wand V_Ed,i = tau_t,i t_ef,i z_i", "Gl. (6.27)")
_add("torsion_laengs", _EC2, "6.3.2 (3)",
     "Torsionslaengsbewehrung sum(Asl) fyd / u_k = T_Ed/(2 A_k) cot(theta)",
     "Gl. (6.28)")
_add("torsion_interaktion", _EC2, "6.3.2 (4)",
     "Interaktion Torsion/Querkraft: T_Ed/T_Rd,max + V_Ed/V_Rd,max <= 1,0",
     "Gl. (6.29)")
_add("TRdmax", _EC2, "6.3.2 (4)",
     "Torsionsdruckstrebe T_Rd,max = 2 nu alpha_cw fcd A_k t_ef,i sin(theta) cos(theta)",
     "Gl. (6.30)", din1045="10.4.2")
_add("nu_torsion", _NA, "NDP zu 6.3.2 (4)",
     "nu = 0,525 nu_2 (Vollquerschnitt) / 0,75 nu_2 (Kastenquerschnitt)")
_add("torsion_ohne", _NA, "NDP zu 6.3.2 (5), Gl. (6.31aDE)/(6.31bDE)",
     "Torsionsbewehrung darf entfallen wenn T_Ed <= V_Ed bw/4,5 und "
     "V_Ed [1 + 4,5 T_Ed/(V_Ed bw)] <= V_Rd,c")
_add("TRdc", _EC2, "6.3.2 (5)",
     "Risstorsionsmoment T_Rd,c = fctd t_ef 2 A_k", "Gl. (6.31)")
_add("torsion_konstr", _EC2, "9.2.3",
     "Konstruktive Durchbildung der Torsionsbewehrung: geschlossene Buegel, "
     "s <= u_k/8 und <= kleinste Querschnittsabmessung; Laengsstab in jeder Ecke")

# ---------------------------------------------------------------------------
# 7. GRENZZUSTAENDE DER GEBRAUCHSTAUGLICHKEIT (GZG)
# ---------------------------------------------------------------------------
_add("w_max", _NA, "NDP zu 7.3.1 (5), Tab. 7.1DE",
     "Zulaessige Rissbreite w_max je Expositionsklasse")
_add("As_min_riss", _EC2, "7.3.2 (2)",
     "Mindestbewehrung fuer die Rissbreitenbegrenzung As,min sigma_s = kc k fct,eff Act",
     "Gl. (7.1)")
_add("wk", _EC2, "7.3.4 (1)", "Rechnerische Rissbreite wk = sr,max (eps_sm - eps_cm)",
     "Gl. (7.8)", din1045="11.2.4")
_add("eps_sm", _EC2, "7.3.4 (2)", "Mittlere Dehnungsdifferenz Stahl - Beton", "Gl. (7.9)")
_add("sr_max", _NA, "NDP zu 7.3.4 (3), Gl. (7.11DE)",
     "sr,max = phi/(3,6 rho_p,eff) <= sigma_s phi/(3,6 fct,eff)")
_add("hc_eff", _EC2, "7.3.2 (3)",
     "Wirkungszone hc,ef = min(2,5(h-d) ; (h-x)/3 ; h/2)")
_add("durchbiegung_ld", _EC2, "7.4.2 (2)",
     "Begrenzung der Verformung ueber die Biegeschlankheit l/d",
     "Gl. (7.16a)/(7.16b)", din1045="11.3.2")
_add("durchbiegung_NA", _NA, "NDP zu 7.4.2",
     "Vereinfachte Regel l/d <= 35 (bzw. l/d <= 35 (7/l) bei erhoehten Anforderungen)")
_add("durchbiegung_rech", _EC2, "7.4.3 (3)",
     "Verformungsberechnung: Interpolation Zustand I / Zustand II",
     "Gl. (7.18)/(7.19)")
_add("kriechen", _EC2, "3.1.4, Bild 3.1", "Kriechzahl phi(inf,t0)")

# ---------------------------------------------------------------------------
# 8. VERBUND UND VERANKERUNG
# ---------------------------------------------------------------------------
_add("fbd", _EC2, "8.4.2 (2)", "Verbundspannung fbd = 2,25 eta1 eta2 fctd",
     "Gl. (8.2)", din1045="12.5")
_add("lb_rqd", _EC2, "8.4.3 (2)",
     "Grundwert der Verankerungslaenge lb,rqd = (phi/4)(sigma_sd/fbd)", "Gl. (8.3)")
_add("lbd", _EC2, "8.4.4 (1)",
     "Verankerungslaenge lbd = alpha1..alpha5 lb,rqd >= lb,min", "Gl. (8.4)/(8.6)")
_add("l0", _EC2, "8.7.3",
     "Uebergreifungslaenge l0 = alpha1..alpha6 lb,rqd >= l0,min", "Gl. (8.10)")

# ---------------------------------------------------------------------------
# 9. BEWEHRUNGS- UND KONSTRUKTIONSREGELN
# ---------------------------------------------------------------------------
_add("As_min", _EC2, "9.2.1.1 (1)",
     "As,min = 0,26 (fctm/fyk) bt d >= 0,0013 bt d", "Gl. (9.1N)", din1045="13.1.1")
_add("As_max", _NA, "NDP zu 9.2.1.1 (3)",
     "As,max = 0,04 Ac ausserhalb von Uebergreifungsstoessen")
_add("robustheit", _NA, "NDP zu 9.2.1.1 (1)",
     "Mindestbewehrung fuer duktiles Bauteilverhalten (Robustheitsbewehrung) aus Mcr")
_add("rho_w_min", _NA, "NDP zu 9.2.2 (5)",
     "Mindestquerkraftbewehrungsgrad rho_w,min = 0,16 fctm/fyk (Balken)")
_add("s_max", _NA, "NDP zu 9.2.2 (6)/(8), Tab. NA.9.1",
     "Groesste Laengs- und Querabstaende der Buegel")
_add("stababstand", _EC2, "8.2 (2)",
     "Lichter Stababstand >= max(phi ; d_g + 5 mm ; 20 mm)")


def ref(id_):
    """Gibt die unter `id_` hinterlegte Normstelle zurueck."""
    if id_ not in _N:
        raise KeyError("Normstelle nicht hinterlegt: " + repr(id_))
    return _N[id_]


def alle():
    """Vollstaendige Liste der vom Programm verwendeten Normstellen."""
    return list(_N.values())


def normentabelle():
    """Texttabelle mit allen verwendeten DIN-Normstellen."""
    zeilen = ["VERWENDETE NORMSTELLEN", "=" * 96]
    aktuell = None
    for n in sorted(alle(), key=lambda x: (x.norm, x.abschnitt)):
        if n.norm != aktuell:
            aktuell = n.norm
            zeilen.append("")
            zeilen.append("--- " + aktuell + " ---")
        gl = " " + n.gleichung if n.gleichung else ""
        zeilen.append("  {:<32s}{:<22s} {}".format(n.abschnitt, gl, n.titel))
    return "\n".join(zeilen)
