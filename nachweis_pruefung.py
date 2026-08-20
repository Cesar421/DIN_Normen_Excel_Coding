# -*- coding: utf-8 -*-
"""
Pruefung der Rechenkerne `din_balken` und `din_pfahl`.

Die Ergebnisse werden gegen geschlossene analytische Loesungen und gegen
bekannte Referenzwerte der deutschen Praxis geprueft. Aufruf:

    python nachweis_pruefung.py
"""

import math
import sys

import numpy as np

from din_balken.baustoffe import Beton, Betonstahl, stabflaeche_n
from din_balken.querschnitt import Querschnitt
from din_balken.schnittgroessen import Durchlauftraeger, Auflager, GELENKIG, EINGESPANNT
from din_balken.biegung import (bemessung_biegung, momententragfaehigkeit,
                                xi_grenz, bemessungsdiagramm)
from din_balken.querkraft import V_Rd_c, cot_theta_NA, V_Rd_max, asw_mindest
from din_balken.torsion import (ersatzhohlquerschnitt, T_Rd_c, T_Rd_max,
                                asw_torsion, asl_torsion, bemessung_torsion,
                                bewehrung_entbehrlich)
from din_balken.gebrauchstauglichkeit import zustand_II
from din_balken.bemessung import EingabeBalken, bemessung_balken

from din_pfahl.kreisquerschnitt import (Kreisquerschnitt, interaktionsdiagramm,
                                        M_Rd_bei_N, mindestbewehrung_pfahl,
                                        betondeckung_pfahl)
from din_pfahl.bettung import pfahl_horizontal, Bodenschicht, bettungsmodul
from din_pfahl.tragfaehigkeit import axiale_tragfaehigkeit, GAMMA_R
from din_pfahl.bemessung_pfahl import EingabePfahl, bemessung_pfahl

FEHLER = []


def pruef(name, ist, soll, tol=1e-3):
    d = abs(ist - soll)
    ok = d / max(abs(soll), 1e-12) <= tol
    print("  [{}] {:<56s} {:>13.5f}  (Soll {:.5f})"
          .format("OK  " if ok else "NEIN", name, ist, soll))
    if not ok:
        FEHLER.append(name)
    return ok


print("=" * 100)
print("PRUEFUNG DER RECHENKERNE  -  DIN EN 1992-1-1 + NA | DIN EN 1536 / "
      "DIN 1054 / EA-Pfaehle")
print("=" * 100)

C30 = Beton("C30/37")
B500 = Betonstahl("B500B")

# ---------------------------------------------------------------------------
print("\n1. BAUSTOFFE  [EC2 3.1 / 3.2 + NA]")
pruef("fcd = 0,85*30/1,5", C30.fcd, 17.0)
pruef("fctm = 0,30*30^(2/3)", C30.fctm, 0.30 * 30 ** (2 / 3))
pruef("Ecm = 22000*(38/10)^0,3", C30.Ecm, 22000 * 3.8 ** 0.3)
pruef("fyd = 500/1,15", B500.fyd, 500 / 1.15)
pruef("eps_ud (NA NDP zu 3.2.7(2))", B500.eps_ud, 25.0)

print("\n2. DRUCKZONE PARABEL-RECHTECK  [EC2 3.1.7 (1)]")
pruef("alpha_R(3,5 permil) = 1 - 2/(3*3,5)", C30.alpha_R(3.5), 1 - 2 / 10.5)
_m = ((5 / 12) * 2.0 ** 2 + 0.5 * (3.5 ** 2 - 2.0 ** 2)) / 3.5 ** 2
pruef("k_a(3,5 permil) analytisch", C30.k_a(3.5), 1 - _m / (1 - 2 / 10.5), 1e-6)
pruef("alpha_R(2,0 permil) = 2/3", C30.alpha_R(2.0), 2 / 3)
pruef("k_a(2,0 permil) = 3/8", C30.k_a(2.0), 0.375)

print("\n3. DUKTILITAETSGRENZE UND mu_lim  [NA NDP zu 5.5 (4)]")
xl = xi_grenz(C30, 1.0)["xi_lim"]
aR, ka = C30.alpha_R(3.5), C30.k_a(3.5)
pruef("xi_lim (delta = 1,0)", xl, 0.45)
pruef("mu_lim = alpha_R xi (1 - k_a xi)", aR * xl * (1 - ka * xl), 0.29610, 1e-4)
pruef("omega_lim = alpha_R xi_lim", aR * xl, 0.36429, 1e-4)
pruef("zeta_lim = 1 - k_a xi_lim", 1 - ka * xl, 0.81281, 1e-4)
pruef("xi_lim C55/67", xi_grenz(Beton("C55/67"), 1.0)["xi_lim"], 0.35)

print("\n4. BIEGUNG: M_Ed -> As -> M_Rd  [EC2 6.1]")
qs = Querschnitt(b=300., h=600., d1=50., d2=50.)
for M in (50., 150., 300., 456.65, 700.):
    r = bemessung_biegung(qs, C30, B500, M)
    mr = momententragfaehigkeit(qs, C30, B500, r.As1, r.As2)
    pruef("M_Rd(As) fuer M_Ed = %.1f kNm" % M, mr["M_Rd"], M, 2e-3)
pruef("xi bei mu_Eds = 0,296", bemessung_biegung(qs, C30, B500, 456.65).xi, 0.45, 3e-3)
r = bemessung_biegung(qs, C30, B500, 200., N_Ed=-400.)
pruef("M_Rd mit N_Ed = -400 kN",
      momententragfaehigkeit(qs, C30, B500, r.As1, r.As2, N_Ed=-400.)["M_Rd"],
      200., 2e-3)
qsT = Querschnitt(b=300., h=600., d1=50., d2=50., typ="plattenbalken",
                  b_eff=1500., hf=120.)
for M in (300., 900., 1500.):
    r = bemessung_biegung(qsT, C30, B500, M)
    pruef("M_Rd Plattenbalken fuer M_Ed = %.0f kNm" % M,
          momententragfaehigkeit(qsT, C30, B500, r.As1, r.As2)["M_Rd"], M, 2e-3)

print("\n5. SCHNITTGROESSEN  (geschlossene Loesungen)")
t = Durchlauftraeger(6.0, [Auflager(0.0), Auflager(6.0)])
t.strecke(0, 6, 10.0)
e = t.berechnen()
pruef("Einfeldtraeger: M_max = qL^2/8", e.M.max(), 45.0)
pruef("Einfeldtraeger: V_max = qL/2", e.V.max(), 30.0)
pruef("Einfeldtraeger: w = 5qL^4/384EI [mm]", e.w.max(),
      5 * 10 * 6.0 ** 4 / (384 * t.EI) * 1000, 2e-3)
t = Durchlauftraeger(3.0, [Auflager(0.0, EINGESPANNT)])
t.strecke(0, 3, 10.0)
e = t.berechnen()
pruef("Kragarm: M(0) = -qL^2/2", e.M[0], -45.0)
pruef("Kragarm: w = qL^4/8EI [mm]", e.w[-1], 10 * 3.0 ** 4 / (8 * t.EI) * 1000, 3e-3)
t = Durchlauftraeger(12.0, [Auflager(0.0), Auflager(6.0), Auflager(12.0)])
t.strecke(0, 12, 10.0)
e = t.berechnen()
pruef("Zweifeldtraeger: M_Stuetze = -qL^2/8", e.M.min(), -45.0)
pruef("Zweifeldtraeger: M_Feld = 0,0703 qL^2", e.M.max(), 0.0703125 * 10 * 36, 1e-3)
t = Durchlauftraeger(18.0, [Auflager(0.0), Auflager(6.0), Auflager(12.0),
                            Auflager(18.0)])
t.strecke(0, 18, 10.0, art="Q")
e = t.berechnen(1.0, 1.0, lambda x: x <= 12.0 + 1e-9)
pruef("3 Felder, Last in 1+2: M_B = -0,11667 wL^2",
      float(np.interp(6.0, e.x, e.M)), -0.116667 * 10 * 36, 2e-3)
pruef("3 Felder, Last in 1+2: M_C = -0,03333 wL^2",
      float(np.interp(12.0, e.x, e.M)), -0.033333 * 10 * 36, 3e-3)

print("\n6. QUERKRAFT  [EC2 6.2 + NA]")
bw, d = 300., 550.
As_l = stabflaeche_n(4, 20.)
vc = V_Rd_c(C30, bw, d, As_l)
k_h = 1 + math.sqrt(200 / d)
v_h = 0.10 * k_h * (100 * As_l / (bw * d) * 30) ** (1 / 3)
pruef("k = 1 + sqrt(200/d)", vc["k"], k_h)
pruef("C_Rd,c = 0,15/1,5 (NA)", vc["C_Rdc"], 0.10)
pruef("V_Rd,c [kN] Gl. (6.2a)", vc["V_Rdc"], v_h * bw * d / 1e3, 2e-3)
pruef("kappa_1 = 0,0525 (d <= 600)", vc["kappa1"], 0.0525)
pruef("kappa_1 = 0,0375 (d >= 800)", V_Rd_c(C30, bw, 850., As_l)["kappa1"], 0.0375)
z = 495.
ct = cot_theta_NA(C30, bw, z, 300.)
vcc_h = 0.5 * 0.48 * 30 ** (1 / 3) * bw * z / 1e3
pruef("V_Rd,cc [kN] Gl. (6.7bDE)", ct["V_Rdcc"], vcc_h)
pruef("cot(theta) Gl. (6.7aDE)", ct["cot_theta"], 1.2 / (1 - vcc_h / 300.), 1e-4)
pruef("cot(theta) begrenzt auf 3,0", cot_theta_NA(C30, bw, z, 130.)["cot_theta"], 3.0)
Ac_ = 300. * 600.
N_c = -0.20 * C30.fcd * Ac_ / 1e3
pruef("cot(theta) begrenzt auf 1,0 (N_Ed = %.0f kN)" % N_c,
      cot_theta_NA(C30, bw, z, 5000., N_c, Ac_)["cot_theta"], 1.0)
vm = V_Rd_max(C30, bw, z, 2.0)
pruef("nu_1 = 0,75 (<= C50/60)", vm["nu1"], 0.75)
pruef("V_Rd,max [kN] Gl. (6.9)", vm["V_Rdmax"], bw * z * 0.75 * 17.0 / 2.5 / 1e3)
pruef("rho_w,min = 0,16 fctm/fyk", asw_mindest(C30, B500, bw)["rho_w_min"],
      0.16 * C30.fctm / 500.)

print("\n7. TORSION  [EC2 6.3 + NA NDP zu 6.3.2]")
qs_t = Querschnitt(b=400., h=700., d1=40., d2=40.)
ehq = ersatzhohlquerschnitt(qs_t)
A_h, u_h = 400. * 700., 2 * (400. + 700.)
t_h = A_h / u_h
pruef("t_ef = A/u", ehq["t_ef"], t_h)
pruef("A_k = (b-t_ef)(h-t_ef)", ehq["A_k"], (400 - t_h) * (700 - t_h))
pruef("u_k = 2[(b-t_ef)+(h-t_ef)]", ehq["u_k"], 2 * ((400 - t_h) + (700 - t_h)))
pruef("T_Rd,c = fctd t_ef 2 A_k [kNm]", T_Rd_c(C30, ehq)["T_Rdc"],
      C30.fctd * t_h * 2 * (400 - t_h) * (700 - t_h) / 1e6)
cot3 = 3.0
th = math.atan(1 / cot3)
trm = T_Rd_max(C30, ehq, cot3, kasten=False)
pruef("nu = 0,525 nu_2 (Vollquerschnitt)", trm["nu"], 0.525)
pruef("T_Rd,max [kNm] Gl. (6.30)", trm["T_Rdmax"],
      2 * 0.525 * 17.0 * ehq["A_k"] * t_h * math.sin(th) * math.cos(th) / 1e6)
pruef("nu = 0,75 nu_2 (Kastenquerschnitt)",
      T_Rd_max(C30, ehq, cot3, kasten=True)["nu"], 0.75)
T_Ed = 45.0
pruef("asw,T = T_Ed/(2 A_k fywd cot) [mm2/mm]",
      asw_torsion(T_Ed, ehq, B500.fyd, cot3),
      T_Ed * 1e6 / (2 * ehq["A_k"] * B500.fyd * cot3))
pruef("sum(Asl) = T_Ed cot u_k/(2 A_k fyd) [mm2]",
      asl_torsion(T_Ed, ehq, B500.fyd, cot3),
      T_Ed * 1e6 * cot3 * ehq["u_k"] / (2 * ehq["A_k"] * B500.fyd))
# Gl. (6.31aDE): T_Ed <= V_Ed bw / 4,5
entb = bewehrung_entbehrlich(400. * 0.4 / 4.5 / 1e3, 400., 400., 1e6)
pruef("Gl. (6.31aDE) Grenzwert V_Ed bw/4,5 [kNm]", entb["grenze_a"],
      400. * 400. / 4.5 / 1e3)
tor = bemessung_torsion(qs_t, C30, B500, 45.0, 200.0, cot3, 900.0, 300.0)
pruef("Interaktion T/T_Rd,max + V/V_Rd,max", tor.interaktion,
      45.0 / trm["T_Rdmax"] + 200.0 / 900.0, 1e-4)
pruef("s_max Torsion = min(u_k/8 ; min(b,h))", tor.s_max,
      min(ehq["u_k"] / 8.0, 400.0))

print("\n8. GEBRAUCHSTAUGLICHKEIT  [EC2 7.3 / 7.4]")
As, ae = 1885., 15.
q2 = Querschnitt(b=300., h=600., d1=50., d2=50.)
z2 = zustand_II(q2, As, 0.0, ae)
kk = ae * As
x_h = (-kk + math.sqrt(kk ** 2 + 2 * 300. * kk * 550.)) / 300.
pruef("x_II (geschlossene Formel)", z2["x"], x_h, 1e-4)
pruef("I_II (geschlossene Formel)", z2["I"],
      300. * x_h ** 3 / 3 + ae * As * (550. - x_h) ** 2, 1e-3)

print("\n9. PFAHL: KREISQUERSCHNITT  [EC2 6.1 / DIN EN 1536 7.6]")
kq = Kreisquerschnitt(D=900., c_nom=60., phi_l=20., n_l=10, phi_w=10.)
pruef("Ac = pi D^2/4", kq.Ac, math.pi * 900. ** 2 / 4)
pruef("D_s = D - 2(c+phi_w) - phi_l", kq.D_s, 900 - 2 * 70 - 20)
pruef("d_eff = D/2 + D_s/pi", kq.d_eff, 450 + 740 / math.pi)
dg = interaktionsdiagramm(kq, C30, B500, n_punkte=300)
N_druck_h = -(kq.Ac * C30.fcd + kq.As_ges * (B500.sigma_s(C30.eps_c2) - C30.fcd)) / 1e3
pruef("N_Rd zentrischer Druck [kN]", dg["N_druck_max"], N_druck_h, 2e-3)
pruef("N_Rd zentrischer Zug [kN]", dg["N_zug"], kq.As_ges * B500.fyd / 1e3)
pruef("M_Rd(N = N_Zug) = 0", M_Rd_bei_N(dg, dg["N_zug"]), 0.0, 1.0)
pruef("c_nom Bohrpfahl D >= 600 mm", betondeckung_pfahl(900.)["c_nom"], 60.0)
pruef("c_nom Bohrpfahl D < 600 mm", betondeckung_pfahl(500.)["c_nom"], 50.0)
pruef("c_nom unter Stuetzfluessigkeit", betondeckung_pfahl(900., True)["c_nom"], 75.0)
pruef("As,min EN 1536 (0,5 < Ac <= 1,0 m2) [mm2]",
      mindestbewehrung_pfahl(kq)["As_min"], 2500.0)
pruef("As,min EN 1536 (Ac <= 0,5 m2) = 0,5 % Ac",
      mindestbewehrung_pfahl(Kreisquerschnitt(D=600.))["As_min"],
      0.005 * math.pi * 600. ** 2 / 4)
pruef("As,min EN 1536 (Ac > 1,0 m2) = 0,25 % Ac",
      mindestbewehrung_pfahl(Kreisquerschnitt(D=1500.))["As_min"],
      0.0025 * math.pi * 1500. ** 2 / 4)

print("\n10. PFAHL: BETTUNGSMODULVERFAHREN  [EA-Pfaehle 6.3]")
D_m, L_p, EI_p, Es = 0.9, 20.0, 3.0e6, 20000.
ks = bettungsmodul(Es, D_m)
k_b = ks * D_m
lam = (k_b / (4 * EI_p)) ** 0.25
pruef("k_s = E_s/D (D <= 1 m)", ks, Es / D_m)
pruef("k_s = E_s/1 m (D > 1 m)", bettungsmodul(Es, 1.5), Es / 1.0)
sch = [Bodenschicht(0.0, L_p, E_s=Es)]
rh = pfahl_horizontal(L_p, D_m, EI_p, sch, H=200., kopf="frei", n_elem=400)
pruef("w_Kopf = 2 H lambda/k [mm]", rh.w_kopf, 2 * 200 * lam / k_b * 1000, 2e-3)
pruef("M_max = 0,3224 H/lambda [kNm]", rh.M_max, 0.3224 * 200 / lam, 2e-3)
pruef("z(M_max) = pi/(4 lambda) [m]", rh.z_Mmax, math.pi / (4 * lam), 5e-3)
rm = pfahl_horizontal(L_p, D_m, EI_p, sch, M_kopf=300., kopf="frei", n_elem=400)
pruef("w_Kopf bei M0 = 2 M0 lambda^2/k [mm]", rm.w_kopf,
      2 * 300 * lam ** 2 / k_b * 1000, 3e-3)
pruef("M(0) = M_Kopf", rm.M[0], 300.0, 1e-3)
rf = pfahl_horizontal(L_p, D_m, EI_p, sch, H=200., kopf="eingespannt", n_elem=400)
pruef("w_Kopf eingespannt = H lambda/k [mm]", rf.w_kopf,
      200 * lam / k_b * 1000, 2e-3)
pruef("M(0) eingespannt = -H/(2 lambda) [kNm]", rf.M[0], -200 / (2 * lam), 3e-3)

print("\n11. PFAHL: AXIALE TRAGFAEHIGKEIT  [DIN EN 1997-1 7.6 + DIN 1054]")
sch = [Bodenschicht(0.0, 4.0, q_s_k=30.), Bodenschicht(4.0, 10.0, q_s_k=70.),
       Bodenschicht(10.0, 15.0, q_s_k=120.)]
tr = axiale_tragfaehigkeit(0.9, sch, 1800., 2500.)
U, Ab = math.pi * 0.9, math.pi * 0.81 / 4
Rs_h = 30 * U * 4 + 70 * U * 6 + 120 * U * 5
pruef("R_s,k = sum(q_s,k U dl) [kN]", tr.R_s_k, Rs_h)
pruef("R_b,k = q_b,k A_b [kN]", tr.R_b_k, 1800 * Ab)
pruef("R_c,d = R_b,k/1,10 + R_s,k/1,10 [kN]", tr.R_c_d,
      1800 * Ab / 1.10 + Rs_h / 1.10)
pruef("gamma_b (BS-P) = 1,10", GAMMA_R["BS-P"]["gamma_b"], 1.10)
pruef("gamma_s,t Zug (BS-P) = 1,15", GAMMA_R["BS-P"]["gamma_s_t"], 1.15)
pruef("gamma_b (BS-A) = 1,00", GAMMA_R["BS-A"]["gamma_b"], 1.00)
tr_a = axiale_tragfaehigkeit(0.9, sch, 1800., 2500., situation="BS-A")
pruef("R_c,d (BS-A) = R_c,k", tr_a.R_c_d, tr_a.R_c_k)

print("\n12. GESAMTBEMESSUNG (innere Konsistenz)")
e = EingabeBalken(L=7.0, auflager=[(0.0, GELENKIG, .3), (7.0, GELENKIG, .3)],
                  b=300., h=600., g_k=15., q_k=20., phi_laengs=20.,
                  phi_laengs_oben=12., phi_buegel=8.)
b = bemessung_balken(e)
q_Ed = 1.35 * (15. + 25 * 0.18) + 1.5 * 20.
pruef("Balken: M_Ed,max = q_Ed L^2/8", float(np.max(b["einhuellende"]["Mmax"])),
      q_Ed * 49 / 8, 1e-3)
pruef("Balken: V_Ed,max = q_Ed L/2",
      float(np.max(np.abs(b["einhuellende"]["Vmax"]))), q_Ed * 7 / 2, 1e-3)
pruef("Balken: Nachweise erfuellt", 1.0 if b["ok_gesamt"] else 0.0, 1.0)

e2 = EingabeBalken(L=7.0, auflager=[(0.0, GELENKIG, .3), (7.0, GELENKIG, .3)],
                   b=400., h=700., g_k=18., q_k=22., T_Ed=45.0,
                   phi_laengs=20., phi_laengs_oben=16., phi_buegel=10.)
b2 = bemessung_balken(e2)
t2 = b2["torsion"]
pruef("Balken mit Torsion: A_k", t2.A_k, ehq["A_k"])
pruef("Balken mit Torsion: Interaktion <= 1", 1.0 if t2.interaktion <= 1.0 else 0.0, 1.0)
pruef("Balken mit Torsion: Nachweise erfuellt", 1.0 if b2["ok_gesamt"] else 0.0, 1.0)

ep = EingabePfahl(D=900., L=15.0, schichten=sch, q_b_k=1800., N_Ed=-2500.,
                  N_k=-1850., H_Ed=150.)
bp = bemessung_pfahl(ep)
pruef("Pfahl: R_c,d [kN]", bp["tragfaehigkeit"].R_c_d,
      1800 * Ab / 1.10 + Rs_h / 1.10, 1e-3)
pruef("Pfahl: M_Rd >= M_Ed",
      1.0 if bp["M_Rd"] >= bp["M_Ed_max"] else 0.0, 1.0)
pruef("Pfahl: Nachweise erfuellt", 1.0 if bp["ok_gesamt"] else 0.0, 1.0)

print("\n" + "=" * 100)
if FEHLER:
    print("FEHLGESCHLAGEN (%d):" % len(FEHLER))
    for f in FEHLER:
        print("   -", f)
    sys.exit(1)
print("ALLE PRUEFUNGEN BESTANDEN")
print("=" * 100)
