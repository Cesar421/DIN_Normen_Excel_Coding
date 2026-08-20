# -*- coding: utf-8 -*-
"""
GRENZZUSTAENDE DER GEBRAUCHSTAUGLICHKEIT (GZG).

Norm: DIN EN 1992-1-1:2011-01, Abschnitt 7  + DIN EN 1992-1-1/NA:2013-04

  7.3.2  Mindestbewehrung zur Rissbreitenbegrenzung   Gl. (7.1)
  7.3.4  Rechnerische Rissbreite wk                   Gl. (7.8)/(7.9) + NA (7.11DE)
  7.4.2  Begrenzung ueber die Biegeschlankheit l/d    Gl. (7.16a)/(7.16b) + NA
  7.4.3  Verformungsberechnung (Zustand I/II)         Gl. (7.18)/(7.19)

Einwirkungskombination: quasi-staendig  G + psi_2 Q  (DIN EN 1990, Gl. 6.16b).
"""

import math

import numpy as np

from .normen import ref


# ---------------------------------------------------------------------------
# Zustand I (ungerissen) und Zustand II (gerissen)
# ---------------------------------------------------------------------------
def _integrale(querschnitt, x, n=200):
    """Flaeche, statisches Moment und Traegheitsmoment der Zone 0..x."""
    ys = set(np.linspace(0.0, max(x, 1e-9), n + 1).tolist())
    if 0.0 < querschnitt.hf < x:
        ys.update((querschnitt.hf - 1e-9, querschnitt.hf + 1e-9))
    ys = np.array(sorted(ys))
    bs = np.array([querschnitt.breite(y) for y in ys])
    return (np.trapz(bs, ys), np.trapz(bs * (x - ys), ys),
            np.trapz(bs * (x - ys) ** 2, ys))


def zustand_I(querschnitt, As1, As2, alpha_e):
    """
    Ideeller ungerissener Querschnitt (Zustand I).
    Liefert x_I (ab oberem Rand), I_I und das Widerstandsmoment unten.
    """
    d, d2, h = querschnitt.d, querschnitt.d2, querschnitt.h
    A_ges, S_ges, _ = _integrale(querschnitt, h, 400)
    Sc = A_ges * h - S_ges                       # statisches Moment um oben
    A = A_ges + (alpha_e - 1.0) * (As1 + As2)
    S = Sc + (alpha_e - 1.0) * (As1 * d + As2 * d2)
    x = S / A
    ys = set(np.linspace(0.0, h, 401).tolist())
    if 0.0 < querschnitt.hf < h:
        ys.update((querschnitt.hf - 1e-9, querschnitt.hf + 1e-9))
    ys = np.array(sorted(ys))
    bs = np.array([querschnitt.breite(y) for y in ys])
    I = np.trapz(bs * (ys - x) ** 2, ys)
    I += (alpha_e - 1.0) * (As1 * (d - x) ** 2 + As2 * (x - d2) ** 2)
    return dict(x=x, I=I, A=A, W_unten=I / max(h - x, 1e-9))


def zustand_II(querschnitt, As1, As2, alpha_e, tol=1e-9):
    """
    Gerissener ideeller Querschnitt (Zustand II).
    Gleichgewicht der statischen Momente um die Nulllinie:

        int_0^x b(y)(x-y) dy + (alpha_e-1) As2 (x-d2) - alpha_e As1 (d-x) = 0
    """
    d, d2, h = querschnitt.d, querschnitt.d2, querschnitt.h

    def f(x):
        _, S, _ = _integrale(querschnitt, x)
        return S + (alpha_e - 1.0) * As2 * (x - d2) - alpha_e * As1 * (d - x)

    lo, hi = 1e-6, h
    if f(hi) < 0:
        return dict(x=h, I=float("nan"), ok=False)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol * h:
            break
    x = 0.5 * (lo + hi)
    _, _, I = _integrale(querschnitt, x, 400)
    I += (alpha_e - 1.0) * As2 * (x - d2) ** 2 + alpha_e * As1 * (d - x) ** 2
    return dict(x=x, I=I, ok=True)


def stahlspannung(querschnitt, beton, stahl, M_kNm, As1, As2=0.0, phi_kriech=0.0):
    """
    Stahlspannung im Zustand II unter dem Moment M [kNm]:
        sigma_s = alpha_e M (d - x) / I_II
    """
    Ece = beton.Ec_eff(phi_kriech)
    alpha_e = stahl.Es / Ece
    z2 = zustand_II(querschnitt, As1, As2, alpha_e)
    M = abs(M_kNm) * 1.0e6
    sigma_s = alpha_e * M * (querschnitt.d - z2["x"]) / z2["I"] if z2["I"] > 0 else 0.0
    return dict(sigma_s=sigma_s, sigma_c=M * z2["x"] / z2["I"] if z2["I"] > 0 else 0.0,
                x_II=z2["x"], I_II=z2["I"], alpha_e=alpha_e, Ec_eff=Ece)


# ---------------------------------------------------------------------------
# 7.3.4  Rissbreite
# ---------------------------------------------------------------------------
def rissbreite(querschnitt, beton, stahl, M_qs_kNm, As1, phi_stab,
               As2=0.0, phi_kriech=2.0, w_max=0.30, langzeit=True,
               fct_eff=None):
    """
    Rechnerische Rissbreite.
    DIN EN 1992-1-1, 7.3.4, Gl. (7.8)/(7.9) und DIN EN 1992-1-1/NA Gl. (7.11DE).

        wk            = sr,max (eps_sm - eps_cm)                       Gl. (7.8)
        eps_sm-eps_cm = [sigma_s - kt fct,eff/rho_p,eff (1+alpha_e rho_p,eff)]/Es
                        >= 0,6 sigma_s/Es                              Gl. (7.9)
        sr,max        = phi/(3,6 rho_p,eff) <= sigma_s phi/(3,6 fct,eff)  (7.11DE)
        Ac,eff        = b hc,ef ; hc,ef = min(2,5(h-d) ; (h-x)/3 ; h/2)  7.3.2 (3)

    kt = 0,4 (Langzeiteinwirkung) / 0,6 (Kurzzeiteinwirkung).
    """
    st = stahlspannung(querschnitt, beton, stahl, M_qs_kNm, As1, As2, phi_kriech)
    sigma_s, x = st["sigma_s"], st["x_II"]
    h, d, b = querschnitt.h, querschnitt.d, querschnitt.bw
    fcte = beton.fctm if fct_eff is None else fct_eff
    kt = 0.4 if langzeit else 0.6

    hc_ef = min(2.5 * (h - d), (h - x) / 3.0, h / 2.0)
    Ac_eff = b * hc_ef
    rho_eff = As1 / Ac_eff if Ac_eff > 0 else 0.0
    alpha_e = st["alpha_e"]

    if rho_eff <= 0 or sigma_s <= 0:
        return dict(wk=0.0, w_max=w_max, ok=True, sigma_s=sigma_s,
                    hinweis="Keine massgebende Zugbeanspruchung")

    de = (sigma_s - kt * fcte / rho_eff * (1.0 + alpha_e * rho_eff)) / stahl.Es
    de_min = 0.6 * sigma_s / stahl.Es
    min_massgebend = de < de_min
    de = max(de, de_min)

    sr_a = phi_stab / (3.6 * rho_eff)
    sr_b = sigma_s * phi_stab / (3.6 * fcte)
    sr_max = min(sr_a, sr_b)
    wk = sr_max * de

    return dict(wk=wk, w_max=w_max, ok=wk <= w_max + 1e-9,
                ausnutzung=wk / w_max if w_max > 0 else 0.0,
                sigma_s=sigma_s, sigma_c=st["sigma_c"], x_II=x, I_II=st["I_II"],
                alpha_e=alpha_e, hc_ef=hc_ef, Ac_eff=Ac_eff, rho_p_eff=rho_eff,
                eps_sm_cm=de, min_massgebend=min_massgebend, kt=kt,
                fct_eff=fcte, sr_max=sr_max, sr_a=sr_a, sr_b=sr_b, phi=phi_stab,
                M_qs=M_qs_kNm,
                normen=[ref("wk"), ref("eps_sm"), ref("sr_max"), ref("hc_eff"),
                        ref("w_max"), ref("komb_QS")])


def mindestbewehrung_riss(querschnitt, beton, stahl, sigma_s=None,
                          frueher_zwang=False, biegung=True):
    """
    Mindestbewehrung zur Rissbreitenbegrenzung.
    DIN EN 1992-1-1, 7.3.2 (2), Gl. (7.1):

        As,min sigma_s = kc k fct,eff Act

    kc = 0,4 (reine Biegung, Rechteck) ; 1,0 (zentrischer Zug)
    k  = 1,0 (h <= 300 mm) ... 0,65 (h >= 800 mm), linear interpoliert
    Act = gezogene Betonflaeche unmittelbar vor der Rissbildung

    Ohne Angabe von sigma_s wird fyk angesetzt (auf der sicheren Seite). Ein
    kleinerer Wert nach Tab. 7.2DE/7.3DE des NA (Begrenzung des Stabdurchmessers
    bzw. des Stababstandes) fuehrt zu groesserer Bewehrung, erlaubt aber den
    Nachweis konkreter Durchmesser; dieser Weg ist hier NICHT implementiert.
    """
    h = querschnitt.h
    fcte = beton.fct_eff(frueher_zwang)
    if h <= 300.0:
        k = 1.0
    elif h >= 800.0:
        k = 0.65
    else:
        k = 1.0 + (0.65 - 1.0) * (h - 300.0) / 500.0
    kc = 0.4 if biegung else 1.0
    ss = stahl.fyk if sigma_s is None else sigma_s
    Act = querschnitt.bw * (h - querschnitt.y_schwerpunkt()) if biegung \
        else querschnitt.Ac
    return dict(As_min=kc * k * fcte * Act / ss, kc=kc, k=k, fct_eff=fcte,
                Act=Act, sigma_s=ss, normen=[ref("As_min_riss")])


# ---------------------------------------------------------------------------
# 7.4  Verformungen
# ---------------------------------------------------------------------------
def zulaessige_schlankheit(beton, rho, rho_druck=0.0, K=1.0, l_m=5.0,
                           sigma_s=310.0, empfindlich=False):
    """
    Begrenzung der Verformung ueber die Biegeschlankheit.
    DIN EN 1992-1-1, 7.4.2 (2), Gl. (7.16a)/(7.16b):

        rho <= rho0:  l/d = K [11 + 1,5 sqrt(fck) rho0/rho
                                 + 3,2 sqrt(fck) (rho0/rho - 1)^1,5]
        rho >  rho0:  l/d = K [11 + 1,5 sqrt(fck) rho0/(rho-rho')
                                 + (1/12) sqrt(fck) sqrt(rho'/rho0)]
        rho0 = 1e-3 sqrt(fck)

    Korrektur ueber die Stahlspannung: Faktor 310/sigma_s  (7.4.2 (2)).
    Zusaetzlich die vereinfachte Regel des NA (NDP zu 7.4.2):
        l/d <= 35   ;  bei erhoehten Anforderungen: l/d <= 35 * 7/l_eff
    """
    fck = beton.fck
    rho0 = 1.0e-3 * math.sqrt(fck)
    rho = max(rho, 1e-6)
    sf = math.sqrt(fck)
    if rho <= rho0:
        ld = K * (11.0 + 1.5 * sf * rho0 / rho
                  + 3.2 * sf * (rho0 / rho - 1.0) ** 1.5)
        gl = "Gl. (7.16a)"
    else:
        ld = K * (11.0 + 1.5 * sf * rho0 / max(rho - rho_druck, 1e-6)
                  + (1.0 / 12.0) * sf * math.sqrt(max(rho_druck, 0.0) / rho0))
        gl = "Gl. (7.16b)"
    ld_korr = min(ld * (310.0 / sigma_s) if sigma_s > 0 else ld, 1.5 * ld)

    ld_NA = 35.0
    if empfindlich and l_m > 7.0:
        ld_NA = 35.0 * 7.0 / l_m
    return dict(ld_7_16=ld, ld_korrigiert=ld_korr, ld_NA=ld_NA, rho0=rho0,
                rho=rho, K=K, gleichung=gl,
                normen=[ref("durchbiegung_ld"), ref("durchbiegung_NA")])


def kruemmung(querschnitt, beton, stahl, M_kNm, As1, As2, phi_kriech,
              beta=0.5, eps_cs=0.0):
    """
    Mittlere Kruemmung 1/r [1/mm] nach DIN EN 1992-1-1, 7.4.3 (3), Gl. (7.18):

        alpha = zeta alpha_II + (1 - zeta) alpha_I
        zeta  = 1 - beta (sigma_sr/sigma_s)^2 = 1 - beta (Mcr/M)^2   Gl. (7.19)
        beta = 0,5 bei Langzeit- oder wiederholter Einwirkung

    Enthaelt die Schwindkruemmung (7.4.3 (6)):  1/r_cs = eps_cs alpha_e S / I
    """
    M = abs(M_kNm) * 1.0e6
    Ece = beton.Ec_eff(phi_kriech)
    ae = stahl.Es / Ece
    z1 = zustand_I(querschnitt, As1, As2, ae)
    z2 = zustand_II(querschnitt, As1, As2, ae)
    Mcr = beton.fctm * z1["W_unten"]

    k1 = M / (Ece * z1["I"]) if z1["I"] > 0 else 0.0
    k2 = M / (Ece * z2["I"]) if z2["I"] > 0 else 0.0
    zeta = 0.0 if (M <= Mcr or M <= 0) else 1.0 - beta * (Mcr / M) ** 2
    zeta = min(max(zeta, 0.0), 1.0)

    def _S(z):
        return As1 * (querschnitt.d - z["x"]) - As2 * (z["x"] - querschnitt.d2)

    kcs1 = eps_cs * ae * _S(z1) / z1["I"] if z1["I"] > 0 else 0.0
    kcs2 = eps_cs * ae * _S(z2) / z2["I"] if z2["I"] > 0 else 0.0

    return dict(kappa=zeta * (k2 + kcs2) + (1.0 - zeta) * (k1 + kcs1),
                zeta=zeta, Mcr=Mcr / 1.0e6, kappa_I=k1, kappa_II=k2,
                x_I=z1["x"], I_I=z1["I"], x_II=z2["x"], I_II=z2["I"],
                alpha_e=ae, Ec_eff=Ece, kappa_cs=zeta * kcs2 + (1 - zeta) * kcs1,
                normen=[ref("durchbiegung_rech"), ref("kriechen")])


def durchbiegung_feld(x_m, M_kNm, querschnitt, beton, stahl, As1, As2,
                      phi_kriech, x_anf, x_end, beta=0.5, eps_cs=0.0, n=101,
                      kragarm=False):
    """
    Durchbiegung eines Abschnitts durch Integration der wirklichen Kruemmung
    (DIN EN 1992-1-1, 7.4.3).

    Feld zwischen zwei Auflagern (Durchbiegung bezogen auf die Sehne):
        w(x) = (x/L) int_0^L (L-t) k(t) dt - int_0^x (x-t) k(t) dt

    Kragarm mit Einspannung bei x_anf (bezogen auf die Tangente dort):
        w(x) = int_0^x (x-t) k(t) dt
    """
    L = (x_end - x_anf) * 1000.0
    if L <= 0:
        return dict(w_max=0.0, x=np.array([]), w=np.array([]))
    xs = np.linspace(x_anf, x_end, n)
    Ms = np.interp(xs, x_m, M_kNm)
    ks = np.array([kruemmung(querschnitt, beton, stahl, m, As1, As2, phi_kriech,
                             beta, eps_cs)["kappa"] * (1.0 if m >= 0 else -1.0)
                   for m in Ms])
    t = (xs - x_anf) * 1000.0
    I1 = np.trapz((L - t) * ks, t)
    w = np.zeros(n)
    for i in range(n):
        anteil = np.trapz((t[i] - t[:i + 1]) * ks[:i + 1], t[:i + 1])
        w[i] = -anteil if kragarm else (t[i] / L) * I1 - anteil
    j = int(np.argmax(np.abs(w)))
    return dict(w_max=float(w[j]), x_max=float(xs[j]), x=xs, w=w, L=L,
                normen=[ref("durchbiegung_rech")])


def nachweis_durchbiegung(x_m, M_qs_kNm, querschnitt, beton, stahl, As1, As2,
                          felder, phi_kriech=2.0, eps_cs=0.0, grenze=250.0,
                          grenze_2=500.0, K=1.0, sigma_s=310.0, kragarme=()):
    """
    Vollstaendiger Verformungsnachweis.

    Uebliche Grenzwerte (DIN EN 1992-1-1, 7.4.1 (4)/(5)):
        w <= l/250   (Gesamtdurchbiegung, quasi-staendige Kombination)
        w <= l/500   (Durchbiegung nach Einbau verformungsempfindlicher Bauteile)

    `kragarme` : Indizes in `felder`, die Kragarme mit Einspannung am Anfang
    sind. Der Grenzwert wird mit der tatsaechlichen Kraglaenge gebildet
    (konservativ; in der Literatur wird teils 2 l angesetzt).
    """
    erg = []
    for iv, (a, b) in enumerate(felder):
        krag = iv in kragarme
        f = durchbiegung_feld(x_m, M_qs_kNm, querschnitt, beton, stahl, As1, As2,
                              phi_kriech, a, b, eps_cs=eps_cs, kragarm=krag)
        L = (b - a) * 1000.0
        rho = As1 / (querschnitt.bw * querschnitt.d)
        sch = zulaessige_schlankheit(beton, rho,
                                     As2 / (querschnitt.bw * querschnitt.d),
                                     K=K, l_m=b - a, sigma_s=sigma_s)
        erg.append(dict(feld=(a, b), L=L, w_max=abs(f["w_max"]),
                        x_max=f.get("x_max", 0.0), w_grenz=L / grenze,
                        ok=abs(f["w_max"]) <= L / grenze,
                        ausnutzung=abs(f["w_max"]) / (L / grenze),
                        w_grenz2=L / grenze_2, x=f["x"], w=f["w"],
                        ld_vorh=L / querschnitt.d, ld_zul=sch["ld_korrigiert"],
                        ld_NA=sch["ld_NA"], schlankheit=sch, kragarm=krag,
                        normen=[ref("durchbiegung_rech"), ref("durchbiegung_ld"),
                                ref("durchbiegung_NA"), ref("komb_QS")]))
    return erg
