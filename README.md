# BetonPlännen — Bemessung von Stahlbetonbalken und Bohrpfählen nach DIN

Zwei Programme mit **grafischer Oberfläche** und zwei **Jupyter-Notebooks** Schritt für Schritt:

| Bauteil | Nachweise | Programm | Notebook |
|---|---|---|---|
| **Stahlbetonbalken** | Biegung · **Querkraft** · **Torsion** (+ Interaktion), Rissbreite, Verformung | `gui_stahlbetonbalken.py` | `Bemessung_Stahlbetonbalken_DIN.ipynb` |
| **Bohrpfahl** | axiale und horizontale Tragfähigkeit, **Längs- und Querbewehrung** (M-N-Interaktion, Wendel) | `gui_bohrpfahl.py` | `Bemessung_Bohrpfahl_DIN.ipynb` |

Jeder Berichtsblock nennt den **angewendeten DIN-Abschnitt und die Gleichung**.

---

## Normgrundlage

### Balken

| Kurz | Norm | Inhalt |
|---|---|---|
| **EC2** | `DIN EN 1992-1-1:2011-01` | *Bemessung und Konstruktion von Stahlbeton- und Spannbetontragwerken — Teil 1-1* |
| **NA** | `DIN EN 1992-1-1/NA:2013-04` (+A1:2015-12) | *Nationaler Anhang* — national festgelegte Parameter (**NDP**) |
| **EC0** | `DIN EN 1990:2010-12` + `/NA` | Einwirkungskombinationen |
| **EC1** | `DIN EN 1991-1-1:2010-12` + `/NA` | Wichten, Eigengewicht |
| **DIN 488** | `DIN 488-1:2009-08` | Betonstahl B500A / B500B |
| **DIN 1045-2** | `DIN EN 206-1 / DIN 1045-2` | Expositionsklassen |

### Bohrpfahl (zusätzlich)

| Kurz | Norm | Inhalt |
|---|---|---|
| **EN 1536** | `DIN EN 1536:2015-10` | Bohrpfähle: Betondeckung, Mindestbewehrung, Bewehrungskorb |
| **EC7** | `DIN EN 1997-1:2009-09` + `/NA` | Eurocode 7 — Geotechnik |
| **DIN 1054** | `DIN 1054:2010-12` (+A1) | Teilsicherheitsbeiwerte, Nachweis GEO-2 |
| **EA-Pfähle** | Empfehlungen des AK „Pfähle“ (DGGT), 2. Auflage | $q_{b,k}$/$q_{s,k}$, WSL, Bettungsmodulverfahren |

> `DIN 1045-1:2008-08` wurde 2010 zurückgezogen. Der entsprechende Abschnitt ist im Code
> als historische Referenz vermerkt.

---

## Installation und Start

```bash
pip install -r requirements.txt

python gui_stahlbetonbalken.py     # Balken: Biegung, Querkraft, Torsion
python gui_bohrpfahl.py            # Bohrpfahl
python nachweis_pruefung.py        # 96 Prüfungen des Rechenkerns

jupyter notebook Bemessung_Stahlbetonbalken_DIN.ipynb
jupyter notebook Bemessung_Bohrpfahl_DIN.ipynb
```

Voraussetzungen: Python ≥ 3.8, `numpy`, `matplotlib`, `tkinter` (im Standardumfang enthalten).
Für die Notebooks zusätzlich `jupyter`.

---

## Inhalt

```
BetonPlännen/
├── gui_stahlbetonbalken.py                 Oberfläche Balken (Tkinter + matplotlib)
├── gui_bohrpfahl.py                        Oberfläche Bohrpfahl
├── Bemessung_Stahlbetonbalken_DIN.ipynb    Notebook Balken (ausgeführt)
├── Bemessung_Bohrpfahl_DIN.ipynb           Notebook Bohrpfahl (ausgeführt)
├── nachweis_pruefung.py                    96 Prüfungen gegen geschlossene Lösungen
├── notebook_balken_erzeugen.py             erzeugt das Balken-Notebook neu
├── notebook_pfahl_erzeugen.py              erzeugt das Pfahl-Notebook neu
├── requirements.txt
│
├── din_balken/                             Rechenkern Balken
│   ├── normen.py                Verzeichnis der ~70 verwendeten DIN-Normstellen
│   ├── baustoffe.py             EC2 3.1 / 3.2, DIN 488-1, Betondeckung EC2 4.4
│   ├── querschnitt.py           Rechteck und Plattenbalken (b_eff, EC2 5.3.2.1)
│   ├── schnittgroessen.py       Durchlaufträger (Weggrößenverfahren) + Laststellungen
│   ├── biegung.py               GZT Biegung, EC2 6.1
│   ├── querkraft.py             GZT Querkraft, EC2 6.2 + deutsche NDP
│   ├── torsion.py               GZT Torsion, EC2 6.3 + NDP zu 6.3.2
│   ├── gebrauchstauglichkeit.py GZG Rissbreite und Verformung, EC2 7.3 / 7.4
│   ├── konstruktion.py          Verankerung und Mindestbewehrung, EC2 8 / 9.2
│   ├── bemessung.py             Ablaufsteuerung und Bericht
│   └── grafiken.py              10 Abbildungen
│
└── din_pfahl/                              Rechenkern Bohrpfahl
    ├── normen_pfahl.py          Normstellen EN 1536 / EC7 / DIN 1054 / EA-Pfähle
    ├── kreisquerschnitt.py      Kreisquerschnitt, M-N-Interaktion (Punkte A/B/C)
    ├── bettung.py               Bettungsmodulverfahren (Winkler), Knicklast
    ├── tragfaehigkeit.py        axialer Widerstand, Widerstands-Setzungs-Linie
    ├── bemessung_pfahl.py       Ablaufsteuerung, Querkraft am Kreisquerschnitt
    └── grafiken_pfahl.py        7 Abbildungen
```

---

## Umfang — Balken

| Schritt | Inhalt | DIN-Abschnitt |
|---|---|---|
| 1 | Baustoffe: `f_cd = α_cc·f_ck/γ_C`, Parabel-Rechteck | EC2 3.1.6, 3.1.7 · NA NDP zu 3.1.6 (1)P |
| 2 | Dauerhaftigkeit: `c_nom = c_min + Δc_dev` | EC2 4.4.1 Gl. (4.1)/(4.2) · NA Tab. 4.4DE |
| 3 | Einwirkungen: `1,35 G + 1,50 Q` / `G + ψ₂ Q` | EC0 Gl. (6.10) und (6.16b) |
| 4 | Schnittgrößen: Durchlaufträger, feldweise Laststellung | EC2 5.1.3 |
| 5 | **Biegung GZT**: `μ_Eds → ξ, ζ, ω → A_s1`, Druckbewehrung | EC2 6.1 · NA NDP zu 5.5 (4) |
| 6 | Mindest-/Höchstbewehrung, Robustheitsbewehrung | EC2 9.2.1.1 Gl. (9.1N) · NA |
| 7 | **Querkraft GZT**: `V_Rd,c`, `cot θ`, `V_Rd,max`, `A_sw/s`, `s_max` | EC2 6.2 · NA Gl. (6.3aDE), (6.7aDE), (6.7bDE), Tab. NA.9.1 |
| 8 | **Torsion GZT**: Ersatzhohlquerschnitt, `T_Rd,c`, `T_Rd,max`, Interaktion mit V, Längs- und Bügelbewehrung | EC2 6.3 Gl. (6.26)–(6.31) · NA NDP zu 6.3.2 · EC2 9.2.3 |
| 9 | **Rissbreite GZG**: `w_k = s_r,max(ε_sm−ε_cm)` | EC2 7.3.4 Gl. (7.8)/(7.9) · NA Gl. (7.11DE) |
| 10 | **Verformung GZG**: Interpolation Zustand I/II | EC2 7.4.2 Gl. (7.16) und 7.4.3 Gl. (7.18)/(7.19) |
| 11 | Verankerung, Übergreifung, *Versatzmaß* `a_l = z·cot θ/2` | EC2 8.4, 9.2.1.3 Gl. (9.2), 9.2.1.4 Gl. (9.3) |

**Querschnitte**: Rechteck und Plattenbalken. **Systeme**: Einfeld-, Kragarm-, Durchlaufträger
mit *n* Feldern, Einspannungen. **Beanspruchung**: Biegung mit Längskraft (`N_Ed`, Druck negativ),
Querkraft und Torsion einschließlich Interaktion.

## Umfang — Bohrpfahl

| Schritt | Inhalt | Normabschnitt |
|---|---|---|
| 1 | Baustoffe, Betondeckung `c_nom ≥ 60/50/75 mm` | DIN EN 1536 6.3 / 7.6.2 |
| 2 | **Axiale Tragfähigkeit** `R_c,d = R_b,k/γ_b + R_s,k/γ_s` | DIN EN 1997-1 7.6.2.3 Gl. (7.8) · DIN 1054 Tab. A 2.3 |
| 3 | Widerstands-Setzungs-Linie, `s_sg` und `s_g = 0,10 D` | EA-Pfähle 5.4.5 |
| 4 | **Horizontallast**: Bettungsmodulverfahren, `k_s = E_s/D` | EA-Pfähle 6.3 |
| 5 | **Längsbewehrung**: M-N-Interaktion des Kreisquerschnitts (Punkte A/B/C) | EC2 6.1 (2)P Bild 6.1, 6.1 (5) |
| 6 | Mindestbewehrung `0,5 % / 25 cm² / 0,25 % A_c` und `0,15|N_Ed|/f_yd` | DIN EN 1536 7.6.3 Tab. 4 · NA zu EC2 NDP zu 9.5.2 (2) |
| 7 | **Querbewehrung (Wendel)**: `V_Rd,c`, `cot θ`, Ganghöhe | EC2 6.2 + 9.5.3 · DIN EN 1536 7.6.4 |
| 8 | Knicknachweis `N_ki = 2√(EI·k)` | EA-Pfähle 4.7 · EC2 5.8 |

---

## Wesentliche Abweichungen des deutschen NA vom Eurocode

Sie sind der Grund, warum eine „allgemeine Eurocode-2-Rechnung“ **keine** deutsche Bemessung ist:

| Größe | EC2 (empfohlen) | **DIN EN 1992-1-1/NA** | Abschnitt |
|---|---|---|---|
| `α_cc` | 1,00 | **0,85** | NDP zu 3.1.6 (1)P |
| `ε_ud` | 0,9·ε_uk | **25 ‰** | NDP zu 3.2.7 (2) |
| `k₁`, `k₂` (Umlagerung) | 0,44 / 1,25 | **0,64 / 0,80** → `x_u/d ≤ 0,45` | NDP zu 5.5 (4) |
| `C_Rd,c` | 0,18/γ_C = 0,12 | **0,15/γ_C = 0,10** | NDP zu 6.2.2 (1) |
| `v_min` | 0,035·k^1,5·√f_ck | **Gl. (6.3aDE)** | NDP zu 6.2.2 (1) |
| `cot θ` | 1,0…2,5 frei | **Formel Gl. (6.7aDE)**, 1,0…3,0 | NDP zu 6.2.3 (2) |
| `ν₁` (Querkraft) | 0,6(1−f_ck/250) | **0,75·ν₂** | NDP zu 6.2.3 (3) |
| **`ν` (Torsion)** | 0,6(1−f_ck/250) | **0,525·ν₂ (voll) / 0,75·ν₂ (Kasten)** | NDP zu 6.3.2 (4) |
| **Verzicht Torsionsbew.** | Gl. (6.31) | **Gl. (6.31aDE)/(6.31bDE)** | NDP zu 6.3.2 (5) |
| `z` (Querkraft) | 0,9d | **min{0,9d; d−2c_v,l; d−c_v,l−30}** | NDP zu 6.2.3 (1) |
| `s_r,max` | k₃c + k₁k₂k₄⌀/ρ_p,eff | **⌀/(3,6·ρ_p,eff)**, Gl. (7.11DE) | NDP zu 7.3.4 (3) |
| `ρ_w,min` | 0,08√f_ck/f_yk | **0,16·f_ctm/f_yk** | NDP zu 9.2.2 (5) |
| `Δc_dev` | 10 mm | **15 mm** (10 mm bei XC1) | NDP zu 4.4.1.3 (1)P |

---

## Abbildungen

**Balken** (10): Bewehrter Querschnitt · σ-ε-Linien der Baustoffe · Statisches System und
Einhüllende von `M`/`V` · Dehnungen und Spannungen im Bruchzustand · Allgemeines
Bemessungsdiagramm `μ–ω–ξ–ζ` · Querkraft mit Bügelbereichen · **Torsion: Ersatzhohlquerschnitt
und Interaktionsdiagramm T–V** · Zugkraftdeckungslinie · Durchbiegung · Ausnutzungsgrade

**Bohrpfahl** (7): Längsschnitt mit Baugrund und Widerstandsanteilen · Querschnitt mit
Bewehrungskorb · Anteile am Pfahlwiderstand · Widerstands-Setzungs-Linie · `w(z)`, `M(z)`, `V(z)`
aus dem Bettungsmodulverfahren · **M-N-Interaktionsdiagramm** · Ausnutzungsgrade

Beide Oberflächen exportieren den Bericht als `.txt` und alle Abbildungen als `.pdf`.

---

## Prüfung des Rechenkerns

`python nachweis_pruefung.py` vergleicht **96 Größen** mit geschlossenen Lösungen:

* `α_R = 0,80952` und `k_a = 0,41597` (exakte Integration der Parabel-Rechteck-Druckzone)
* `μ_lim = 0,296`, `ω_lim = 0,3643`, `ζ_lim = 0,8128` bei `x_u/d = 0,45`
* Hin- und Rückrechnung `M_Ed → A_s → M_Rd` für Rechteck, Plattenbalken und mit Längskraft
* Durchlaufträger gegen die **Dreimomentengleichung**, einschließlich feldweiser Laststellung
* `V_Rd,c`, `V_Rd,cc`, `cot θ`, `V_Rd,max`, `ρ_w,min` gegen die NA-Formeln
* **Torsion**: `t_ef`, `A_k`, `u_k`, `T_Rd,c`, `T_Rd,max`, `a_sw,T`, `ΣA_sl`, Interaktion
* Nulllinie und Trägheitsmoment im Zustand II gegen die geschlossene Formel
* **Pfahl**: zentrischer Druck/Zug des Kreisquerschnitts, `c_nom` und `A_s,min` nach EN 1536
* **Bettungsmodulverfahren** gegen die analytische Lösung des elastisch gebetteten Balkens
  (`w(0) = 2Hλ/k`, `M_max = 0,3224 H/λ` bei `z = π/(4λ)`, eingespannter Kopf)
* `R_c,d` und die Teilsicherheitsbeiwerte nach DIN 1054, Tab. A 2.3

---

## ⚠️ Hinweis

Berechnungs- und **Lehrwerkzeug**. Alle Ergebnisse sind von einem verantwortlichen Ingenieur
gegen die **gültige Ausgabe** der jeweiligen Norm zu prüfen.

**Balken** — nicht erfasst: Durchstanzen (EC2 6.4), Kippen schlanker Träger (5.9),
Ermüdung (6.8), Vorspannung (5.10), außergewöhnliche und seismische Bemessungssituationen,
Brandfall (EN 1992-1-2), Leichtbeton (Kap. 11).

**Bohrpfahl** — `q_b,k` und `q_s,k` sind **Eingabewerte des Anwenders** (EA-Pfähle
Tab. 5.12–5.15 oder Probebelastung); das Programm enthält bewusst keine Erfahrungswerte-Tabellen.
Nicht erfasst: Pfahlgruppenwirkung (EA-Pfähle 8), negative Mantelreibung, zyklische und
dynamische Einwirkungen, Erdbeben. Die Bettung ist **linear** angesetzt — bei großen
Kopfverschiebungen ist ein nichtlineares p-y-Verfahren zu verwenden (EA-Pfähle 6.4).

Die Zahlenwerte der NA-Tabellen (Betondeckung Tab. 4.4DE, `w_max` Tab. 7.1DE, Bügelabstände
Tab. NA.9.1) sowie die NDP zur Torsion (`ν = 0,525/0,75`, Gl. (6.31aDE)/(6.31bDE)) sind gegen
die für Sie geltende Normausgabe zu prüfen.
