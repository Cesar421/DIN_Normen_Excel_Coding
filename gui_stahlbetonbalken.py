# -*- coding: utf-8 -*-
"""
===============================================================================
 BEMESSUNG VON STAHLBETONBALKEN NACH DIN
 Grafische Benutzeroberflaeche (Tkinter + matplotlib)
===============================================================================
 Normgrundlage
   DIN EN 1992-1-1:2011-01        Eurocode 2, Teil 1-1
   DIN EN 1992-1-1/NA:2013-04     Nationaler Anhang (deutsche NDP)
   DIN EN 1990 / DIN EN 1991-1-1  Einwirkungen und Kombinationen (+ NA)
   DIN 488-1:2009-08              Betonstahl B500A / B500B

 Nachgewiesen werden BIEGUNG (6.1), QUERKRAFT (6.2) und TORSION (6.3)
 einschliesslich der Interaktion, sowie Rissbreite (7.3) und Durchbiegung (7.4).

 Jeder Berichtsblock nennt den angewendeten DIN-ABSCHNITT und die GLEICHUNG.

 Start:   python gui_stahlbetonbalken.py
===============================================================================
"""

import os
import re
import sys
import traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from din_balken.baustoffe import BETONKLASSEN, DURCHMESSER, EXPOSITION
from din_balken.schnittgroessen import GELENKIG, EINGESPANNT
from din_balken.bemessung import EingabeBalken, bemessung_balken, bericht_text
from din_balken import grafiken as G
from din_balken.normen import normentabelle

# ---------------------------------------------------------------------------
BG = "#f4f4f1"
CARD = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
ACC = "#2a78d6"
OKC = "#0ca30c"
NOKC = "#d03b3b"

SYSTEME = {
    "Einfeldtraeger":            ("0; L", "g; g"),
    "Kragarm":                   ("0", "e"),
    "Zweifeldtraeger":           ("0; L/2; L", "g; g; g"),
    "Dreifeldtraeger":           ("0; L/3; 2L/3; L", "g; g; g; g"),
    "Eingespannt - gelenkig":    ("0; L", "e; g"),
    "Beidseitig eingespannt":    ("0; L", "e; e"),
    "Benutzerdefiniert":         ("0; L", "g; g"),
}


def _hilfe(widget, text):
    tip = {"win": None}

    def rein(_e):
        if tip["win"]:
            return
        w = tk.Toplevel(widget)
        w.wm_overrideredirect(True)
        w.wm_geometry("+%d+%d" % (widget.winfo_rootx() + 18,
                                  widget.winfo_rooty() + 18))
        tk.Label(w, text=text, justify="left", bg="#ffffe0", fg=INK,
                 relief="solid", bd=1, font=("Segoe UI", 8),
                 wraplength=400).pack()
        tip["win"] = w

    def raus(_e):
        if tip["win"]:
            tip["win"].destroy()
            tip["win"] = None

    widget.bind("<Enter>", rein)
    widget.bind("<Leave>", raus)


class Feld:
    """Zeile aus Beschriftung und Eingabeelement."""

    def __init__(self, parent, zeile, text, wert="", breite=14, einheit="",
                 auswahl=None, hilfe="", typ="entry"):
        ttk.Label(parent, text=text, style="Feld.TLabel").grid(
            row=zeile, column=0, sticky="w", padx=(8, 4), pady=2)
        self.var = tk.StringVar(value=str(wert))
        if typ == "combo":
            self.w = ttk.Combobox(parent, textvariable=self.var, width=breite - 2,
                                  values=auswahl, state="readonly")
        elif typ == "check":
            self.var = tk.BooleanVar(value=bool(wert))
            self.w = ttk.Checkbutton(parent, variable=self.var)
        else:
            self.w = ttk.Entry(parent, textvariable=self.var, width=breite)
        self.w.grid(row=zeile, column=1, sticky="w", pady=2)
        if einheit:
            ttk.Label(parent, text=einheit, style="Einheit.TLabel").grid(
                row=zeile, column=2, sticky="w", padx=(4, 8))
        if hilfe:
            lab = ttk.Label(parent, text="?", style="Hilfe.TLabel", cursor="hand2")
            lab.grid(row=zeile, column=3, sticky="w")
            _hilfe(lab, hilfe)

    def get(self):
        return self.var.get()

    def f(self, standard=0.0):
        try:
            return float(str(self.var.get()).replace(",", "."))
        except (TypeError, ValueError):
            return standard

    def i(self, standard=0):
        try:
            return int(float(str(self.var.get()).replace(",", ".")))
        except (TypeError, ValueError):
            return standard

    def b(self):
        return bool(self.var.get())


class App(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("Bemessung von Stahlbetonbalken - DIN EN 1992-1-1 + "
                   "Nationaler Anhang  |  Biegung, Querkraft, Torsion")
        self.geometry("1580x960")
        self.minsize(1200, 740)
        self.configure(bg=BG)
        self.bericht = None
        self.canvases = []
        self._stile()
        self._aufbau()
        self.after(250, self.berechnen)

    # -- Stile -----------------------------------------------------------
    def _stile(self):
        s = ttk.Style(self)
        try:
            s.theme_use("clam")
        except tk.TclError:
            pass
        s.configure(".", background=BG, foreground=INK, font=("Segoe UI", 9))
        s.configure("TFrame", background=BG)
        s.configure("TLabelframe", background=BG, foreground=ACC, borderwidth=1,
                    relief="solid")
        s.configure("TLabelframe.Label", background=BG, foreground=ACC,
                    font=("Segoe UI", 9, "bold"))
        s.configure("Feld.TLabel", background=BG, foreground=INK2)
        s.configure("Einheit.TLabel", background=BG, foreground="#898781",
                    font=("Segoe UI", 8))
        s.configure("Hilfe.TLabel", background=BG, foreground=ACC,
                    font=("Segoe UI", 8, "bold"))
        s.configure("Titel.TLabel", background=BG, foreground=INK,
                    font=("Segoe UI", 13, "bold"))
        s.configure("Unter.TLabel", background=BG, foreground=INK2,
                    font=("Segoe UI", 8))
        s.configure("Aktion.TButton", font=("Segoe UI", 10, "bold"),
                    padding=(10, 7))
        s.configure("TNotebook", background=BG)
        s.configure("TNotebook.Tab", padding=(11, 6), font=("Segoe UI", 9))

    # -- Aufbau ----------------------------------------------------------
    def _aufbau(self):
        kopf = ttk.Frame(self)
        kopf.pack(fill="x", padx=10, pady=(8, 4))
        ttk.Label(kopf, text="Bemessung von Stahlbetonbalken",
                  style="Titel.TLabel").pack(anchor="w")
        ttk.Label(kopf, text="DIN EN 1992-1-1:2011-01 (Eurocode 2) + "
                             "DIN EN 1992-1-1/NA:2013-04  |  DIN EN 1990 / "
                             "DIN EN 1991-1-1  |  DIN 488-1 (B500)  -  "
                             "Biegung, Querkraft und Torsion",
                  style="Unter.TLabel").pack(anchor="w")

        koerper = ttk.Frame(self)
        koerper.pack(fill="both", expand=True, padx=10, pady=6)

        links = ttk.Frame(koerper, width=440)
        links.pack(side="left", fill="y")
        links.pack_propagate(False)
        cnv = tk.Canvas(links, bg=BG, highlightthickness=0, width=420)
        scr = ttk.Scrollbar(links, orient="vertical", command=cnv.yview)
        self.form = ttk.Frame(cnv)
        self.form.bind("<Configure>",
                       lambda e: cnv.configure(scrollregion=cnv.bbox("all")))
        cnv.create_window((0, 0), window=self.form, anchor="nw")
        cnv.configure(yscrollcommand=scr.set)
        cnv.pack(side="left", fill="both", expand=True)
        scr.pack(side="right", fill="y")
        cnv.bind_all("<MouseWheel>",
                     lambda e: cnv.yview_scroll(int(-e.delta / 120), "units"))
        self._formular(self.form)

        rechts = ttk.Frame(koerper)
        rechts.pack(side="left", fill="both", expand=True, padx=(10, 0))
        self.nb = ttk.Notebook(rechts)
        self.nb.pack(fill="both", expand=True)

        self.tab_txt = ttk.Frame(self.nb)
        self.nb.add(self.tab_txt, text="  Bemessungsbericht  ")
        rahmen = ttk.Frame(self.tab_txt)
        rahmen.pack(fill="both", expand=True)
        self.txt = tk.Text(rahmen, wrap="none", font=("Consolas", 9), bg=CARD,
                           fg=INK, bd=0, padx=10, pady=8)
        sy = ttk.Scrollbar(rahmen, orient="vertical", command=self.txt.yview)
        sx = ttk.Scrollbar(self.tab_txt, orient="horizontal", command=self.txt.xview)
        self.txt.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        self.txt.pack(side="left", fill="both", expand=True)
        sy.pack(side="right", fill="y")
        sx.pack(fill="x")
        self.txt.tag_configure("h", foreground=ACC, font=("Consolas", 9, "bold"))
        self.txt.tag_configure("ok", foreground=OKC, font=("Consolas", 9, "bold"))
        self.txt.tag_configure("no", foreground=NOKC, font=("Consolas", 9, "bold"))
        self.txt.tag_configure("din", foreground="#4a3aa7")

        fuss = ttk.Frame(self)
        fuss.pack(fill="x", padx=10, pady=(0, 8))
        self.lbl_status = ttk.Label(fuss, text="Bereit", style="Unter.TLabel")
        self.lbl_status.pack(side="left")
        ttk.Label(fuss, text="Berechnungs- und Lehrwerkzeug: Ergebnisse gegen die "
                             "gueltige Normausgabe pruefen.",
                  style="Unter.TLabel").pack(side="right")

    def _formular(self, f):
        r = 0
        g = ttk.LabelFrame(f, text=" 1. STATISCHES SYSTEM   [EC2 5.3.2 / 5.1.3] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=(6, 4)); r += 1
        self.c_system = Feld(g, 0, "Systemtyp", "Einfeldtraeger", breite=30,
                             typ="combo", auswahl=list(SYSTEME))
        self.c_system.w.bind("<<ComboboxSelected>>", self._system_geaendert)
        self.c_L = Feld(g, 1, "Gesamtlaenge L", 7.0, einheit="m",
                        hilfe="Gesamtlaenge des Traegers (Summe aller Felder).")
        self.c_auflager = Feld(g, 2, "Auflagerlagen", "0; 7.0", breite=26,
                               einheit="m",
                               hilfe="Abszissen der Auflager, durch ';' getrennt. "
                                     "Ausdruecke mit L sind zulaessig, z.B. "
                                     "'0; L/2; L'.")
        self.c_typen = Feld(g, 3, "Typ (g=gelenkig, e=eingesp.)", "g; g",
                            breite=26,
                            hilfe="Ein Zeichen je Auflager: 'g' gelenkig, "
                                  "'e' eingespannt.")
        self.c_breite_a = Feld(g, 4, "Auflagerbreite", 0.30, einheit="m",
                               hilfe="Fuer die Bemessungsquerkraft im Abstand d "
                                     "vom Auflagerrand (EC2 6.2.1 (8)).")

        g = ttk.LabelFrame(f, text=" 2. QUERSCHNITT   [EC2 5.3.2.1] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=4); r += 1
        self.c_typ = Feld(g, 0, "Querschnittstyp", "rechteck", breite=26,
                          typ="combo", auswahl=["rechteck", "plattenbalken"],
                          hilfe="rechteck = Rechteckquerschnitt ; "
                                "plattenbalken = Platte in der Druckzone (Feld).")
        self.c_b = Feld(g, 1, "Stegbreite  b_w", 300, einheit="mm")
        self.c_h = Feld(g, 2, "Gesamthoehe  h", 600, einheit="mm")
        self.c_beff = Feld(g, 3, "Mitwirkende Breite  b_eff", 1200, einheit="mm",
                           hilfe="Nur Plattenbalken. DIN EN 1992-1-1, 5.3.2.1 (3), "
                                 "Gl. (5.7): b_eff = SUM b_eff,i + b_w.")
        self.c_hf = Feld(g, 4, "Plattendicke  h_f", 150, einheit="mm")

        g = ttk.LabelFrame(f, text=" 3. BAUSTOFFE   [EC2 3.1 / 3.2, DIN 488-1] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=4); r += 1
        self.c_beton = Feld(g, 0, "Betonfestigkeitsklasse", "C30/37", breite=26,
                            typ="combo", auswahl=BETONKLASSEN,
                            hilfe="EC2 3.1.2, Tab. 3.1. fcd = 0,85 fck/1,5 "
                                  "(NA NDP zu 3.1.6 (1)P).")
        self.c_stahl = Feld(g, 1, "Betonstahl", "B500B", breite=26, typ="combo",
                            auswahl=["B500A", "B500B"],
                            hilfe="DIN 488-1. fyd = 500/1,15 = 434,8 N/mm2.")
        self.c_expo = Feld(g, 2, "Expositionsklasse", "XC1", breite=26,
                           typ="combo", auswahl=list(EXPOSITION),
                           hilfe="EC2 4.2, Tab. 4.1. Bestimmt c_min,dur "
                                 "(NA Tab. 4.4DE) und w_max (NA Tab. 7.1DE).")
        self.c_dg = Feld(g, 3, "Groesstkorn  d_g", 16, einheit="mm",
                         hilfe="EC2 8.2 (2): lichter Abstand >= d_g + 5 mm.")

        g = ttk.LabelFrame(f, text=" 4. EINWIRKUNGEN   [EC0 6.4.3.2 / EC1] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=4); r += 1
        self.c_gk = Feld(g, 0, "Staendige Last  g_k", 15.0, einheit="kN/m",
                         hilfe="OHNE Eigengewicht, wenn das folgende Feld "
                               "aktiviert ist.")
        self.c_eg = Feld(g, 1, "Eigengewicht ansetzen", True, typ="check",
                         hilfe="25,0 kN/m3 x b_w h  (DIN EN 1991-1-1, Tab. A.1).")
        self.c_qk = Feld(g, 2, "Veraenderliche Last  q_k", 20.0, einheit="kN/m")
        self.c_einzel = Feld(g, 3, "Einzellasten", "", breite=26,
                             hilfe="Format:  x,P,Art ; x,P,Art\n"
                                   "x in m, P in kN (nach unten), Art G oder Q.\n"
                                   "Beispiel: 3.0,50,Q ; 5.0,30,G")
        self.c_N = Feld(g, 4, "Laengskraft  N_Ed", 0.0, einheit="kN",
                        hilfe="DRUCK NEGATIV. Umrechnung des Momentes: "
                              "M_Eds = M_Ed - N_Ed z_s1.")
        self.c_gG = Feld(g, 5, "gamma_G", 1.35,
                         hilfe="DIN EN 1990/NA, Tab. NA.A.1.2(B).")
        self.c_gQ = Feld(g, 6, "gamma_Q", 1.50)
        self.c_psi2 = Feld(g, 7, "psi_2 (quasi-staendig)", 0.30,
                           hilfe="DIN EN 1990/NA, Tab. NA.A.1.1. Wohnraeume 0,3 | "
                                 "Bueros 0,3 | Verkaufsraeume 0,6 | Lager 0,8 | "
                                 "Schnee 0,0.")
        self.c_last = Feld(g, 8, "Feldweise Laststellung", True, typ="check",
                           hilfe="Durchlaeuft die 2^n Laststellungen der "
                                 "veraenderlichen Last (DIN EN 1992-1-1, 5.1.3).")
        self.c_delta = Feld(g, 9, "delta (Umlagerung)", 1.00,
                            hilfe="NA NDP zu 5.5 (4): xu/d <= (delta-0,64)/0,80. "
                                  "delta = 1,0 -> xu/d <= 0,45.")

        g = ttk.LabelFrame(f, text=" 5. TORSION   [EC2 6.3 + NA NDP zu 6.3.2] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=4); r += 1
        self.c_T = Feld(g, 0, "Torsionsmoment  T_Ed", 0.0, einheit="kNm",
                        hilfe="Bemessungswert des Torsionsmomentes. 0 = kein "
                              "Torsionsnachweis.")
        self.c_gleichgew = Feld(g, 1, "Gleichgewichtstorsion", True, typ="check",
                                hilfe="EC2 6.3.1 (2): Ist das Torsionsmoment "
                                      "fuer das Gleichgewicht NICHT erforderlich "
                                      "(Vertraeglichkeitstorsion), darf im GZT "
                                      "auf den Nachweis verzichtet werden.")
        self.c_kasten = Feld(g, 2, "Kastenquerschnitt", False, typ="check",
                             hilfe="NA NDP zu 6.3.2 (4): nu = 0,75 nu_2 beim "
                                   "Kastenquerschnitt, sonst 0,525 nu_2.")

        g = ttk.LabelFrame(f, text=" 6. BEWEHRUNG   [EC2 8.2 / 9.2] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=4); r += 1
        self.c_phi = Feld(g, 0, "Durchmesser unten", 20, breite=26, typ="combo",
                          auswahl=[str(int(p)) for p in DURCHMESSER], einheit="mm")
        self.c_phi_o = Feld(g, 1, "Durchmesser oben", 12, breite=26, typ="combo",
                            auswahl=[str(int(p)) for p in DURCHMESSER], einheit="mm")
        self.c_phi_w = Feld(g, 2, "Buegeldurchmesser", 8, breite=26, typ="combo",
                            auswahl=["6", "8", "10", "12", "14"], einheit="mm")
        self.c_schenkel = Feld(g, 3, "Buegelschenkel", 2, breite=26, typ="combo",
                               auswahl=["2", "4", "6"],
                               hilfe="Anzahl der lotrechten Schenkel je Buegel. "
                                     "Bei Torsion sind GESCHLOSSENE Buegel "
                                     "erforderlich (EC2 9.2.3).")

        g = ttk.LabelFrame(f, text=" 7. GEBRAUCHSTAUGLICHKEIT   [EC2 7] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=4); r += 1
        self.c_kriech = Feld(g, 0, "Kriechzahl  phi(inf,t0)", 2.0,
                             hilfe="EC2 3.1.4, Bild 3.1. Innenraeume ~2,0-2,5 ; "
                                   "im Freien ~1,5-2,0.")
        self.c_schwind = Feld(g, 1, "Schwinddehnung  eps_cs", 0.0, einheit="[-]",
                              hilfe="EC2 3.1.4 (6). z.B. 5e-4 = 0,5 Promille. "
                                    "0 = nicht beruecksichtigt.")
        self.c_grenz = Feld(g, 2, "Durchbiegungsgrenze  L/", 250,
                            hilfe="EC2 7.4.1 (4): L/250 unter quasi-staendiger "
                                  "Kombination.")

        bg = ttk.Frame(f)
        bg.grid(row=r, column=0, sticky="ew", padx=6, pady=(10, 14)); r += 1
        ttk.Button(bg, text="BERECHNEN", style="Aktion.TButton",
                   command=self.berechnen).pack(fill="x", pady=(0, 6))
        b2 = ttk.Frame(bg); b2.pack(fill="x")
        ttk.Button(b2, text="Bericht .txt", command=self.export_txt)\
            .pack(side="left", expand=True, fill="x", padx=(0, 3))
        ttk.Button(b2, text="Grafiken .pdf", command=self.export_pdf)\
            .pack(side="left", expand=True, fill="x", padx=3)
        ttk.Button(b2, text="Normen", command=self.zeige_normen)\
            .pack(side="left", expand=True, fill="x", padx=(3, 0))

    # -- Ereignisse ------------------------------------------------------
    def _system_geaendert(self, _e=None):
        sys_ = self.c_system.get()
        if sys_ in SYSTEME and sys_ != "Benutzerdefiniert":
            a, t = SYSTEME[sys_]
            self.c_auflager.var.set(self._auflager_text(a, self.c_L.f(7.0)))
            self.c_typen.var.set(t)

    @staticmethod
    def _zahl(txt, L):
        """Wertet einen einfachen arithmetischen Ausdruck mit L aus."""
        t = str(txt).strip().replace(",", ".")
        if not t:
            raise ValueError("Leerer Wert")
        if not set(t) <= set("0123456789.+-*/() L"):
            raise ValueError("Ungueltiger Ausdruck: %r" % txt)
        t = re.sub(r"([\d\)])\s*L", r"\1*L", t)      # '2L/3' -> '2*L/3'
        return float(eval(t.replace("L", "(%r)" % float(L)),
                          {"__builtins__": {}}, {}))

    @classmethod
    def _auflager_text(cls, ausdruck, L):
        out = []
        for t in ausdruck.split(";"):
            if not t.strip():
                continue
            try:
                out.append("%g" % cls._zahl(t, L))
            except Exception:
                out.append(t.strip())
        return "; ".join(out)

    # -- Eingaben lesen --------------------------------------------------
    def _lesen(self):
        L = self.c_L.f(7.0)
        xs = [self._zahl(t, L) for t in self.c_auflager.get().split(";") if t.strip()]
        typen = [t.strip().lower()[:1] or "g"
                 for t in self.c_typen.get().split(";") if t.strip()]
        while len(typen) < len(xs):
            typen.append("g")
        br = self.c_breite_a.f(0.30)
        auflager = [(x, EINGESPANNT if t == "e" else GELENKIG, br)
                    for x, t in zip(xs, typen)]
        if not auflager:
            raise ValueError("Es muss mindestens ein Auflager definiert werden.")

        einzel = []
        for teil in self.c_einzel.get().split(";"):
            teil = teil.strip()
            if not teil:
                continue
            p = [q.strip() for q in teil.split(",")]
            if len(p) < 2:
                raise ValueError("Fehlerhafte Einzellast: '%s'" % teil)
            einzel.append((float(p[0]), float(p[1]),
                           (p[2].upper() if len(p) > 2 else "G")))

        return EingabeBalken(
            L=L, auflager=auflager,
            querschnittstyp=self.c_typ.get(), b=self.c_b.f(300), h=self.c_h.f(600),
            b_eff=self.c_beff.f(1200), hf=self.c_hf.f(150),
            betonklasse=self.c_beton.get(), stahlsorte=self.c_stahl.get(),
            expositionsklasse=self.c_expo.get(), d_g=self.c_dg.f(16),
            g_k=self.c_gk.f(0), q_k=self.c_qk.f(0), eigengewicht=self.c_eg.b(),
            einzellasten=einzel, psi_2=self.c_psi2.f(0.3),
            gamma_G=self.c_gG.f(1.35), gamma_Q=self.c_gQ.f(1.50),
            laststellungen=self.c_last.b(), N_Ed=self.c_N.f(0.0),
            T_Ed=self.c_T.f(0.0),
            gleichgewichtstorsion=self.c_gleichgew.b(),
            kastenquerschnitt=self.c_kasten.b(),
            phi_laengs=self.c_phi.f(20), phi_laengs_oben=self.c_phi_o.f(12),
            phi_buegel=self.c_phi_w.f(8), n_schenkel=self.c_schenkel.i(2),
            delta=self.c_delta.f(1.0), phi_kriech=self.c_kriech.f(2.0),
            eps_cs=self.c_schwind.f(0.0),
            grenze_durchbiegung=self.c_grenz.f(250),
        )

    # -- Berechnung ------------------------------------------------------
    def berechnen(self):
        self.lbl_status.config(text="Berechnung laeuft ...", foreground=INK2)
        self.update_idletasks()
        try:
            e = self._lesen()
            b = bemessung_balken(e)
        except Exception as exc:
            messagebox.showerror("Fehler in der Berechnung",
                                 "%s\n\n%s" % (exc, traceback.format_exc()))
            self.lbl_status.config(text="Fehler: %s" % exc, foreground=NOKC)
            return
        self.bericht = b
        self._bericht_anzeigen(b)
        self._bilder_anzeigen(b)
        ok = b["ok_gesamt"]
        t_txt = "  |  T_Ed = %.0f kNm" % e.T_Ed if abs(e.T_Ed) > 1e-9 else ""
        self.lbl_status.config(
            text=("NACHWEISE ERFUELLT  |  " if ok else "NACHWEISE NICHT ERFUELLT  |  ")
                 + "unten %d fi%.0f  |  oben %d fi%.0f  |  Buegel fi%.0f%s"
                 % (b["n_unten"], e.phi_laengs, b["n_oben"], e.phi_laengs_oben,
                    e.phi_buegel, t_txt),
            foreground=(OKC if ok else NOKC))

    def _bericht_anzeigen(self, b):
        self.txt.config(state="normal")
        self.txt.delete("1.0", "end")
        for zeile in bericht_text(b, mit_normen=False).split("\n"):
            tag = ""
            zs = zeile.strip()
            if zs.startswith("[OK]"):
                tag = "ok"
            elif zs.startswith("[NEIN]") or zs.startswith("!!") or "!!" in zeile:
                tag = "no"
            elif zs.startswith("- DIN") or zs.startswith("- EA") \
                    or zs.startswith("Normgrundlage"):
                tag = "din"
            elif (zeile and zeile[0].isdigit() and ". " in zeile[:5]) \
                    or zs.startswith("---") or zs.startswith("==="):
                tag = "h"
            elif "GESAMTERGEBNIS" in zeile:
                tag = "ok" if b["ok_gesamt"] else "no"
            self.txt.insert("end", zeile + "\n", tag)
        self.txt.config(state="disabled")

    def _bilder_anzeigen(self, b):
        for tab, canvas, fig in self.canvases:
            canvas.get_tk_widget().destroy()
            plt.close(fig)
            self.nb.forget(tab)
        self.canvases = []
        for name, fig in G.alle_bilder(b):
            tab = ttk.Frame(self.nb)
            self.nb.add(tab, text="  %s  " % name)
            c = FigureCanvasTkAgg(fig, master=tab)
            c.draw()
            NavigationToolbar2Tk(c, tab).update()
            c.get_tk_widget().pack(fill="both", expand=True)
            self.canvases.append((tab, c, fig))

    # -- Export ----------------------------------------------------------
    def export_txt(self):
        if not self.bericht:
            return messagebox.showinfo("Bericht", "Bitte zuerst berechnen.")
        pfad = filedialog.asksaveasfilename(
            defaultextension=".txt", initialfile="Bemessungsbericht_Balken.txt",
            filetypes=[("Text", "*.txt")])
        if not pfad:
            return
        with open(pfad, "w", encoding="utf-8") as fh:
            fh.write(bericht_text(self.bericht, mit_normen=True))
        messagebox.showinfo("Bericht", "Gespeichert:\n%s" % pfad)

    def export_pdf(self):
        if not self.bericht:
            return messagebox.showinfo("Grafiken", "Bitte zuerst berechnen.")
        pfad = filedialog.asksaveasfilename(
            defaultextension=".pdf", initialfile="Grafiken_Balken.pdf",
            filetypes=[("PDF", "*.pdf")])
        if not pfad:
            return
        G.pdf_export(self.bericht, pfad)
        messagebox.showinfo("Grafiken", "Gespeichert:\n%s" % pfad)

    def zeige_normen(self):
        w = tk.Toplevel(self)
        w.title("Verwendete Normstellen")
        w.geometry("1100x700")
        t = tk.Text(w, wrap="none", font=("Consolas", 9), bg=CARD, fg=INK,
                    padx=10, pady=8)
        sb = ttk.Scrollbar(w, orient="vertical", command=t.yview)
        t.configure(yscrollcommand=sb.set)
        t.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        t.insert("1.0", normentabelle())
        t.config(state="disabled")


if __name__ == "__main__":
    App().mainloop()
