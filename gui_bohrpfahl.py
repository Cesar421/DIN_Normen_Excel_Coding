# -*- coding: utf-8 -*-
"""
===============================================================================
 BEMESSUNG VON BOHRPFAEHLEN NACH DIN
 Grafische Benutzeroberflaeche (Tkinter + matplotlib)
===============================================================================
 Normgrundlage
   DIN EN 1536:2015-10        Bohrpfaehle - Ausfuehrung, Bewehrung, Deckung
   DIN EN 1997-1:2009-09      Eurocode 7 - Geotechnik
   DIN 1054:2010-12           Teilsicherheitsbeiwerte, Nachweis GEO-2
   EA-Pfaehle (DGGT)          Erfahrungswerte, WSL, Bettungsmodulverfahren
   DIN EN 1992-1-1 + NA       Querschnittsbemessung (M-N, Querkraft)

 Bemessen werden LAENGSBEWEHRUNG (M-N-Interaktion) und QUERBEWEHRUNG (Wendel)
 sowie die axiale und horizontale Tragfaehigkeit.

 Start:   python gui_bohrpfahl.py
===============================================================================
"""

import os
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

from din_balken.baustoffe import BETONKLASSEN, DURCHMESSER
from din_pfahl.bettung import Bodenschicht
from din_pfahl.bemessung_pfahl import (EingabePfahl, bemessung_pfahl,
                                       bericht_text)
from din_pfahl import grafiken_pfahl as GP
from din_pfahl.normen_pfahl import normentabelle

BG = "#f4f4f1"
CARD = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
ACC = "#2a78d6"
OKC = "#0ca30c"
NOKC = "#d03b3b"

# Beispielprofil (Startwerte); q_s,k und q_b,k sind vom Anwender einzugeben
BEISPIEL_SCHICHTEN = [
    ("Auffuellung / Weichschicht", 0.0, 4.0, 8000.0, 30.0, 25.0),
    ("Sand, mitteldicht", 4.0, 10.0, 25000.0, 70.0, 0.0),
    ("Sand, dicht", 10.0, 15.0, 60000.0, 120.0, 0.0),
]


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
        self.title("Bemessung von Bohrpfaehlen - DIN EN 1536 | DIN EN 1997-1 / "
                   "DIN 1054 | EA-Pfaehle | DIN EN 1992-1-1 + NA")
        self.geometry("1600x980")
        self.minsize(1220, 760)
        self.configure(bg=BG)
        self.bericht = None
        self.canvases = []
        self._stile()
        self._aufbau()
        self.after(250, self.berechnen)

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
        s.configure("Treeview", background=CARD, fieldbackground=CARD,
                    font=("Consolas", 8.5), rowheight=21)
        s.configure("Treeview.Heading", font=("Segoe UI", 8, "bold"))
        s.configure("TNotebook.Tab", padding=(11, 6), font=("Segoe UI", 9))

    def _aufbau(self):
        kopf = ttk.Frame(self)
        kopf.pack(fill="x", padx=10, pady=(8, 4))
        ttk.Label(kopf, text="Bemessung von Bohrpfaehlen",
                  style="Titel.TLabel").pack(anchor="w")
        ttk.Label(kopf, text="DIN EN 1536:2015-10  |  DIN EN 1997-1 + "
                             "DIN 1054:2010-12 (GEO-2)  |  EA-Pfaehle (DGGT)  |  "
                             "DIN EN 1992-1-1 + NA  -  Laengs- und Querbewehrung",
                  style="Unter.TLabel").pack(anchor="w")

        koerper = ttk.Frame(self)
        koerper.pack(fill="both", expand=True, padx=10, pady=6)

        links = ttk.Frame(koerper, width=470)
        links.pack(side="left", fill="y")
        links.pack_propagate(False)
        cnv = tk.Canvas(links, bg=BG, highlightthickness=0, width=450)
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
        ttk.Label(fuss, text="q_b,k und q_s,k sind Eingabewerte nach EA-Pfaehle "
                             "Tab. 5.12-5.15 oder Probebelastung.",
                  style="Unter.TLabel").pack(side="right")

    def _formular(self, f):
        r = 0
        g = ttk.LabelFrame(f, text=" 1. PFAHLGEOMETRIE   [DIN EN 1536, 7.6] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=(6, 4)); r += 1
        self.c_D = Feld(g, 0, "Pfahldurchmesser  D", 900, einheit="mm",
                        hilfe="Uebliche Bohrpfahldurchmesser: 300 bis 3000 mm.")
        self.c_L = Feld(g, 1, "Pfahllaenge  L", 15.0, einheit="m")
        self.c_korb = Feld(g, 2, "Laenge Bewehrungskorb", 0.0, einheit="m",
                           hilfe="0 = volle Pfahllaenge. Der Korb muss "
                                 "mindestens bis unter den Bereich mit "
                                 "Zugspannungen reichen.")
        self.c_stuetz = Feld(g, 3, "Betonage unter Stuetzfluessigkeit", False,
                             typ="check",
                             hilfe="DIN EN 1536, 7.6.2: dann c_nom >= 75 mm.")

        g = ttk.LabelFrame(f, text=" 2. BAUSTOFFE   [EC2 3.1 / 3.2, DIN EN 1536 6.3] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=4); r += 1
        self.c_beton = Feld(g, 0, "Betonfestigkeitsklasse", "C25/30", breite=26,
                            typ="combo", auswahl=BETONKLASSEN,
                            hilfe="DIN EN 1536, 6.3: Pfahlbeton i.d.R. "
                                  ">= C25/30, weiche bis fliessfaehige "
                                  "Konsistenz.")
        self.c_stahl = Feld(g, 1, "Betonstahl", "B500B", breite=26, typ="combo",
                            auswahl=["B500A", "B500B"])
        self.c_dg = Feld(g, 2, "Groesstkorn  d_g", 16, einheit="mm",
                         hilfe="DIN EN 1536, 7.6.3: lichter Stababstand "
                               ">= 100 mm; >= 80 mm bei d_g <= 20 mm.")

        g = ttk.LabelFrame(f, text=" 3. BAUGRUND   [EA-Pfaehle, Tab. 5.12-5.15] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=4); r += 1
        ttk.Label(g, text="Schichten (Doppelklick zum Aendern):",
                  style="Feld.TLabel").grid(row=0, column=0, columnspan=4,
                                            sticky="w", padx=8, pady=(4, 2))
        spalten = ("name", "z_o", "z_u", "E_s", "q_s_k", "c_u")
        self.tree = ttk.Treeview(g, columns=spalten, show="headings", height=6)
        for c, t, w in zip(spalten,
                           ("Schicht", "z_o [m]", "z_u [m]", "E_s [kN/m2]",
                            "q_s,k [kN/m2]", "c_u [kN/m2]"),
                           (150, 55, 55, 85, 85, 75)):
            self.tree.heading(c, text=t)
            self.tree.column(c, width=w, anchor="center" if c != "name" else "w")
        self.tree.grid(row=1, column=0, columnspan=4, sticky="ew", padx=8)
        for s in BEISPIEL_SCHICHTEN:
            self.tree.insert("", "end", values=(s[0], s[1], s[2], s[3], s[4], s[5]))
        self.tree.bind("<Double-1>", self._schicht_bearbeiten)
        bb = ttk.Frame(g)
        bb.grid(row=2, column=0, columnspan=4, sticky="ew", padx=8, pady=4)
        ttk.Button(bb, text="+ Schicht", width=11,
                   command=self._schicht_neu).pack(side="left")
        ttk.Button(bb, text="- Schicht", width=11,
                   command=self._schicht_loeschen).pack(side="left", padx=4)
        self.c_qb = Feld(g, 3, "Spitzendruck  q_b,k", 1800.0, einheit="kN/m2",
                         hilfe="Charakteristischer Wert nach EA-Pfaehle, "
                               "Tab. 5.12/5.14 (abhaengig von q_c bzw. c_u) "
                               "oder aus einer Probebelastung.")
        self.c_mantel_ab = Feld(g, 4, "Mantelreibung ab Tiefe", 0.0, einheit="m",
                                hilfe="Tiefe, ab der Mantelreibung angesetzt "
                                      "wird (z.B. Auffuellung ohne Ansatz).")
        self.c_situation = Feld(g, 5, "Bemessungssituation", "BS-P", breite=26,
                                typ="combo", auswahl=["BS-P", "BS-T", "BS-A"],
                                hilfe="DIN 1054, Tab. A 2.3: BS-P staendig, "
                                      "BS-T voruebergehend, BS-A "
                                      "aussergewoehnlich.")

        g = ttk.LabelFrame(f, text=" 4. EINWIRKUNGEN AM PFAHLKOPF   [EC0 6.10] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=4); r += 1
        self.c_N = Feld(g, 0, "Laengskraft  N_Ed", -2500.0, einheit="kN",
                        hilfe="Bemessungswert, DRUCK NEGATIV.")
        self.c_Nk = Feld(g, 1, "N_k (charakteristisch)", -1850.0, einheit="kN",
                         hilfe="Fuer die Setzungsermittlung aus der "
                               "Widerstands-Setzungs-Linie.")
        self.c_H = Feld(g, 2, "Horizontalkraft  H_Ed", 150.0, einheit="kN")
        self.c_M = Feld(g, 3, "Kopfmoment  M_Ed", 0.0, einheit="kNm")
        self.c_kopf = Feld(g, 4, "Kopfausbildung", "frei", breite=26, typ="combo",
                           auswahl=["frei", "eingespannt"],
                           hilfe="frei = frei drehbar; eingespannt = drehstarr "
                                 "in die Kopfplatte eingebunden.")
        self.c_Nzug = Feld(g, 5, "Zugkraft  N_Ed,Zug", 0.0, einheit="kN",
                           hilfe="Positiver Wert, falls der Pfahl auf Zug "
                                 "beansprucht wird (DIN EN 1997-1, 7.6.3).")

        g = ttk.LabelFrame(f, text=" 5. BEWEHRUNG   [DIN EN 1536, 7.6.3/7.6.4] ")
        g.grid(row=r, column=0, sticky="ew", padx=6, pady=4); r += 1
        self.c_phi = Feld(g, 0, "Laengsstabdurchmesser", 20, breite=26,
                          typ="combo", auswahl=[str(int(p)) for p in DURCHMESSER
                                                if p >= 16], einheit="mm",
                          hilfe="DIN EN 1536, 7.6.3: phi >= 16 mm.")
        self.c_n = Feld(g, 1, "Stabanzahl (0 = automatisch)", 0, breite=26,
                        hilfe="0 = die kleinste erforderliche Anzahl wird "
                              "ermittelt. Mindestens 6 Staebe.")
        self.c_phi_w = Feld(g, 2, "Wendeldurchmesser", 10, breite=26,
                            typ="combo", auswahl=["6", "8", "10", "12", "14", "16"],
                            einheit="mm",
                            hilfe="DIN EN 1536, 7.6.4: phi >= 6 mm; "
                                  "EC2 9.5.3: phi_w >= phi_l/4.")

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

    # -- Schichtentabelle ------------------------------------------------
    def _schicht_neu(self):
        kinder = self.tree.get_children()
        z_o = float(self.tree.item(kinder[-1])["values"][2]) if kinder else 0.0
        self.tree.insert("", "end", values=("Neue Schicht", z_o, z_o + 3.0,
                                            20000.0, 60.0, 0.0))

    def _schicht_loeschen(self):
        for it in self.tree.selection():
            self.tree.delete(it)

    def _schicht_bearbeiten(self, ereignis):
        item = self.tree.identify_row(ereignis.y)
        spalte = self.tree.identify_column(ereignis.x)
        if not item or not spalte:
            return
        idx = int(spalte[1:]) - 1
        x, y, w, h = self.tree.bbox(item, spalte)
        wert = self.tree.item(item)["values"][idx]
        ed = ttk.Entry(self.tree)
        ed.place(x=x, y=y, width=w, height=h)
        ed.insert(0, str(wert))
        ed.focus_set()

        def fertig(_e=None):
            werte = list(self.tree.item(item)["values"])
            werte[idx] = ed.get()
            self.tree.item(item, values=werte)
            ed.destroy()

        ed.bind("<Return>", fertig)
        ed.bind("<FocusOut>", fertig)

    def _schichten_lesen(self):
        schichten = []
        for it in self.tree.get_children():
            v = self.tree.item(it)["values"]
            schichten.append(Bodenschicht(
                z_o=float(str(v[1]).replace(",", ".")),
                z_u=float(str(v[2]).replace(",", ".")),
                E_s=float(str(v[3]).replace(",", ".")),
                q_s_k=float(str(v[4]).replace(",", ".")),
                c_u_k=float(str(v[5]).replace(",", ".")),
                name=str(v[0])))
        if not schichten:
            raise ValueError("Es muss mindestens eine Bodenschicht "
                             "definiert werden.")
        return schichten

    # -- Berechnung ------------------------------------------------------
    def _lesen(self):
        return EingabePfahl(
            D=self.c_D.f(900), L=self.c_L.f(15.0),
            unter_stuetzfluessigkeit=self.c_stuetz.b(),
            betonklasse=self.c_beton.get(), stahlsorte=self.c_stahl.get(),
            d_g=self.c_dg.f(16), schichten=self._schichten_lesen(),
            q_b_k=self.c_qb.f(1800.0), mantel_ab_tiefe=self.c_mantel_ab.f(0.0),
            situation=self.c_situation.get(),
            N_Ed=self.c_N.f(-2500.0), N_k=self.c_Nk.f(-1850.0),
            H_Ed=self.c_H.f(0.0), M_Ed_kopf=self.c_M.f(0.0),
            kopf=self.c_kopf.get(), N_Ed_zug=self.c_Nzug.f(0.0),
            phi_l=self.c_phi.f(20), n_l=self.c_n.i(0),
            phi_w=self.c_phi_w.f(10),
            laenge_bewehrungskorb=self.c_korb.f(0.0))

    def berechnen(self):
        self.lbl_status.config(text="Berechnung laeuft ...", foreground=INK2)
        self.update_idletasks()
        try:
            e = self._lesen()
            b = bemessung_pfahl(e)
        except Exception as exc:
            messagebox.showerror("Fehler in der Berechnung",
                                 "%s\n\n%s" % (exc, traceback.format_exc()))
            self.lbl_status.config(text="Fehler: %s" % exc, foreground=NOKC)
            return
        self.bericht = b
        self._bericht_anzeigen(b)
        self._bilder_anzeigen(b)
        qs, qk = b["querschnitt"], b["querkraft"]
        ok = b["ok_gesamt"]
        self.lbl_status.config(
            text=("NACHWEISE ERFUELLT  |  " if ok else "NACHWEISE NICHT ERFUELLT  |  ")
                 + "%d fi%.0f (rho = %.2f %%)  |  Wendel fi%.0f/%.0f mm  |  "
                   "R_c,d = %.0f kN"
                 % (qs.n_l, qs.phi_l, 100 * qs.rho_l, qs.phi_w,
                    qk["s_gewaehlt"], b["tragfaehigkeit"].R_c_d),
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
        for name, fig in GP.alle_bilder(b):
            tab = ttk.Frame(self.nb)
            self.nb.add(tab, text="  %s  " % name)
            c = FigureCanvasTkAgg(fig, master=tab)
            c.draw()
            NavigationToolbar2Tk(c, tab).update()
            c.get_tk_widget().pack(fill="both", expand=True)
            self.canvases.append((tab, c, fig))

    def export_txt(self):
        if not self.bericht:
            return messagebox.showinfo("Bericht", "Bitte zuerst berechnen.")
        pfad = filedialog.asksaveasfilename(
            defaultextension=".txt", initialfile="Bemessungsbericht_Bohrpfahl.txt",
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
            defaultextension=".pdf", initialfile="Grafiken_Bohrpfahl.pdf",
            filetypes=[("PDF", "*.pdf")])
        if not pfad:
            return
        GP.pdf_export(self.bericht, pfad)
        messagebox.showinfo("Grafiken", "Gespeichert:\n%s" % pfad)

    def zeige_normen(self):
        w = tk.Toplevel(self)
        w.title("Verwendete Normstellen - Pfahlbemessung")
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
