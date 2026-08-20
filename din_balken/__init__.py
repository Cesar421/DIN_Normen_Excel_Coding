# -*- coding: utf-8 -*-
"""
din_balken - Bemessung von Stahlbetonbalken nach deutschen Normen.

Normgrundlage
-------------
DIN EN 1992-1-1:2011-01  (Eurocode 2, Teil 1-1)
DIN EN 1992-1-1/NA:2013-04 (+ A1:2015-12)  Nationaler Anhang
DIN EN 1990 / DIN EN 1991-1-1 (+ NA)       Einwirkungen und Kombinationen
DIN 488-1:2009-08                          Betonstahl B500A / B500B

Nachgewiesen werden BIEGUNG (EC2 6.1), QUERKRAFT (6.2) und TORSION (6.3)
einschliesslich ihrer Interaktion, sowie die Grenzzustaende der
Gebrauchstauglichkeit (Rissbreite 7.3, Verformung 7.4).

Jedes Ergebnis fuehrt die angewendete Normstelle mit (Modul `normen`).

HINWEIS: Berechnungs- und Lehrwerkzeug. Die Ergebnisse sind von einem
verantwortlichen Ingenieur gegen die gueltige Ausgabe der Norm zu pruefen.
"""

__version__ = "2.0.0"

from . import normen                          # noqa: F401
from .baustoffe import (                      # noqa: F401
    Beton, Betonstahl, BETONKLASSEN, DURCHMESSER, EXPOSITION,
    stabflaeche, stabflaeche_n, betondeckung, stabwahl,
)
from .querschnitt import Querschnitt, mitwirkende_plattenbreite   # noqa: F401
from .schnittgroessen import (                # noqa: F401
    Durchlauftraeger, Auflager, Streckenlast, Einzellast,
    GELENKIG, EINGESPANNT, bemessungsquerkraft, wert_bei,
)
from .biegung import (                        # noqa: F401
    bemessung_biegung, momententragfaehigkeit, bemessungsdiagramm, xi_grenz,
    dehnungszustand, betondruckkraft,
)
from .querkraft import (                      # noqa: F401
    bemessung_querkraft, innerer_hebelarm, V_Rd_c, cot_theta_NA, V_Rd_max,
    asw_erforderlich, asw_mindest, groesster_buegelabstand, versatzmass,
)
from .torsion import (                        # noqa: F401
    bemessung_torsion, ersatzhohlquerschnitt, T_Rd_c, T_Rd_max,
    asw_torsion, asl_torsion, bewehrung_entbehrlich, s_max_torsion,
)
from .gebrauchstauglichkeit import (          # noqa: F401
    rissbreite, mindestbewehrung_riss, nachweis_durchbiegung,
    zustand_I, zustand_II, stahlspannung, kruemmung, zulaessige_schlankheit,
)
from .konstruktion import (                   # noqa: F401
    mindestbewehrung_biegung, hoechstbewehrung, robustheitsbewehrung,
    verankerungslaenge, uebergreifungslaenge, platznachweis, d1_schaetzung,
    zugkraft,
)
from .bemessung import (                      # noqa: F401
    EingabeBalken, bemessung_balken, bericht_text,
)
