# -*- coding: utf-8 -*-
"""
din_pfahl - Bemessung von Bohrpfaehlen nach deutschen Normen.

Normgrundlage
-------------
DIN EN 1536:2015-10        Bohrpfaehle (Ausfuehrung, Bewehrung, Betondeckung)
DIN EN 1997-1:2009-09      Eurocode 7 - Geotechnik
DIN 1054:2010-12           Ergaenzende Regelungen (Teilsicherheitsbeiwerte)
EA-Pfaehle (DGGT)          Erfahrungswerte, WSL, Bettungsmodulverfahren
DIN EN 1992-1-1 + NA       Querschnittsbemessung (M-N, Querkraft)
DIN EN 1990 + NA           Einwirkungskombinationen

Umfang
------
  * axiale Tragfaehigkeit  R_c,d = R_b,k/gamma_b + R_s,k/gamma_s
  * Widerstands-Setzungs-Linie
  * horizontal belasteter Pfahl (Bettungsmodulverfahren, Winkler)
  * M-N-Interaktionsdiagramm des Kreisquerschnitts
  * Laengsbewehrung und Mindestbewehrung (DIN EN 1536, Tab. 4)
  * Querbewehrung (Wendel) einschliesslich Querkraftnachweis
  * Knicknachweis (nur bei sehr weichen Boeden erforderlich)

HINWEIS: Berechnungs- und Lehrwerkzeug. q_b,k und q_s,k sind Eingabewerte des
Anwenders. Die Ergebnisse sind von einem verantwortlichen Ingenieur zu pruefen.
"""

__version__ = "1.0.0"

from . import normen_pfahl                     # noqa: F401
from .kreisquerschnitt import (                # noqa: F401
    Kreisquerschnitt, interaktionsdiagramm, M_Rd_bei_N, schnittgroessen_bei_x,
    dehnungsebene, erforderliche_bewehrung, mindestbewehrung_pfahl,
    mindestbewehrung_druckglied, betondeckung_pfahl, konstruktive_pruefung,
)
from .bettung import (                         # noqa: F401
    Bodenschicht, pfahl_horizontal, bettungsmodul, knicklast_gebettet,
)
from .tragfaehigkeit import (                  # noqa: F401
    axiale_tragfaehigkeit, widerstands_setzungs_linie, pfahlkopfsetzung,
    GAMMA_R,
)
from .bemessung_pfahl import (                 # noqa: F401
    EingabePfahl, bemessung_pfahl, bericht_text, querkraft_kreis,
)
