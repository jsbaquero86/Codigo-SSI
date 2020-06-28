# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 16:57:33 2020

@author: sebas
"""

"""
    Script with all unit transformation to kN-m-C --> kPa, kN-m, rad
"""

kN = 1.0            # [kN] basic force units.
m = 1.0             # [m] basic length units.
sec = 1.0           # [sec] basic time units.

minut = 60*sec      # minuts time unit.
hour = 60*minut     # hour time unit.

cm = 1e-2*m         # centimeters length units.
mm = 1e-3*m         # milimeters length units.
km = 1e3*m          # kilometer length units.
inch = 0.0254*m     # inches length units.
ft = 12*inch        # feet legth units.
mile = 1.6*km       # miles length units.

cm2 = cm**2         # area in cm.
cm3 = cm**3         # volume in cm.
cm4 = cm**4         # forth in cm.
m2 = m**2           # area in meters.
m3 = m**3           # volume in meters.
m4 = m**4           # forth in meters.
inch2 = inch**2     # area in inches.
inch3 = inch**3     # volume in inches.
inch4 = inch**4     # forth in inches.
ft2 = ft**2         # area in feet.
ft3 = ft**3         # volume in feet.
ft4 = ft**4         # forth in feet.

N = 1e-3*kN         # Newtons force units.
MN = 1e3*kN         # MegaNewtons force untis.
GN = 1e6*kN         # GigaNewtons force units.
kgf = 9.81*N        # kilogram-force units.
tonf = 9.81*kN      # Ton-force untis.
lbf = 1/2.204*kgf   # pound-force unit.
kip = 1e3*lbf       # kilo-pund-force units.

mps = m/sec         # velocity m/s.
kmph = km/hour      # velocity km/h.
mph = mile/hour     # velocity miles/h.

MPa = MN/m2         # Mega-Pascals pressure units.
Pa = N/m2           # Pascals pressure units.
GPa = GN/m2         # Giga-Pascals pressure units.
kgfcm2 = kgf/cm2    # kgf/cm2 pressure units.
psi = lbf/inch2     # pressure units - US.
ksi = kip/inch2     # pressure units - US.
