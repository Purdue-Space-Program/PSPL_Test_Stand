from press_master import *
from pint import UnitRegistry
import numpy as np

S_g_gn2 = Q_(0.9669,'dimensionless')

# Choked flow regulator Cv calculations
C_v_fu = (scfm_gn2_fu.to('ft^3/min').m * 2 * np.sqrt(S_g_gn2.m)) / P_BOTTLE.to('psi').m  # Note:units in psia and SCFM
C_v_lox = (scfm_gn2_lox.to('ft^3/min').m * 2 * np.sqrt(S_g_gn2.m)) / P_BOTTLE.to('psi').m  # Note:units in psia and SCFM

print("Required Cv for Methane: ", f"{C_v_fu: .4f}")
print("Required Cv for LOX: ", f"{C_v_lox: .4f}")