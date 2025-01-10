import numpy as np
import math
from saving_data import *
import json


def scalarproduct(ket1, ket2):  # norm. Skalarprodukt, vergleicht nur ob zwei Kets gleich sind oder nicht und gibt
    # entsprechend 1 oder 0 aus
    ket1 = np.asarray(ket1)
    ket2 = np.asarray(ket2)
    if np.array_equal(ket1, ket2):
        return 1  # gibt 1 aus wenn beide kets gleich sind (Bra = ket)
    else:
        return 0

























