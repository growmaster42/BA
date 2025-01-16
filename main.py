from Masochismus.old_versions.saving_data import *
from expectation_values import *
from saving_data import *
from concurrent.futures import ProcessPoolExecutor
import time
from plotting import degeneracy_plot

for i in range(-1, 2):
    degeneracy_plot(0.5, i, 0)
    degeneracy_plot(1, i, 0)
    degeneracy_plot(1.5, i, 0)
    degeneracy_plot(2, i, 0)
    degeneracy_plot(2.5, i, 0)