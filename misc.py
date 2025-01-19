from Masochismus.old_versions.saving_data import *
from saving_data import *

def scalarproduct(ket1, ket2):  # norm. Skalarprodukt, vergleicht nur ob zwei Kets gleich sind oder nicht und gibt
    # entsprechend 1 oder 0 aus
    ket1 = np.asarray(ket1)
    ket2 = np.asarray(ket2)
    if np.array_equal(ket1, ket2):
        return 1  # gibt 1 aus wenn beide kets gleich sind (Bra = ket)
    else:
        return 0


def show_data(spin, min_spin, max_spin):
    for num_spins in range(min_spin, max_spin + 1):
        j_ij = 1
        system = load_system(spin, num_spins, j_ij)
        print(f"--------SYSTEM: s={spin}, num_spins = {num_spins}, j_ij = {j_ij}--------")
        print(f"eigenvalue ground state: \n{system['eigenvalues'][0]:.9f}")
        print("degeneracy ground state:\n ", system['degeneracies'][0])






















