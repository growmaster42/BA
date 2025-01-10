import numpy as np


def basvec(spin, number_of_spins):
    number_of_states = int(2 * spin + 1)
    dim = number_of_states ** number_of_spins
    # generates an array with zeros, dim x and number_of_spins entries per ket
    ketname = np.zeros((dim, number_of_spins), dtype=int)
    # makes rev_ketname globally available
    global rev_ketname
    # for loop that creates the parameter a that will be divided by number_of_states using divmod
    for i in range(dim - 1, -1, -1):
        # e.g. 5 spins 1/2 will result in dim = 32 basis vectors
        a = dim - 1 - i
        # for loop that runs through every entry of the ket
        for j in range(number_of_spins - 1, -1, -1):
            # creating a tuple that corresponds to modulus, x is a tuple that which contains the result of the
            # division and the rest (a/number_of_states, rest)
            x = divmod(a, number_of_states)
            # this puts the second value of the tuple (the rest) in the previously created ketname,
            ketname[i, j] = x[1]
            # this overwrites the previously generated a with the result of the division in divmod
            a = int(x[0])
            # ketname is being flipped for no reason just to have it nicely arranged
    rev_ketname = ketname[::-1]
    # this returns the ascending basis vectors as an array and terminates the function
    return rev_ketname


# function that returns the vector number of any given basis vector
def vector_number(vector, spin):
    m = 0
    for i, entry in enumerate(vector[::-1]):
        m += entry * (2 * spin + 1) ** i
    return int(m)


