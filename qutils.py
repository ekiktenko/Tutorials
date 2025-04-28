import numpy as np
from numpy import pi as PI
from numpy import ndarray, sqrt, sin, cos, exp, log2, log, real, imag
from scipy.linalg import expm
from itertools import count

# Some basic 'constants'
SIGMA = np.array(
    [np.eye(2), [[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]]
)  # x-, y-, z-Pauli matrices

PSI_PLUS = np.array([1, 1]) / sqrt(2)  # |+>
PSI_MINUS = np.array([1, -1]) / sqrt(2)  # |->

SPIN_CR = np.array([[0, 1], [0, 0]])  # spin raising
SPIN_AN = np.transpose(SPIN_CR)  # spin lowering


# Some special fuctions
def put_paulis_on_pos(pauli_str: list[int], ind_list: list[int], q_num: int):
    """Generates 'q_num'-qubit Pauli string with Pauli operators taken from
    'pauli_str' and put on corresponding positions from 'ind_list'.
    For expample for (pauli_str=[3,2,1], ind_list=[0,3,4], q_num=5)
    the function returns ZIIYX operator.
    """
    assert len(pauli_str) == len(ind_list)
    assert np.max(ind_list) < q_num and np.min(ind_list) >= 0
    M = 1.0
    for i in range(q_num):
        if i in ind_list:
            M = np.kron(M, SIGMA[pauli_str[ind_list.index(i)]])
        else:
            M = np.kron(M, SIGMA[0])
    return M


def put_matr_on_pos(matr: ndarray, pos: int, q_num: int, dim: int = 2):
    """Put matr on position pos within tensor product of q_num identities"""
    assert 0 <= pos < q_num
    M = 1.0
    for i in range(q_num):
        if i == pos:
            M = np.kron(M, matr)
        else:
            M = np.kron(M, np.eye(dim))
    return M


def hermitian_check(H: ndarray):
    """Checks whether a given matrix 'rho' is Hermitian"""
    assert H.ndim == 2 and H.shape[0] == H.shape[1]
    assert np.isclose(np.linalg.norm(H - H.conj().transpose()), 0)


def densmatr_check(rho: ndarray, need_normalization: bool = True):
    """Checks whether a given matrix 'rho' is a fair dansity matrix"""
    hermitian_check(rho)
    eigvals = np.linalg.eigh(rho)[0]
    assert np.isclose(np.min(eigvals), 0) or np.min(eigvals) > 0
    if need_normalization:
        assert np.isclose(np.sum(eigvals), 1)


def plogp(p: float, base: int = 2):
    """Makes calculation of p * log_base(p) for entropy functions
    avoiding possible numerical errors related to p near 0"""
    if np.isclose(p, 0):
        return 0
    elif p > 0:
        return p * np.log2(p) / np.log2(base)
    else:
        raise (ValueError)


def vnentropy(rho: ndarray, base: int = 2):
    """Calculates von Neuman entropy of a given density matrix 'rho'"""
    densmatr_check(rho)
    eigs, _ = np.linalg.eigh(rho)
    entropy = 0
    for eigval in eigs:
        entropy -= plogp(eigval, base=base)
    return entropy


def sqrtm_(M: ndarray):
    """Self-made matrix sqrt function to avoid erros"""
    hermitian_check(M)
    E, U = np.linalg.eigh(M)
    return U @ np.diag(np.sqrt(np.abs(E))) @ U.conj().transpose()


def fidelity(rho_psi_1: ndarray, rho_psi_2: ndarray):
    """Computes fidelity between normalized density matrix/state and another density matrix/state"""
    if rho_psi_1.ndim == 1:
        assert np.isclose(np.linalg.norm(rho_psi_1), 1)
    if rho_psi_2.ndim == 1:
        assert np.isclose(np.linalg.norm(rho_psi_2), 1)
    if rho_psi_1.ndim == 2:
        assert np.isclose(np.linalg.trace(rho_psi_1), 1)
    if rho_psi_2.ndim == 2:
        assert np.isclose(np.linalg.trace(rho_psi_2), 1)

    if rho_psi_1.ndim == 1 and rho_psi_2.ndim == 1:  # <psi1|rho2|psi1>
        return np.abs(rho_psi_1 @ rho_psi_2) ** 2
    elif rho_psi_1.ndim == 1 and rho_psi_2.ndim == 2:  # <psi1|rho2|psi1>
        return np.real(np.einsum("i,ij,j->", rho_psi_1.conj(), rho_psi_2, rho_psi_1))
    elif rho_psi_1.ndim == 2 and rho_psi_2.ndim == 1:  # <psi2|rho1|psi2>
        return np.real(np.einsum("i,ij,j->", rho_psi_2.conj(), rho_psi_1, rho_psi_2))
    elif rho_psi_1.ndim == 2 and rho_psi_2.ndim == 2:  # <psi2|rho1|psi2>
        return (np.real(np.trace(sqrtm_(rho_psi_1) @ sqrtm_(rho_psi_2)))) ** 2


def get_reduced_densmatr(psi_rho: ndarray, ind_list: list[int], n: int, dim: int = 2):
    """Make a reduced density matrix from pure state 'psi' for qubits from ind_list"""
    assert np.min(ind_list) >= 0 and np.max(ind_list) < n
    m = len(ind_list)
    if psi_rho.ndim == 2:
        density_matrix_case = True
        densmatr_check(psi_rho)
        rho_ = np.reshape(psi_rho, [dim] * n * 2)
    elif psi_rho.ndim == 1:
        density_matrix_case = False
        psi_ = np.reshape(psi_rho, [dim] * n)
    else:
        raise ValueError
    ind_for_psi = list(range(n))
    ind_for_psi_conj = list(range(n))
    cnt = count(n)
    for pos in ind_list:
        ind_for_psi_conj[pos] = next(cnt)
    ind_for_rho = ind_list + list(range(n, n + m))

    if density_matrix_case:
        rho_out = np.einsum(rho_, ind_for_psi + ind_for_psi_conj, ind_for_rho)
    else:
        rho_out = np.einsum(
            psi_, ind_for_psi, psi_.conj(), ind_for_psi_conj, ind_for_rho
        )

    rho_out = np.reshape(rho_out, (dim**m, dim**m))
    densmatr_check(rho_out)
    return rho_out


def mutinfo(psi_rho, i, j, n, dim=2):
    rho_AB = get_reduced_densmatr(psi_rho, [i, j], n)
    rho_A = get_reduced_densmatr(rho_AB, [0], 2)
    rho_B = get_reduced_densmatr(rho_AB, [1], 2)
    return (
        vnentropy(rho_A, base=dim)
        + vnentropy(rho_B, base=dim)
        - vnentropy(rho_AB, base=dim)
    )


def get_mutinfo_matrix(psi_rho, n):
    info_matr = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            info_cur = mutinfo(psi_rho, i, j, n)
            info_matr[i, j] = info_cur
            info_matr[j, i] = info_cur
    return info_matr


def rho_thermal(hmlt: ndarray, beta: float = 1.0):
    """Makes Gibbs density matrix for Hamiltonian 'hmlt' and inverse temperature 'beta'"""
    hermitian_check(hmlt)
    rho = expm(-hmlt * beta)
    return rho / np.trace(rho)


def anticom(A: ndarray, B: ndarray):
    """Standard anicommutator of A and B"""
    assert A.shape == B.shape and A.shape[0] == A.shape[1], A.ndim == 2
    return A @ B + B @ A