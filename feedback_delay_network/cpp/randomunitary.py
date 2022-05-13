import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt

def haar_measure(n, seed=0):
    rng = np.random.default_rng(seed)

    z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    q, r = linalg.qr(z)
    d = r.diagonal()
    q *= d / abs(d)
    return q

def random_ortho(dim, seed=0):
    """Comments from `scipy.stats.ortho_group` are preserved."""
    rng = np.random.default_rng(seed)

    H = np.eye(dim)
    for n in range(dim):
        x = rng.normal(size=(dim - n,))
        norm2 = np.dot(x, x)
        x0 = x[0].item()

        # random sign, 50/50, but chosen carefully to avoid roundoff error
        D = np.sign(x[0]) if x[0] != 0 else 1
        x[0] += D * np.sqrt(norm2)
        x /= np.sqrt((norm2 - x0**2 + x[0]**2) / 2.)

        # Householder transformation
        H[:, n:] = -D * (H[:, n:] - np.outer(np.dot(H[:, n:], x), x))
        # print(H)
    return H

def random_special_ortho(dim, seed=0):
    """Comments from `scipy.stats.special_ortho_group` are preserved."""
    rng = np.random.default_rng(seed)

    H = np.eye(dim)
    D = np.empty((dim,))
    for n in range(dim - 1):
        x = rng.normal(size=(dim - n,))
        norm2 = np.dot(x, x)
        x0 = x[0].item()
        D[n] = np.sign(x[0]) if x[0] != 0 else 1
        x[0] += D[n] * np.sqrt(norm2)
        x /= np.sqrt((norm2 - x0**2 + x[0]**2) / 2.)

        # Householder transformation
        H[:, n:] -= np.outer(np.dot(H[:, n:], x), x)

    D[-1] = (-1)**(dim - 1) * D[:-1].prod()

    # Equivalent to np.dot(np.diag(D), H) but faster, apparently
    H = (D * H.T).T
    return H

def random_ortho_close_to_identity(dim=4, identityAmount=0.1, seed=0):
    rng = np.random.default_rng(seed)

    H = np.eye(dim)
    for n in range(dim):
        x = rng.normal(size=(dim - n,)) * identityAmount
        x[0] = 1
        norm2 = np.dot(x, x)
        x0 = x[0].item()

        D = np.sign(x[0]) if x[0] != 0 else 1
        x[0] += D * np.sqrt(norm2)
        x /= np.sqrt((norm2 - x0**2 + x[0]**2) / 2.)

        H[:, n:] = -D * (H[:, n:] - np.outer(np.dot(H[:, n:], x), x))
    return H

def random_ortho_circulant(dim=4, seed=0):
    """
    Eq. (24) and (25) in the paper below.

    - Rocchesso, Davide, and Julius O. Smith. "Circulant and elliptic feedback delay networks for artificial reverberation." IEEE Transactions on Speech and Audio Processing 5.1 (1997): 51-63.
    """
    rng = np.random.default_rng(seed)
    source = rng.uniform(0, 1, dim)
    scale = 2 / np.sum(source)
    sqrt = np.sqrt(source)
    circ = np.repeat(sqrt, dim).reshape((dim, dim))
    return circ * sqrt * scale - np.eye(dim)

def random_triangular(dim=4, seed=0):
    """
    For feedback delay network, the advantage of this method is the ability to
    provide control for reverb time. However, the paper referenced in
    `random_ortho_circulant` argues that reverb in real room doesn't decay uniformly
    across all frequencies.

    - Jot, Jean-Marc, and Antoine Chaigne. "Digital delay networks for designing artificial reverberators." Audio Engineering Society Convention 90. Audio Engineering Society, 1991.

    ## Construct diagonal components fron seconds.
    ```
    sampleRate = 48000  # Set as constant for testing.
    reverbSeconds = 2
    delaySeconds = rng.uniform(0, 1, dim)
    diag = np.power(10, -3 * delaySeconds * reverbSeconds / sampleRate)
    ```
    """
    rng = np.random.default_rng(seed)
    mat = []
    summed = 0
    for idx in range(0, dim):
        row = np.zeros(dim)
        randomized = rng.uniform(0, 1, dim - idx - 1)
        row[idx + 1:] = randomized
        summed += np.sum(randomized)
        mat.append(row)
    scale = 1 / summed
    for idx in range(0, dim):
        mat[idx][idx] = scale - 1
        mat[idx][idx + 1:] *= scale
    return np.vstack(mat)

def random_schroeder(dim=4, seed=0):
    """
    Section IV. A. in the paper below.

    - Schlecht, Sebastian J., and Emanuel AP Habets. "On lossless feedback delay networks." IEEE Transactions on Signal Processing 65.6 (2016): 1554-1564.
    """
    if dim < 2:
        print("Error: dim must be greater than or equals to 2.")
        exit()
    rng = np.random.default_rng(seed)
    diag = rng.uniform(0, 1, dim)
    mat = np.zeros((dim, dim))
    np.fill_diagonal(mat, diag)
    mat[dim - 1][:-2] = np.full(dim - 2, -diag[-2])
    mat[dim - 1][-2] = 1 - diag[-2] * diag[-2]
    return mat

def random_absorbent(dim=4, seed=0):
    """
    Eq. (10) in the paper below.

    This one looks like equivalent to lattice all-pass topology.

    - Schlecht, Sebastian J., and Emanuël AP Habets. "Time-varying feedback matrices in feedback delay networks and their application in artificial reverberation." The Journal of the Acoustical Society of America 138.3 (2015): 1389-1398.
    """
    if dim < 2:
        print("Error: dim must be greater than or equals to 2.")
        exit()
    if dim % 2 != 0:
        print("Error: dim must be even integer.")
        exit()
    rng = np.random.default_rng(seed)
    half = dim // 2
    diag = rng.uniform(0, 1, half)

    G = np.zeros((half, half))  # Allpass factors.
    np.fill_diagonal(G, diag)

    G2 = np.zeros((half, half))
    np.fill_diagonal(G2, diag * diag)

    Q = random_ortho(half, seed)  # Q can be any unitary matrix.

    return np.vstack((
        np.hstack((-Q * G, Q)),
        np.hstack((np.eye(half) - G2, G)),
    ))

def hadamard_sylvester(dim=4):
    """
    Sylvester's construction of Hadamard matrix using following initial condition:

    ```
    [[1,  1],
     [1, -1]]
    ```
    """
    mat = np.empty((dim, dim), dtype=int)

    mat[0][0] = 1  #/ np.sqrt(dim)

    start = 1
    end = 2
    while start < dim:
        for row in range(start, end):
            for col in range(start, end):
                value = mat[row - start][col - start]
                mat[row - start][col] = value  # Upper right.
                mat[row][col - start] = value  # Lower left.
                mat[row][col] = -value  # Lower right.
                print(
                    f"({col:2d}, {row:2d}): {value:2d}, {mat[row - start][col]:2d}, {mat[row][col - start]:2d}, {mat[row][col]:2d}"
                )
        start *= 2
        end *= 2
        print(start, end)
    print(mat)
    exit()
    return mat

def random_conference(dim=4):
    """
    Following conditions are required. Reference: https://oeis.org/A000952

    - `dim mod 4 == 2`.
    - `dim - 1` is sum of 2 squared integer.
    """
    def isConference(n):
        if n % 4 != 2:
            print("non conference")
            return
        dim - 1

    isConference(dim)
    modulo = dim - 1
    value = 1 / np.sqrt(modulo)

    quadraticResidue = set([np.mod(i * i, modulo) for i in range(1, modulo)])
    if 0 in quadraticResidue:
        quadraticResidue.remove(0)

    legendreSymbol = [
        0 if i == 0 else value if i in quadraticResidue else -value for i in range(modulo)
    ]

    # print(quadraticResidue)
    # print(legendreSymbol)

    mat = np.empty((dim, dim))
    mat[0][0] = 0
    for i in range(1, dim):
        mat[0][i] = value
        mat[i][0] = value

    for shift in range(modulo):
        mat[shift + 1][1:] = np.roll(legendreSymbol, shift)
    return mat

def weighing_matrix(dim=4):
    """
    W(n, n - 1) is conference matrix.

    https://en.wikipedia.org/wiki/Weighing_matrix
    """
    pass

def testMatrix():
    # mat = haar_measure(4)
    # mat = random_ortho(4, None)
    # mat = random_special_ortho(4)
    # mat = random_ortho_close_to_identity()
    # mat = random_ortho_circulant(4, None)
    # mat = random_triangular(4, None)
    # mat = random_schroeder(4, None)
    # mat = random_absorbent(4, None)
    mat = hadamard_sylvester(4)
    # mat = random_conference(62)

    det_mat = np.linalg.det(mat)
    sum_eig_mat = np.sum(linalg.eig(mat)[0])

    inv = np.linalg.inv(mat)
    det_inv = np.linalg.det(inv)
    sum_eig_inv = np.sum(linalg.eig(inv)[0])

    dotted = mat.dot(mat.conj().T)

    print(
        f"--- matrix\ndet: {det_mat}\nsum(eig): {sum_eig_mat}\ndet(inv): {det_inv}\neig(inv): {sum_eig_inv}\n",
        mat,
        "--- dotted",
        dotted,
        sep="\n\n",
    )

    cmap = plt.get_cmap("magma")
    for idx, row in enumerate(dotted):
        plt.plot(np.abs(row), color=cmap(idx / row.shape[0]), label=str(idx))
    plt.show()

def testRotation():
    # mat = random_ortho_close_to_identity()
    # mat /= np.linalg.norm(mat)

    mat = random_ortho_circulant()
    eye = np.ones(4)

    print(mat, end="\n\n\n\n")

    # dotted = mat.dot(mat.conj().T)
    # cmap = plt.get_cmap("magma")
    # for idx, row in enumerate(dotted):
    #     plt.plot(np.abs(row), color=cmap(idx / row.shape[0]), label=str(idx))
    # plt.show()

    for _ in range(16):
        print(eye, end="\n\n")
        eye = eye * mat
        # eye /= np.sqrt(np.sum(np.multiply(eye, eye)))
        eye /= np.linalg.norm(eye)

        dotted = eye.dot(eye.conj().T)
        cmap = plt.get_cmap("magma")
        for idx, row in enumerate(dotted):
            plt.plot(np.abs(row), color=cmap(idx / row.shape[0]), label=str(idx))
        plt.show()

def testSort():
    mat = random_ortho(4)

    axis = 0
    indices = np.argmax(np.abs(mat), axis=axis)
    st = np.take(mat, indices, axis)

    print(mat, indices, st, sep="\n\n")

    dotted = st.dot(st.conj().T)
    cmap = plt.get_cmap("magma")
    for idx, row in enumerate(dotted):
        plt.plot(np.abs(row), color=cmap(idx / row.shape[0]), label=str(idx))
    plt.show()

if __name__ == "__main__":
    testMatrix()
    # testRotation()
    # testSort()
