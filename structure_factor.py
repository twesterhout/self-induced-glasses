import numpy as np
import h5py

for β in [0.1, 0.2, 0.3, 0.5]:
    with h5py.File("data/structure_n=20_λ=2.5_β={}_seed=569485.h5".format(β), "r") as f:
        S = np.asarray(f["data"])
    print("[*] Processing β={} ...".format(β))
    n = S.shape[0]
    L = int(np.sqrt(n))
    assert L * L == n
    x = np.arange(n) % L
    y = np.arange(n) // L
    delta_x = x.reshape(-1, 1) - x.reshape(1, -1)
    delta_x = np.minimum(delta_x, L + delta_x)
    delta_y = y.reshape(-1, 1) - y.reshape(1, -1)
    delta_y = np.minimum(delta_y, L + delta_y)
    qs = np.linspace(-2 * np.pi / L, 2 * np.pi / L, 31)

    out = np.zeros((len(qs), len(qs)))
    for q_i, q_x in enumerate(qs):
        for q_j, q_y in enumerate(qs):
            out[q_i, q_j] = np.dot(S.flatten(), np.cos(q_x * delta_x + q_y * delta_y).flatten())

    np.savetxt("S(R)_β={}.dat".format(β), S[0].reshape(L, L))
    np.savetxt("S(Q)_β={}.dat".format(β), out)
