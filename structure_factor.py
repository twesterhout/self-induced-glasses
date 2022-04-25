import numpy as np
import h5py
import glob

def fourier_transform_structure_factor(S, number_points: int):
    number_bits = S.shape[0]
    length = int(np.sqrt(number_bits))
    assert length * length == number_bits
    x = np.arange(number_bits) % length
    y = np.arange(number_bits) // length
    delta_x = x.reshape(-1, 1) - x.reshape(1, -1)
    delta_x = np.minimum(delta_x, length + delta_x)
    delta_y = y.reshape(-1, 1) - y.reshape(1, -1)
    delta_y = np.minimum(delta_y, length + delta_y)

    qs = np.linspace(-4 * np.pi, 4 * np.pi, number_points)

    out = np.zeros((len(qs), len(qs)))
    for q_i, q_x in enumerate(qs):
        for q_j, q_y in enumerate(qs):
            out[q_i, q_j] = np.dot(S.reshape(-1), np.cos(q_x * delta_x + q_y * delta_y).reshape(-1))
    return out

for filename in glob.glob("data/mean/structure_n=25_λ=2.5_β=*.h5"):
    print("[*] Processing {} ...".format(filename))
    with h5py.File(filename, "r") as f:
        Sr = np.asarray(f["data"])
    Sq = fourier_transform_structure_factor(Sr, number_points=51)
    np.savetxt(filename.replace(".h5", ".dat").replace("structure_", "structure_fourier_"), Sq)
