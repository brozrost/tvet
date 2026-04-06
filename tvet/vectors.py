import numpy as np

def normalize_vectors(vectors: np.ndarray) -> np.ndarray:
    vectors = np.asarray(vectors, dtype=np.double)

    if vectors.ndim == 1 and vectors.shape[0] == 3:
        n = np.linalg.norm(vectors)
        return vectors / (n if n > 0.0 else 1.0)

    if vectors.ndim != 2 or vectors.shape[1] != 3:
        raise ValueError(f"Expected shape (N,3) or (3,), got {vectors.shape}")
    
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    norms = np.where(norms > 0.0, norms, 1.0)

    return vectors / norms

def rotate_x(a, phi):
    x, y, z = a

    x_ = x
    y_ = y * np.cos(phi) - z * np.sin(phi)
    z_ = y * np.sin(phi) + z * np.cos(phi)

    a_ = np.array([x_, y_, z_], dtype=np.double)

    return a_

def rotate_y(a, phi):
    x, y, z = a

    x_ = x * np.cos(phi) + z * np.sin(phi)
    y_ = y
    z_ = - x * np.sin(phi) + z * np.cos(phi)

    a_ = np.array([x_, y_, z_], dtype=np.double)

    return a_

def rotate_z(a, phi):
    x, y, z = a

    x_ = x * np.cos(phi) - y * np.sin(phi)
    y_ = x * np.sin(phi) + y * np.cos(phi)
    z_ = z

    a_ = np.array([x_, y_, z_], dtype=np.double)

    return a_