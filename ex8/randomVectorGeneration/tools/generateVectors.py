import numpy as np

PI = np.pi

def generate_orthonormal_vectors(num_vectors):
    """
    Generates `num_vectors` orthonormal vectors using the Gram-Schmidt process.

    :param num_vectors: The number of orthonormal vectors to generate.
    :return: A list of orthonormal vectors.
    """
    vectors = []
    for _ in range(num_vectors):
        v = np.random.randn(3)
        if len(vectors) > 0:
            for existing_vector in vectors:
                v -= np.dot(v, existing_vector) * np.array(existing_vector)
        v /= np.linalg.norm(v)
        vectors.append(v.tolist())
    return vectors


def generate_orthonormal_vectors_from_angles(num_vectors):
    """
    Generates `num_vectors` orthonormal vectors from random angles.
    see: Annexe A of Approche micromécanique du comportement du combustible
    dioxyde d'uranium, Julian Soulacroix, 2014

    :param num_vectors: The number of orthonormal vectors to generate.
    :return: A list of orthonormal vectors.
    """
    vectors = []
    for _ in range(num_vectors):
        phi1, phi2 = np.random.uniform(0, PI, 2)
        psi = (lambda  x : np.acos(1-2*x))(np.random.uniform(0, 1))
        c1 = np.cos(phi1)
        c2 = np.cos(psi)
        c3 = np.cos(phi2)
        s1 = np.sin(phi1)
        s2 = np.sin(psi)
        s3 = np.sin(phi2)
        v = [c1*c3-c2*s1*s3, c3*s1+c1*c2*s3, s2*s3]
        if len(vectors) > 0:
            for existing_vector in vectors:
                v -= np.dot(v, existing_vector) * np.array(existing_vector)
        v /= np.linalg.norm(v)
        vectors.append(v.tolist())
    return vectors

if __name__ == '__main__':
    vectors = generate_orthonormal_vectors_from_angles(2)
    print(list(x for x in vectors))
    vectors = generate_orthonormal_vectors(2)
    print(list(x for x in vectors))
