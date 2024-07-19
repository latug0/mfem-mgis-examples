import numpy as np


def check_orthogonality(vectors):
    """
        Checks if the given vectors are orthogonal to each other.

        The function calculates the dot product between each pair of vectors and
        verifies that it is within a small tolerance to zero, indicating orthogonality.

        :param vectors: A list of vectors to be checked for orthogonality.
        :return: True if all vectors are orthogonal to each other, False otherwise.

        Test:
        >>> vectors = [
        ...     [1, 0, 0],
        ...     [0, 1, 0],
        ...     [0, 0, 1]
        ... ]
        >>> check_orthogonality(vectors)
        True

        >>> non_orthogonal_vectors = [
        ...     [1, 1, 0],
        ...     [2, 2, 0]
        ... ]
        >>> check_orthogonality(non_orthogonal_vectors)
        non-orthogonal vectors:  4
        False
        """
    passed = True
    for i in range(len(vectors)):
        for j in range(i + 1, len(vectors)):
            dot_product = np.dot(vectors[i], vectors[j])
            if abs(dot_product) > 1e-16:
                passed = False
                print("non-orthogonal vectors: ", abs(dot_product))
                break
    return passed


def check_norm(vectors):
    """
        Checks if the given vectors have a norm of 1 (i.e., are unit vectors).

        The function calculates the norm of each vector and verifies that it is
        close to 1 within a small tolerance. If any vector does not meet this
        criterion, the function prints the norm of the non-unit vector.

        :param vectors: A list of vectors to be checked for unit norm.
        :return: True if all vectors have a norm of 1, False otherwise.

        Test:
        >>> vectors = [
        ...     [1, 0, 0],
        ...     [0, 1, 0],
        ...     [0, 0, 1]
        ... ]
        >>> check_norm(vectors)
        True

        >>> non_unit_vectors = [
        ...     [1, 1, 0],  # Norm is sqrt(2)
        ...     [0, 1, 0]
        ... ]
        >>> check_norm(non_unit_vectors)
        non-unit vector norm :  1.4142135623730951
        False
        """
    passed = True
    for v in vectors:
        norm_v = np.linalg.norm([float(x) for x in v])
        if abs(norm_v - 1.0) > 1e-16:
            passed = False
            print("non-unit vector norm : ", norm_v)
            break
    return passed


if __name__ == '__main__':
    vectors = [[0, 1, 0], [0, 1, 0], [0, 0, 2]]
    check_norm(vectors)
    check_orthogonality(vectors)
