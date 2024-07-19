"""
Author: maxence-wz
Date: 2023-10-06
Description: This module contains functions to generate random vectors.
             The main function generates a specified number of vectors,
             orthonormal 2 by 2 and exports them to a text file `vectors.txt`.
"""

from tools.generateVectors import generate_orthonormal_vectors, generate_orthonormal_vectors_from_angles
from tools.tools import check_norm, check_orthogonality

def main(nbVec:int, algoKey:str='default'):
    """
    Generates `nbVec` orthonormal vectors 2 by 2 and exports them to a text file `vectors.txt`.

    :param nbVec: The number of vectors to generate.
    :param algo: The algorithm to use for generating vectors ('angles' or 'default').
    """
    vectors = []
    dictAlgo = {
        'angles': generate_orthonormal_vectors_from_angles,
        'default': generate_orthonormal_vectors
    }
    funcAlgo = dictAlgo.get(algoKey)
    # Generates vectors
    while len(vectors) < nbVec:
        orthonormal_vectors = funcAlgo(2)
        if check_norm(orthonormal_vectors):
            if check_orthogonality(orthonormal_vectors):
                vectors.append(orthonormal_vectors[0])
                vectors.append(orthonormal_vectors[1])
    # Exports vectors
    with open("vectors.txt", "w") as file:
        for i, v in enumerate(vectors, start=1):
            line = f"{' '.join(format(x, '.32f') for x in v)}"
            file.write(line + "\n")

    # Display
    print(f"{nbVec} vectors have been exported to the file vectors.txt.")


if __name__ == '__main__':
    main(nbVec=16,algoKey='angles')

    vec1 = [0.71287845123819126857966921306797,
            0.40877656662388073272040855954401,
            0.56982982752698174699901301210048]
    vec2 = [-0.62079290640125961431294854264706,
            0.74581447521998256444675234888564,
            0.24161319482639823097436249099701]
    check_norm([vec1,vec2])
    check_orthogonality([vec1,vec2])