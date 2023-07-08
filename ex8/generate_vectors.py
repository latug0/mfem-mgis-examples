import numpy as np

def generate_orthonormal_vectors(num_vectors):
    vectors = []
    for _ in range(num_vectors):
        v = np.random.randn(3)
        if len(vectors) > 0:
            for existing_vector in vectors:
                v -= np.dot(v, existing_vector) * np.array(existing_vector)
        v /= np.linalg.norm(v)
        vectors.append(v.tolist())
    return vectors

def check_orthonormality(vectors):
    passed = True
    for i in range(len(vectors)):
        for j in range(i+1, len(vectors)):
            dot_product = np.dot(vectors[i], vectors[j])
            if abs(dot_product) > 1e-16:
                passed = False
                print("non-orthogonal vectors: ", abs(dot_product))
                break
    return passed

def check_norm(vectors):
    passed = True
    for v in vectors:
        norm_v = np.linalg.norm([float(x) for x in v])
        if abs(norm_v - 1.0) > 1e-16:
            passed = False
            print("non-unit vector norm : ", norm_v)
            break
    return passed

vectors = []
while (len(vectors) < 16):
    # Utilisation du code pour générer 16 vecteurs orthonormaux avec une précision de 16 chiffres significatifs
    orthonormal_vectors = generate_orthonormal_vectors(2)
    if (check_norm(orthonormal_vectors)):
        if check_orthonormality(orthonormal_vectors):
            vectors.append(orthonormal_vectors[0])
            vectors.append(orthonormal_vectors[1])

# Affichage des vecteurs
with open("vecteurs.txt", "w") as file:
    for i, v in enumerate(vectors, start=1):
        line = f"{' '.join(format(x, '.16f') for x in v)}"
        file.write(line + "\n")

# Afficher un message de confirmation
print("Les vecteurs ont été exportés dans le fichier vecteurs.txt.")