import numpy as np
import matplotlib.pyplot as plt


def remplir(sequence1, sequence2):
    # initialisation des scores
    addition = -2
    identique = 2
    substitution = -1
    # remplissage de la grille de 0
    grille = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
    # initialisation de la 2e ligne et la 2e colonne
    for i in range(len(sequence1) + 1):
        grille[i][0] = i * addition

    for j in range(len(sequence2) + 1):
        grille[0][j] = j * addition

    for i in range(1, len(sequence1) + 1):
        for j in range(1, len(sequence2) + 1):
            if sequence1[i - 1] == sequence2[j - 1]:
                nv_match = grille[i - 1][j - 1] + identique
            else:
                nv_match = grille[i - 1][j - 1] + substitution
            nv_deletion = grille[i - 1][j] + addition
            nv_insertion = grille[i][j - 1] + addition
            grille[i][j] = max(nv_match, nv_insertion, nv_deletion)
    return grille


def affiche(grille, sequence1, sequence2):
    fig, ax = plt.subplots()
    ax.set_axis_off()
    ax.table(cellText=grille, rowLabels=sequence1, colLabels=sequence2, loc="center")
    plt.show()


if __name__ == "__main__":
    sequencea = "ACGGCTAT"
    sequenceb = "ACTGTAG"
    grille1 = remplir(sequencea, sequenceb)
    print(np.shape(grille1))
    print(len(sequenceb))
    affiche(grille1, sequencea, sequenceb)
