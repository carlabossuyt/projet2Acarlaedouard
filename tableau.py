from pylab import *
import pandas as pd
import matplotlib.pyplot as plt


def remplir(sequencea, sequenceb):
    gap_score = -2
    match = 2
    mismatch = -1
    grille = np.zeros((len(sequencea) + 1, len(sequenceb) + 1))

    for i in range(len(sequencea) + 1):
        grille[i][0] = i * gap_score

    for j in range(len(sequenceb) + 1):
        grille[0][j] = j * gap_score


def affiche(sequencea, sequenceb):
    plt.axis('tight')
    plt.axis('off')
    plt.table(cellText=grille, colLabels=sequencea, rowLabels=sequenceb, loc="center")
    plt.show()


if __name__ == "__main__":
    sequencea = ["A", "C", "G", "G", "C", "T", "A", "T"]
    sequenceb = ["A", "C", "T", "G", "T", "A", "G"]
    remplir(sequencea, sequenceb)

