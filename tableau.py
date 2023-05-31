from pylab import *
import numpy as np
import matplotlib.pyplot as plt


def remplir(sequencea, sequenceb):
      #initialisation des scores
      addition = -2
      identique = 2
      substitution = -1
      #remplissage de la grille de 0
      grille = np.zeros((len(sequencea) + 1, len(sequenceb) + 1))
      #initialisation de la 2e ligne et la 2e colonne
      for i in range(len(sequencea)+1):
            grille[i][0] = i * addition

      for j in range(len(sequenceb)+1):
            grille[0][j] = j * addition

      for i in range(1, len(sequencea) + 1):
            for j in range(1, len(sequenceb) + 1):
                  if sequencea[i-1] == sequenceb[j-1]:
                        return grille[i-1][j-1] + identique
                  else:
                        return grille[i-1][j-1] + substitution
                  return grille[i-1][j] + addition
                  return grille[i][j-1] + addition
      return grille


def affiche(grille,sequencea, sequenceb):
      fig, ax = plt.subplots()
      ax.set_axis_off()
      ax.table(cellText=grille, rowLabels=sequenceb, colLabels=sequencea, loc="center")
      plt.show()


if __name__ == "__main__":
      sequencea = "A C G G C T A T"
      sequenceb = "A C T G T A G"
      grille1 = remplir(sequencea, sequenceb)
      affiche(grille1,sequencea,sequenceb)
