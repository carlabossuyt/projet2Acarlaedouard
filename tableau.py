import numpy as np
import matplotlib.pyplot as plt
import operator

def remplir(sequence1, sequence2):
    """
    Remplissage de la grille
    :param sequence1: str
    :param sequence2: str
    :return: grille
    """
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


def affiche(grille, sequence1, sequence2, cheminoptimal):
    """
    Affichage de la grille
    :param grille: np.array
    :param sequence1: str
    :param sequence2: str
    """
    #creation de la figure et des axes
    fig, ax = plt.subplots()
    #suppression des axes
    ax.set_axis_off()
    #creation d'une table avec la grille d'alignement , les séquences en tant que labels de lignes et de colonnes
    table = ax.table(cellText=grille, rowLabels="-" + sequence1, colLabels="-" + sequence2, loc="center")
    for x,y in cheminoptimal:
        cellule = table.get_celld()[x+1, y]
        cellule.set_facecolor("pink")
    #affichage de la figure
    plt.show()
    print(grille)

def remonter(grille, sequence1, sequence2):
    def verif():
        i, j= parcours[-1] #recupere les coordonnées de la dernière case
        score=grille[x][y]
        if position == 1 and ((max - 2) == score): #position haut
            return True
        elif position == 3 and ((max-2)==score): #position gauche
            return True
        elif position == 2 and sequence1[x-1] == sequence2[y-1] and max +2 == score: #match
            return True
        elif position == 2 and sequence1[x-1] != sequence2[y-1] and max -1 == score: #dismatch
            return True
        else:
            return False
    x= len(grille)-1
    y = len(grille[0])-1 #on parcours les colonnes dans une ligne fixe
    parcours=[(x,y)]
    while x!=1 and y!=1:
        antecedents=[((x,y-1), grille[x][y-1], 1), #antecedant haut = position 1
                     ((x-1, y-1), grille[x-1][y-1], 2), #antecedant diago = position 2
                     ((x-1, y), grille[x-1][y], 3)] ##antecedant gauche = position 3
        antecedents.sort(key=operator.itemgetter(1), reverse=True) #permet de trier et sélectionner le + grand antecedents
        for cord_ant_max, score_ant_max, position in antecedents:
            max = score_ant_max
            if verif():
                parcours.append(cord_ant_max) #ajoute les coordonnées à la liste parcours
                break #retire le maximum si il le score n'est pas bon
        x, y =parcours[-1]
    return parcours



if __name__ == "__main__":
    sequencea = "ACTGTAG"
    sequenceb = "ACGGCTAT"
    grille1 = remplir(sequencea, sequenceb) #remplissage de la grille d'alignement
    cheminoptimal = remonter(grille1, sequencea, sequenceb)
    affiche(grille1, sequencea, sequenceb, cheminoptimal) #affichage de la grille d'alignement


    print(remonter(grille1, sequencea, sequenceb))
