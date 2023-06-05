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
    :param cheminoptimal: list
    """
    #creation de la figure et des axes
    fig, ax = plt.subplots()
    #suppression des axes
    ax.set_axis_off()
    #creation d'une table avec la grille d'alignement , les séquences en tant que labels de lignes et de colonnes
    table = ax.table(cellText=grille, rowLabels="-" + sequence1, colLabels="-" + sequence2, loc="center")
    for x,y in cheminoptimal:
        cellule = table.get_celld()[x+1, y] #table.cell pour accéder à une cellule spécifique dans un tableau
        cellule.set_facecolor("pink")
    #affichage de la figure
    plt.show()
    print(grille)

def remonter(grille, sequence1, sequence2):
    """"
    Remonte le chemin optimal danns la grille
    :param grille: np.array
    :param sequence1: str
    :param sequence2: str
    :return: list"""
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
    y = len(grille[0])-1 #on parcourt les colonnes dans une ligne fixe
    parcours=[(x,y)]
    while x!=1 and y!=1:
        antecedents=[((x,y-1), grille[x][y-1], 1), #antecedant haut = position 1
                     ((x-1, y-1), grille[x-1][y-1], 2), #antecedant diago = position 2
                     ((x-1, y), grille[x-1][y], 3)] ##antecedant gauche = position 3
        antecedents.sort(key=operator.itemgetter(1), reverse=True) #permet de trier et sélectionner le plus grand antécédent
        for cord_ant_max, score_ant_max, position in antecedents:
            max = score_ant_max
            if verif():
                parcours.append(cord_ant_max) #ajoute les coordonnées à la liste parcours
                break #retire le maximum si le score n'est pas bon
        x, y =parcours[-1]
    return parcours

def blosum(sequence1, sequence2):
    """
    Calcul les scores Blosum entre 2 séquences
    :param sequence1: str
    :param sequence2: str
    :return:
    """
    li_bl = []
    f = open("blosum62.txt", "r")
    first = None
    for ligne in f.readlines():
        if "#" not in ligne:
            if first is None: #permet de récupérer la 1ère ligne
                first = ligne.split()
            else:
                decoupe = ligne.split()
                li_bl.append(decoupe[1:]) #permet d'enlever le 1er caractère de chaque ligne
    print(first)
    for l in li_bl:
        print(l) #afficher les lignes

    indexes = {}
    # {"A": 0, "R": 1, etc}
    i = 0
    for lettre in first:
        indexes[lettre] = i
        i += 1

    blossum_values = {}
    # {("R", "Q")): 1, ...}
    for l1 in first:
        for l2 in first:
            blossum_values[(l1, l2)] = li_bl[indexes[l1]][indexes[l2]]

    print(blossum_values[("A", "D")])

    #n_colonne = first.index(sequence1[5])
    #n_ligne = first.index(sequence2[4])
    #print(str("le score blosum de la colonne"), n_colonne, str("et de la ligne"), n_ligne, str("est"), li_bl[n_colonne][n_ligne])

def adn2arn(sequence):
    """
    Convertit une séqeunce d'ADN en ARN
    :param sequence: str
    :return: str
    """
    #initialiser le résultat (ARN)
    arn = ""
    #parcours de la chaine d'ADN
    for b in sequence:
        if b == "T":
            arn = arn + "U"
        else:
            arn += b
    return arn


def arn2protein(arn):
    """
    Convertit une séquence d'ARN en séquence d'acides aminés/ en protéine
    :param arn: str
    :return: str
    """
    table = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I','GUU': 'V', 'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V', 'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V', 'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V','UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A', 'UCC': 'S', 'CCC': 'P', 'ACC':'T', 'GCC': 'A', 'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'UCG': 'S','CCG': 'P', 'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU':'D', 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 'UAA': 'Stop', 'CAA': 'Q','AAA': 'K', 'GAA': 'E', 'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 'UGU':'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R', 'AGC': 'S','GGC': 'G', 'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', 'UGG': 'W', 'CGG':'R', 'AGG': 'R', 'GGG': 'G'}
    res= ''
    for i in range(0, len(arn), 3):
        code = arn[i:i+3]
        if code in table:
            if table[code] == "Stop":
                return res
            else:
                res += table[code]
        else:
            return res
    return res



if __name__ == "__main__":
    sequencea = "ACTGTAG"
    sequenceb = "ACGGCTAT"
    grille1 = remplir(sequencea, sequenceb) #remplissage de la grille d'alignement
    cheminoptimal = remonter(grille1, sequencea, sequenceb)
    affiche(grille1, sequencea, sequenceb, cheminoptimal) #affichage de la grille d'alignement
    print(remonter(grille1, sequencea, sequenceb))
    blosum(sequencea, sequenceb)

    print(sequencea, str("devient l'arn1: "), adn2arn(sequencea))
    print(sequenceb, str("devient l'arn2: "), adn2arn(sequenceb))

    protein1=adn2arn(sequencea)
    protein2=adn2arn(sequenceb)
    print(str("l'arn1 devient la protéine1: "), arn2protein(protein1))
    print(str("l'arn2 devient la protéine2: "), arn2protein(protein2))
