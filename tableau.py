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
    Remonte le chemin optimal dans la grille
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


def align(parcours, sequence1, sequence2):
    alignement1 = ""
    nv_sequence1 = ""
    nv_sequence2 = ""
    parcours.reverse()
    print(parcours)
    for i in range(1, len(parcours)):
        if parcours[i][0] == parcours[i-1][0] + 1 and parcours[i][1] == parcours[i-1][1] + 1:
            alignement1 += "|"
            nv_sequence1 += sequence1[parcours[i][0]-1] #0 : abscisse du mot 1
            nv_sequence2 += sequence2[parcours[i][1]-1] #-1 = pour le décalage
        elif parcours[i][0] == parcours[i][0] + 1 and parcours[i][1] == parcours[i-1][1]:
            nv_sequence1 += "-"
            alignement1 += " "
            nv_sequence2 += sequence2[parcours[i][1] - 1]
        else:
            alignement1 += " "
            nv_sequence2 += "-"
            nv_sequence1 += sequence1[parcours[i][0] - 1]
    return nv_sequence1, alignement1, nv_sequence2


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
            blossum_values[(l1, l2)] = int(li_bl[indexes[l1]][indexes[l2]])

    print(blossum_values)
    return blossum_values
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

def remplir_blosum(blosum, sequence1, sequence2):
    """
    Remplissage de la grille
    :param sequence1: str
    :param sequence2: str
    :return: grille
    """
    # remplissage de la grille de 0
    grille2 = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
    # initialisation de la 2e ligne et la 2e colonne
    print(grille2)
    for i in range(1,len(sequence1)+1):
        grille2[i][0] = grille2[i-1][0] -4

    for j in range(1,len(sequence2)+1):
        grille2[0][j] = grille2[0][j-1] -4
    for i in range(1, len(sequence1) + 1):
        for j in range(1, len(sequence2) + 1):
            # print(i, j, grille[i - 1][j - 1], sequence1[i - 1], sequence2[j - 1], blosum[sequence1[i - 1], sequence2[j - 1]])
            nv_match = grille2[i - 1][j - 1] + blosum[sequence1[i - 1], sequence2[j - 1]]
            nv_deletion = grille2[i - 1][j] -4
            nv_insertion = grille2[i][j - 1] -4
            grille2[i][j] = max(nv_match, nv_insertion, nv_deletion)
    return grille2



def remonter_blosum(blosum_dict, grille, sequence1, sequence2):
    """"
    Remonte le chemin optimal dans la grille
    :param grille: np.array
    :param sequence1: str
    :param sequence2: str
    :return: list"""
    def verif():
        i, j= parcours[-1] #recupere les coordonnées de la dernière case
        score=grille[x][y]
        #a, b, c = sequence1[i - 1], sequence2[j - 1],  blosum_dict[sequence1[i - 1], sequence2[j - 1]]
        if position == 1 and ((max - 4) == score): #position haut
            return True
        elif position == 3 and ((max-4) == score): #position gauche
            return True
        elif position == 2 and max + blosum_dict[sequence1[i - 1], sequence2[j - 1]] == score:
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
        print(x, y, antecedents)
        for cord_ant_max, score_ant_max, position in antecedents:
            max = score_ant_max
            if verif():
                parcours.append(cord_ant_max) #ajoute les coordonnées à la liste parcours
                break #retire le maximum si le score n'est pas bon
        x, y =parcours[-1]
    return parcours




if __name__ == "__main__":
    sequencea = "ACTGTAG"
    sequenceb = "ACGGCTAT"
    grille1 = remplir(sequencea, sequenceb) #remplissage de la grille d'alignement
    cheminoptimal = remonter(grille1, sequencea, sequenceb)
    affiche(grille1, sequencea, sequenceb, cheminoptimal) #affichage de la grille d'alignement
    parcours = remonter(grille1, sequencea, sequenceb)
    alignement = align(parcours, sequencea, sequenceb)
    print(alignement[0])
    print(alignement[1])
    print(alignement[2])


    print(sequencea, str("devient l'arn1: "), adn2arn(sequencea))
    print(sequenceb, str("devient l'arn2: "), adn2arn(sequenceb))

    protein1=adn2arn(sequencea)
    protein2=adn2arn(sequenceb)
    print(str("l'arn1 devient la protéine1: "), arn2protein(protein1))
    print(str("l'arn2 devient la protéine2: "), arn2protein(protein2))



    sequencea2="MKTKIAEYLKALLKNTEKYL"
    sequenceb2="VIENEIAYIKDPVFGIVRNRVSA"
    blosum_dict = blosum(sequencea2, sequenceb2)
    grille2 = remplir_blosum(blosum_dict, sequencea2, sequenceb2)
    print("1", grille2)
    cheminoptimal2 = remonter_blosum(blosum_dict, grille2, sequencea2, sequenceb2)
    print(2, cheminoptimal2)
    affiche(grille2, sequencea2, sequenceb2, cheminoptimal2)
    print(3)
    parcours2 = remonter_blosum(blosum_dict, grille2, sequencea2, sequenceb2)
