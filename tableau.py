import numpy as np
import matplotlib.pyplot as plt
import operator


def affiche(grille, sequence1, sequence2, cheminoptimal):
    """
    Affichage de la grille
    :param grille: np.array, grille d'alignement
    :param sequence1: str, 1ère séquence
    :param sequence2: str, 2ème séquence
    :param cheminoptimal: list, chemin optimal de la grille
    """
    fig, ax = plt.subplots() #creation de la figure et des axes
    ax.set_axis_off() #suppression des axes
    #creation d'une table avec la grille d'alignement , séquences = labels de lignes et de colonnes
    table = ax.table(cellText=grille, rowLabels="-" + sequence1, colLabels="-" + sequence2, loc="center")
    for x,y in cheminoptimal:
        cellule = table.get_celld()[x+1, y] #table.cell pour accéder à une cellule spécifique dans un tableau
        cellule.set_facecolor("pink") #cellule du chemin optimal en couleur
    plt.show() #affichage de la figure
    #print(grille)


def remplir(sequence1, sequence2):
    """
    Remplissage de la grille
    :param sequence1: str, 1ère séquence
    :param sequence2: str, 2ème séquence
    :return: grille, grille remplie
    """
    # initialisation des scores
    addition = -2
    identique = 2
    substitution = -1
    # remplissage de la grille de 0
    grille = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
    # initialisation de la 2ème ligne et la 2ème colonne
    for i in range(len(sequence1) + 1):
        grille[i][0] = i * addition
    for j in range(len(sequence2) + 1):
        grille[0][j] = j * addition
    #remplissage de la grille avec les scores maximaux
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


def remonter(grille, sequence1, sequence2):
    """"
    Remonte le chemin optimal dans la grille
    :param grille: np.array, grille remplie
    :param sequence1: str, 1ère séquence
    :param sequence2: str, 2ème séquence
    :return: list, chemin optimal de la grille
    """
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
    parcours.append((0,0))
    return parcours


def align(parcours, sequence1, sequence2):
    """
    Alignement des 2 séquences d'ADN grâce à la liste parcours
    :param parcours: list, liste des positions des acides aminés
    :param sequence1: str, 1ère séquence
    :param sequence2: str, 2ème séquence
    :return: tuple, tuple contenant la nouvelle séquence 1, la nouvelle séquence 2 et l'alignement
    """
    alignement1 = ""
    nv_sequence1 = ""
    nv_sequence2 = ""
    parcours.reverse() #inversion de l'ordre des positions parcourues pour avoir l'alignement à l'endroit
    #print(parcours)
    for i in range(1, len(parcours)): #parcours des positions
        if parcours[i][0] == parcours[i-1][0] + 1 and parcours[i][1] == parcours[i-1][1] + 1: #si les positions se succèdent horizontalement et verticalement
            if sequence1[parcours[i][0]-1] == sequence2[parcours[i][1]-1]: # si les 2 acides aminés sont identiques
                alignement1 += "|"
            else:
                alignement1 += " "
            nv_sequence1 += sequence1[parcours[i][0]-1] #ajout de le l'acide aminé à la nouvelle séquence 1, 0 : abscisse du mot 1
            nv_sequence2 += sequence2[parcours[i][1]-1] #ajout de le l'acide aminé à la nouvelle séquence 2, -1 = pour le décalage
        elif parcours[i][0] == parcours[i-1][0] + 1 and parcours[i][1] == parcours[i-1][1]: #si les positions se succèdent horizontalement
            nv_sequence2 += "-" #ajouter d'un tiret à la nouvelle séquence 2
            alignement1 += " "
            nv_sequence1 += sequence1[parcours[i][0] - 1] #ajout de l'acide aminé à la nouvelle séquence1
        else: #si les positions se succèdent verticalement
            alignement1 += " "
            nv_sequence1 += "-"
            nv_sequence2 += sequence2[parcours[i][1] - 1]
    return nv_sequence1, alignement1, nv_sequence2 #tuple retournant les 2 nouvelles séquences et l'alignement


def blosum(sequence1, sequence2):
    """
    Calcul les scores Blosum entre 2 séquences
    :param sequence1: str, la 1ère séquence
    :param sequence2: str, la 2ème séquence
    :return:dict, dictionnaire contenant les scores Blosum entre les acides aminés des 2 séquences
    """
    li_bl = [] #création d'une liste pour stocker les lignes du texte blosum62
    f = open("blosum62.txt", "r") #ouverture du fichier blosum62 en mode lecture
    first = None #variable qui stocke la 1ère ligne du fichier
    for ligne in f.readlines(): #parcours des lignes du fichier
        if "#" not in ligne: #exclure les lignes contenant les lignes avec "#"
            if first is None: #permet de récupérer la 1ère ligne
                first = ligne.split() #split = séparer les nombres
            else: #ligne contenant les scores blosum
                decoupe = ligne.split()
                li_bl.append(decoupe[1:]) #permet d'enlever le 1er caractère de chaque ligne
    #print(first)
    #for l in li_bl:
        #print(l) #afficher les lignes
    indexes = {} #dictionnaire qui stocke les index des acides aminés
    # {"A": 0, "R": 1, etc}
    i = 0
    for lettre in first:
        indexes[lettre] = i
        i += 1
    blossum_values = {} #dictionnaire qui stocke les scores blosum entre les acides aminés
    # {("R", "Q")): 1, ...}
    for l1 in first:
        for l2 in first:
            blossum_values[(l1, l2)] = int(li_bl[indexes[l1]][indexes[l2]])
    print(blossum_values)
    return blossum_values #retourne le dictionnaire des scores Blosum
    #n_colonne = first.index(sequence1[5])
    #n_ligne = first.index(sequence2[4])
    #print(str("le score blosum de la colonne"), n_colonne, str("et de la ligne"), n_ligne, str("est"), li_bl[n_colonne][n_ligne])


def remplir_blosum(blosum, sequence1, sequence2):
    """
    Remplissage de la grille
    :param sequence1: str, 1ère séquence
    :param sequence2: str, 2ème séquence
    :return: np.array, grille des scores remplie
    """
    # remplissage de la grille de 0
    grille2 = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
    # initialisation de la 2e ligne et la 2e colonne
    #print(grille2)
    #remplissage de la grille en calculant des scores différents pour chaque antécédants
    for i in range(1,len(sequence1)+1):
        grille2[i][0] = grille2[i-1][0] -4
    for j in range(1,len(sequence2)+1):
        grille2[0][j] = grille2[0][j-1] -4
    for i in range(1, len(sequence1) + 1):
        for j in range(1, len(sequence2) + 1):
            # print(i, j, grille[i - 1][j - 1], sequence1[i - 1], sequence2[j - 1], blosum[sequence1[i - 1], sequence2[j - 1]])
            nv_match = grille2[i - 1][j - 1] + blosum[sequence1[i - 1], sequence2[j - 1]] #score de l'anétécédent diagonale
            nv_deletion = grille2[i - 1][j] -4 #score de l'antécédant de gauche
            nv_insertion = grille2[i][j - 1] -4 #score de l'antécédant du haut
            grille2[i][j] = max(nv_match, nv_insertion, nv_deletion) #sélection du score maximal
    return grille2 #retours de la grille de score


def remonter_blosum(blosum_dict, grille, sequence1, sequence2):
    """"
    Remonte le chemin optimal dans la grille en partant du bas à droite
    :param blosum_dict: dict, matrice blosum contenant les scores pour chaque paires d'acides aminés
    :param grille: np.array, grille des scores remplie
    :param sequence1: str, 1ère séquence
    :param sequence2: str, 2ème séquence
    :return: list, liste des coordonnées parcourues dans le chemin optimal
    """
    def verif():
        """
        Verifie si la case actuelle est en accord avec la case antérieur
        :return: bool, True si le score de la case antécédente correspond à la case antérieure, False sinon
        """
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
        antecedents.sort(key=operator.itemgetter(1), reverse=True) #permet de trier les 3 scores des antécédents par score décroissant pour sélectionner le plus grand
        #print(x, y, antecedents)
        for cord_ant_max, score_ant_max, position in antecedents:
            max = score_ant_max
            if verif():
                parcours.append(cord_ant_max) #ajoute les coordonnées à la liste parcours
                break #retire le maximum si le score n'est pas bon
        x, y =parcours[-1] #met à jour les coordonnées de la case
    parcours.append((0,0))
    return parcours #retourne la liste des coordonnées parcourues dans le chemin optimal


def align_blosum(parcours, sequence1, sequence2):
    """
    Alignement des 2 séquences de protéines
    :param parcours: list, liste des positions des acides aminés
    :param sequence1: 1ère séquence
    :param sequence2: 2ème séquence
    :return: tuple: tuple contenant la nouvelle séquence 1, la nouvelle séquence 2 et l'aligneme
    """
    alignement1 = ""
    nv_sequence1 = ""
    nv_sequence2 = ""
    blos = blosum(sequence1, sequence2)
    parcours.reverse() #inversion de parcours
    print(parcours)
    for i in range(1, len(parcours)):
        x, y = parcours[i][0], parcours[i][1]
        if x == parcours[i - 1][0] + 1 and y == parcours[i - 1][1] + 1: #diagonale
            key = (sequence1[x-1], sequence2[y-1])
            score = blos[key]
            if sequence1[x - 1] == sequence2[y - 1]: #match
                alignement1 += "|"
            elif score > 0: #mismatch avec un score positif
                alignement1 += ":"
            else: #mismatch avec un score négatif
                alignement1 += "."
            nv_sequence1 += sequence1[x - 1]  # Ajout du caractère de la séquence 1 à la nouvelle séquence alignée
            nv_sequence2 += sequence2[y - 1] # Ajout du caractère de la séquence 2 à la nouvelle séquence alignée
        elif x == parcours[i - 1][0] + 1 and y == parcours[i - 1][1]: #gauche
            alignement1 += " "
            nv_sequence1 += sequence1[x - 1] # ajouter caractère de la séquence 1 à la nouvelle séquence 1
            nv_sequence2 += "-" # Ajout d'un gap "-" à la nouvelle séquence 2
        else: #haut
            alignement1 += " "
            nv_sequence1 += "-" # Ajout d'un gap "-" à la nouvelle séquence 1
            nv_sequence2 += sequence2[y - 1] # ajout du carcatère de la séquence 2 à la nouvelle séquence2
    return nv_sequence1, alignement1, nv_sequence2 #retourne les 2 nouvelles séquences alignées et l'alignement


def adn2arn(sequence):
    """
    Convertit une séqeunce d'ADN en ARN
    :param sequence: str, la séquence d'ADN à convertir en ARN
    :return: str, la séquence d'ADN résultante
    """
    #initialiser le résultat (ARN)
    arn = ""
    #parcours de la chaine d'ADN
    for b in sequence:
        if b == "T":
            #si la base est un "T", la remplacer par un "U" dans l'ARN
            arn = arn + "U"
        else:
            #si la base est différente de T, l'ajouter telle qu'elle est dans l'ARN
            arn += b
    #retourner la séquence d'ARN créée
    return arn


def arn2protein(arn):
    """
    Convertit une séquence d'ARN en protéine
    :param arn: str, séquence d'ARN à convetir en protéine
    :return: str, la protéine résultante
    """
    #table de correspondance entre les codons et les acides aminés
    table = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I','GUU': 'V', 'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V', 'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V', 'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V','UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A', 'UCC': 'S', 'CCC': 'P', 'ACC':'T', 'GCC': 'A', 'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'UCG': 'S','CCG': 'P', 'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU':'D', 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 'UAA': 'Stop', 'CAA': 'Q','AAA': 'K', 'GAA': 'E', 'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 'UGU':'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R', 'AGC': 'S','GGC': 'G', 'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', 'UGG': 'W', 'CGG':'R', 'AGG': 'R', 'GGG': 'G'}
    res= '' #permet de stocker la protéine crée
    #parcourir la séquence d'ARN en prenant 3 caractères = codon
    for i in range(0, len(arn), 3):
        code = arn[i:i+3] # extraire le codon sélectionné
        if code in table: #vérifier si codon est dans la table
            if table[code] == "Stop": #si codon est un Stop, arrêter et retourner la protéine
                return res
            else:
                res += table[code]
        else:
            return res #ajoute l'acide aminé correspondant au codon dans la protéine
    return res


