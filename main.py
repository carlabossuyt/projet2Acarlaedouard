from tableau import *

if __name__ == "__main__":
    #ADN
    sequencea = "ACTGTAGTCAGT"
    sequenceb = "ACGGCTATGCCTAGC"
    grille1 = remplir(sequencea, sequenceb) #remplissage de la grille d'alignement
    cheminoptimal = remonter(grille1, sequencea, sequenceb)
    affiche(grille1, sequencea, sequenceb, cheminoptimal) #affichage de la grille d'alignement
    parcours = remonter(grille1, sequencea, sequenceb)
    alignement = align(parcours, sequencea, sequenceb)

    #alignement ADN:
    print(alignement[2])
    print(alignement[1])
    print(alignement[0])

    #protéines
    sequencea2="MKTKIAEYLKALLKNTEKYL"
    sequenceb2="VIENEIAYIKDPVFGIVRNRVSA"
    blosum_dict = blosum(sequencea2, sequenceb2)
    grille2 = remplir_blosum(blosum_dict, sequencea2, sequenceb2)
    #print("1", grille2)
    cheminoptimal2 = remonter_blosum(blosum_dict, grille2, sequencea2, sequenceb2)
    #print(2, cheminoptimal2)
    affiche(grille2, sequencea2, sequenceb2, cheminoptimal2)
    #print(3)
    parcours2 = remonter_blosum(blosum_dict, grille2, sequencea2, sequenceb2)

    #alignement protéine
    alignement2 = align_blosum(parcours2, sequencea2, sequenceb2)
    print(alignement2[0])
    print(alignement2[1])
    print(alignement2[2])

    #ADN en ARN et ARN en protéine
    print(str("L'ADN1: "), sequencea, str("devient l'arn1: "), adn2arn(sequencea))
    print(str("L'ADN2: "), sequenceb, str("devient l'arn2: "), adn2arn(sequenceb))
    protein1=adn2arn(sequencea)
    protein2=adn2arn(sequenceb)
    print(str("l'arn1 devient la protéine1: "), arn2protein(protein1))
    print(str("l'arn2 devient la protéine2: "), arn2protein(protein2))
