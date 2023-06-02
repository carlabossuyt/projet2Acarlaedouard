
import tableau

if __name__ == "__main__":
    sequencea = "ACTGTAG"
    sequenceb = "ACGGCTAT"
    grille1 = remplir(sequencea, sequenceb) #remplissage de la grille d'alignement
    affiche(grille1, sequencea, sequenceb) #affichage de la grille d'alignement
