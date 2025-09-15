"""
Travail pratique 1 : Conversion de fichiers en génétique
Auteurs: Kentz Dorcelus
Date: 15 Septembre 2025
Description: Programme pour convertir des fichiers génétiques du format complexe (.cmp) vers le format simplifié (.sim)
"""
import sys, re, os, tkinter as tk, turtle
from tkinter import ttk


class ConvertisseurFichierGenetique:
    """Classe principale pour la conversion de fichiers génétiques."""
    COMPLEMENTS = {'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    def __init__(self):
        self.nom_organisme = self.sequence_complete = ""
        self.genes, self.bases_totales, self.lignes_lues = [], 0, 0
    
    def complement_inverse(self, sequence):
        """Crée la séquence complémentaire inverse."""
        return ''.join(self.COMPLEMENTS.get(base, base) for base in sequence)[::-1]
    
    def parser_fichier_cmp(self, nom_fichier):
        """Parse un fichier .cmp et extrait les informations génétiques."""
        if not os.path.exists(nom_fichier):
            raise FileNotFoundError(f"Le fichier {nom_fichier} n'existe pas.")
        
        with open(nom_fichier, 'r') as fichier:
            lignes = fichier.readlines()
        
        self.lignes_lues = len(lignes)
        
        # Extraire le nom de l'organisme et le nombre de bases
        for ligne in lignes:
            ligne = ligne.strip()
            
            # Extraire le nombre de bases de la ligne LOCUS
            if ligne.startswith('LOCUS'):
                parties = ligne.split()
                for i, partie in enumerate(parties):
                    if partie.isdigit():
                        self.bases_totales = int(partie)
                        # Vérifier si le prochain élément est "bp"
                        if i + 1 < len(parties) and parties[i + 1].lower() == "bp":
                            break
            
            # Extraire le nom de l'organisme
            elif ligne.startswith('SOURCE'):
                # Le nom de l'organisme suit "SOURCE"
                correspondance_organisme = re.search(r'SOURCE\s+(.+?)(?:\s*\(|$)', ligne)
                if correspondance_organisme:
                    self.nom_organisme = correspondance_organisme.group(1).strip()
        
        # Construire la séquence complète à partir de la section ORIGIN
        origine_commencee = False
        sequence_temp = ""
        
        for ligne in lignes:
            ligne = ligne.strip()
            if ligne.startswith('ORIGIN'):
                origine_commencee = True
                continue
            elif ligne.startswith('//'):
                break
            elif origine_commencee and ligne:  # Vérifier que la ligne n'est pas vide
                # Enlever les numéros de ligne et espaces
                ligne_sequence = re.sub(r'^\s*\d+\s*', '', ligne)
                ligne_sequence = re.sub(r'\s+', '', ligne_sequence)
                sequence_temp += ligne_sequence.lower()
        
        # Vérifier que la séquence extraite correspond à la taille attendue
        if len(sequence_temp) > 0:
            self.sequence_complete = sequence_temp
            
            # Si la taille de la séquence ne correspond pas au nombre de bases indiqué
            if self.bases_totales > 0 and len(self.sequence_complete) != self.bases_totales:
                print(f"Attention: La taille de la séquence extraite ({len(self.sequence_complete)} bases) "
                      f"ne correspond pas à la taille indiquée ({self.bases_totales} bases)")
        
        # Extraire les gènes
        self.extraire_genes(lignes)
    
    def extraire_genes(self, lignes):
        """Extrait les informations des gènes du fichier."""
        self.genes = []
        compteur_genes = 1
        gene_courant = None
        gene_nom = None
        genes_vus = set()  # Pour suivre les gènes déjà traités
        
        for i, ligne in enumerate(lignes):
            ligne = ligne.strip()
            
            # Début de la section FEATURES
            if ligne.startswith('FEATURES'):
                continue
                
            # Chercher les lignes "gene" avec des positions
            if re.match(r'\s*gene\s+', ligne):
                # Vérifier si la ligne suivante contient un nom de gène
                if i + 1 < len(lignes):
                    ligne_suivante = lignes[i + 1].strip()
                    # Chercher d'abord un nom de gène
                    correspondance_nom = re.search(r'/gene="([^"]+)"', ligne_suivante)
                    if not correspondance_nom:
                        # Si pas de nom de gène, chercher un locus_tag
                        correspondance_nom = re.search(r'/locus_tag="([^"]+)"', ligne_suivante)
                    
                    if correspondance_nom:
                        gene_nom = correspondance_nom.group(1)
                        
                        # Si on a déjà traité ce gène, on passe au suivant
                        if gene_nom in genes_vus:
                            continue
                        
                        # Extraire la position
                        if 'complement(' in ligne:
                            # Gène complémentaire
                            correspondance = re.search(r'complement\(<?(\d+)\.\.>?(\d+)\)', ligne)
                            if correspondance:
                                debut = int(correspondance.group(1))
                                fin = int(correspondance.group(2))
                                type_gene = "complémentaire"
                                
                                # Extraire la séquence et créer le complément inverse
                                sequence_brute = self.sequence_complete[debut-1:fin]
                                sequence_traitee = self.complement_inverse(sequence_brute)
                                
                                gene_courant = {
                                    'number': compteur_genes,
                                    'start': debut,
                                    'end': fin,
                                    'type': type_gene,
                                    'sequence': sequence_traitee
                                }
                                
                        else:
                            # Gène codant
                            correspondance = re.search(r'<?(\d+)\.\.>?(\d+)', ligne)
                            if correspondance:
                                debut = int(correspondance.group(1))
                                fin = int(correspondance.group(2))
                                type_gene = "codant"
                                
                                # Extraire la séquence directement
                                sequence = self.sequence_complete[debut-1:fin]
                                
                                gene_courant = {
                                    'number': compteur_genes,
                                    'start': debut,
                                    'end': fin,
                                    'type': type_gene,
                                    'sequence': sequence
                                }
                
                        # Si on a trouvé un gène valide avec un nom, l'ajouter à la liste
                        if gene_courant:
                            self.genes.append(gene_courant)
                            genes_vus.add(gene_nom)  # Marquer le gène comme traité
                            compteur_genes += 1
                            gene_courant = None
                            gene_nom = None
    
    def formater_sequence(self, sequence, longueur_max=80):
        """Formate une séquence avec un maximum de caractères par ligne."""
        lignes_formatees = []
        for i in range(0, len(sequence), longueur_max):
            lignes_formatees.append(sequence[i:i+longueur_max])
        return '\n'.join(lignes_formatees)
    
    def ecrire_fichiers_sim(self, nom_fichier_base):
        """Écrit les fichiers .sim pour chaque gène."""
        nom_base = os.path.splitext(nom_fichier_base)[0]
        messages_sauvegarde = []
        
        for gene in self.genes:
            nom_fichier_sortie = f"{nom_base}_gene_{gene['number']}.sim"
            
            with open(nom_fichier_sortie, 'w') as fichier:
                # Ligne d'en-tête
                en_tete = f">{self.nom_organisme}:{gene['number']}:{gene['start']}:{gene['end']}:{gene['type']}"
                fichier.write(en_tete + '\n')
                
                # Séquence formatée
                sequence_formatee = self.formater_sequence(gene['sequence'])
                fichier.write(sequence_formatee + '\n')
            
            messages_sauvegarde.append(f"Gène {gene['number']} écrit dans le fichier {nom_fichier_sortie}")
        
        return messages_sauvegarde
    
    def afficher_resume(self, nom_fichier, messages_sauvegarde=None):
        """Affiche le résumé de l'analyse."""
        nombre_codants = sum(1 for gene in self.genes if gene['type'] == 'codant')
        nombre_complementaires = sum(1 for gene in self.genes if gene['type'] == 'complémentaire')
        
        # Afficher d'abord le résumé
        print(f"Nombre de lignes lues dans le fichier {nom_fichier} : {self.lignes_lues}")
        print(f"Organisme : {self.nom_organisme}")
        print(f"{len(self.genes)} gènes trouvés : {nombre_codants} codant(s) / {nombre_complementaires} complémentaire(s)")
        print(f"Séquence complète : {self.bases_totales} bases")
        
        # Afficher ensuite les messages de sauvegarde si fournis
        if messages_sauvegarde:
            for message in messages_sauvegarde:
                print(message)


class InterfaceGraphique:
    """Interface graphique pour le convertisseur génétique."""
    
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("420-555 Travail pratique 1")
        self.root.geometry("800x600")
        self.convertisseur = ConvertisseurFichierGenetique()
        self.configurer_interface()
    
    def configurer_interface(self):
        """Configure l'interface graphique."""
        cadre_principal = tk.Frame(self.root, padx=5, pady=5)
        cadre_principal.pack(fill=tk.BOTH, expand=True)
        
        # Zone de contrôles (gauche)
        self.cadre_controles = tk.Frame(cadre_principal, bg='lightblue', relief='raised', bd=2, width=250)
        self.cadre_controles.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        self.cadre_controles.pack_propagate(False)
        
        # Conteneur pour les contrôles
        conteneur = tk.Frame(self.cadre_controles, bg='lightblue')
        conteneur.pack(expand=True)
        
        # Champ de saisie et bouton
        self.variable_nom_fichier = tk.StringVar(value="complexe_2.cmp")
        tk.Entry(conteneur, textvariable=self.variable_nom_fichier, width=25, justify='center').pack(pady=(0, 10))
        tk.Button(conteneur, text="Fichier à traiter", command=self.traiter_fichier, bg='white', relief='raised', bd=1).pack()
        
        # Label "Frame tkinter"
        tk.Label(self.cadre_controles, text="Frame tkinter", bg='lightblue', font=('Arial', 8), fg='blue').pack(side=tk.BOTTOM, anchor=tk.W, padx=5, pady=5)
        
        # Zone d'animation (droite)
        self.configurer_canvas_turtle(cadre_principal)
    
    def configurer_canvas_turtle(self, parent):
        self.cadre_animation = tk.Frame(parent, bg='white', relief='sunken', bd=1)
        self.cadre_animation.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(5, 0))
        self.canvas = tk.Canvas(self.cadre_animation, bg="white", width=600, height=500)
        self.canvas.pack(expand=True, padx=2, pady=2)
        self.ecran = turtle.TurtleScreen(self.canvas)
        self.tortue_animatrice = turtle.RawTurtle(self.ecran)
        self.tortue_animatrice.speed(0)
        self.tortue_animatrice.hideturtle()
    
    def dessiner_texte_resume(self, texte_resume):
        """Dessine le texte de résumé."""
        self.tortue_animatrice.clear()
        self.tortue_animatrice.penup()
        self.tortue_animatrice.color('black')
        
        # Position fixe pour le texte
        y_start = 150
        espacement = 30
        
        # Afficher chaque ligne du résumé
        y_position = y_start
        for ligne in texte_resume.split('\n'):
            self.tortue_animatrice.goto(0, y_position)
            self.tortue_animatrice.write(ligne, align="center", font=("Arial", 14, "bold"))
            y_position -= espacement

    def animer_asterisque(self):
        tortue_anim = turtle.RawTurtle(self.ecran)
        tortue_anim.hideturtle()
        tortue_anim.speed(0)
        
        for phase in [range(9), range(8, -1, -1)]:
            for i in phase:
                tortue_anim.clear()
                tortue_anim.penup()
                tortue_anim.goto(0, -50)
                tortue_anim.pendown()
                
                for j in range(i):
                    tortue_anim.setheading(-j * 45)
                    tortue_anim.forward(20)
                    tortue_anim.backward(20)
                
                self.ecran.update()
                self.root.after(50)
        
        tortue_anim.clear()
    
    def redessiner_texte_resume(self):
        if hasattr(self, 'texte_resume_actuel'):
            self.tortue_animatrice.penup()
            y_position = 80
            for ligne in self.texte_resume_actuel.split('\n'):
                if ligne.strip():
                    self.tortue_animatrice.goto(0, y_position)
                    self.tortue_animatrice.write(ligne, align="center", font=("Arial", 12, "normal"))
                    y_position -= 25
    
    def afficher_texte_resultat(self, texte_resume):
        self.texte_resume_actuel = texte_resume
        self.tortue_animatrice.clear()
        self.tortue_animatrice.penup()
        y_position = 80
        for ligne in texte_resume.split('\n'):
            if ligne.strip():
                self.tortue_animatrice.goto(0, y_position)
                self.tortue_animatrice.write(ligne, align="center", font=("Arial", 12, "normal"))
                y_position -= 25
    
    def traiter_fichier(self):
        """Traite le fichier sélectionné."""
        nom_fichier = self.variable_nom_fichier.get().strip()
        if not nom_fichier:
            self.afficher_erreur("Veuillez entrer un nom de fichier.")
            return
        
        try:
            self.convertisseur = ConvertisseurFichierGenetique()
            self.convertisseur.parser_fichier_cmp(nom_fichier)
            self.convertisseur.ecrire_fichiers_sim(nom_fichier)
            
            # Afficher le résumé et l'animation
            self.tortue_animatrice.clear()
            self.dessiner_texte_resume(self.obtenir_texte_resume(nom_fichier))
            self.animer_asterisque()
            
            # Message final
            self.tortue_animatrice.penup()
            self.tortue_animatrice.goto(0, -100)
            self.tortue_animatrice.write("Sauvegarde terminée.", align="center", font=("Arial", 14, "bold"))
        
        except (FileNotFoundError, Exception) as e:
            self.afficher_erreur(str(e) if isinstance(e, FileNotFoundError) else f"Erreur: {str(e)}")
    
    def afficher_erreur(self, message):
        """Affiche un message d'erreur."""
        self.tortue_animatrice.clear()
        self.tortue_animatrice.penup()
        self.tortue_animatrice.goto(0, 0)
        self.tortue_animatrice.write(message, align="center", font=("Arial", 12, "normal"))
    
    def obtenir_texte_resume(self, nom_fichier):
        """Génère le texte de résumé."""
        nb_codants = sum(1 for g in self.convertisseur.genes if g['type'] == 'codant')
        nb_compl = len(self.convertisseur.genes) - nb_codants
        return "\n".join([
            f"Nombre de lignes lues dans le fichier {nom_fichier} : {self.convertisseur.lignes_lues}",
            f"Organisme : {self.convertisseur.nom_organisme}",
            f"{len(self.convertisseur.genes)} gènes trouvés : {nb_codants} codant(s) / {nb_compl} complémentaire(s)",
            f"Séquence complète : {self.convertisseur.bases_totales} bases",
            f"Sauvegarde des {len(self.convertisseur.genes)} gènes en cours ..."
        ])
    
    def executer(self):
        self.root.mainloop()


def main():
    """Fonction principale du programme."""
    if len(sys.argv) < 2 or sys.argv[1] not in ["-t", "-g"]:
        print("Usage:\n  python tp1.py -t nom_du_fichier.cmp  # Mode textuel\n  python tp1.py -g  # Mode graphique")
        return
    
    if sys.argv[1] == "-g":
        InterfaceGraphique().executer()
        return
        
    if len(sys.argv) < 3:
        print("Erreur: Nom de fichier requis pour le mode textuel\nUsage: python tp1.py -t nom_du_fichier.cmp")
        return
        
    try:
        conv = ConvertisseurFichierGenetique()
        conv.parser_fichier_cmp(sys.argv[2])
        messages = conv.ecrire_fichiers_sim(sys.argv[2])
        conv.afficher_resume(sys.argv[2], messages)
    except Exception as e:
        print(f"Erreur: {str(e)}")


if __name__ == "__main__":
    main()