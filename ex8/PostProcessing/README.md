# Post-Traitement (PostProcessing)

Ce dossier regroupe l'ensemble des outils, scripts (Bash, Python) et données utilisés pour analyser les performances informatiques et les résultats physiques des simulations.

## Structure du répertoire

Le dossier est divisé en trois sous-dossiers thématiques :

### HPC_scaling
Dédié à l'analyse des performances et de la scalabilité du code sur supercalculateur.
- **Contenu :** Scripts d'agrégation des temps de calcul, calculs de speedup, et tracés des courbes de scalabilité forte (Strong Scaling) et faible (Weak Scaling).

### SolverPreconditioner
Dédié au profilage et à l'optimisation de l'algèbre linéaire.
- **Contenu :** Scripts d'extraction des logs (ex: `aggregation_SvPc_mech.sh`), fichiers CSV de résultats, et scripts Python pour comparer l'impact des différents couples solveurs / préconditionneurs (HyprePCG, MUMPS, etc).

### Physics_postprocessing
Dédié à l'analyse purement physique et thermomécanique des résultats.
- **Contenu :** Traitement des champs physiques (température, déformation), extraction des variables d'état MFront (gonflement de l'U3Si2, plasticité de l'ALFENI), et génération des courbes d'évolution temporelle.

---

## Prérequis
Pour exécuter les scripts contenus dans ces dossiers, l'environnement suivant est généralement nécessaire :
- Python 3
- Bibliothèques : `pandas`, `matplotlib`, `numpy`
- Bash standard (outils `grep`, `sed`, `awk`)