# Projet ECMA - MPRO 2022

Ce projet consiste à trouver le plus court chemin robuste avec différents algorithmes.

Les instances considérées sont dans le répertoire instances/

Pour executer le code, il faut se rendre dans le répertoire src/

Puis via la commande : julia main.jl

Il faut ensuite préciser l'algorithme avec --algo [nom_algo]  avec nom_algo parmi :
-s  (pour l'algorithme statique )
-bc (pour l'algorithme du branch&cut)
-d (pour l'algorithme dual)
-pc (pour l'algorithme des plans coupants)
-pck (pour l'algorithme des plans coupants knackpack)

Vous pouvez aussi rentrer d'autres paramètre avec :
--file [nom_instance]

et : --time [temps]  (le temps maximum de résolution en secondes)
