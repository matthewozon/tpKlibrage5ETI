tpKlibrage5ETI
==============

Code partiel du TP de calibrage

Maintenant, c'est plus ou moins la correction.

Ce petit projet a pour but de simuler l'acquisition d'une scene par une camera dont on donne les parametres (int, ext).
On se propose aussi de faire l'estimation des parametres utilises (calibrage) en utilisant les points de la mire dans 
le repere objet et dans le repere image. On fait l'estimation des parametres physique par la methode de Tsai.

Pour utiliser ce projet, il est possible le compiler en utilisant le cmake file puis make.

cmake CMakeLists.txt ; make

Seul les library standard sont utilisees, donc il devrait etre possible de le compiler sous toutes les platformes... je 
ne l'ai teste que sous Ubuntu 14.04 sur un HP ProBook 6560b (intel core i3-2350M CPU @ 2.30GHz × 4)
