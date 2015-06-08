# code_geo_szeliski_tripleint
Décomposition géométrique - affinité par Szeliski (avec translations) - homographie particulière par image triple intégrale

Pour l'executer il suffit de placer tous les fichiers dans un fichier de départ en effaçant l'ancien viho



Il subsiste un petit décalage, probablement causé par homo_box.h (? a voir)


Changement par rapport à la master version :

1 viho modifie pour mon ordi : seul le .pgm marche

2 regle le pb du blanc sur les bord dans homo_box.h

3 on rogne l'image symetrisé au bord à la fin pour limiter l'effet Gibbs

4 Un peu plus de commentaire

5 quelque modification dans homo_box.h : on a supprimer le "pix" pour juste prendre assez de marge dans tous les cas
Pas mal de temps à debugger : finalement un version un peu plus propre et commentée



clang -I/usr/X11/include viho.c iio.c ftr.c -L/usr/X11/lib -lX11 -lfftw3
