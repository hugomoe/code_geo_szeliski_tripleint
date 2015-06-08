# code_geo_szeliski_tripleint
Décomposition géométrique - affinité par Szeliski (avec translations) - homographie particulière par image triple intégrale


Il subsiste un petit décalage, probablement causé par homo_box.h (? a voir)



clang -I/usr/X11/include viho.c iio.c ftr.c -L/usr/X11/lib -lX11 -lfftw3
