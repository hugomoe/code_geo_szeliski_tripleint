#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

//Ce programme homobow utilise l'image 4-int (modifiée) pour interpoler l'image par une fonction affine par morceaux (interp par spline d'orde 1)
//et pour calculer la convoléde cette fonction par le filtre g3 (trois porte de largeur d convolées entre elles)
//Le programme fonctionne une erreur de typage à été corrigé par rapport aux versions précédentes 

float absf(float a){if(a>0){return a;}else{return -a;}}

void build_quatre_int(float *img,double *Img,int wh){
/**
  * @param
  *     img : colonne/ligne à intégrer pour une couleur fixée
  *     wh : taille de img
  *     Img : colonne/ligne contenant les image 1-int 2-int 3-int 4-int de img aux point 0...wh (lorsque img est interpoler par une fonction constante par morceaux
  */
Img[0] = 0; // on initialise les valeurs des image 1-int 2-int 3-int 4-int à 0 en 0
Img[wh+1] = 0;
Img[2*wh+2] = 0;
Img[3*wh+3] = 0;

int u;
double q,l;
for(u=1;u<wh+1;u++){
    q = (double) img[u-1];
    Img[u] = Img[u-1] +  q;                //première intégrale 0...wh
    Img[u+wh+1] = Img[u+wh] + Img[u-1] + q/2;      //seconde intégrale  wh+1...2*wh+1
    Img[u+2*wh+2] = Img[u+2*wh+1] + Img[wh+u] + Img[u-1]/2 + q/6;//troisième intégrale 2*wh+2...3*wh+2
    Img[u+3*wh+3] = Img[u+3*wh+2] + Img[2*wh+1+u] + Img[wh+u]/2 + Img[u-1]/6 + q/24;//quatrième intégrale 3*wh+3...4*wh+3
}
}


double eval_quatre_int(float *img,double *Img,double xy,int wh){
/**
  * @param
  *     img : colonne/ligne à intégrer pour une couleur fixée
  *     wh : taille de img
  *     Img : colonne/ligne contenant les image 1-int 2-int 3-int 4-int de img aux point 0...wh (lorsque img est interpoler par une fonction constante par morceaux
  *     xy : point où évaluer Img
  *     xys : partie entière de xy
  */

    double s,r,Q;

    if(xy>=wh){
    // on calcule le polynome de degrès 3 si on sort de l'img par au dessus
       r = (double) (xy-wh);
       s =  Img[4*wh+3] + r*(Img[3*wh+2] + r*(Img[2*wh+1]/2 + r * Img[wh]/6));
    }
    else if(xy<0){
    // zeros si on sort par en dessous
       s = 0;
    }
    else{
    // On calcule F4(xys) + polynome d'interp entre xys et xys+1 si (0<=xys<wh)

      int xys = floor(xy);
      r = (double) (xy-xys);
      Q = (double) img[xys];

      s = Img[3*wh+3+xys] + r*(Img[2*wh+2+xys] + r * (Img[wh+1+xys]/2 + r *(Img[xys]/6 + r*Q/24)));
    }
    return s;
}

float convol_img(float *img,double *Img,double xy,float d,int wh){
//Ce programme permet de calculer en xy grâce à la 4_int la convolé par g3 de img interpolée initialement par une fonctions affines par morceaux
//l'écart type de g3 est défini en fonction de d
/**
  * @param
  *     img : colonne/ligne à intégrer pour une couleur fixée
  *     wh : taille de img
  *     Img : colonne/ligne contenant les image 1-int 2-int 3-int 4-int de img aux point 0...wh (lorsque img est interpoler par une fonction constante par morceaux
  *     xy : point où évaluer Img
  *     d : distance de zoom
  *     D :  taille des fonctions porte utilisées pour la convolution
  *     xys : partie entière de xy
  */


//On limite l'écart type de la triple porte sinon il y a des problèmes numériques
double d_aux = 0.64*pow(d,2)-0.49;
if(d_aux<0.001){d_aux=0.001;}
double D = 2*sqrt(d_aux);

float s;

if(D>=wh/2 && xy>=0 && xy<=wh){s = (float) Img[wh]/wh;} // moyenne de img

else {

  // points on l'on doit évaluer l'image 4-int pour réaliser la convolution
  double xy1,xy2,xy3,xy4,xy5,xy6,xy7,xy8;

  xy1 = (double) xy + (3*D+1)/2;
  xy2 = (double) xy + (3*D-1)/2;
  xy3 = (double) xy + (D+1)/2;
  xy4 = (double) xy + (D-1)/2;
  xy5 = (double) xy + (1-D)/2;
  xy6 = (double) xy - (1+D)/2;
  xy7 = (double) xy + (1-3*D)/2;
  xy8 = (double) xy - (1+3*D)/2;

  // a. = valeur de l'image 4-int en .
  // b. = valeur de la derivé discrète de l'image 4-int
  // c. = valeur de la derivé discrète seconde de l'image 4-int
  // d. = valeur de la derivé discrète troisième de l'image 4-int
  // e. = valeur de la derivé discrète quatrième de l'image 4-int = valeur de la convolution

  double a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,c1,c2,c3,d1,d2,e1;

  a1 = eval_quatre_int(img,Img,xy1,wh);
  a2 = eval_quatre_int(img,Img,xy2,wh);
  a3 = eval_quatre_int(img,Img,xy3,wh);
  a4 = eval_quatre_int(img,Img,xy4,wh);
  a5 = eval_quatre_int(img,Img,xy5,wh);
  a6 = eval_quatre_int(img,Img,xy6,wh);
  a7 = eval_quatre_int(img,Img,xy7,wh);
  a8 = eval_quatre_int(img,Img,xy8,wh);

  // dérivé discrète de taille 1 (la convolution est faite avec une fonction porte de taille 1)
  //pour passer d'une interpolation constante par morceaux à linéaire par morceaux 

  b1 = a1 - a2;
  b2 = a3 - a4;
  b3 = a5 - a6;
  b4=  a7 - a8;

  // dérivés pour  convoler

  c1 = (b1 - b2)/D;
  c2 = (b2 - b3)/D;
  c3 = (b3 - b4)/D;

  d1 = (c1 - c2)/D;
  d2 = (c2 - c3)/D;

  e1 = (d1 - d2)/D;

  s = (float) e1;
   }

return s;
}


int apply_homo(float *img,float *img_f,int w,int h,int w_f,int h_f,int mu,int nu,int mu_f,int nu_f,double H[9]){
/*
  * @param
  *     img, img_f : les images d'entrÈe et de sortie
  *     w,h, w_f,h_f : les dimensions des images
  *     mu,nu, mu_f,nu_f : les coordonnÈes du pixel en haut a gauche des images final et initiale
  *     H : homographie telle que b=c=s=0
  * Un pixel ayant une epaisseur de 1, on considere que son antÈcÈdent est d'epaisseur d
  * (d la valeur absolue de la dÈrivÈe de l'homographie en ce point)
  * Dans le code, x et y reprÈsentent les coordonnÈes reelles, float, avec decentrage
  * alors que i et j reprÈsentent les indexes dans le tableau, int, centrÈs en haut ‡ gauche
  * On pourrait Èviter certains dÈcentrage (-mu, -nu) qui seront compensÈs dans linear_int,
  * mais cela permet d'Ítre cohÈrent dans les notations
  */
	int i,j,l;

    //w_aux,h_aux, mu_aux,nu_aux pour l'image intermÈdiaire img_aux

    int w_aux = w_f; //la 2nde Ètape laisse inchangÈe x, donc w_f=w_aux
    int h_aux = h; //la 1ere Ètape laisse inchangÈe y, donc c'est noir en dehors de cette Èpaisseur
    int mu_aux = mu_f;

    int nu_aux = nu;

    float *imgw = malloc(w*sizeof(float));
	double *Img = malloc(4*(w+1)*sizeof(double));
	float *img_aux = malloc(w_aux*h_aux*sizeof(float));

	float *img_auxh = malloc(h_aux*sizeof(float));
    double *Img_aux = malloc(4*(h_aux+1)*sizeof(double));
    float *img_aux2 = malloc(w_f*h_f*sizeof(float));

    float flmu = (float) mu, flnu_aux = (float) nu_aux;

	#pragma omp parallel for
	for(l=0;l<3;l++){

		//operations colonnes par colonnes :
		#pragma omp parallel for
		for(j=0;j<h_aux;j++){
            for(i=0;i<w;i++){imgw[i] = img[3*(i+j*w)+l];} // on extrait la colonne
            build_quatre_int(imgw,Img,w); //on construit l'image4 int

			for(i=0;i<w_aux;i++){

				float x = (float) (i+mu_aux);
                float d = absf((H[0]*H[8]-H[6]*H[2])/pow(H[6]*x+H[8],2)); //derivee selon x

				x = (H[0]*x+H[2])/(H[6]*x+H[8]) - flmu;		//on applique l'homographie
				img_aux[i+j*w_aux] = convol_img(imgw,Img,x,d,w);
			}
		}


		//opÈrations lignes par lignes, similaire a la precedente :
		#pragma omp parallel for
		for(i=0;i<w_f;i++){

            for(j=0;j<h_aux;j++){img_auxh[j] = img_aux[i+j*w_aux];} //on extrait une ligne
            build_quatre_int(img_auxh,Img_aux,h_aux); // on construit l'image 4 int

			float x =(float) (i+mu_f);
			float d = absf(H[4]/(H[6]*x+H[8])); //dérivé selon y, qui ne dépend de y

			for(j=0;j<h_f;j++){

				float y = (float) (j+nu_f);
				y = (H[4]*y+H[5])/(H[6]*x+H[8]) - flnu_aux; //on applique l'hographie
				img_aux2[i+j*w_f] = convol_img(img_auxh,Img_aux,y,d,h_aux);

			}
		}
	
		#pragma omp parallel for
		for(i=0;i<w_f*h_f;i++){img_f[3*i+l]=img_aux2[i];}

	}

	return 0;
}
