#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "umax_vmax.h"

/**
  * affinité :
  * 		0 1 2
  * 		3 4 5
  *
  * L'image a trois niveau de couleurs
  *
  * toutes les opérations sont centrées sur le milieu de l'image
  * étant donné l'affinité A (donnée en argument) telle que X=AX'
  * on appliquera en vérité les opérations à B telle que (X-Xc)=B(X'-Xc)
  * (Xc les coordonnées du centre). On en déduit B car pour tout X',
  * AX' = BX'-BXc+Xc
  * donc avec Y = X'-Xc,
  * AY + AXc - Xc = BY
  */

//le filtre comporte 2*TAPS+1 valeurs non nulle.
#define TAPS 4

#define FILTER_TYPE 0 //the index to change the interpolation filter
//list of possible index for interpolation filters
#define RAISED_COSINE 0
#define GAUSSIAN 1

#define VARIANCE 0.36 
#define PERIOD 1.
#define BETA 0.25

//Prec indique la précision des des calculs
#define PREC 10
#define PREC2 20
#define PI 3.14159265358979323

//on a redéfinit la valeur absolue
double absd(double a){if(a>0){return a;}{return -a;}}

//une égalité pour les doubles, qui utilise PREC2
bool eq(double a,double b){if(a<b+pow(2,-PREC2)&& a>b-pow(2,-PREC2)){return true;}{return false;}}


//cette fonction tanspose si cela est utile (cf article Szeliski)
void transpo_opt(float *img,double *a,int wh[2]){
	double c1,c2;
	c1 = sqrt(pow(a[0],2)+pow(a[1],2));
	c2 = sqrt(pow(a[3],2)+pow(a[4],2));
	double a0,a1,a3,a4;
	a1=absd(a[1])/c1; a0=absd(a[0])/c1; a3=absd(a[3])/c2; a4=absd(a[4])/c2;
	if(a0+a4<a1+a3){
		int w = wh[0];
		int h = wh[1];
	//transposition de l'affinite
		double aa[6];
		aa[0]=a[3]; aa[1]=a[4]; aa[2]=a[5]; aa[3]=a[0]; aa[4]=a[1]; aa[5]=a[2];
		a[0]=aa[0]; a[1]=aa[1]; a[2]=aa[2]; a[3]=aa[3]; a[4]=aa[4]; a[5]=aa[5];

	//transposition de l'image
		float *img_t = malloc(3*w*h*sizeof(float));
		int i,j,l;
		for(l=0;l<3;l++){
			for(i=0;i<w;i++){
				for(j=0;j<h;j++){
					img_t[(j+h*i)*3+l]=img[(i+j*w)*3+l];
				}
			}
		}
		for(i=0;i<w*h*3;i++){img[i]=img_t[i];}

	//transposition de w et h
		wh[0]=h; wh[1]=w;
	}
}

//La fonction qui sert de base pour definir le filtre
float filter_fun(float x){
    switch(FILTER_TYPE){
    case GAUSSIAN:
        return exp(-pow(x,2)/(2*VARIANCE));
        break;
    case RAISED_COSINE:
    default:
        if(eq(x,0)){return 1.;}
        else if(eq(fabs(x),PERIOD/2./BETA)){return BETA/2.*sin(PI/2./BETA);}
        else{return sin(PI*x/PERIOD)/(PI*x/PERIOD)*cos(PI*BETA*x/PERIOD)/(1.-pow(2.*BETA*x/PERIOD,2));}
        break;
    }
}



//Cette fonction réalise l interpolation horizontale.
//Les xc_i, yc_f, ... sont les coordonnées du centre de l'image initiale / final
float filter_h(float *img,int w,int h,double xc_i,double yc_i,double xc_f,double yc_f,float *H,int N,int prec,double a0,double a1,int i,int j,int l){
	double x = i;
	double y = j;
	double wf = (double) w, hf = (double) h;
	double M = a0*(x+xc_f-wf/2.) + a1*(y+yc_f-hf/2.) - xc_i+wf/2.;  //Indique ou on doit centrer la gaussienne
	int k = floor(M);
	int p = floor(prec*(M-k));
	float tot = 0;

	int u = k-N;
	if(u<0){u=0;}
	//on realise la convolution
	for(;u<=k+N && u<w;u++){
	       tot += H[(k-u)*prec + p + N*prec]*img[(u+j*w)*3+l];
	}

    return tot;
}

//Fonction similaire mais verticale
float filter_v(float *img,int w,int h,double xc_i,double yc_i,double xc_f,double yc_f,float *H,int N,int prec,double a0,double a1,int i,int j,int l){
	double x = i;
	double y = j;
	double wf = (double) w, hf = (double) h;
	double M = a1*(x+xc_f-wf/2.) + a0*(y+yc_f-hf/2.) - yc_i+hf/2.;
	int k = floor(M);
	int p = floor(prec*(M-k));
	float tot = 0;

	int u = k-N;
	if(u<0){u=0;}

	for(;u<=k+N && u<h;u++){
        tot += H[(k-u)*prec + p + N*prec]*img[(i+u*w)*3+l];
	}

    return tot;
}


//Cette fonction applique un shear horizontale, de paramêtre a0 et a1 
//Le s est la largeur du filtre
//Les xc_i, yc_f, ... sont les coordonnées du centre de l'image initiale / final
int apply_rh(float *img1,float *img2,int w,int h,double xc_i,double yc_i,double xc_f,double yc_f,double s,double a0,double a1){

	//si c'est l'identité on ne veut pas convoler, on renvoit alors l'image initiale
	if(eq(a0,1.) && eq(a1,0.) && eq(xc_i,xc_f)){
		for(int i=0;i<3*w*h;i++){img2[i]=img1[i];}
		return 0;
	}

//printf("apply_rh(%f,%f,%f)\n",s,a0,a1);

	//on déclare H qui est un precalcul de coefficient du filtre
	int N = ceil(abs(s*TAPS));
	int prec = pow(2,PREC);
	float precf = (float) prec;
	//float sf = (float) s;
	float *H = malloc((2*N+1)*prec*sizeof(float));
	if(H==NULL){printf("apply_rh : H n'a pas ete cree");return 1;} //gestion d'erreur du malloc

	int k,p,i,j;

	//on remplit le H a l'aide de la fonction du filtre filter_fun
	for(p=0;p<prec;p++){
		float Htot = 0;
		float pf = (float) p;
		for(k=-N;k<=N;k++){
                float kf = (float) k;
			Htot += H[k*prec + p + N*prec] = filter_fun((kf+pf/precf)/s);  //s est la largeur du flitre
		}
		for(k=-N;k<=N;k++){H[k*prec + p + N*prec] = H[k*prec + p + N*prec]/Htot;} //On normalise le filtre, pour chaque p
	}

	//On fait la convolution
	for(int l=0;l<3;l++){
		for(j=0;j<h;j++){
			for(i=0;i<w;i++){
				img2[(i+j*w)*3 + l] = filter_h(img1,w,h,xc_i,yc_i,xc_f,yc_f,H,N,prec,a0,a1,i,j,l);
			}
		}
	}

	return 0;
}


//Fonction similaire mais verticale
int apply_rv(float *img1,float *img2,int w,int h,double xc_i,double yc_i,double xc_f,double yc_f,double s,double a0,double a1){

	if(eq(a0,1.) && eq(a1,0.) && eq(yc_i,yc_f)){
		for(int i=0;i<3*w*h;i++){img2[i]=img1[i];}
		return 0;
	}

//printf("apply_rv(%f,%f,%f)\n",s,a0,a1);

	int N = ceil(abs(s*TAPS));
	int prec = pow(2,PREC);
	float sf = (float) s;
	float precf = (float) prec;
	float *H = malloc((2*N+1)*prec*sizeof(float));
	if(H==NULL){printf("apply_rv : H n'a pas ete cree");return 1;}

	int k,p,i,j;

	for(p=0;p<prec;p++){
		float Htot = 0;
		float pf = (float) p;
		for(k=-N;k<=N;k++){
            float kf = (float) k;
			Htot += H[k*prec + p + N*prec] = filter_fun((kf+pf/precf)/s);
		}
		for(k=-N;k<=N;k++){H[k*prec + p + N*prec] = H[k*prec + p + N*prec]/Htot;}
	}

	for(int l=0;l<3;l++){
		for(i=0;i<w;i++){
			for(j=0;j<h;j++){
				img2[(i+j*w)*3+l] = filter_v(img1,w,h,xc_i,yc_i,xc_f,yc_f,H,N,prec,a0,a1,i,j,l);
			}
		}
	}
	return 0;
}



//Fonction principale, qui prend deux pointeurs vers des images, leurs dimensions, ainsi qu'une affinite.
//Elle copie le resultat de la transformation dans img_f, l'image finale
int apply_affinite(float *img,float *img_f,int w,int h,int w_f,int h_f,double *affinity){
	int i,j,l;
	
	//copie des arguments
	float *img_i = malloc(3*w*h*sizeof(float));
	for(i=0;i<3*w*h;i++){img_i[i]=img[i];}
	double a[6];
	for (i=0;i<6;i++){a[i]=affinity[i];} //utile lors des tests.

	//on transpose eventuellement l'image et l'affinite
	int wh[2] = {w,h};
	transpo_opt(img_i,a,wh);
	w=wh[0];
	h=wh[1];

    	//On calcul umax et vmax, cf. Szeliski
	double umax,vmax;
	double A[2][2] = {a[0],a[1],a[3],a[4]};
    	int test = umax_vmax(&umax,&vmax,A);  //dans le fichier "umax_vmax.h"
    	if(test==0){
        	printf("@apply_affinite : erreur dans umax_vmax\n");
        exit(1);
	}

	//Divers valeurs utiles pour calculer les tailles d'images
	//ww et hh sont les dimensions de l'image qu'on traitera (suffisamment grande pour ne pas sortir de l'image)
	int ww,hh;
	if(w<w_f){ww = w_f;} else {ww = w;}
	if(h<h_f){hh = h_f;} else {hh = h;}
	//dw et dh : différence de taille entre img_i et img_f
	int dw = (w_f-w)/2, dh = (h_f-h)/2; //x_f et x ont même parité
	//dwp et dhp leur partie positive (xp=max(0,x))
	int dwp = (dw<0) ? 0 : dw;
	int dhp = (dh<0) ? 0 : dh;

	//declaration d'image intermediaire, avec code d'erreur du malloc
	float *img1 = malloc(3*9*ww*hh*sizeof(float));
	float *img2 = malloc(3*9*ww*hh*sizeof(float));
	if(img1==NULL || img2==NULL){
        printf("apply_affinite : img1 et img2 prennent trop de place. Abandon.");
        exit(0);
    }


	//Ici il y a le choix des conditions aux bord pour éviter l'effet Gibbs. Ici on a opte pour symétrique
    ///** condition aux bords : symétrisation
	int i_sym,j_sym;
	for(l=0;l<3;l++){
		for(i=0;i<3*ww;i++){
			for(j=0;j<3*hh;j++){
                i_sym = i-ww-dwp;
                while(i_sym<0 || i_sym>w-1){i_sym = (i_sym<0) ? -1-i_sym : 2*w-1-i_sym;}
                j_sym = j-hh-dhp;
                while(j_sym<0 || j_sym>h-1){j_sym = (j_sym<0) ? -1-j_sym : 2*h-1-j_sym;}
				img1[(i+3*ww*j)*3+l]=img_i[(i_sym+j_sym*w)*3+l];
			}
		}
	}
    //*/
    /** condition aux bords : périodisation
	int i_sym,j_sym;
	for(l=0;l<3;l++){
		for(i=0;i<3*ww;i++){
			for(j=0;j<3*hh;j++){
                i_sym = i-ww-dwp;
                while(i_sym<0 || i_sym>w-1){i_sym = (i_sym<0) ? i_sym+w : i_sym-w;}
                j_sym = j-hh-dhp;
                while(j_sym<0 || j_sym>h-1){j_sym = (j_sym<0) ? j_sym+h : j_sym-h;}
				img1[(i+3*ww*j)*3+l]=img_i[(i_sym+j_sym*w)*3+l];
			}
		}
	}
    //*/
	/** condition aux bords : tronquature
	float tot[3] = {0,0,0};
	for(l=0;l<3;l++){
		for(i=0;i<w*h;i++){tot[l] += img_i[3*i+l];}
		tot[l] = tot[l]/(ww*hh);
		for(i=0;i<9*ww*hh;i++){
            //img1[3*i+l]=tot[l]; //si on veut mettre la moyenne
            img1[3*i+l]=0; //si on veut mettre du noir
        }
		for(i=0;i<w;i++){
			for(j=0;j<h;j++){
				img1[(i+ww+dwp+3*ww*(j+hh+dhp))*3+l]=img_i[(i+j*w)*3+l];
			}
		}
	}
	//*/


	//Calcul de valeurs auxiliaire pour la decomposition en shear (cf. Szeliski, on a garde les noms de l'article)
	double b0 = a[0] - a[1]*a[3]/a[4];
	double t2 = a[2] - a[1]*a[5]/a[4];
	//rv
	double aux_rv = absd(a[4])*vmax;
	if(aux_rv > 1){aux_rv = 1;}
	double rv = absd(a[1])*umax + aux_rv;
	if(rv>3){rv = 3;}else if(rv<1){rv = 1;}
	//rh
	double aux_rh = absd(b0)*umax;
	if(aux_rh > 1){aux_rh = 1;}
	double rh = absd(a[3]/a[4])*rv*vmax + aux_rh;
	if(rh>3){rh = 3;}else{if(rh<1){rh = 1;}}

	//Calcul du centre des images
	double wf = (double) w, hf = (double) h;
	double wwf = (double) ww, hhf = (double) hh;
	double xc_i, yc_i, xc_f, yc_f; //coordonnées absolues du centre des images initiale et finale
	xc_i = wf/2.;
	yc_i = hf/2.;
	xc_f = xc_i;
	yc_f = rv/a[4] * yc_i;
	
	//application du shear centre
	apply_rv(img1,img2,3*ww,3*hh,xc_i,yc_i,xc_f,yc_f,1/vmax,a[4]/rv,0.);
	
	//translation de -t2 selon x
	xc_i = xc_f - t2;
	yc_i = yc_f;
	xc_f = rh/b0 * xc_i - a[1]/rv*rh/b0 * yc_i;
	yc_f = yc_i;
	
	//deuxieme shear
    	apply_rh(img2,img1,3*ww,3*hh,xc_i,yc_i,xc_f,yc_f,1/umax,b0/rh,a[1]/rv);	
    	
    	//translation de -a[5]*rv/a[4] selon y
	xc_i = xc_f;
	yc_i = yc_f - a[5]*rv/a[4];
	xc_f = xc_i;
	yc_f = hhf/2.;
	
	//troisieme shear
	apply_rv(img1,img2,3*ww,3*hh,xc_i,yc_i,xc_f,yc_f,rv,rv,a[3]*rv/a[4]/rh);
	
	//on recentre
	xc_i = xc_f;
	yc_i = yc_f;
	xc_f = wwf/2.;
	yc_f = hhf/2.;
	
	//quatrieme shear
	apply_rh(img2,img1,3*ww,3*hh,xc_i,yc_i,xc_f,yc_f,rh,rh,0.);


	//On recopie finalement le résultat dans l'image final
	for(l=0;l<3;l++){
		for(i=0;i<w_f;i++){
			for(j=0;j<h_f;j++){
				img_f[(i+j*w_f)*3+l] = img1[(i+ww+(j+hh)*ww*3)*3+l];
			}
		}
	}

	return 0;
}
