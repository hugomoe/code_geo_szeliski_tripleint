#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define PIX 12
//Il aparait un point mysterieux si on fait l'identitŽ, on dŽcale.


//prend en entrée img l'image, l'image ou la mettre, z le zoom et n coeff et w,h dimensions

//on fait l'homographie separément sur les lignes et les colonnes grace à l'intégration
//par convention iimg designe l'image intégrale

//fait la somme dans un nouveau vecteur


float absf(float a){if(a>0){return a;}else{return -a;}}

int build_int_h(float *img1,float *img2,int i,int w,int h){

/**
  * @param
  *     img1 : image d'entrée
  *     img2,j,l : image-ligne intégrale de la ligne i de img1, pour la couleur l
  *     w,h : taille de l'image
  */
	float tot=0;
	for(int u=0;u<h;u++){img2[u] = tot += img1[i+u*w];}
	return 0;
}

int build_int_v(float *img1,float *img2,int j,int w,int h,int l){
/**
  * @param
  *     img1 : image d'entrée
  *     img2,j,l : image-colonne intégrale de la colonne j de img1, pour la couleur l
  *     w,h : taille de l'image
  */
	float tot=0;
	for(int u=0;u<w;u++){img2[u] = tot += img1[3*(u+j*w)+l];}
	return 0;
}


int build_triple(float *iimg,double* iimg3,int wh){

	for(int i=0;i<wh+PIX;i++){
		//iimg3[i]=(double)iimg[i];
		if(i<wh){iimg3[i]=(double)iimg[i];}else{iimg3[i]=iimg[wh-1];} //pour le pb prabole
	}
	for(int i=0;i<2;i++){
		double tot=0;
		for(int u=0;u<wh+PIX;u++){iimg3[u] = tot += iimg3[u];}
	}
	return 0;
}




//version basique, a ameliorer, on ne periodise pas grossierment
double eval_int(double *iimg,int j,int d,int wh){

/**
  * @param
  *     iimg : l'image intégrale (1D) à évaluer
  *     j : l'origine du segment en lequel évaluer iimg
  *     d : l'épaisseur du segment en lequel évaluer iimg (d>0)
  *     wh : la taille (w ou h) de iimg
  */


    //périodisé
/*	int jj = good_modulus(j,wh);
	int jd = good_modulus(j+d,wh);
	//si on essaies d'accéder en dehors de l'image
    float iimg_jd = (0<=jd-1) ? iimg[jd-1] : iimg[wh-1];
    float iimg_jj = (0<=jj-1) ? iimg[jj-1] : iimg[wh-1];
	float dd; //le nombre de point sur lesquels on va prendre la moyenne
	if(jd>jj){dd=(float) (jd-jj);}else{dd = wh-jj+jd;}
	if(jd>jj){return (iimg_jd - iimg_jj)/(float)dd;}
	else{return (iimg[wh-1] - iimg_jj + iimg_jd)/(float)dd;}
//	*/

    //non périodisé
	int jj = j-1;
	int jd = j+d-1;
	if(jd<0){return 0;} //noir hors de l'image
	if(jj>=wh-1){return iimg[wh-1];} //noir hors de l'image

    double iimg_jd = (jd<wh) ? iimg[jd] : iimg[wh-1];
    double iimg_jj = (0<=jj) ? iimg[jj] : 0;

	double dd; //le nombre de point sur lesquels on va prendre la moyenne
	if(jd>jj){
        dd=(double) (jd-jj);
        return (iimg_jd - iimg_jj)/(double)dd;
    }else{return 0;} //si en augmentant jj ou en réduisant jd on a inversé leur sens, c'est qu'on est hors de l'image
}


double eval_triple_int(double *iimg,int j,int d,int wh){
    ///** Version Hugo
	double a,b,c;
	if(j>=wh-PIX-1){return 0;} //ce n'est pas exact, c'est une parabole en dehors de l'image on rogne un peu...
	a=eval_int(iimg,j-d,d,wh);
	b=eval_int(iimg,j,d,wh);
	c=eval_int(iimg,j+d,d,wh);
	return (a-2*b+c)/(double)pow(d,2);
	//*/
	/** Version Simon
	//on suppose iimg de taille infinie
	double eval1 = iimg[j-d];
	double eval2 = iimg[j];
	double eval3 = iimg[j+d];
	double eval4 = iimg[j+2*d];
	return (eval2 - eval1 + eval4 - eval3)/((float) d);
	//*/
}




//interpolation des segments possibles
float triple_int(double *iimg3,float j,float d,int wh,float moyenne){
/**
  * @param
  *     iimg : l'image intégrale 1D à évaluer
  *     j : la coordonnée (float) du début du segment
  *     d : la taille (float) du segment
  *     wh : la taille de l'image
  */
	float a,b,c,dd;
    //float D = d;
    float D = (d>1) ? 2*sqrt(0.64*pow(d,2)-0.49) : d; //par la formule de Morel-Yu pour respecter l'écart-type
	j = j - D/2; //on recale et on recentre
	int id=floor(D);
	int ij=floor(j);
	float x = D-id;
	float y = j-ij;

	if(D<=1){ //à l'intérieur d'un seul pixel, interpolation bilinéaire
		a=eval_triple_int(iimg3,ij,1,wh+PIX);
		b=eval_triple_int(iimg3,ij+1,1,wh+PIX);
		return a*(1-y)+b*y;
	}

	if(D>=wh){return moyenne/(float)wh;} //pour mettre l'horizon en gris

	a=eval_triple_int(iimg3,ij,id,wh+PIX);
	b=eval_triple_int(iimg3,ij+1,id,wh+PIX);
	c=eval_triple_int(iimg3,ij,id+1,wh+PIX);
	dd=eval_triple_int(iimg3,ij+1,id+1,wh+PIX);

	return (1-x)*(1-y)*a + (1-x)*y*b + x*(1-y)*c + x*y*dd;
}

//apply_homo pour H tel que H[1]=H[4]=H[7]=0

int apply_homo(float *img,float *img_f,int w,int h,int w_f,int h_f,int mu,int nu,int mu_f,int nu_f,double H[9]){
/*
  * @param
  *     img, img_f : les images d'entrée et de sortie
  *     w,h, w_f,h_f : les dimensions des images
  *     mu,nu, mu_f,nu_f : les coordonnées du premier pixel
  *     H : homographie telle que b=c=s=0
  * Un pixel ayant une épaisseur de 1, on considère que son antécédent est d'épaisseur d
  * (d la dérivée de l'homographie en ce point)
  * Dans le code, x et y représentent les coordonnées réelles, float, avec décentrage
  * alors que i et j représentent les indexes dans le tableau, int, centrés en haut à gauche
  * On pourrait éviter certains décentrage (-mu, -nu) qui seront compensés dans linear_int,
  * mais cela permet d'être cohérent dans les notations
  */
	int i,j,l;

    //w_aux,h_aux, mu_aux,nu_aux pour l'image intermédiaire img_aux
    int w_aux = w_f; //la 2nde étape laisse inchangée x, donc w_f=w_aux
    int h_aux = h; //la 1ere étape laisse inchangée y, donc c'est noir en dehors de cette épaisseur
        //SAUF SI ON PERIODISE. DANS CE CAS IL FAUT CHOISIR LA TAILLE EN FONCTION DE h_f
        //ET DE LA SURFACE NECESSAIRE
    int mu_aux = mu_f;

    int nu_aux = nu;
	float *img_aux = malloc(w_aux*h_aux*sizeof(float));
	float *img_aux2 = malloc(w_f*h_f*sizeof(float));

    float flmu = (float) mu, flnu_aux = (float) nu_aux;

	for(l=0;l<3;l++){
		float iimgw[w]; //une colonne vide

		//opérations colonnes par colonnes :
		for(j=0;j<h_aux;j++){
			build_int_v(img,iimgw,j,w,h,l);
			double iimgw3[w+PIX];
			build_triple(iimgw,iimgw3,w);

			for(i=0;i<w_aux;i++){
				float x = (float) (i+mu_aux);
				float d = absf((H[0]*H[8]-H[6]*H[2])/pow(H[6]*x+H[8],2)); //dérivée selon x
				x = (H[0]*x+H[2])/(H[6]*x+H[8]) - flmu;
				img_aux[i+j*w_aux] = triple_int(iimgw3,x,d,w,iimgw[w-1]);
			}
		}


		float iimgh[h_aux]; //une ligne vide

		//opérations lignes par lignes :
		for(i=0;i<w_f;i++){
			build_int_h(img_aux,iimgh,i,w_aux,h_aux);
			double iimgh3[h_aux+PIX];
			build_triple(iimgh,iimgh3,h_aux);

			float x =(float) (i+mu_f);
			float d = absf(H[4]/(H[6]*x+H[8])); //dérivée selon y
			for(j=0;j<h_f;j++){
				float y = (float) (j+nu_f);
				y = (H[4]*y+H[5])/(H[6]*x+H[8]) - flnu_aux;
				img_aux2[i+j*w_f] = triple_int(iimgh3,y,d,h_aux,iimgh[h_aux-1]);
			}
		}
		for(i=0;i<w_f*h_f;i++){img_f[3*i+l]=img_aux2[i];}
	}

	return 0;
}
