#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//prend en entrÈe img l'image, l'image ou la mettre, z le zoom et n coeff et w,h dimensions

//on fait l'homographie separÈment sur les lignes et les colonnes grace ‡ l'intÈgration
//par convention iimg designe l'image intÈgrale

//fait la somme dans un nouveau vecteur

//En calculant l'écart type de la gaussienne déjà la si on veut un ecart de type de 0.8, on a

// D=2*sqrt(0.64*d^2-0.49);



float absf(float a){if(a>0){return a;}else{return -a;}}

int build_int_h(float *img1,float *img2,int i,int w,int h){

/**
  * @param
  *     img1 : image d'entrÈe
  *     img2,j,l : image-ligne intÈgrale de la ligne i de img1, pour la couleur l
  *     w,h : taille de l'image
  */
	float tot=0;
	for(int u=0;u<h;u++){img2[u] = tot += img1[i+u*w];}
	return 0;
}

int build_int_v(float *img1,float *img2,int j,int w,int h,int l){
/**
  * @param
  *     img1 : image d'entrÈe
  *     img2,j,l : image-colonne intÈgrale de la colonne j de img1, pour la couleur l
  *     w,h : taille de l'image
  */
	float tot=0;
	for(int u=0;u<w;u++){img2[u] = tot += img1[3*(u+j*w)+l];}
	return 0;
}


int build_triple(float *iimg,double* iimg3,int wh){
/**
  * @param
  *		iimg : ligne/colonne intégré une fois de l'image
  *		iimg3 : image triple intégrale d'une ligne/colonne de l'image
  * 	wh : taille de la ligne/colonne
  */
  
  //on remplit l'image intégrale qui est deux fois plus grande, en effet l'image intégrale n'est pas nul en dehors de l'image
  //trois fois plus grand garantit qu'on n'est jamais en dehors (car d<1.6*wh)
  //on passe en double sinon il y a des probleme de precision en integrant trois fois
	for(int i=0;i<3*wh;i++){
		if(i<wh){iimg3[i]=(double)iimg[i];}else{iimg3[i]=iimg[wh-1];}
	}
	
	//On integre deux fois
	for(int i=0;i<2;i++){
		double tot=0;
		for(int u=0;u<3*wh;u++){iimg3[u] = tot += iimg3[u];}
	}
	return 0;
}




//Evalue la moyenne sur un segment d'une ligne ou d'une colonne 
double eval_int(double *iimg,int j,int d,int wh){

/**
  * @param
  *     iimg : l'image integrale (en une dimension) a evaluer
  *     j : l'origine du segment en lequel Èvaluer iimg
  *     d : longueur du segment en lequel Èvaluer iimg (d>0)
  *     wh : la taille (w ou h) de iimg
  */


    //tentative de periodisation (pas efficace pour l'instant a cause des rotations avant)
/*	int jj = good_modulus(j,wh);
	int jd = good_modulus(j+d,wh);
	//si on essaies d'accÈder en dehors de l'image
    float iimg_jd = (0<=jd-1) ? iimg[jd-1] : iimg[wh-1];
    float iimg_jj = (0<=jj-1) ? iimg[jj-1] : iimg[wh-1];
	float dd; //le nombre de point sur lesquels on va prendre la moyenne
	if(jd>jj){dd=(float) (jd-jj);}else{dd = wh-jj+jd;}
	if(jd>jj){return (iimg_jd - iimg_jj)/(float)dd;}
	else{return (iimg[wh-1] - iimg_jj + iimg_jd)/(float)dd;}
//	*/

    //version sans periodisation
    
    //borne du segment
	int jj = j-1;
	int jd = j+d-1;
	if(jd<0){return 0;} //noir hors de l'image
	if(jj>=wh-1){return 0;}//noir hors de l'image
	
    double iimg_jd = iimg[jd];
    double iimg_jj = (0<=jj) ? iimg[jj] : 0; //il y a des zero avant le début de l'image

	double dd; //le nombre de point sur lesquels on va prendre la moyenne
	dd=(double) (jd-jj);
    return (iimg_jd - iimg_jj)/(double)dd;
}


//Fonction qui renvoie la convolution de l'image par la convolution de trois portes
double eval_triple_int(double *iimg,int j,int d,int wh){

/**
  * @param
  *     iimg : l'image integrale (en une dimension) a evaluer
  *     j : l'origine de la porte convolé en lequel evaluer iimg
  *     d : longueur de la porte convolé en lequel evaluer iimg (d>0)
  *     wh : la taille (w ou h) de iimg
  *		pix : marge possible en dehors de l'image
  */

	double a,b,c;
	if(j>=wh){return 0;} //dans ce cas on est centré en dehors de l'image
	a=eval_int(iimg,j-d,d,wh);
	b=eval_int(iimg,j,d,wh);
	c=eval_int(iimg,j+d,d,wh); //c'est ici qu'il est utile d'avoir iimg plus grande que l'image
	return (a-2*b+c)/(double)pow(d,2);		//cette formule provient de l'integration par partie de la convolution
}




double cubic_interpolation_dbl(double v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}


//interpolation des segments possibles
float triple_int(double *iimg3,float j,float d,int wh,float moyenne){
/**
  * @param
  *     iimg3 : l'image intÈgrale 1D ‡ Èvaluer
  *     j : la coordonnÈe (float) du dÈbut du segment
  *     d : la taille (float) du segment
  *     wh : la taille de l'image
  */
	float a,b,c,dd;
	float d_aux = 0.64*pow(d,2)-0.49;  //cf. formalisation papier, un peu flou...
	if(d_aux <=0){d_aux=0;}
	float D = 2*sqrt(d_aux); //ou d si on est naif
	//D=d;  //TOut ce passage est à revoir un peu, on a beaucoup de flou...
	
	j = j - D/2; //on decale le centre
	int id=floor(D);
	int ij=floor(j);
	float x = D-id;
	float y = j-ij;
	
	//a l'interieur d'un seul pixel, interpolation bilinaire entre les deux pixels voisin 
	//(il existe mieux que l'interpolation lineaire...))
	if(id<=1){
		double Y[4];
		for(int i=0;i<4;i++){Y[i]=eval_triple_int(iimg3,ij-1+i,1,3*wh);}
		return cubic_interpolation_dbl(Y,y);
	}
	
	//on est oblige de garantir id petit sinon on sort de iimg3
	if(id>=wh-1){id=wh-2;}
	
	double Y1[4];
	for(int i=0;i<4;i++){Y1[i]=eval_triple_int(iimg3,ij-1+i,id,3*wh);}
	double Y2[4];
	for(int i=0;i<4;i++){Y2[i]=eval_triple_int(iimg3,ij-1+i,id+1,3*wh);}

	a = cubic_interpolation_dbl(Y1,y);
	b = cubic_interpolation_dbl(Y2,y);

	return (1-x)*a + b*x;
}



//apply_homo est concu pour H une homographie tel que H[1]=H[4]=H[7]=0
//elle separe la transformation selon les lignes et les colonnes

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
	float *img_aux = malloc(w_aux*h_aux*sizeof(float));
	float *img_aux2 = malloc(w_f*h_f*sizeof(float));

    float flmu = (float) mu, flnu_aux = (float) nu_aux;

	for(l=0;l<3;l++){
		float iimgw[w]; //une colonne vide

		//operations colonnes par colonnes :
		for(j=0;j<h_aux;j++){
			build_int_v(img,iimgw,j,w,h,l);		//on extrait le bonne colonne
			double iimgw3[3*w];   		//On met trois parce que le D peut croitre plus vite que d.
			build_triple(iimgw,iimgw3,w); 		//on construit la triple integrale, dans une image plus grande

			for(i=0;i<w_aux;i++){
				float x = (float) (i+mu_aux);
				float d = absf((H[0]*H[8]-H[6]*H[2])/pow(H[6]*x+H[8],2)); //derivee selon x
				x = (H[0]*x+H[2])/(H[6]*x+H[8]) - flmu;		//on applique l'homographie
				img_aux[i+j*w_aux] = triple_int(iimgw3,x,d,w,iimgw[w-1]);		//on realise la convolution
			}
		}


		float iimgh[h_aux]; //une ligne vide
		
		//opÈrations lignes par lignes, similaire a la precedente :
		for(i=0;i<w_f;i++){
			build_int_h(img_aux,iimgh,i,w_aux,h_aux);
			double iimgh3[3*h_aux];
			build_triple(iimgh,iimgh3,h_aux);
			
			float x =(float) (i+mu_f);
			float d = absf(H[4]/(H[6]*x+H[8])); //dÈrivÈe selon y, qui ne dépend de y
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
