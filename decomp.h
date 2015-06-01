#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "affine.h"
#include "homo_box_triple_integrale.h"

//Pour les homographie
// 0 1 2
// 3 4 5
// 6 7 8

//Pour les affinité
// 0 1 2
// 3 4 5


void decomp(double H[9],double A[6],double H0[9],double B[6]){
/**
  * @param
  *     H,A,H0,B : H = A H0 B
  */
    double a=H[0], b=H[1], p=H[2], c=H[3], d=H[4], q=H[5], r=H[6], s=H[7], t=H[8];
    double t0, t1;
	//on suppose r,s != 0,0



    //il existe une infinité de couples (t0,t1) viables. On en choisit un simple
    if(fabs(a*s-b*r)<fabs(c*s-d*r)){
        t0 = 0.;
        t1 = -((a*s-b*r)*(a*r+b*s)+(c*s-d*r)*(c*r+d*s))/(pow(r,2.)+pow(s,2.))/(c*s-d*r);
    }else{
        t0 = -((a*s-b*r)*(a*r+b*s)+(c*s-d*r)*(c*r+d*s))/(pow(r,2.)+pow(s,2.))/(a*s-b*r);;
        t1 = 0.;
    }



    //on change les notations : on translate les coefficients (H devient T_(t0,t1)*H)
    a += t0*r;
    b += t0*s;
    p += t0*t;
    c += t1*r;
    d += t1*s;
    q += t1*t;



    double Nphi = sqrt(pow(r,2.)+pow(s,2.));
    B[0] = r/Nphi;
	B[1] = s/Nphi;
	B[2] = 0.;
	B[3] = -s/Nphi;
	B[4] = r/Nphi;
	B[5] = 0.;

    double Npsi = sqrt(pow(a*s-b*r,2.)+pow(c*s-d*r,2.));
	A[0] = (c*s-d*r)/Npsi;
	A[1] = (a*s-b*r)/Npsi;
	A[2] = -t0;
	A[3] = -(a*s-b*r)/Npsi;
	A[4] = (c*s-d*r)/Npsi;
	A[5] = -t1;

    H0[0] = -(a*d-b*c)*Nphi/Npsi;
    H0[1] = 0.;
    H0[2] = (p*(c*s-d*r)-q*(a*s-b*r))/Npsi;
    H0[3] = 0.;
    H0[4] = -Npsi/Nphi;
    H0[5] = (p*(a*s-b*r)+q*(c*s-d*r))/Npsi;
    H0[6] = Nphi;
    H0[7] = 0.;
    H0[8] = t;




//double theta = 3.14159265358979323*60./180.;
//A[0] = cos(theta); A[1] = sin(theta); A[2] = 0; A[3] = -sin(theta); A[4] = cos(theta); A[5] = 0;
//B[0] = cos(theta); B[1] = -sin(theta); B[2] = 0; B[3] = sin(theta); B[4] = cos(theta); B[5] = 0;

//printf("\nA :\n%f %f %f\n%f %f %f\n",A[0],A[1],A[2],A[3],A[4],A[5]);
//printf("H0 :\n%f %f %f\n%f %f %f\n%f %f %f\n",H0[0],H0[1],H0[2],H0[3],H0[4],H0[5],H0[6],H0[7],H0[8]);
//printf("B :\n%f %f %f\n%f %f %f\n",B[0],B[1],B[2],B[3],B[4],B[5]);

}

void apply_homo_final(float *img,float *img_f,int w,int h,int w_f,int h_f,double H[3][3]){

/**
  * @param
  *     img : image d'entrée
  *     img_f : image de sortie
  *     w,h : dimension de l'image d'entrée
  *     H : matrice de l'homographie img_f(x,y)=img(H(x,y))
  */
//LE BLOC DOCUMENTATION PLUS BAS NE SERVIRA QUE QUAND TOUT FONCTIONNERA
//(w1,h1,mu1,nu1)

//printf("w x h : %d x %d\n",w,h);printf("w_f x h_f : %d x %d\n",w_f,h_f);
	double *a = *H;

printf("------------------- debut de apply_homo_final -------------------");
printf("\nH :\n");
printf("%f %f %f\n",H[0][0],H[0][1],H[0][2]);
printf("%f %f %f\n",H[1][0],H[1][1],H[1][2]);
printf("%f %f %f\n",H[2][0],H[2][1],H[2][2]);
	if(a[6]==0 && a[7]==0){ //cas d'une affinité
        for(int i=0;i<8;i++){a[i]=a[i]/a[8];}
        a[8]=1.;
		double A[6];
		for(int i=0;i<6;i++){A[i]=a[i];}
		apply_affinite(img,img_f,w,h,w_f,h_f,A);
	}
	else{
		double A[6];
		double H0[9];
		double B[6];

		decomp(a,A,H0,B);










	double x4[4] = {0,w_f,w_f,0};
        double y4[4] = {0,0,h_f,h_f};
        //x1,y1 = H(x4,y4)
        double x1[4];
        double y1[4];
/**
        for(int k=0;k<4;k++){
            //si (x4[k]*a[6]+y4[k]*a[7]+a[8]==0?), !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!prévoir un cas pour simuler l'infini
            x1[k]=(x4[k]*a[0]+y4[k]*a[1]+a[2])/(x4[k]*a[6]+y4[k]*a[7]+a[8]);
            y1[k]=(x4[k]*a[3]+y4[k]*a[4]+a[5])/(x4[k]*a[6]+y4[k]*a[7]+a[8]);
        }
        double temp_min_x1 = x1[0];
        double temp_min_y1 = y1[0];
        double temp_max_x1 = x1[0];
        double temp_max_y1 = y1[0];
        for(int k=1;k<4;k++){
            if(x1[k]<temp_min_x1){temp_min_x1 = x1[k];}
            if(y1[k]<temp_min_y1){temp_min_y1 = y1[k];}
            if(x1[k]>temp_max_x1){temp_max_x1 = x1[k];}
            if(y1[k]>temp_max_y1){temp_max_y1 = y1[k];}
        }
        int min_x1 = floor(temp_min_x1);
        int min_y1 = floor(temp_min_y1);
        int max_x1 = ceil(temp_max_x1);
        int max_y1 = ceil(temp_max_y1);
        int w1 = max_x1-min_x1, h1 = max_y1-min_y1, mu1 = min_x1, nu1 = min_y1;
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!à faire
        //tronquer l'image si on ne périodise pas
        //périodiser puis tronquer la périodisée en ne gardant au plus qu'une période si on périodise
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pour l'instant, on oublie ça
        //x1,y1 sont les sommets du quadrilatère utile
        //si ils dépassaient de l'image, peut-être vaut-il mieux leur donner maintenant la valeur "coins de img1"
*/
        int w1=w;
        int h1=h;
        int mu1=0;
        int nu1=0;
printf("\nw1 x h1 : %d x %d\n",w1,h1);
printf("mu1 x nu1 : %d x %d\n",mu1,nu1);
		x1[0]=0; x1[1]=w; x1[2]=w; x1[3]=0;
		y1[0]=0; y1[1]=0; y1[2]=h; y1[3]=h;



		//x2,y2 = AA(x1,y1) où AA = A^(-1)
		double det_A = A[0]*A[4]-A[3]*A[1];
        double AA[6] = {
            A[4]/det_A,
            -A[1]/det_A,
            (A[1]*A[5]-A[4]*A[2])/det_A,
            -A[3]/det_A,
            A[0]/det_A,
            (-A[0]*A[5]+A[3]*A[2])/det_A
        };

        double x2[4];
        double y2[4];
        for(int k=0;k<4;k++){
            x2[k]=x1[k]*AA[0]+y1[k]*AA[1]+AA[2];
            y2[k]=x1[k]*AA[3]+y1[k]*AA[4]+AA[5];
        }
        double temp_min_x2 = x2[0];
        double temp_min_y2 = y2[0];
        double temp_max_x2 = x2[0];
        double temp_max_y2 = y2[0];
        for(int k=1;k<4;k++){
            if(x2[k]<temp_min_x2){temp_min_x2 = x2[k];}
            if(y2[k]<temp_min_y2){temp_min_y2 = y2[k];}
            if(x2[k]>temp_max_x2){temp_max_x2 = x2[k];}
            if(y2[k]>temp_max_y2){temp_max_y2 = y2[k];}
        }

        int min_x2 = floor(temp_min_x2);
        int min_y2 = floor(temp_min_y2);
        int max_x2 = ceil(temp_max_x2);
        int max_y2 = ceil(temp_max_y2);
        int w2 = max_x2-min_x2, h2 = max_y2-min_y2, mu2 = min_x2, nu2 = min_y2;
        //augmenter w2, h2 pour avoir même parité que w1, h1
        if((w2-w1)%2!=0){w2++;}
        if((h2-h1)%2!=0){h2++;}
/* pour contraindre la taille de img2
int ns2 = 512; //ns la new size
mu2 = mu2+(w2-ns2)/2; nu2 = nu2+(h2-ns2)/2;
w2 = ns2; h2 = ns2;*/
printf("w2 x h2 : %d x %d\n",w2,h2);
printf("mu2 x nu2 : %d x %d\n",mu2,nu2);



        //x3,y3 = B(x4,y4)
        double x3[4];
        double y3[4];
        for(int k=0;k<4;k++){
            x3[k]=x4[k]*B[0]+y4[k]*B[1]+B[2];
            y3[k]=x4[k]*B[3]+y4[k]*B[4]+B[5];
        }
        double temp_min_x3 = x3[0];
        double temp_min_y3 = y3[0];
        double temp_max_x3 = x3[0];
        double temp_max_y3 = y3[0];
        for(int k=1;k<4;k++){
            if(x3[k]<temp_min_x3){temp_min_x3 = x3[k];}
            if(y3[k]<temp_min_y3){temp_min_y3 = y3[k];}
            if(x3[k]>temp_max_x3){temp_max_x3 = x3[k];}
            if(y3[k]>temp_max_y3){temp_max_y3 = y3[k];}
        }
        int min_x3 = floor(temp_min_x3);
        int min_y3 = floor(temp_min_y3);
        int max_x3 = ceil(temp_max_x3);
        int max_y3 = ceil(temp_max_y3);
        int w3 = max_x3-min_x3, h3 = max_y3-min_y3, mu3 = min_x3, nu3 = min_y3;
        //augmenter w3, h3 pour avoir même parité que w_f, h_f
        if((w3-w_f)%2!=0){w3++;}
        if((h3-h_f)%2!=0){h3++;}
/* pour contraindre la taille de img3
int ns3 = 512; //ns la new size
mu3 = mu3+(w3-ns3)/2; nu3 = nu3+(h3-ns3)/2;
w3 = ns3; h3 = ns3;
*/
printf("w3 x h3 : %d x %d\n",w3,h3);
printf("mu3 x nu3 : %d x %d\n",mu3,nu3);



        int w4 = w_f, h4 = h_f;
        int mu4 = 0, nu4 = 0;
printf("w4 x h4 : %d x %d\n",w4,h4);
printf("mu4 x nu4 : %d x %d\n\n",mu4,nu4);



        //à ce stade, on a pris en compte une grande partie des translations de A et , via les mu,nu
        //on modifie donc A et B en conséquence

        A[2] = mu2*A[0] + nu2*A[1] + A[2] - mu1;
        A[5] = mu2*A[3] + nu2*A[4] + A[5] - nu1;
        B[2] = mu4*B[0] + nu4*B[1] + B[2] - mu3;
        B[5] = mu4*B[3] + nu4*B[4] + B[5] - nu3;



        /** allers et retours de rotation (A, A^(-1), A, A^(-1), ...)
        float *img1 = malloc(3*w1*h1*sizeof(float));
        for(int i=0;i<3*w1*h1;i++){img1[i]=img[i];}

        float *img2 = malloc(3*w2*h2*sizeof(float));
        float *img4 = malloc(3*w4*h4*sizeof(float));


printf("------------------------ 0/8\n");
        apply_affinite(img1,img2,w1,h1,w2,h2,A);
printf("+++--------------------- 1/8\n");
        apply_affinite(img2,img4,w2,h2,w4,h4,B);
printf("++++++------------------ 2/8\n");
        apply_affinite(img4,img2,w4,h4,w2,h2,A);
printf("+++++++++--------------- 3/8\n");
        apply_affinite(img2,img4,w2,h2,w4,h4,B);
printf("++++++++++++------------ 4/8\n");
        apply_affinite(img4,img2,w4,h4,w2,h2,A);
printf("+++++++++++++++--------- 5/8\n");
        apply_affinite(img2,img4,w2,h2,w4,h4,B);
printf("++++++++++++++++++------ 6/8\n");
        apply_affinite(img4,img2,w4,h4,w2,h2,A);
printf("+++++++++++++++++++++--- 7/8\n");
        apply_affinite(img2,img4,w2,h2,w4,h4,B);
printf("++++++++++++++++++++++++ 8/8\n");

        for(int l=0;l<3;l++){
            for(int i=0;i<w_f;i++){
                for(int j=0;j<h_f;j++){
                    img_f[(i+w_f*j)*3+l]=img4[(i+w4*j)*3+l];
                }
            }
        }
        //*/
        /** 6 rotations identiques
        w2 = 2*w1; h2 = 2*h1;
        w3 = w2; h3 = h2;
        float *img1 = malloc(3*w1*h1*sizeof(float));
        for(int i=0;i<3*w1*h1;i++){img1[i]=img[i];}

        float *img2 = malloc(3*w2*h2*sizeof(float));
        float *img4 = malloc(3*w4*h4*sizeof(float));
        float *img3 = malloc(3*w3*h3*sizeof(float));

printf("------------------------------ 0/6\n");
        apply_affinite(img1,img2,w1,h1,w2,h2,A);
printf("+++++------------------------- 1/6\n");
        apply_affinite(img2,img3,w2,h2,w3,h3,A);
printf("++++++++++-------------------- 2/6\n");
        apply_affinite(img3,img2,w3,h3,w2,h2,A);
printf("+++++++++++++++--------------- 3/6\n");
        apply_affinite(img2,img3,w2,h2,w3,h3,A);
printf("++++++++++++++++++++---------- 4/6\n");
        apply_affinite(img3,img2,w3,h3,w2,h2,A);
printf("+++++++++++++++++++++++++----- 5/6\n");
        apply_affinite(img2,img4,w2,h2,w4,h4,A);
printf("++++++++++++++++++++++++++++++ 6/6\n");

        for(int l=0;l<3;l++){
            for(int i=0;i<w_f;i++){
                for(int j=0;j<h_f;j++){
                    img_f[(i+w_f*j)*3+l]=img4[(i+w4*j)*3+l];
                }
            }
        }
        //*/
        ///** code original (application de la décomposition)
printf("\ntransformation de l'image en cours\n");
        float *img1 = malloc(3*w1*h1*sizeof(float));
		for(int i=0;i<3*w1*h1;i++){img1[i]=img[i];}

        float *img2 = malloc(3*w2*h2*sizeof(float));
		apply_affinite(img1,img2,w1,h1,w2,h2,A);
		//free(img1);
printf("application de A  : done\n");

		float *img3 = malloc(3*w3*h3*sizeof(float));
		apply_homo(img2,img3,w2,h2,w3,h3,mu2,nu2,mu3,nu3,H0);
		//free(img2);
printf("application de H0 : done\n");

		float *img4 = malloc(3*w_f*h_f*sizeof(float));
		apply_affinite(img3,img4,w3,h3,w4,h4,B);
		//free(img3);
printf("application de B  : done\n");

		for(int l=0;l<3;l++)
            for(int i=0;i<w_f;i++)
                for(int j=0;j<h_f;j++)
                    {img_f[(i+w_f*j)*3+l]=img4[(i+w4*j)*3+l];}
        //free(img4);
        //*/
	}
printf("-------------------- fin de apply_homo_final --------------------\n");
}










//decomp pour la decompo H_N,lambda
//void decomp(double H[3][3],double A[6],double H0[9],double B[6]){
/**
  * @param
  *     H,A,H0,B : H = A H0 B
  */
/*
	double *a = *H;
	for(int i=0;i<8;i++){a[i]=a[i]/a[8];}
	a[8]=1.;

	double N = sqrt( pow(a[6],2) + pow(a[7],2) );
	double delta1 = a[1]*a[6] - a[0]*a[7];
	double delta2 = a[4]*a[6] - a[3]*a[7];
	double gamma1 = a[0]*a[6] + a[1]*a[7];
	double gamma2 = a[3]*a[6] + a[4]*a[7];
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6] - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	double s = -( delta1*(gamma1/N - a[2]*N) + delta2*(gamma2/N - a[5]*N) )/( pow(delta1,2) + pow(delta2,2) );
	double b = pow(s*delta1 - a[2]*N + gamma1/N,2) + pow(s*delta2 - a[5]*N + gamma2/N,2) + pow(delta1/N,2) + pow(delta2/N,2);
	double lambda = sqrt(( b + sqrt( absd(pow(b,2) - 4*pow(det,2)) ) )/2);

printf("s:%f lambda:%f N:%f\n",s,lambda,N);

	B[0] = a[6]/N;
	B[1] = a[7]/N;
	B[2] = 0;
	B[3] = -a[7]/N;
	B[4] = a[6]/N;
	B[5] = s;
	A[5] = a[5] - s*delta2/N;
	A[2] = a[2] - s*delta1/N;
	//if(lambda<=1){lambda = 1;}
    A[0] = (s*delta1 - a[2]*N + gamma1/N)/lambda;
    A[1] = (delta1/N)/lambda;
    A[3] = (s*delta2 - a[5]*N + gamma2/N)/lambda;
    A[4] = (delta2/N)/lambda;

    H0[0] = lambda;
    H0[1] = 0.;
    H0[2] = 0.;
    H0[3] = 0.;
    H0[4] = lambda;
    H0[5] = 0.;
    H0[6] = N;
    H0[7] = 0.;
    H0[8] = 1.;

//printf("A :\n%f %f %f\n%f %f %f\n",A[0],A[1],A[2],A[3],A[4],A[5]);
//printf("B :\n%f %f %f\n%f %f %f\n",B[0],B[1],B[2],B[3],B[4],B[5]);

}
*/
