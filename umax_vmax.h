#include <stdio.h>
#include <math.h>
#include <stdlib.h>

bool equals(double a,double b){
    //if(a<b+pow(2,-PREC)&& a>b-pow(2,-PREC)){return true;}else{return false;}
    if(a==b){return true;}else{return false;}
}

void aiv(double M[2][10], int *n, double u1, double v1, double u2, double v2){
/**
  * @name
  *     AjouterIntersectionsVerticales
  * @param
  *     M : tableau de taille >=n+1 rempli de 0 à n-1
  *     (u1,v1),(u2,v2) : extrémités du segment à intersecter avec les verticales
  *         (du carré [-1,1]^2) ; equals(u1,u2)=false
  */
    //on définit A et B les points d'intersections des droites
    double ua = 1.;
    double ub = -1.;
    double va = v1+(ua-u1)*(v2-v1)/(u2-u1);
    double vb = v1+(ub-u1)*(v2-v1)/(u2-u1);
    //on ajoute A et B à M si ils sont sur les segments
    if(fabs(va)<=1. && ((u1<=ua && ua<=u2) || (u2<=ua && ua<=u1))){
        M[0][*n] = ua;
        M[1][*n] = va;
        *n = *n+1;
    }
    if(fabs(vb)<=1. && ((u1<=ub && ub<=u2) || (u2<=ub && ub<=u1))){
        M[0][*n] = ub;
        M[1][*n] = vb;
        *n = *n+1;
    }
}

void aih(double M[2][10], int *n, double u1, double v1, double u2, double v2){
/**
  * @name
  *     AjouterIntersectionsHorizontales
  * @param
  *     M : tableau de taille >=n+1 rempli de 0 à n-1
  *     (u1,v1),(u2,v2) : extrémités du segment à intersecter avec les horizontales
  *         (du carré [-1,1]^2) ; equals(v1,v2)=false
  */
    //on définit A et B les points d'intersections des droites
    double va = 1.;
    double vb = -1.;
    double ua = u1+(va-v1)*(u2-u1)/(v2-v1);
    double ub = u1+(vb-v1)*(u2-u1)/(v2-v1);
    //on ajoute A et B à M si ils sont sur les segments
    if(fabs(ua)<=1. && ((v1<=va && va<=v2) || (v2<=va && va<=v1))){
        M[0][*n] = ua;
        M[1][*n] = va;
        *n = *n+1;
    }
    if(fabs(ub)<=1. && ((v1<=vb && vb<=v2) || (v2<=vb && vb<=v1))){
        M[0][*n] = ub;
        M[1][*n] = vb;
        *n = *n+1;
    }
}

int umax_vmax(double *u, double *v, double A[2][2]){
/**
  * @param
  *     u, v : contiendront u_max,v_max
  *     A : matrice 2x2 telle que X=AX'
  * renvoie 0 si il y a eu une erreur, 1 sinon
  */
    double detA = A[0][0]*A[1][1]-A[0][1]*A[1][0];
    //a = A^(-T)
    double a[2][2] = {A[0][0]/detA, -A[1][0]/detA, -A[0][1]/detA, A[1][1]/detA};
    // (u1,v1) = a(1,1)
    double u1 = a[0][0]+a[0][1], v1 = a[1][0]+a[1][1];
    // (u2,v2) = a(1,-1)
    double u2 = a[0][0]-a[0][1], v2 = a[1][0]-a[1][1];

    int n = 0; //taille de la liste M
    double M[2][10];

    if(fabs(u1)<=1 && fabs(v1)<=1){
        M[0][n]=u1;
        M[1][n]=v1;
        n++;
    }
    if(fabs(u2)<=1 && fabs(v2)<=1){
        M[0][n]=u2;
        M[1][n]=v2;
        n++;
    }

    if(equals(u1,u2)){
        if(equals(v1,v2)){
            printf("@umax_vmax : affinité non inversible\n");
            return 0;
        }else if(equals(v1,-v2)){
            *u = fmin(fabs(u1),1.);
            *v = fmin(fabs(v1),1.);
            return 1;
        }else{
            aih(M,&n,u1,v1,u2,v2);
            aih(M,&n,u1,v1,-u2,-v2);
            aiv(M,&n,u1,v1,-u2,-u2);
        }
    }else if(equals(u1,-u2)){
        if(equals(v1,-v2)){
            printf("@umax_vmax : affinité non inversible\n");
            return 0;
        }else if(equals(v1,v2)){
            *u = fmin(fabs(u1),1.);
            *v = fmin(fabs(v1),1.);
            return 1;
        }else{
            aih(M,&n,u1,v1,-u2,-v2);
            aih(M,&n,u1,v1,u2,v2);
            aiv(M,&n,u1,v1,u2,u2);
        }
    }else{
        aiv(M,&n,u1,v1,u2,v2);
        aiv(M,&n,u1,v1,-u2,-v2);
        if(equals(v1,v2)){
            aih(M,&n,u1,v1,-u2,-v2);
        }else if(equals(v1,-v2)){
            aih(M,&n,u1,v1,u2,v2);
        }else{
            aih(M,&n,u1,v1,u2,v1);
            aih(M,&n,u1,v1,-u2,-v2);
        }
    }
    if(n==0){
        *u = 1;
        *v = 1;
        return 1;
    }else{
        *u=fabs(M[0][0]);
        *v=fabs(M[1][0]);
        int i;
        for(i = 0;i<n;i++){
            *u=fmax(*u,fabs(M[0][i]));
            *v=fmax(*v,fabs(M[1][i]));
        }
        return 1;
    }
}
