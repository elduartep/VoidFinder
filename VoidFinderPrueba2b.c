#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parameters.h"
#include "nr3.h"
#include "ran.h"
#include "eigen_sym.h"



  char catalogo_halos[500];		//	halos, are considered for setting randoms and improving small voids
  char catalogo_esferas[500];		//	spherical voids defined as the greatest sphere in each family
  char catalogo_union[500];		//	non-spherical voids catalog
  char catalogo_voids[500];		//	individual spherical members of the non-spherical voids
  char catalogo_huecos[500];		//	spheres for making the voids
  char catalogo_huecos_respaldo[500];	//	same as last one, but without HD updating
  char catalogo_particulas[500];	//	dark matter catalog
  char catalogo_randoms[500];		//	random catalog for improving small non-spherical voids
  char particulas_voids[500];		//	particles belonging to voids
  char densidad_radios[500];		//	radoius for spheres centered in the initial grid
  char randoms_radios[500];		//	radious for spheres centered in the random catalog
  char semi_ejes[500];			//	semi-axes for non-spherical voids



const int   nc=np;			//	number of cells in each direction for making a ~KDtree
const int   nh=np;			//	numero of cells in each direction for making fase 1
const int   NumPart=np*np*np;		//	number of particles
const int   NumCel=nc*nc*nc;		//	number of cells for storing the particles in the ~KDtree
const int   NumVoid=10000000;		//	maximum number of voids we expect to find (for allocation purposes)
const int   nc2=nc*nc;
const int   nc24=nc2/4;

const float rho_cr=2.77526627;		//	[10¹¹ h² Ms / Mpc³]
const float Mp=pow(Lbox,3.)*Om*rho_cr/NumPart;	// particles mass [10¹¹ Ms/h]
const float pi=4.*atan(1.);		//	number pi
const float limite=0.2*Om*rho_cr;	//	voids overdensity [10¹¹ Ms h² / Mpc³]
const float halo=200.*Om*rho_cr;	//	halos overdensity [10¹¹ Ms h² / Mpc³]

//	en http://arxiv.org/pdf/1403.5499v2.pdf solo toman en cuenta voids con radio > 2 rp
//	donde rp es la distancia media entre particulas
//	y a la hora de calcular el perfil de densidad solamente muestran el resultado para r > rp
const float rad_min_dos=0.2;		//	minimum HD radious in cell units, for randoms
const float rad_min_uno=2.;		//	minimum radious in cell units, for cell spheres

//	para el crecimiento de las esferas centradas en la red
//  radio maximo de busqueda en unidades de nc, equivale a 32Mpc     /4=8Mpc
const int   incremento=(int)(8.*nc/Lbox);//	increment the radious of search when a void is sphere is bigger than our initial guess,
const int   otro=4*incremento;		//	radious of search in cell units
const int   semi=6;			//	=2, then rad-accuracy is half the cell size
const int   bin=semi*otro;		//	maximum number of log bins
const int   semi_otro_max=768;		//	maximum value for semi*otro considering nc(128) and Lbox(256)



//	para refinar los radios estimados en centros de red y centrados en randoms
const float despl_max=0.05;		//	maximum steep for improving the rad valus, in rad units
const float despl_min=0.01;		//	minimum steep for improving the rad value, in rad units



int   *Ocu;			//	number of particles in a given cell
int  **Id;			//	index of all particles in a given cell

float *Rho;			//	cell density
float *Rad;			//	cell rad

float *X,*Y,*Z;			//	dark matter particles position

float *Radv;			//	void rad
float *Xv,*Yv,*Zv;		//	void center

float *Xr,*Yr,*Zr;		//	random positions for HD refinement
float *Rr;			//	radom rad
int    NumRan;			//	number of randoms
int    r1=0;			//	==1 -> r=1,     ==0 -> r=r

float *Xh,*Yh,*Zh;		//	haloes position
float *Rh;			//	haloes rad
int    NumHalos;		//	number of haloes

int   *ID_random;		//	<0:random >0:pixel
int   *IDV;			//	void counter
int    NumHuecos;		//	number of cell spheres with right density

float *cont;	//	cell spheres radial density, ref: number of particles in each shell
float *vol;	//	Mp/vol for each shell


float *randomxin,*randomyin,*randomzin;		//	randoms in void for elipticity computation
MatDoub tensor(3,3);				//	stress tensor
float valor[3];					//	stress tensor eigen-values
int NumVoids;					//	fumber of voids



float rad_rad[9];	//	aux variables for step 2
float rho_rho[9];
float x_x[9];
float y_y[9];
float z_z[9];



int   id_void=-1;		//	voids counter
float rad_min_actual;		//	minimum radious for each fase of the process





void archivos(void){
  sprintf(catalogo_halos,                  "halos_%s",dark_matter_file); // halos, I use them for the sampling of randoms
  sprintf(catalogo_esferas,	         "esferas_%s",dark_matter_file); // spherical voids
  sprintf(catalogo_union,	           "union_%s",dark_matter_file); // all the members of the non-spherical voids
  sprintf(catalogo_voids,  	           "voids_%s",dark_matter_file); // non-spherical voids
  sprintf(catalogo_huecos,   		  "huecos_%s",dark_matter_file); // spheres
  sprintf(catalogo_huecos_respaldo,   "huecos_NHD_%s",dark_matter_file); // spheres before the HD fase
  sprintf(catalogo_particulas,	                 "%s",dark_matter_file); // dark matter
  sprintf(catalogo_randoms,	         "randoms_%s",dark_matter_file); // catalog of randoms for improvement of small voids
  sprintf(particulas_voids,	"particulas_voids_%s",dark_matter_file); // particles belonging to voids
  sprintf(densidad_radios, 	          "radios_%s",dark_matter_file); // cell spheres radious
  sprintf(randoms_radios,  	  "randoms_radios_%s",dark_matter_file); // random spheres radious
  sprintf(semi_ejes,	  	     "elipticidad_%s",dark_matter_file); // principal axes for non-spherical voids

}







// read spheres radious
void lee_huecos(void){
  printf("\nlee_huecos_NHD-respaldo...\n");fflush(stdout);
  float x,y,z,r;
  id_void=-1;
  FILE * LC;
  LC=fopen(catalogo_huecos_respaldo,"r");
  while(fscanf(LC,"%f %f %f %f\n",&x,&y,&z,&r)!=EOF){
    id_void++;
    Xv[id_void]=x*nc/Lbox;
    Yv[id_void]=y*nc/Lbox;
    Zv[id_void]=z*nc/Lbox;
    Radv[id_void]=r*nc/Lbox;}
  fclose(LC);
  printf("lee_huecos_NHD-respaldo...hecho\n");fflush(stdout);
}





// count haloes, for allocating purposes
void cuenta_halos(void){
  printf("\ncuenta_halos\n");fflush(stdout);
  float x,y,z,r,masa;
  NumHalos=0;
  FILE * LC;
  LC=fopen(catalogo_halos,"r");
  while(fscanf(LC,"%f %f %f %f %f/n",&x,&y,&z,&masa,&r)!=EOF){
    NumHalos++;}
  fclose(LC);
  printf("cuenta_halos...hecho\n");fflush(stdout);
}


// read haloes
void lee_halos(void){
  printf("\nlee_halos\n");fflush(stdout);
  float x,y,z,r,masa;
  NumHalos=0;
  FILE * LC;
  LC=fopen(catalogo_halos,"r");
  while(fscanf(LC,"%f %f %f %f %f\n",&x,&y,&z,&masa,&r)!=EOF){
    Xh[NumHalos]=x*nc/Lbox;
    Yh[NumHalos]=y*nc/Lbox;
    Zh[NumHalos]=z*nc/Lbox;
    Rh[NumHalos]=r*nc/Lbox;
    NumHalos++;}
  fclose(LC);
  printf("lee_halos...hecho\n");fflush(stdout);
}



// allocate halo variables
void alloca_halos(void){
  printf("\nalloca_halos\n");fflush(stdout);
  printf("NumHalos=%d\n",NumHalos);fflush(stdout);
  if(!(Rh = (float*) calloc (NumHalos , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Rh.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Xh = (float*) malloc (NumHalos * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Xh.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Yh = (float*) malloc (NumHalos * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Yh.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Zh = (float*) malloc (NumHalos * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Zh.\n");
      fflush(stderr);
      exit(0);
    }
  printf("alloca halos...hecho\n");
  fflush(stdout);
}







// allocate void variables
void allocar(void){
  printf("\nallocar...\n");fflush(stdout);
  if(!(Rho = (float*) calloc (NumCel , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Rho.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Ocu = (int*) calloc (NumCel , sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory Ocu.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Rad = (float*) calloc (NumCel , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Rad.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(X = (float*) malloc (NumPart * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory X.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Y = (float*) malloc (NumPart * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Y.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Z = (float*) malloc (NumPart * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Z.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Radv = (float*) calloc (NumVoid , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Radv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Xv = (float*) malloc (NumVoid * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Xv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Yv = (float*) malloc (NumVoid * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Yv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Zv = (float*) malloc (NumVoid * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Zv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(cont = (float*) malloc (bin * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory cont.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(vol = (float*) malloc (bin * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory vol.\n");
      fflush(stderr);
      exit(0);
    }
  printf("allocar...hecho\n");
  fflush(stdout);
}







// allocate ~kdtree
void allocaId(void){
  printf("\nalloca Id...\n");
  fflush(stdout);
  Id = (int **) malloc(NumCel * sizeof(int *));
  for(int p = 0; p < NumCel; p++) { 
   if(Ocu[p]>0)
    Id[p] = (int *) malloc(Ocu[p] * sizeof(int));
   else
    Id[p] = (int *) malloc( sizeof(int));
   }
  printf("alloca Id...hecho\n");
  fflush(stdout);
}






void free_pixelizado(void){
  printf("\nfree_pixelizado...\n");fflush(stdout);
  free(Rho);			//	densidad CIC de las células
  free(Rad);			//	radio de la esfera subdensa más grande centrada en cada célula
  printf("free_pixelizado...hecho\n");fflush(stdout);
}



void free_halos(void){
  printf("\nfree_halos...\n");fflush(stdout);
  free(Xh);
  free(Yh);
  free(Zh);
  free(Rh);
  printf("free_halos...hecho\n");fflush(stdout);
}









// generate some random positions for centering more spheres around small voids
// those spheres are going to be added to those voids for improving the non-sphericity characteristic
void genera_catalogo_random(void){
  printf("\nllena_random...\n");fflush(stdout);
  Ran   myran(1050);		//	random: pagina 343 nr3
  float nr1,nr2,nr3,aux1,aux2,aux3;
  float x,y,z;
  float aux_numero;
  int   i,j,k,numero,discriminante;


  FILE * FI;
  FI=fopen(catalogo_randoms,"w+");

  NumRan=0;

  for(i=0;i<=id_void;i++){

    //	esta es la modificación más reciente
    // solo genero randoms para huecos de radio menor a 5*rad_min_uno
    if(Radv[i]<5*rad_min_uno){

      // numero de randoms para este void
      // 20 * (2 R 1.5)^3 /(4pi/3 (1.5^3-0.9^3)r^3)
      aux_numero = ceil(20*pow(Radv[i]/rad_min_uno,3));
      numero = (int) aux_numero;

      for(k=0;k<numero;k++){

//  	  ahora va desde 0.95*Radv hasta 1.5*Radv (para poder refinar hasta segundo orden)
          nr1=(float) myran.doub();
          nr1*=0.6;
          nr1+=0.9;
          nr1*=Radv[i];			//	r\in(0.9,1.5)R
          nr2=(float) myran.doub();
          nr2*=pi;			//	theta\in(0,pi)
          nr3=(float) myran.doub();
          nr3*=2.*pi;			//	phi\in(0,2pi)

          aux1=nr1*sin(nr2)*cos(nr3);	//	coordenadas esfericas
          aux2=nr1*sin(nr2)*sin(nr3);
          aux3=nr1*cos(nr2);

 	  //	condiciones de frontera periodicas
          x=fmod(Xv[i]+aux1+1.*nc,1.*nc);
          y=fmod(Yv[i]+aux2+1.*nc,1.*nc);
          z=fmod(Zv[i]+aux3+1.*nc,1.*nc);

          discriminante=0;			//	fuera de voids y halos

          //	comparo con los otros voids
          #pragma omp parallel
          {
          int discr_priv=0;
          float dist_priv;
          float dx,dy,dz;
          #pragma omp for
          for(j=0;j<=id_void;j++){
            dx=x-Xv[j];
            dy=y-Yv[j];
            dz=z-Zv[j];
            if(dx>nc*0.5)              dx-=nc*1.;
            else if(dx<-nc*0.5)        dx+=nc*1.;
            if(dy>nc*0.5)              dy-=nc*1.;
            else if(dy<-nc*0.5)        dy+=nc*1.;
            if(dz>nc*0.5)              dz-=nc*1.;
            else if(dz<-nc*0.5)        dz+=nc*1.;
            dist_priv=dx*dx+dy*dy+dz*dz;
            if( (dist_priv<(Radv[j]*Radv[j]*0.9))){	//	dentro no
              discr_priv=1;
          }}
          #pragma omp critical
          {
          if(discr_priv>0)
            discriminante=discr_priv;
          }
          }

          if(discriminante==0){		//	continua, ahora busca dentro de los halos
            #pragma omp parallel
            {
            int discr_priv=0;
            float dist_priv;
            float dx,dy,dz;
            #pragma omp for
            for(j=0;j<NumHalos;j++){
              dx=x-Xh[j];
              dy=y-Yh[j];
              dz=z-Zh[j];
              if(dx>nc*0.5)                dx-=nc*1.;
              else if(dx<-nc*0.5)          dx+=nc*1.;
              if(dy>nc*0.5)                dy-=nc*1.;
              else if(dy<-nc*0.5)          dy+=nc*1.;
              if(dz>nc*0.5)                dz-=nc*1.;
              else if(dz<-nc*0.5)          dz+=nc*1.;
              dist_priv=dx*dx+dy*dy+dz*dz;
              if(dist_priv<=Rh[j]*Rh[j]){
                discr_priv=1;
            }}
            #pragma omp critical
            {
            if(discr_priv>0)
              discriminante=discr_priv;
            }
            }
          }		//	cierra el if que da paso a la busqueda en halos

          if(discriminante==0){		//	fuera de voids
            fprintf(FI,"%f %f %f\n",x,y,z);
            NumRan++;		//	número de randoms totales
          }
      }		//	cierra el for sobre el numero k de randoms que le corresponden al void i
//printf("%d ",numero);fflush(stdout);
//if(i%100==0)
//printf("\n");fflush(stdout);
    }		//	cierra el else que filtra los huecos con radio menor que 5*rad-min_uno
  }		//	cierra el for sobre el void i, alrededor del que quiero colocar los randoms
  fclose(FI);
  printf("\nllena_random...hecho\n");fflush(stdout);
}






// count randoms if they were already generated 
void cuenta_catalogo_random(void){
  printf("\ncuenta_catalogo_random...\n");fflush(stdout);
  FILE * FI;
  FI=fopen(catalogo_randoms,"r");
  NumRan=0;
  float x,y,z;
  while(fscanf(FI,"%f %f %f\n",&x,&y,&z)!=EOF)
    NumRan++;
  fclose(FI);
  printf("cuenta_catalogo_random...hecho\n");fflush(stdout);
}






void alloca_random(void){
  printf("\nid_void=%d\n",id_void);fflush(stdout);
  printf("NumRan=%d\n",NumRan);fflush(stdout);
  printf("alloca_random...");fflush(stdout);

  if(!(Xr = (float*) calloc (NumRan , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Xr.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Yr = (float*) calloc (NumRan , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Yr.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Zr = (float*) calloc (NumRan , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Zr.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Rr = (float*) calloc (NumRan , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Rr.\n");
      fflush(stderr);
      exit(0);
    }
  printf("hecho\n");fflush(stdout);
}







// read random catalog
void lee_catalogo_random(void){
  printf("\nlee_catalogo_random...\n");fflush(stdout);
  FILE * FI;
  FI=fopen(catalogo_randoms,"r");
  for(int t=0;t<NumRan;t++)
    fscanf(FI,"%f %f %f\n",&Xr[t],&Yr[t],&Zr[t]);
  fclose(FI);
  printf("lee_catalogo_random...hecho\n");fflush(stdout);
}




// print file with the radious of spheres around the randoms
void imprime_random(void){
  printf("\nimprime_random...\n");fflush(stdout);
FILE * FI;
FI=fopen(randoms_radios,"w+");
for(int t=0;t<NumRan;t++)
  fprintf(FI,"%f %f %f %f\n",Xr[t],Yr[t],Zr[t],Rr[t]);
fclose(FI);
  printf("imprime_random...hecho\n");fflush(stdout);
}



// read the file with the radious around randoms
void lee_random(void){
  printf("\nlee_random...\n");fflush(stdout);
FILE * FI;
FI=fopen(randoms_radios,"r");
for(int t=0;t<NumRan;t++)
  fscanf(FI,"%f %f %f %f\n",&Xr[t],&Yr[t],&Zr[t],&Rr[t]);
fclose(FI);
  printf("lee_random...hecho\n");fflush(stdout);
}










//	set rad=0 to spheres inside the last one that was refined
void anula(float x,float y,float z,float radio){
  // Dados x, y, z y rad_aux, anula Rad[ijk] para todas las células ijk in void

  int i,j,k,ii,jj,kk,iii,jjj,kkk;
  ii=floor(x-radio+0.5);
  jj=floor(y-radio+0.5);
  kk=floor(z-radio+0.5);
  iii=floor(x+radio+0.5);
  jjj=floor(y+radio+0.5);
  kkk=floor(z+radio+0.5);
  float dy,dz;

  for(k=kk;k<=kkk;k++){
    dz=z-1.*k;
    for(j=jj;j<=jjj;j++){
      dy=y-1.*j;
//      #pragma omp parallel
      {
      float dist_priv,dx_priv;
      int ijk_priv;
//      #pragma omp for
      for(i=ii;i<=iii;i++){
        dx_priv=x-1.*i;
	dist_priv=dx_priv*dx_priv+dy*dy+dz*dz;
	if(dist_priv<(radio*radio)){		//	deja una franja sin anular cerca al borde?
	  ijk_priv=((i+nc)%nc)+nc*((j+nc)%nc)+nc2*((k+nc)%nc);
	  Rad[ijk_priv]=-1.;
	  Rho[ijk_priv]=halo;
	}
      }
      }
    }
  }

}








// set rad=0 to all random spheres inside the last one that was refined
void anula_random(float x,float y,float z,float radio){
  // Dados x, y, z y rad_aux, anula Rr[i] para todas los randoms i \in void
  int i;
  double radio2m=radio*radio*0.9025;	//	vecindad de radio 0.95
  #pragma omp parallel
  {
  float dist_priv;
  float dx,dy,dz;
  #pragma omp for
  for(i=0;i<NumRan;i++){
    if(Rr[i]>0){
      dx=x-Xr[i];
      dy=y-Yr[i];
      dz=z-Zr[i];
      if(dx>nc*0.5)              dx-=nc*1.;
      else if(dx<-nc*0.5)        dx+=nc*1.;
      if(dy>nc*0.5)              dy-=nc*1.;
      else if(dy<-nc*0.5)        dy+=nc*1.;
      if(dz>nc*0.5)              dz-=nc*1.;
      else if(dz<-nc*0.5)        dz+=nc*1.;
      dist_priv=dx*dx+dy*dy+dz*dz;
      if(dist_priv<(radio2m)){		//	deja una franja sin anular cerca al borde?
        Rr[i]=-1.;
      }
    }
  }
  }
}










// save all the particles inside voids in a file
void imprime(float x,float y,float z,float radio){
// busca la esfera subdensa más grande alrededor de x,y,z, y guarda su rádio y densidad
  int i,j,k,ii,jj,kk,iii,jjj,kkk,ijk,l;
  float fx,fy,fz;
  float dista;

  FILE * WRI;
  WRI=fopen(particulas_voids,"a");

  float Rmax2=radio*radio;

  float aux=radio+sqrt(3.);

  ii=floor(x-aux+0.5);	//	el (3/2)^0.5 se debe a que en la diagonal de 45º
  jj=floor(y-aux+0.5);	//	las células a esta distancia corregida
  kk=floor(z-aux+0.5);	//	tienen partículas a la distancia rad+2
  iii=floor(x+aux+0.5);
  jjj=floor(y+aux+0.5);
  kkk=floor(z+aux+0.5);

  float aux_aux=Rmax2+3./4+radio*sqrt(3.);

  for(k=kk;k<=kkk;k++){
    if(k<0) fz=-1.*nc;
    else if(k>=nc) fz=1.*nc;
    else fz=0.;
    for(j=jj;j<=jjj;j++){
      if(j<0) fy=-1*nc;
      else if(j>=nc) fy=1.*nc;
      else fy=0.;
      for(i=ii;i<=iii;i++){
        if(i<0) fx=-1.*nc;
        else if(i>=nc) fx=1.*nc;
        else fx=0.;
        ijk=(i+nc)%nc+((j+nc)%nc)*nc+((k+nc)%nc)*nc2;

	if(Ocu[ijk]>0){

          dista=pow(x-1.*i,2.)+pow(y-1.*j,2.)+pow(z-1.*k,2.);

          if(dista<aux_aux){  //	celula que tiene partículas dentro del rádio requerido

//	    #pragma omp parallel
//	    {
	    int id_priv;
//	    int int_priv;
	    float dist_priv;
//	    float aux2_priv;
//	    int cont_priv[bin]={0};
//	    #pragma omp for
            for(l=0;l<Ocu[ijk];l++){
              id_priv=Id[ijk][l];
              dist_priv=pow(x-X[id_priv]-fx,2.)+pow(y-Y[id_priv]-fy,2.)+pow(z-Z[id_priv]-fz,2.);
              if (dist_priv<Rmax2){
		fprintf(WRI,"%f %f %f %f\n",X[id_priv]*Lbox/nc,Y[id_priv]*Lbox/nc,Z[id_priv]*Lbox/nc,radio*Lbox/nc);
	      }
	    }
//	    #pragma omp critical
//	    {
//	    for(l=0;l<bin;l++)
//	    cont[l]+=1.*cont_priv[l];
//	    }
//	    }
	  }	//	if de la distancia al centro
	}	//	if del número de ocupación
      }
    }
  }
fclose(WRI);
}


















//	given a coordinate x,y,z
//	this routine tells me if it is inside any sphere
int criterio(float x, float y, float z){	//	1: OK, 0: está dentro de un hueco, no calcular
  //	el centro no debe estar dentro de los huecos ya identificados

  int entero=1;	//	es la variable que voy a retornar
  int aux_void;

  #pragma omp parallel
  {
  float tx,ty,tz;
  float dist_priv;
  int   entero_priv=1;
  #pragma omp for
  for(aux_void=0;aux_void<=id_void;aux_void++){
    tx = Xv[aux_void]-x;
         if (tx<-nc*0.5)	tx+=1.*nc;
    else if (tx>nc*0.5)		tx-=1.*nc;

    ty = Yv[aux_void]-y;
         if (ty<-nc*0.5)	ty+=1.*nc;
    else if (ty>nc*0.5)		ty-=1.*nc;

    tz = Zv[aux_void]-z;
         if (tz<-nc*0.5)	tz+=1.*nc;
    else if (tz>nc*0.5)		tz-=1.*nc;

    dist_priv=sqrt(tx*tx+ty*ty+tz*tz);
    //	si este centro está dentro de un hueco
    if(dist_priv<Radv[aux_void]*0.95){
      entero_priv=0;
      aux_void=id_void+200;	//	para que salga del pragma for
      ///////////// no se si esto funcione, pero es para ahorrar tiempo
    }
  }
  #pragma omp critical
  {
  if(entero_priv==0)
    entero=0;
  }
  }
//printf("entero=%d\n",entero);fflush(stdout);
return entero;

}













// refine rad and position for a given spheres, given initial guesses for x y z r
void crece_refina(float x,float y,float z,float &rho,float &radio){
  int i,j,k,ii,jj,kk,iii,jjj,kkk,ijk,l;
  float dx,dy,dz,r,Rmax;
  float dista,c;//,radio_minimo;
  int permanencia=1;

  // busca entre r=rad_min_dos y radio_min_uno*2 cuando r1=true
  // es la parte donde calcula los radios preliminares de los randoms
  if(r1==1){
    Rmax=rad_min_uno*4.;
    r=rad_min_dos;}

  // es la parte donde refina los radios de los randoms
  // busca entre rad_min_actual y radio+1
  else if(r1==2){
    Rmax=radio+3.;
    r=rad_min_dos;
  }

//	marzo 14
//	parece injusto que se busque en un rango fijo de 5 Mpc pata voids de radio 20 y para voids de radio 2
//	para voids grandes parece bien buscar entre 20% más y 20% menos
//	pero pava voids pequeños, r=2, dicha franja está entre 1.6 y 2.4, lo cual parece muy pequeño
//	y dado que la cantidad de voids pequeños es grande, no es bueno que estos tamaños no estén optimizados
//	voy a probar con 20% para ver que pasa, si el tiemp de calculo no se incrementa mucho lo dejo así

  // el la parte donde refina las celdas
  // busca entre radio-2 y radio+3
/*  else{
    Rmax=radio+3.;
    r=Rmax-5.;
    if(r<1.)
      r=1.;
  }*/
  else{
    Rmax=radio*1.2;
    r=radio*0.833334;
  }


  //radio_minimo=rad_min_actual;


  if (Rmax>rad_min_actual){


while(permanencia>0){

  float Rmax2=Rmax*Rmax;
  float r2   =r*r;
  c    	     =pow(Rmax/r,1./bin);	//	constante de proporcionalidad log
  float aux1 =0.5/log(c);

  // volumen y radio de cada bin
  #pragma omp parallel for
  for(i=0;i<bin;i++){
    vol[i]=Mp*3/(4*pi*pow(r*Lbox/nc,3.)*pow(c,3.*i));	//	en [Ms h² / Mpc³]
    cont[i]=0.;						//	# de partículas en cada bin
  }

  float aux=Rmax+sqrt(3.);

  ii=floor(x-aux+0.5);	//	el (3/2)^0.5 se debe a que en la diagonal de 45º
  jj=floor(y-aux+0.5);	//	las células a esta distancia corregida
  kk=floor(z-aux+0.5);	//	tienen partículas a la distancia rad+2
  iii=floor(x+aux+0.5);
  jjj=floor(y+aux+0.5);
  kkk=floor(z+aux+0.5);

  float aux_aux=Rmax2+3./4+Rmax*sqrt(3.);

  for(k=kk;k<=kkk;k++){
    dz=z-1.*k;
    for(j=jj;j<=jjj;j++){
      dy=y-1.*j;
      for(i=ii;i<=iii;i++){
        ijk=((i+nc)%nc)+((j+nc)%nc)*nc+((k+nc)%nc)*nc2;

	if(Ocu[ijk]>0){
          dx=x-1.*i;
           
          dista=dx*dx+dy*dy+dz*dz;

          if(dista<aux_aux){  //	celula que tiene partículas dentro del rádio requerido

//	    #pragma omp parallel
//	    {
	    int id_priv;
	    int int_priv;
	    float dist_priv;
	    float aux2_priv;
            float dx_priv;
            float dy_priv;
            float dz_priv;
//	    int cont_priv[bin]={0};
//	    #pragma omp for
            for(l=0;l<Ocu[ijk];l++){
              id_priv=Id[ijk][l];
              dx_priv=x-X[id_priv];
                   if(dx_priv<-nc*0.5) dx_priv+=1.*nc;
              else if(dx_priv>nc*0.5)  dx_priv-=1.*nc;

              dy_priv=y-Y[id_priv];
                   if(dy_priv<-nc*0.5) dy_priv+=1.*nc;
              else if(dy_priv>nc*0.5)  dy_priv-=1.*nc;

              dz_priv=z-Z[id_priv];
                   if(dz_priv<-nc*0.5) dz_priv+=1.*nc;
              else if(dz_priv>nc*0.5)  dz_priv-=1.*nc;

              dist_priv=dx_priv*dx_priv + dy_priv*dy_priv + dz_priv*dz_priv;

              if (dist_priv<Rmax2){
                if(dist_priv<=r2)
                  aux2_priv=0.1;
                else
                  aux2_priv=aux1*log(dist_priv/r2);
                int_priv=floor(aux2_priv);	//	es el bin radial de esta particula
//                cont_priv[int_priv]+=1;	//	contador entero del número de partículas
                cont[int_priv]+=1.;	//	contador entero del número de partículas
	      }
	    }
//	    #pragma omp critical
//	    {
//	    for(l=0;l<bin;l++)
//	    cont[l]+=1.*cont_priv[l];
//	    }
//	    }
	  }	//	if de la distancia al centro
	}	//	if del número de ocupación
      }
    }
  }

  for(l=1;l<bin;l++){
    cont[l]+=cont[l-1];}

  // hacer las cosas desde fuera hacia dentro nos permite detectar voids con halos dentro
  l=bin-1;
  rho=cont[l]*vol[l];	//	densidad en unidades [Ms h² / Mpc³]

  if(rho<limite){
    Rmax*=1.2;
    r*=1.2;
  }
  else{
    while((rho>=limite)&&(l>-1)){
      rho=cont[l]*vol[l];
      l--;
    }
    permanencia=-1;
  }

}	//	cierra el while

  l++;
  radio=r*pow(c,1.*l);	//	en unidades de nc

}
else{
radio=0.;
rho=halo;
}

//printf("final crece_refina %f %f %f %f\n",x,y,z,radio);fflush(stdout);
}













// grows spheres centered in the random catalog
void crece_random(void){
  printf("crece_random...\n");fflush(stdout);
  r1=1;						//	crece_refina busca desde r=rad_min_actual
  float rho_r;
  for(int i=0;i<NumRan;i++){
//    if(criterio(Xr[i],Yr[i],Zr[i])==1)		//	me dice siestá o nó dentro de un hueco
      crece_refina(Xr[i],Yr[i],Zr[i],rho_r,Rr[i]);//	modifica rho_r y Rr[i]
//    else Rr[i]=0.;
  }
  r1=0;						//	crece_refina busca desde r=r
  printf("crece_random...hecho\n");fflush(stdout);
}










// compares spheres centered in neigbouring places to maximize their size
void refino(int l,float &x,float &y,float &z, float &radio){
  int i,j,k,p,q;

  float escala=radio*despl_max;		//	mayor desplazamiento del centro durante la busqueda
  float presicion=radio*despl_min;	//	error al determinar la posicion del centro

  if (escala>0.5) escala=0.5;
  if (presicion>0.125) presicion=0.125;

  if(r1==2){			//	refino HD
    x_x[0]=Xr[l];
    y_y[0]=Yr[l];
    z_z[0]=Zr[l];
  }
  else{				//	refino celdas
    i=l%nc;
    j=(l-i)%nc2;
    k=(l-i-j)/nc2;
    j/=nc;
    x_x[0]=i*1.;
    y_y[0]=j*1.;
    z_z[0]=k*1.;
  }


  float sx,sy,sz;
  float aux_rad;
//  float aux_rho;

  rad_rad[0]=radio;
//  rho_rho[0]=rho;

if(criterio(x_x[0],y_y[0],z_z[0])==0){		//  me dice si está o nó dentro de un hueco
  rho_rho[0]=1.;rad_rad[0]=0.;}
else{
  crece_refina(x_x[0],y_y[0],z_z[0],rho_rho[0],rad_rad[0]);    //	modifica rho_rho[0] y rad_rad[0]

//int muy_auxiliar=0;
  while(escala>=presicion){
    q=1;
    for(sz=z_z[0]-escala;sz<=z_z[0]+escala*1.5;sz+=2.*escala){
      for(sy=y_y[0]-escala;sy<=y_y[0]+escala*1.5;sy+=2.*escala){
        for(sx=x_x[0]-escala;sx<=x_x[0]+escala*1.5;sx+=2.*escala){
          x_x[q]=(float) fmod(1.*sx+1.*nc,1.*nc);
          y_y[q]=(float) fmod(1.*sy+1.*nc,1.*nc);
          z_z[q]=(float) fmod(1.*sz+1.*nc,1.*nc);
          rad_rad[q]=radio;
          if(criterio(x_x[q],y_y[q],z_z[q])==1)
            crece_refina(x_x[q],y_y[q],z_z[q],rho_rho[q],rad_rad[q]);
          else{rho_rho[q]=1.;	rad_rad[q]=0.;}
          q++;
	}
      }
    }

    aux_rad=0.;
    for(q=8;q>=0;q--){
      if(aux_rad<=rad_rad[q]){
        aux_rad=rad_rad[q];
        p=q;
      }
    }
/*    aux_rho=1.;
    for(q=8;q>=0;q--){
      if(aux_rad==rad_rad[q]){
        if(aux_rho>=rho_rho[q]){
          aux_rho=rho_rho[q];
          p=q;
	}
      }
    }
*/
    if(p==0){   //	si el centro es mejor candidato que sus vecinos reduce el radio de busqueda
      escala*=0.5;
    }
    else{	//	si un vecino es mejor candidato, lo convierto en el centro de busqueda
      x_x[0]=x_x[p];
      y_y[0]=y_y[p];
      z_z[0]=z_z[p];
      rad_rad[0]=rad_rad[p];
//      escala*=0.5;              //	cuando quiera desplazarme debo comentar esta línea
    }


///
///
///
//	si radio<radio_minimo entonces salgo, y anulo todo lo que pueda


  }             //	cierra el while de escala

}		//	cierra el if (si esta fuera de otros huecos)

  radio=rad_rad[0];
  x=x_x[0];
  y=y_y[0];
  z=z_z[0];

}

















inline int fi(int indice,int radio){
  return -radio+indice%(2*radio+1);
  }
inline int fj(int indice,int radio){
  int aux=floor(1.*indice/(2*radio+1));
  return -radio+aux%(2*radio+1);
  }
inline int fk(int indice,int radio){
  return -radio+floor(1.*indice*pow(2.*radio+1.,-2.));
  }













//	read matter cataloge and fill the  ~kdtree
void lee(){
  printf("lee xyz...\n");
  fflush(stdout);

  int i,j,k,ijk,p;
  float x,y,z,u,v,w;

  FILE * LEE;
  LEE = fopen(catalogo_particulas,"r");
  printf("%s\n",catalogo_particulas);

  // llena el vector de posiciones y cuenta las particulas que hay en cada celda
  for (p=0;p<NumPart;p++)
    {
      //    fscanf(LEE,"%f %f %f %f %f %f\n",&x,&y,&z,&u,&v,&w);
    fscanf(LEE,"%f %f %f\n",&x,&y,&z);
    X[p]=fmod(x*nc/Lbox,1.*nc);
    Y[p]=fmod(y*nc/Lbox,1.*nc);
    Z[p]=fmod(z*nc/Lbox,1.*nc);
    i=floor(X[p]+0.5);	//	debo sumar 0.5 porque las células están centradas
    j=floor(Y[p]+0.5);	//	en los vértices enteros del sistema de coordenadas
    k=floor(Z[p]+0.5);	//	Nó en puntos tales como 1.5, 45.5, 178.5
    ijk=(i%nc)+(j%nc)*nc+(k%nc)*nc*nc;
    Ocu[ijk]++;
    }
  fclose (LEE);
  printf("lee xyz...hecho\n");
  fflush(stdout);

  allocaId();

  for(p=0;p<NumCel;p++)
    Ocu[p]=0;

  printf("Llena ID...\n");
  fflush(stdout);
  for(p=0;p<NumPart;p++)//loop sobre las artículas
    {
    i=floor(X[p]+0.5);	//	debo sumar 0.5 porque las células están centradas
    j=floor(Y[p]+0.5);	//	en los vértices enteros del sistema de coordenadas
    k=floor(Z[p]+0.5);	//	Nó en puntos tales como 1.5, 45.5, 178.5
    ijk=(i%nc)+(j%nc)*nc+(k%nc)*nc2;
    Id[ijk][Ocu[ijk]]=p;
    Ocu[ijk]++;
    }

  printf("Llena Id...hecho\n");
  fflush(stdout);

}
















//	computes the cells density
void densidad(){
  printf("\ncalcula la densidad de las celdas...\n");
  fflush(stdout);

  int i,j,k,ii,jj,kk,p;
  float sx,sy,sz;
  float tx,ty,tz;

  for(p=0;p<NumPart;p++)//loop sobre las particulas para calcular la densidad
    {
    i=floor(X[p]);        j=floor(Y[p]);          k=floor(Z[p]);
    sx=X[p]-1.*i;         sy=Y[p]-1.*j;           sz=Z[p]-1.*k;
    tx=1.-sx;             ty=1.-sy;               tz=1.-sz;

    // condiciones de frontera periodicas
    ii=(i+1)%nc;          jj=(j+1)%nc;            kk=(k+1)%nc;

    // incrementando la densidad de las celdas adyacentes
    Rho[i+nc*j +nc2*k] +=Mp*tx*ty*tz;         Rho[ii+nc*j +nc2*k] +=Mp*sx*ty*tz;
    Rho[i+nc*jj+nc2*k] +=Mp*tx*sy*tz;         Rho[ii+nc*jj+nc2*k] +=Mp*sx*sy*tz;
    Rho[i+nc*j +nc2*kk]+=Mp*tx*ty*sz;         Rho[ii+nc*j +nc2*kk]+=Mp*sx*ty*sz;
    Rho[i+nc*jj+nc2*kk]+=Mp*tx*sy*sz;         Rho[ii+nc*jj+nc2*kk]+=Mp*sx*sy*sz;
    // Rho sería la masa física contenida en cada célula
    }

  printf("calcula densidad de las celdas...hecho\n");
  fflush(stdout);
}












// grows a sphere centered in a cell untill it reaches the critical density
void crece(){
  printf("\nradio máximo de búsqueda usando células = %f [Mpc/h]\n",otro*Lbox/nc);
  printf("Encuentra el radio de la mayor esfera subdensa centrada en cada célula...\n");
  fflush(stdout);

  int i,j,k,l,p,ijk;
  float rho;

  FILE * RAD;
  RAD=fopen(densidad_radios,"w+");

  // incremento la densidad de todas las esferas
  int max=incremento;		//	radio máximo en unidades de nc
  int max2=incremento*incremento;

  for(k=0;k<nc;k++){				//	ciclo sobre las células
    for(j=0;j<nc;j++){
      for(i=0;i<nc;i++){

      ijk=i+nc*j+nc2*k;
      //	nuevo: optimiza: solo crezco esferas alrededor de pixeles subdensos
      if(Rho[ijk]>limite)
        Rad[ijk]=0.;

      else{

	int sale=1;
	while(sale>0){	//	buca para un radio, si no funciona lo aumenta,etc

	  // borro la información antigua
	  #pragma omp parallel for
 	  for(p=0;p<semi*max;p++){
	    cont[p]=0.;
	    vol[p]=0.;
	    }

	  #pragma omp parallel
	  {
	  int i_priv;
	  int j_priv;
	  int k_priv;
	  int ijk_priv;
	  int bin_priv;
	  int dist_priv;
	  int vol_priv[semi_otro_max]={0};
	  float mas_priv[semi_otro_max]={0.};
	  #pragma omp for
	  for(p=0;p<(2*max+1)*(2*max+1)*(2*max+1);p++){
	    k_priv=k+fk(p,max);
	    j_priv=j+fj(p,max);
	    i_priv=i+fi(p,max);
            dist_priv=(k_priv-k)*(k_priv-k)+(j_priv-j)*(j_priv-j)+(i_priv-i)*(i_priv-i);
            if(dist_priv<max2){		//	si es vecina de verdad...
	      bin_priv=floor(semi*sqrt(dist_priv));
              ijk_priv=(i_priv+nc)%nc+nc*((j_priv+nc)%nc)+nc2*((k_priv+nc)%nc);
	      mas_priv[bin_priv]+=Rho[ijk_priv];	//	masa física contenida en el bin l
	      vol_priv[bin_priv]++;		//	volumen en unidades de nc
	    }
	  }
	  #pragma omp critical
	  {
	  for(p=0;p<semi*max;p++){
	    cont[p]+=mas_priv[p];
	    vol[p]+=vol_priv[p]*1.;
	  }
	  }
	  }	//	sale del pragma

	  for(l=1;l<semi*max;l++){
	    cont[l]+=cont[l-1];	//	masa física de la esfera l
	    vol[l]+=vol[l-1];	//	volumen de la esfera l en unidades de nc
	  }

	  // hacer las cosas de fuera para dentro nos permite detectar voids con halos dentro
	  l=semi*max-1;
	  rho=cont[l]/vol[l]*pow(nc/Lbox,3);	//	densidad en [Ms h² / Mpc³]

	  if(rho<limite){	//	debo aumentar el radio de busqueda
	    max+=incremento;
	    max2=max*max;
	  }

	  else{		//	si todo esta bien
	    max=incremento;
	    max2=max*max;
	    sale=-1;
	    while((rho>=limite)&&(l>-1)){
	      rho=cont[l]/vol[l]*pow(nc/Lbox,3);	//	densidad en [Ms h² / Mpc³]
	      l--;
	    }
	  }

	  l++;	//	porque l=0 es la esfera de radio 1

	  Rad[ijk]=l*1./semi;		//	radio en unidades de nc
	}	//	cierra el while de sale
      }		//	else: si es subdensa
      fprintf(RAD,"%f\n",Rad[ijk]);
      }		//	i
    }		//	j
  }		//	k

  fclose(RAD);

  printf("ya tenemos los radios de las mayores esferas subdensas centradas en cada célula\n\n");
  fflush(stdout);
}
















// read the cell radious file
void lee_radios(){
  printf("\nlee el archivo con radios y desidades de celulas...\n");
  fflush(stdout);

  float dist;

  FILE * RAD;
  RAD=fopen(densidad_radios,"r");
  for(int p=0;p<NumCel;p++){
    fscanf(RAD,"%f\n",&dist);
    Rad[p]=dist;
  }
  fclose(RAD);
  printf("lee el archivo con radios y desidades de celulas...hecho\n");
  fflush(stdout);
}












//	search for the bigest sphere
int busco(void){
  float radio=rad_min_actual;
//  float rho=halo;
  int l=-1;
  int i;
  #pragma omp parallel
  {
  float radio_priv=radio;	//	=radio_min_actual
//  float rho_priv=rho;
  int   l_priv=l;		//	=-1
  #pragma omp for
  for(i=0;i<NumCel;i++){	//	solo busco el radio más grande
    if(radio_priv<Rad[i]){
      radio_priv=Rad[i];	//	radio en unidades de nc
      l_priv=i;			//	es el índice de la célula con la esfera más grande
    }
//    else if(radio_priv==Rad[i]){
//      if(rho_priv>Den[i]){
//        rho_priv=Den[i];		//	densidad física
//        l_priv=i;			//	es el índice de la célula con la esfera más grande
//      }
//    }
  }
  #pragma omp critical
  {
  if(radio<radio_priv){
    radio=radio_priv;	//	radio en unidades de nc
//    rho=rho_priv;		//	densidad física
    l=l_priv;			//	es el índice de la célula con la esfera más grande
  }
//  else if(radio==radio_priv){
//    if(rho>rho_priv){
//      rho=rho_priv;		//	densidad física
//      l=l_priv;			//	es el índice de la célula con la esfera más grande
//    }
//  }
  }				//	ya tengo el radio máximo
  }
  return l;
}








//  search for the bigest spere around the random catalog
int busco_random(void){
  float radio=rad_min_actual;
//  float rho=1.;
  int l=-1;
  int i;
  #pragma omp parallel
  {
  float radio_priv=radio;
//  float rho_priv=rho;
  int   l_priv=l;
  #pragma omp for
  for(i=0;i<NumRan;i++){	//	solo busco el radio más grande
    if(radio_priv<Rr[i]){
      radio_priv=Rr[i];		//	radio en unidades de nc
      l_priv=i;			//	es el índice de la célula con la esfera más grande
    }
  }
  #pragma omp critical
  {
  if(radio<radio_priv){
    radio=radio_priv;		//	radio en unidades de nc
//    rho=rho_priv;		//	densidad física
    l=l_priv;			//	es el índice de la célula con la esfera más grande
  }
  }				//	ya tengo el radio máximo
  }
  return l;
}









void free_xyz(void){
  free(X);
  free(Y);
  free(Z);
  free(Xv);
  free(Yv);
  free(Zv);
  free(Radv);
  free(Ocu);
  for(int p=0;p<NumCel;p++)
    free(Id[p]);
  free(Id);
}


// count spheres in the file
void cuenta_huecos(void){
  printf("\ncuenta_huecos...\n");fflush(stdout);
  NumHuecos=0;
  float u,v,w,r;
  FILE * HU;
  HU=fopen(catalogo_huecos,"r");
  while(fscanf(HU,"%f %f %f %f\n",&u,&v,&w,&r)!=EOF)
    NumHuecos++;
  fclose(HU);
  printf("NumHuecos=%d\n",NumHuecos);fflush(stdout);
  printf("cuenta_huecos...hecho\n");fflush(stdout);
}






//	realloc void variables
void realloca(void){
  printf("\nrealloca V...\n");
  fflush(stdout);
  if(!(Xv = (float*) calloc (NumHuecos , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Xv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Yv = (float*) calloc (NumHuecos , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Yv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Zv = (float*) calloc (NumHuecos , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Zv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Radv = (float*) calloc (NumHuecos , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Radv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(IDV = (int*) calloc (NumHuecos , sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory IDV.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(ID_random = (int*) calloc (NumHuecos , sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory ID_random.\n");
      fflush(stderr);
      exit(0);
    }
  printf("realloca V...hecho\n");
  fflush(stdout);
}







//	read the spheres catalog, the one was improved by the HD procedure
void lee_huecos_2(void){
  printf("\nlee_huecos2...\n");fflush(stdout);
  float u,v,w,r;

  FILE * PIX;
  // cuenta el número de objetos en el catalogo pixelizado (el de respaldo)
  PIX=fopen(catalogo_huecos_respaldo,"r");
  int NumPix=0;
  while(fscanf(PIX,"%f %f %f %f\n",&u,&v,&w,&r)!=EOF)
    NumPix++;
  fclose(PIX);
  printf("Num Obj NHD=%d\n",NumPix);fflush(stdout);

  FILE * HU;
  // lee el archivo de trabajo (no el de respaldo)
  HU=fopen(catalogo_huecos,"r");
  for(int i=0;i<NumHuecos;i++){
    fscanf(HU,"%f %f %f %f\n",&u,&v,&w,&r);
    Xv[i]=u*nc/Lbox;
    Yv[i]=v*nc/Lbox;
    Zv[i]=w*nc/Lbox;
    Radv[i]=r*nc/Lbox;
    IDV[i]=30000000;		//	es el indice del void al que pertenece
    if(i<NumPix)
      ID_random[i]=30000000;	//	>0 viene de los pixesles
    else			
      ID_random[i]=-3;		//	<0 viene de los random de los random
  }
  fclose(HU);
  printf("NumHuecos=%d\n",NumHuecos);fflush(stdout);
  printf("lee_huecos2...hecho\n");fflush(stdout);
}







// gather spheres into families
void depura(void){
  printf("\ndepura...\n");fflush(stdout);

  FILE * ES;
  ES=fopen(catalogo_esferas,"w+");

  int i;
  float local_max_rad;
  float local_x,local_y,local_z,local_r;
  float dist1,dist2;
  int   local_id_void=-1;	//	indice del void identificado más recientemente
  int   local_id_mayor=0;	//	hueco más gande sin familia
  int   local_id_cerca;		//	el indice del void más cercano

  printf("entra al ciclo de depuración\n");fflush(stdout);
  printf("NumHuecos=%d\n",NumHuecos);fflush(stdout);
 

  while(local_id_mayor>=0){	//	ciclo principal
    local_max_rad=-1.;
    local_id_mayor=-1;


    #pragma omp parallel
    {
    float local_max_rad_priv=local_max_rad;
    int   local_id_mayor_priv=local_id_mayor;
    #pragma omp for
    for(i=0;i<NumHuecos;i++)			//	Busca el hueco más grande
      if((IDV[i]>local_id_void))		//	si aún no tiene familia, y no falló antes
        if(local_max_rad_priv<Radv[i]){
          local_max_rad_priv=Radv[i];
          local_id_mayor_priv=i;}
    #pragma omp critical
    {
    if(local_max_rad<local_max_rad_priv){
      local_max_rad=local_max_rad_priv;
      local_id_mayor=local_id_mayor_priv;}
    }
    }


//  printf("local_id_mayor=%d\n",local_id_mayor);fflush(stdout);


    if(local_id_mayor>=0){		//	esta cercado por voids?
      local_x=Xv[local_id_mayor];	//	pero no lo suficentemente cercano a uno de ellos?
      local_y=Yv[local_id_mayor];
      local_z=Zv[local_id_mayor];
      local_r=Radv[local_id_mayor];

      dist1=1.*nc;
      local_id_cerca=-2;
      #pragma omp parallel
      {
      float dx,dy,dz,dist_priv,dist1_priv=1.*nc;
      int local_id_cerca_priv=-2;
      #pragma omp for
      for(i=0;i<NumHuecos;i++){		//	distancia al vecino más cercano
        if((IDV[i]<=local_id_void)&&(i!=local_id_mayor)&&(IDV[i]>=0)){
          dx=local_x-Xv[i];				//	solo voids ya detectados
               if(dx>nc*0.5)	dx-=1.*nc;		//	pero diferentes del actual candidato
          else if(dx<-nc*0.5)	dx+=1.*nc;		//	y diferentes de los fallidos
          dy=local_y-Yv[i];
               if(dy>nc*0.5)	dy-=1.*nc;
          else if(dy<-nc*0.5)	dy+=1.*nc;
          dz=local_z-Zv[i];
               if(dz>nc*0.5)	dz-=1.*nc;
          else if(dz<-nc*0.5)	dz+=1.*nc;
          dist_priv=sqrt(dx*dx+dy*dy+dz*dz)-Radv[i];
          if(dist1_priv>dist_priv){		//	hallo el void más cercano
            dist1_priv=dist_priv;		//	es la distancia al void más cercano
            local_id_cerca_priv=IDV[i];		//	es el hueco void más cercano
          }
        }
      }
      #pragma omp critical
      {
        if(dist1>dist1_priv){				//	hallo el vecino más cercano
          dist1=dist1_priv;				//	es la distancia al vecino más cercano
          local_id_cerca=local_id_cerca_priv;		//	es el void más cercano
        }
      }
      }

      dist2=1.*nc;
      if(local_id_void>=0){
        #pragma omp parallel
        {
        float dx,dy,dz,dist_priv,dist2_priv=1.*nc;
        #pragma omp for
        for(i=0;i<NumHuecos;i++){		//	distancia al segundo vecino más cercano
          if((IDV[i]<=local_id_void)&&(IDV[i]!=local_id_cerca)&&(i!=local_id_mayor)&&(IDV[i]>=0)){
            dx=local_x-Xv[i];				//	solo voids ya detectados
                 if(dx>nc*0.5)	dx-=1.*nc;		//	excluyendo el más cercano
            else if(dx<-nc*0.5)	dx+=1.*nc;		//	y excluyendo a sí mismo
            dy=local_y-Yv[i];				//	y excluyendo los fallidos
                 if(dy>nc*0.5)	dy-=1.*nc;
            else if(dy<-nc*0.5)	dy+=1.*nc;
            dz=local_z-Zv[i];
                 if(dz>nc*0.5)	dz-=1.*nc;
            else if(dz<-nc*0.5)	dz+=1.*nc;
            dist_priv=sqrt(dx*dx+dy*dy+dz*dz)-Radv[i];
            if(dist2_priv>dist_priv){			//	hallo el segundo vecino más cercano
              dist2_priv=dist_priv;
            }		//	if
          }		//	if
        }		//	for
        #pragma omp critical
        {
          if(dist2>dist2_priv){				//	hallo el segundo void más cercano
            dist2=dist2_priv;				//	es la distancia al segundo void cercano
          }
        }
        }
      }			//	if


      //  si no hay nadie cerca, entonces creo uno nuevo
      if((dist1>weak_overlap*local_r)&&(dist2>weak_overlap*local_r))
        {
        //	a menos que sea muy pequeño, en cuyo caso, lo desecho
        //	o si fue generado por los randoms
        if((local_r<rad_min_uno)||(ID_random[local_id_mayor]<0))
          IDV[local_id_mayor]=-1;

        //	pero si es lo suficientemente grande, entonces creo uno nuevo
        else{
          local_id_void++;
          IDV[local_id_mayor]=local_id_void;
  	  //	catalogo de voids trazados por la esfera central
          fprintf(ES,"%f %f %f %f\n",local_x*Lbox/nc,local_y*Lbox/nc,local_z*Lbox/nc,local_r*Lbox/nc);
        }	//	else: es lo suficientemente grande para ser un void
      }		//	else: no hay nadie cerca


      // si hay solo uno cerca, lo adhiero
      else if((dist1<strong_overlap*local_r)&&(dist2>weak_overlap*local_r)&&(only_spherical==0))
        IDV[local_id_mayor]=local_id_cerca;

      //  si hay dos voids cercanos -> no nos sirve para nada
      else
      //else if((dist1<overlap*local_r)&&(dist2<overlap*local_r))
        IDV[local_id_mayor]=-1;




    }		//	if
  }		//	ciclo principal
  fclose(ES);

  printf("local_id_void=%d\n",local_id_void);fflush(stdout);
  printf("depura...hecho\n");fflush(stdout);

}







// print voids file
void imprime_voids(void){
  printf("\nimprime_voids...\n");fflush(stdout);
  int i;
  FILE * VO;
  VO=fopen(catalogo_voids,"w+");
  for(i=0;i<NumHuecos;i++)
    if(IDV[i]>=0)
      fprintf(VO,"%f %f %f %f %d\n",Xv[i]*Lbox/nc,Yv[i]*Lbox/nc,Zv[i]*Lbox/nc,Radv[i]*Lbox/nc,IDV[i]);
  fclose(VO);
  printf("imprime_voids...hecho\n");fflush(stdout);
}







// allocate random and ellipticity variables
void alloca_ran_elip(int NumRan){
  int b1,i;
  float a1,a2,a3,a4;

  //	cuenta los voids(familias de huecos) y huecos
  FILE * VO;
  VO=fopen(catalogo_voids,"r");
  NumHuecos=0;
  NumVoids=0;
  while(fscanf(VO,"%f %f %f %f %d\n",&a1,&a2,&a3,&a4,&b1)!=EOF){
    Xv[NumHuecos]=a1*nc/Lbox;
    Yv[NumHuecos]=a2*nc/Lbox;
    Zv[NumHuecos]=a3*nc/Lbox;
    Radv[NumHuecos]=a4*nc/Lbox;
    IDV[NumHuecos]=b1;			//	bi es el índice de la familia a la que pertecence
    if(NumVoids<b1)
      NumVoids=b1;
    NumHuecos++;}
  NumVoids++;
  fclose(VO);

  //	alloca el vector que va a contener la cantidad de huecos en cada familia
  if(!(Ocu = (int*) calloc (NumVoids , sizeof(int)))){
      fprintf(stderr, "failed to allocate memory Ocu.\n");
      fflush(stderr);
      exit(0);}

  for(i=0;i<NumVoids;i++)		//	cerando la cantidad de huecos en cada familia
    Ocu[i]=0;

  // asigna los niveles de ocupacion, la cantidad de huecos de cada familia
  for(i=0;i<NumHuecos;i++)
    Ocu[IDV[i]]++;			//	IDV[i] es el índice de la familia a la que pertenece el hueco [i]


  // alloca los Id: van a ser los indices de los huecos de cada familia Id[familias][integrantes]
  if(!(Id = (int **) malloc (NumVoids * sizeof(int *)))){
      fprintf(stderr, "failed to allocate memory Id*.\n");
      fflush(stderr);
      exit(0);}

  for(int p = 0; p < NumVoids; p++) {
    if(!(Id[p] = (int *) malloc (Ocu[p] * sizeof(int)))){
      fprintf(stderr, "failed to allocate memory Id[%d].\n",p);
      fflush(stderr);
      exit(0);}
   }

  for(i=0;i<NumVoids;i++)
    Ocu[i]=0;

  // asigna las ocupaciones, lleno Id, (es más facil llenando de nuevo Ocu)
  for(i=0;i<NumHuecos;i++){
    Id[IDV[i]][Ocu[IDV[i]]]=i;
    Ocu[IDV[i]]++;
  }
  // ya tengo los Id: van a ser los indices de los huecos de cada familia Id[familias][integrantes]

  //	alloca los random que necesito para volumen y elipticidad
  if(!(randomxin = (float*) calloc (NumRan , sizeof(float)))){
      fprintf(stderr, "failed to allocate memory randomxin.\n");
      fflush(stderr);
      exit(0);}
  if(!(randomyin = (float*) calloc (NumRan , sizeof(float)))){
      fprintf(stderr, "failed to allocate memory randomyin.\n");
      fflush(stderr);
      exit(0);}
  if(!(randomzin = (float*) calloc (NumRan , sizeof(float)))){
      fprintf(stderr, "failed to allocate memory randomzin.\n");
      fflush(stderr);
      exit(0);}
}


// computes eigen-values for ellipticity estimation
void calcula_auto_valores(void){
  Symmeig s(tensor,false);		//	false == solo autovalores
  valor[0]=(float) s.d[0];
  valor[1]=(float) s.d[1];
  valor[2]=(float) s.d[2];
  //	ordeno los aouto-valores, el mayor/menor queda guardado en auto[2]/[0]
  float aux_orden;
  int j,k;
  for(j=0;j<2;j++){
    for(k=j+1;k<3;k++){
      if(valor[j]>valor[k]){
        aux_orden=valor[j];
        valor[j]=valor[k];
        valor[k]=aux_orden;
      }
    }
  }
}






//	computes center, vaolume and ellipticity for non-spherical voids
void volumen_elipticidad(void){
  Ran myran(601);                //      random: pagina 343 nr3  (semilla)
  printf("\ncalcula volumen...\n");fflush(stdout);
  int i,j,k;


  // calcula volumen y centro geométrico usando randoms
  float volumen;
  float radio;
  int esta_dentro;
  int NumRan=100000;


  float randomx,randomy,randomz;
  int contador;

  alloca_ran_elip(NumRan);	//	ahi voy a guardan los randoms in void para la elipticidad
  //	tambien leo el archivo union, alloco y lleno las variables que contienen los indices de familias e integrantes

  float inf[3];
  float sup[3];
  float inf_tmp,sup_tmp,VolCaja;
  float x,y,z,dx,dy,dz,dist;
  FILE * UN;
  UN=fopen(catalogo_union,"w+");	//	voids == familias de esferas
  FILE * EL;
  EL=fopen(semi_ejes,"w+");		//	semiejes de las familias


  //	tensor de inercia y autovalores
  float Ixx,Iyy,Izz,Ixy,Iyz,Izx;
  float aux,auy,auz;

  for(i=0;i<NumVoids;i++){

    //	calculando los limites fisicos de la caja más pequeña que contiene a cada void
    if(Ocu[i]>1){
      for(k=0;k<3;k++){
        inf[k]=3.*nc;
        sup[k]=-3.*nc;}
      for(j=0;j<Ocu[i];j++){
        // condiciones de frontera periodicas:
        // llevo cada centro de hueco al lugar más cercano del primero
        dx=Xv[Id[i][0]]-Xv[Id[i][j]];
        dy=Yv[Id[i][0]]-Yv[Id[i][j]];
        dz=Zv[Id[i][0]]-Zv[Id[i][j]];
             if(dx<-nc*0.5)
               Xv[Id[i][j]]-=1.*nc;
        else if(dx>nc*0.5)
               Xv[Id[i][j]]+=1.*nc;

             if(dy<-nc*0.5)
               Yv[Id[i][j]]-=1.*nc;
        else if(dy>nc*0.5)
               Yv[Id[i][j]]+=1.*nc;

             if(dz<-nc*0.5)
               Zv[Id[i][j]]-=1.*nc;
        else if(dz>nc*0.5)
               Zv[Id[i][j]]+=1.*nc;

        // limites de la caja que contiene al void
        inf_tmp=Xv[Id[i][j]]-Radv[Id[i][j]];
        sup_tmp=Xv[Id[i][j]]+Radv[Id[i][j]];
        if(inf[0]>inf_tmp)
          inf[0]=inf_tmp;
        if(sup[0]<sup_tmp)
          sup[0]=sup_tmp;
        
        inf_tmp=Yv[Id[i][j]]-Radv[Id[i][j]];
        sup_tmp=Yv[Id[i][j]]+Radv[Id[i][j]];
        if(inf[1]>inf_tmp)
          inf[1]=inf_tmp;
        if(sup[1]<sup_tmp)
          sup[1]=sup_tmp;
        
        inf_tmp=Zv[Id[i][j]]-Radv[Id[i][j]];
        sup_tmp=Zv[Id[i][j]]+Radv[Id[i][j]];
        if(inf[2]>inf_tmp)
          inf[2]=inf_tmp;
        if(sup[2]<sup_tmp)
          sup[2]=sup_tmp;
      }		//	limites de la caja
 


      // volumen de la caja
      VolCaja=(sup[0]-inf[0])*(sup[1]-inf[1])*(sup[2]-inf[2]);



      //	calcula el centro geométrico y el volumen a partir de los randoms in void
      contador=0;		// volumen del void
      x=0.;			// posiciones del centro geométrico
      y=0.;
      z=0.;
      for(k=0;k<NumRan;k++){
        // llamo numeros aleatorios uniformemente dentro de esta caja:randomx randomy randomz
        randomx=inf[0]+myran.doub()*(sup[0]-inf[0]);
        randomy=inf[1]+myran.doub()*(sup[1]-inf[1]);
        randomz=inf[2]+myran.doub()*(sup[2]-inf[2]);
        esta_dentro=0;
        for(j=0;j<Ocu[i];j++){
          dx=randomx-Xv[Id[i][j]];
          dy=randomy-Yv[Id[i][j]];
          dz=randomz-Zv[Id[i][j]];
          dist=sqrt(dx*dx+dy*dy+dz*dz);
          if(dist<=Radv[Id[i][j]]){		//	random in void
            esta_dentro=1;
          }
        }
        if(esta_dentro==1){
          x+=randomx;
          y+=randomy;
          z+=randomz;
          randomxin[contador]=randomx;
          randomyin[contador]=randomy;
          randomzin[contador]=randomz;
          contador++;
        }
      }		//	k: randoms

      volumen=(float)contador;
      x/=volumen;
      y/=volumen;
      z/=volumen;
      volumen*=VolCaja/NumRan;		//	en unidades de nc
      radio=pow(3./pi*0.25*volumen,1./3.);

      fprintf(UN,"%f %f %f %f\n",
                 fmod(x+1.*nc,1.*nc)*Lbox/nc,
                 fmod(y+1.*nc,1.*nc)*Lbox/nc,
                 fmod(z+1.*nc,1.*nc)*Lbox/nc,
                 radio*Lbox/nc);     


      //	ya está el volumen, sigue la elipticidad
      Ixx=0.;	Iyy=0.;	Izz=0.;	Ixy=0.;	Iyz=0.;	Izx=0.;

      for(k=0;k<contador;k++){		//	número de randoms in void
        aux=randomxin[k]-x;
        auy=randomyin[k]-y;
        auz=randomzin[k]-z;
        Ixx+=auy*auy+auz*auz;
        Iyy+=auz*auz+aux*aux;
        Izz+=aux*aux+auy*auy;
        Ixy-=aux*auy;
        Iyz-=auy*auz;
        Izx-=auz*aux;
      }

      tensor[0][0]=(double) Ixx;
      tensor[1][1]=(double) Iyy;
      tensor[2][2]=(double) Izz;
      tensor[0][1]=(double) Ixy;      tensor[1][0]=tensor[0][1];
      tensor[1][2]=(double) Iyz;      tensor[2][1]=tensor[1][2];
      tensor[2][0]=(double) Izx;      tensor[0][2]=tensor[2][0];

      calcula_auto_valores();

      //	vueno, vamos a imprimir esos tres valores, para despues hacer calculos con ellos
      fprintf(EL,"%f %f %f %f\n",radio*Lbox/nc,valor[2],valor[1],valor[0]);

    }						//	if (no esferico) tiene más de dos huecos
    else{					//	else (esferico)
      fprintf(UN,"%f %f %f %f\n",		//	coordenadas y radio
                 Xv[Id[i][0]]*Lbox/nc,
                 Yv[Id[i][0]]*Lbox/nc,
                 Zv[Id[i][0]]*Lbox/nc,
                 Radv[Id[i][0]]*Lbox/nc);

      fprintf(EL,"%f %f %f %f\n",Radv[Id[i][0]]*Lbox/nc,1.,1.,1.);	//	elipticidad
    }

  }		//	ciclo sobre los voids

  fclose(UN);
  fclose(EL);
  printf("calcula volumen...hecho\n");fflush(stdout);
}














int main(int argc, char **argv){
  printf("comienza\n");
  fflush(stdout);
  omp_set_num_threads(32);


  FILE * ESC;	//	catalogo de trabajo
  FILE * RES;	//	catalogo de respaldo (no HD)
  int l;
  float x,y,z,radio;

  archivos();		//	carga el nombre de los archivos

  if(solo_depura==0){
    allocar();
    lee();

    if(solo_HD==1)
      lee_huecos();
    else{				//	calcula prixelizado y refina
      rad_min_actual=rad_min_uno;
      densidad();
      if(calcula==1)
        crece();
      else
        lee_radios();

      // search for the bigest cell sphere
      l=busco();	//	es el índice de la célula con la mayor esfera subdensa
      printf("ya buscó la célula con el rádio más grande\n");
      fflush(stdout);

      // de aquí en adelante x, y, z son las coordenadas del centro del cada void
      radio=Rad[l];printf("su radio es %f\n",Rad[l]);fflush(stdout);
      Rad[l]=-1.;

      refino(l,x,y,z,radio);		//	improves the sphere l by using particles positions

      printf("radio refinado del void más grande = %f\n",radio);
      printf("centro refinado del void más grande = (%f, %f, %f)\n\n",x,y,z);  fflush(stdout);

      // escribo en los dos catalogos (trabajo y respaldo)
      ESC = fopen(catalogo_huecos,"w+");
      RES = fopen(catalogo_huecos_respaldo,"w+");
      fprintf(ESC,"%f %f %f %f\n",x*Lbox/nc,y*Lbox/nc,z*Lbox/nc,radio*Lbox/nc);
      fprintf(RES,"%f %f %f %f\n",x*Lbox/nc,y*Lbox/nc,z*Lbox/nc,radio*Lbox/nc);
      fclose(ESC);
      fclose(RES);

      id_void++;
      Xv[id_void]=x;	//	en unidades de nc
      Yv[id_void]=y;
      Zv[id_void]=z;
      Radv[id_void]=radio;

      anula(x,y,z,radio);		// set to 0 the radious of spheres inside sphere l

      FILE * WR;
      WR=fopen(particulas_voids,"w+");
      fclose(WR);
      //  imprime(x,y,z,radio);

      printf("ya encontramos, refinamos, anulamos e imprimímos el Void 0\n\n");
      fflush(stdout);








      int condicion=1;		//	keep going
      printf("entra al ciclo principal\n");
      fflush(stdout);

      while(condicion>=0){		//	main cycle

        // search for the next bigest sphere
        condicion=busco();


        if(condicion==-1){
          printf("everything is done for now\n\n");
          fflush(stdout);
        }

        else{
          radio=Rad[condicion];
          Rad[condicion]=-1.;
          refino(condicion,x,y,z,radio);	//	refina la esfera de la célula l=condicion

          if(radio>rad_min_actual){
            // escribo en los dos catalogos (trabajo y respaldo)
            ESC = fopen(catalogo_huecos,"a");
            RES = fopen(catalogo_huecos_respaldo,"a");
	    fprintf(ESC,"%f %f %f %f\n",x*Lbox/nc,y*Lbox/nc,z*Lbox/nc,radio*Lbox/nc);
	    fprintf(RES,"%f %f %f %f\n",x*Lbox/nc,y*Lbox/nc,z*Lbox/nc,radio*Lbox/nc);
            fclose(ESC);
            fclose(RES);
            id_void++;
            Xv[id_void]=x;	//	en unidades de nc
            Yv[id_void]=y;
            Zv[id_void]=z;
            Radv[id_void]=radio;

            anula(x,y,z,radio);		// set to 0 the radious of spheres inside sphere l
          }
          else
            Rad[condicion]=-1.;		
        }	//	cierra el else
      }	//	cierra el ciclo principal


    }	//	cierra el calculo prixelizado + refinado (else solo_HD=0)

















    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    //	DH improvement
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    if(no_HD==0){
      printf("comienza HD\n");fflush(stdout);
      rad_min_actual=rad_min_dos;
      free_pixelizado();
      if(ya_radios_randoms==0){

        if(ya_catalogo_randoms==0){
          cuenta_halos();
          alloca_halos();
          lee_halos();
 
          genera_catalogo_random();
          free_halos();
        }

        cuenta_catalogo_random();
        alloca_random();
        lee_catalogo_random();

        crece_random();
        imprime_random();	//	incluye los radios calculados en crece_random()
      }				//	cierra ya randoms

      else{			//	else: calcula randoms
        cuenta_catalogo_random();
        alloca_random();
        lee_random();		//	ya incluye los radios calculados en crece_random()
      }

      int condicion=1;		//	permanencia dentro del ciclo de alta resolución
      printf("entra al ciclo de alta resolución\n");
      fflush(stdout);




      while(condicion>=0){		//	HD cycle

        condicion=busco_random();    // search for the bigest random sphere


        if(condicion==-1){
          printf("en este punto deberia salir del ciclo de alta resolución\n\n");
          fflush(stdout);
        }
        else{
          radio=Rr[condicion];
//          printf("Rr[%d]=%f\n",condicion,Rr[condicion]);fflush(stdout);
          r1=2;
          refino(condicion,x,y,z,radio);	//	refina la esfera random com índice l=condicion
          r1=0;

          if((radio>rad_min_actual)&&(Rr[condicion]>rad_min_actual)){
            //	solo escribo en el archivo de trabajo (no en el de respaldo = sin HD)
            ESC = fopen(catalogo_huecos,"a");
            fprintf(ESC,"%f %f %f %f\n",x*Lbox/nc,y*Lbox/nc,z*Lbox/nc,radio*Lbox/nc);
            fclose(ESC);
            id_void++;
            Xv[id_void]=x;	//	en unidades de nc
            Yv[id_void]=y;
            Zv[id_void]=z;
            Radv[id_void]=radio;
            anula_random(x,y,z,radio);		//	anula el rádio de los randoms \in void


          }
          else
            Rr[condicion]=-1;			//	ya no es candidato a centro de void

        }	//	cierra el else
      }		//	cierra el ciclo principal
      printf("de pura casualidad llega a este punto?\n");fflush(stdout);
    }		//	if no_HD==0

    printf("de pura casualidad llega a este punto? 2?\n");fflush(stdout);
    free_xyz();	//	libera X Y Z ID

  }	//	acaba solo_depura==0





  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  //    gathering spheres into families
  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  printf("depurando\n");fflush(stdout);

  cuenta_huecos();	//	lee el archivo de coordenadas y radios de los huecos
  realloca();		//	realloca Xv Yv Zv, y llena dos ID
  lee_huecos_2();	//	lee el archivo de coordenadas y radios de los huecos

  if(solo_vol_elip==0){
    depura();
    imprime_voids();
  }

  //  volumen_elipticidad();//	lee el archivo de voids: huecos que se sobrelapan
			//	calcula volumen y centro de cada unión
			//	imprime centro y radio esferico equivalente

  printf("este es el fin del camino\n");fflush(stdout);



  return 0;

}



