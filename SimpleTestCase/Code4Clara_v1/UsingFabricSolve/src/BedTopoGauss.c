#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define min(A,B) (( A < B ) ? A : B)
#define max(A,B) (( A > B ) ? A : B)
#define PI 3.1415926535897932


double hb_fun(double, double),H_fun(double, double );
double Interpolate3D(double, double, int, double *, double *, double *, double, double ,double);

int nsdata;
double *xsdata,*ysdata,*hsdata;
int nbdata;
double *xbdata,*ybdata,*hbdata;

double H=100, Hmin=1.0e-3;
double xBumpb=790500,yBumpb=155006.25,ABumpb=1500.0,SxBumpb=1.5e3,SyBumpb=1e3;
double xBumpH=790500,yBumpH=155006.25,ABumpH=10.0,SxBumpH=1e3,SyBumpH=0.75e3;

int main(int argc,char *argv[]){

  FILE *fichero,*ficheronew;
  int nm,i,particiones,par;
  int basurai,*uno,*dos;
  double basuraf,*x,*y,*z,zz,hb,hs;
  char name[100],namenew[100],command[100];


  /////////////////////////////////////////////////////////////////////////

  if(argc!=2){
    printf("Esta version del programa requiere el numero de particiones. Adios!\n");
    exit(-1);
  }

  particiones=atoi(argv[1]);


  if(particiones>1){

    nm=0;
    for(par=0;par<particiones;par++){
      sprintf(name,"mesh/partitioning.%d/part.%d.nodes",particiones,par+1);
      if((fichero=fopen(name,"r"))==NULL){
	printf("No se puede leer %s (0)\n",name);
	exit(-1);
      }
      while(fscanf(fichero,"%d %d %lg %lg %lg",&basurai,&basurai,&basuraf,&basuraf,&basuraf)!=EOF){
	nm++;
      }
      fclose(fichero);
    }


    uno=(int *)calloc(nm,sizeof(int));
    dos=(int *)calloc(nm,sizeof(int));
    x=(double *)calloc(nm,sizeof(double));
    y=(double *)calloc(nm,sizeof(double));
    z=(double *)calloc(nm,sizeof(double));

    for(i=0,par=0;par<particiones;par++){
      sprintf(name,"mesh/partitioning.%d/part.%d.nodes",particiones,par+1);
      sprintf(namenew,"mesh/partitioning.%d/part.%d.nodes_new",particiones,par+1);
      if((fichero=fopen(name,"r"))==NULL){
	printf("No se puede leer %s\n",name);
	exit(-1);
      }
      if((ficheronew=fopen(namenew,"w"))==NULL){
	printf("No se puede escribir %s\n",namenew);
	exit(-1);
      }
      while(fscanf(fichero,"%d %d %lg %lg %lg",&uno[i],&dos[i],&x[i],&y[i],&z[i])!=EOF){     
	hb=hb_fun(x[i],y[i]);
	hs=hb+max(Hmin,H_fun(x[i],y[i]));
	z[i]=(1.0-z[i]/H)*hb+(z[i]/H)*hs;

	fprintf(ficheronew,"%d %d %.8lg %.8lg %.8lg\n",uno[i],dos[i],x[i],y[i],z[i]);
	i++;
      }
      fclose(fichero);
      fclose(ficheronew);
    
      sprintf(command,"mv %s %s",namenew,name);
      system(command);
    }    
  }
  else{

    if((fichero=fopen("mesh/mesh.nodes","r"))==NULL){
      printf("No se puede leer mesh.nodes\n");
      exit(-1);
    }
 
    nm=0;
    while(fscanf(fichero,"%d %d %lg %lg %lg",&basurai,&basurai,&basuraf,&basuraf,&basuraf)!=EOF){
      nm++;
    }

    fclose(fichero);

    uno=(int *)calloc(nm,sizeof(int));
    dos=(int *)calloc(nm,sizeof(int));
    x=(double *)calloc(nm,sizeof(double));
    y=(double *)calloc(nm,sizeof(double));
    z=(double *)calloc(nm,sizeof(double));

    if((fichero=fopen("mesh/mesh.nodes","r"))==NULL){
      printf("No se puede leer mesh.nodes\n");
      exit(-1);
    }

    for(i=0;i<nm;i++){
      fscanf(fichero,"%d %d %lg %lg %lg",&uno[i],&dos[i],&x[i],&y[i],&z[i]);
    }
    fclose(fichero);


    if((fichero=fopen("mesh/mesh.nodes","w"))==NULL){
      printf("No se puede escribir mesh.nodes\n");
      exit(-1);
    }
    
    for(i=0;i<nm;i++){
      hb=hb_fun(x[i],y[i]);
      hs=hb+max(Hmin,H_fun(x[i],y[i]));      
      z[i]=(1.0-z[i]/H)*hb+(z[i]/H)*hs;
      
      fprintf(fichero,"%d %d %.8lg %.8lg %.8lg\n",uno[i],dos[i],x[i],y[i],z[i]);
  }
  fclose(fichero);
  }

  return 0;
};


double H_fun(double xx, double yy)
{
  double zz;

  zz=ABumpH*exp(-((xx-xBumpH)/SxBumpH*(xx-xBumpH)/SxBumpH+(yy-yBumpH)/SyBumpH*(yy-yBumpH)/SyBumpH));
  return(zz);
};

double hb_fun(double xx, double yy)
{
  double zz;

  zz=ABumpb*exp(-((xx-xBumpb)/SxBumpb*(xx-xBumpb)/SxBumpb+(yy-yBumpb)/SyBumpb*(yy-yBumpb)/SyBumpb));
  return(zz);
};
