#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define min(A,B) (( A < B ) ? A : B)
#define max(A,B) (( A > B ) ? A : B)
#define PI 3.1415926535897932


double hb_fun(double, double),hs_fun(double, double );
double Interpolate3D(double, double, int, double *, double *, double *, double, double ,double);

int nsdata;
double *xsdata,*ysdata,*hsdata;
int nbdata;
double *xbdata,*ybdata,*hbdata;

double H=100, Hmin=1.0;

int main(int argc,char *argv[]){

  FILE *fichero,*ficheronew;
  int nm,i,particiones,par;
  int basurai,*uno,*dos;
  double basuraf,*x,*y,*z,zz,hb,hs;
  char name[100],namenew[100],command[100];

  //
  //First We read the data from testcase.xyz
  //
  ////////////////////////////////////////////////////////////////////////
  if((fichero=fopen("testcase.xyz","r"))==NULL){
    printf("No se puede leertestcase.xyz\n");
    exit(-1);
  }
  nbdata=0;
  while(fscanf(fichero,"%lg %lg %lg",&basuraf,&basuraf,&basuraf)!=EOF){
    nbdata++;
  }
  fclose(fichero);

  xbdata=(double *)calloc(nbdata,sizeof(double));
  ybdata=(double *)calloc(nbdata,sizeof(double));
  hbdata=(double *)calloc(nbdata,sizeof(double));
  
  if((fichero=fopen("testcase.xyz","r"))==NULL){
    printf("No se puede leer testcase.xyz\n");
    exit(-1);
  }
  for(i=0;i<nbdata;i++){
    fscanf(fichero,"%lg %lg %lg",&xbdata[i],&ybdata[i],&hbdata[i]);
    //    printf("%lg %lg %lg\n",xbdata[i],ybdata[i],hbdata[i]);
  }
  fclose(fichero);

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
	hs=hb+Hmin;
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

      hs=hb+Hmin;
      z[i]=(1.0-z[i]/H)*hb+(z[i]/H)*hs;
      
      fprintf(fichero,"%d %d %.8lg %.8lg %.8lg\n",uno[i],dos[i],x[i],y[i],z[i]);
  }
  fclose(fichero);
  }

  return 0;
};


double hs_fun(double xx, double yy)
{
  return(Interpolate3D(xx,yy,nsdata,xsdata,ysdata,hsdata,1.0,200,50));
};

double hb_fun(double xx, double yy)
{
  return(Interpolate3D(xx,yy,nbdata,xbdata,ybdata,hbdata, -1.0,500, 10));
};

double Interpolate3D(double xi,double yi, int nn,double *xx,double *yy, double *zz,double r1, double r2, double rd)
{
  double zi=-6,z1=0.0,z2=0.0,w1=0,w2=0,r;
  int ii;

  for(ii=0;ii<nn;ii++){

    r=sqrt((xx[ii]-xi)*(xx[ii]-xi)+(yy[ii]-yi)*(yy[ii]-yi));
    
    if(r<r1){
      z1+=zz[ii];
      w1+=1.0;
    }
    else if(r<r2){
      z2+=zz[ii]*exp(-r/rd);
      w2+=exp(-r/rd);
    }
  }
  if(w1>0.0){
    zi=z1/w1;
  }
  else if(w2>0.0){
    zi=z2/w2;
  }
  
  return(zi);
};
