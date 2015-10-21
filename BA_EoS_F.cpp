#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "BA_EoS_F.h"
#define zero 1.e-10

#define _ORDER_  'epnm'
#define MN  939.603 
FILE* infile;
EoS_F::EoS_F(void)
{
//	mu_min = 939.603 ; // in case it is not define in the EoS file.
	EOSTABFILEsize= peos;
}

EoS_F::~EoS_F(void)
{
}

double aprfu_z(double x,double *X,double *Y,double *dYdX,int size)
{int c=1;double _a_,_b_;
while (x<X[c]&&c<size) c++;
*dYdX=_a_=(Y[c]-Y[c-1])/(X[c]-X[c-1]);
_b_=Y[c-1]-_a_*X[c-1];
return _a_*x+_b_;
}


void EoS_F::Read_EoS_FILE_(char* EOSTABFILE)
{
	/*Reading Tableised EoSfile*/
	//double* metrMJUB,*metrNB,*metrP;
	
	int contr,sw; 
	double *inoNB, *inoP, *inomju, *temp;
	double dmju_start,
	mju_start=2994.,
	mju_end=931.05,
	Nc,Pc;

	/*metrMJUB=(double *)malloc((peos+2)*sizeof(double));	
	metrNB=(double *)malloc((peos+2)*sizeof(double)); 
	metrP=(double *)malloc((peos+2)*sizeof(double));
	*/
    temp=(double *)malloc(15*sizeof(double));

	if((infile =fopen(EOSTABFILE, "r"))==NULL){
		printf("There is no file %s to open\n",EOSTABFILE);
		exit(1);
		};
	char k;
	int FBW_key = 0; 
	while (((k=fgetc(infile))!=EOF) && (k!='\n'));
	sw=0;
	//	for (sw=0;sw<=EOSTABFILEsize; sw++)
	double sp,sn,sm,se,spi=0;	
do{ 
	switch (_ORDER_){
		case 'mpne' :k = fscanf(infile,"%lg %lg %lg %lg",&sm,&sp,&sn,&se); break;
		default : k = fscanf(infile,"%lg %lg %lg %lg",&se,&sp,&sn,&sm);
	};
	/*if(sw == 3)*/ 
//	printf("%g   %g  \n",sm,sp);
    if(sw == 4&&spi < sp) FBW_key = 1;
	spi = sp;
	sw++;
	} while(k!=EOF&&sw<2000);
	fclose(infile);
	EOSTABFILEsize = sw-1;
	infile =fopen(EOSTABFILE, "r");

    inomju=(double *)malloc((EOSTABFILEsize+2)*sizeof(double));
	inoNB =(double *)malloc((EOSTABFILEsize+2)*sizeof(double)); 
	inoP  =(double *)malloc((EOSTABFILEsize+2)*sizeof(double));
	double eee,xxx;	
	int sww;
	while (((k=fgetc(infile))!=EOF) && (k!='\n'));
	for (sw=0;sw<=EOSTABFILEsize; sw++)
	{
		if(FBW_key) sww = EOSTABFILEsize-sw-1; else sww = sw;
		switch(_ORDER_){
		case 'mpne': fscanf(infile,"%lg %lg %lg %lg",&inomju[sww],&inoP[sww],&inoNB[sww],&eee); break;
		case 'mnep': fscanf(infile,"%lg %lg %lg %lg",&inomju[sww],&inoNB[sww],&eee,&inoP[sww]); break;
		//inomju[sww]=(eee+inoP[sww])/inoNB[sww];
		//case 'mepxn':	fscanf(infile,"%lg %lg %lg %lg %lg",&inomju[sw],&eee,&inoP[sw],&xxx,&inoNB[sw]); break;
		case 'npem': fscanf(infile,"%lg %lg %lg %lg",&inoNB[sw],&inoP[sw],&eee,&inomju[sw]); break;
		default   :	 fscanf(infile,"%lg %lg %lg %lg",&eee,&inoP[sww],&inoNB[sww],&inomju[sww]);
		}
	}
	fclose(infile);
	free(temp);


	mju_start=inomju[0];
	mju_end=inomju[EOSTABFILEsize-1];
	if(MN>mju_end)mju_end = MN;
	if (mju_end>mju_start) {dmju_start=mju_start; mju_start=mju_end; mju_end=dmju_start;}

	dmju_start=(mju_end-mju_start)/peos;
	muf[0]=mju_start;
	/*Tableisation of EoS*/
	for (contr=0;contr<=peos;contr++)
	{double temp;
	nf[contr]=aprfu_z(muf[contr],inomju,inoNB,&temp,EOSTABFILEsize);
	Pf[contr]=aprfu_z(muf[contr],inomju,inoP,&temp,EOSTABFILEsize);
	epsf[contr] = nf[contr]*muf[contr] - Pf[contr]; 
	muf[contr+1]=muf[contr]+dmju_start;
	}
	free(inoNB);free(inoP);free(inomju);

}

double	EoS_F::get_n_B(double e_)
{double nn,a;
nn= aprfu_z(e_,epsf,nf,&a,peos);
return nn;
}

double	EoS_F::get_P(double e_)
{double pp,a;
pp= aprfu_z(e_,epsf,Pf,&a,peos);
return pp;
}

double 	EoS_F::get_mu(double e_)
{double mumu,a;
mumu= aprfu_z(e_,epsf,muf,&a,peos);
return mumu;
}

