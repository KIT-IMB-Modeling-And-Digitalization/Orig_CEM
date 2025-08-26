#include <stdio.h>
#include <math.h>

/* This software was developed at the National Institute of */
/* Standards and Technology by employees of the Federal */
/* Government in the course of their official duties. Pursuant */
/* to title 17 Section 105 of the United States Code this */
/* software is not subject to copyright protection and is in */
/* the public domain. CEMHYD3D is an experimental system. NIST */
/* assumes no responsibility whatsoever for its use by other */
/* parties, and makes no guarantees, expressed or implied, */
/* about its quality, reliability, or any other characteristic. */
/* We would appreciate acknowledgement if the software is used. */
/* This software can be redistributed and/or modified freely */
/* provided that any derivative works bear some notice that */
/* they are derived from it, and any modified versions bear */
/* some notice that they have been modified. */

#define POROSITY 0
#define C3S 1
#define C2S 2
#define C3A 3
#define C4AF 4
#define GYPSUM 5
#define HEMIHYD 6
#define ANHYDRITE 7
#define POZZ 8
#define INERT 9
#define SLAG 10
#define ASG 11		  /* aluminosilicate glass */
#define CAS2 12
#define CH 13
#define CSH 14
#define C3AH6 15
#define ETTR 16
#define ETTRC4AF 17     /* Iron-rich stable ettringite phase */
#define AFM 18
#define FH3 19
#define POZZCSH 20
#define SLAGCSH 21
#define CACL2 22
#define FREIDEL 23	  /* Freidel's salt */
#define STRAT 24	  /* stratlingite (C2ASH8) */
#define GYPSUMS 25	/* Gypsum formed from hemihydrate and anhydrite */
#define CACO3 26
#define AFMC 27
#define INERTAGG 28
#define ABSGYP 29
#define FLYASH 30
#define EMPTYP 45         /* Empty porosity due to self desiccation */

main(){
	FILE *infile,*outfile;
	char filein[80],fileout[80];
	int valin;
	int valout,i1,j1;
	int i,nskip,dx,dy,j,scf,iscale,dxtot,dytot;
	static short int image[1000][1000];
	static short int red[50],green[50],blue[50];

	for(i1=0;i1<50;i1++){
		red[i1]=0;
		green[i1]=0;
		blue[i1]=0;
	}
	for(i1=11;i1<28;i1++){
		red[i1]=192;
		green[i1]=255;
		blue[i1]=128;
	}
	red[CH]=0;
	green[CH]=0;
	blue[CH]=255;
	red[CSH]=255;
	green[CSH]=128;
	blue[CSH]=0;
	red[C3S]=255;
	green[C2S]=255;
	blue[C2S]=255;
	green[C3A]=255;
	red[C4AF]=255;
        green[C4AF]=255;
	red[GYPSUM]=128;
	green[GYPSUM]=128;
	blue[GYPSUM]=128;
	red[HEMIHYD]=64;
	green[HEMIHYD]=64;
	blue[HEMIHYD]=64;
	red[ANHYDRITE]=192;
	green[ANHYDRITE]=192;
	blue[ANHYDRITE]=192;
	red[POZZ]=64;
	green[POZZ]=128;
	blue[POZZ]=128;
	red[INERT]=128;
	green[INERT]=128;
	blue[INERT]=255;
	red[POZZCSH]=255;
	green[POZZCSH]=128;
	blue[POZZCSH]=128;
	red[INERTAGG]=255;
	blue[INERTAGG]=255;
        red[CACO3]=64;
        green[CACO3]=128;
        blue[CACO3]=255;
	red[FLYASH]=255;
	green[FLYASH]=255;
	blue[FLYASH]=255;
        red[SLAGCSH]=255;
        green[SLAGCSH]=0;
        blue[SLAGCSH]=128;
        red[SLAG]=0;
        green[SLAG]=128;
        blue[SLAG]=0;
	red[EMPTYP]=255;
	green[EMPTYP]=255;
	blue[EMPTYP]=255;
	
	printf("Enter name of file with raw (3-D image) data \n");
	scanf("%s",filein);
	printf("%s\n",filein);
	printf("Enter name of image file to create \n");
	scanf("%s",fileout);
	printf("%s\n",fileout);

	printf("Enter number of pixels to skip at start of image file \n");
	scanf("%d",&nskip);
	printf("%d\n",nskip);
	dx=dy=100;

	printf("Enter factor to scale image by\n");
	scanf("%d",&iscale);
	printf("%d\n",iscale);

	infile=fopen(filein,"r");
	outfile=fopen(fileout,"w");

	for(i=0;i<nskip;i++){
		fscanf(infile,"%d",&valin);
	}

	fprintf(outfile,"P3\n");
	fprintf(outfile,"%d %d\n",dx*iscale,dy*iscale);
	fprintf(outfile,"255\n");
	dxtot=dx*iscale;
	dytot=dy*iscale;
	for(j=0;j<dy;j++){
	for(i=0;i<dx;i++){

		fscanf(infile,"%d",&valin);
		valout=valin;
		for(j1=0;j1<iscale;j1++){
		for(i1=0;i1<iscale;i1++){
			image[i*iscale+i1][j*iscale+j1]=valout;
		}
		}
	}
	}
	fclose(infile);

	for(j=0;j<dytot;j++){
	for(i=0;i<dxtot;i++){
		fprintf(outfile,"%d %d %d\n",red[image[i][j]],green[image[i][j]],blue[image[i][j]]);
	}
	}

	fclose(outfile);
}
