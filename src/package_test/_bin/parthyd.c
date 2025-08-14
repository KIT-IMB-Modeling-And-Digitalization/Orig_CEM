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

#include <stdio.h>  // Required for FILE, fopen, fprintf, fclose

// External global variables
// extern int SYSIZE;
extern char mic[SYSIZE][SYSIZE][SYSIZE];
extern char micorig[SYSIZE][SYSIZE][SYSIZE];
extern short int micpart[SYSIZE][SYSIZE][SYSIZE];


// Cement phase definitions
// extern int C3S, C2S, C3A, C4AF;

// Hydration process variables
extern int cyccnt;      // Hydration cycle count
extern double alpha_cur; // Current hydration level

// File handling
extern char phrname[80];  // Filename for hydration output


/* Routine to assess relative particle hydration */
void parthyd(){
	int norig[50000],nleft[50000];
	int ix,iy,iz;
	char valmic,valmicorig;
	int valpart,partmax;
        float alpart;
	FILE *phydfile;

        /* Initialize the particle count arrays */
	for(ix=0;ix<50000;ix++){
		nleft[ix]=norig[ix]=0;
	}
	phydfile=fopen(phrname,"a");
	fprintf(phydfile,"%d %f\n",cyccnt,alpha_cur);

	partmax=0;
        /* Scan the microstructure pixel by pixel and update counts */
	for(ix=0;ix<SYSIZE;ix++){
	for(iy=0;iy<SYSIZE;iy++){
	for(iz=0;iz<SYSIZE;iz++){

		if(micpart[ix][iy][iz]!=0){
			valpart=micpart[ix][iy][iz];
			if(valpart>partmax){partmax=valpart;}
			valmic=mic[ix][iy][iz];
			if((valmic==C3S)||(valmic==C2S)||(valmic==C3A)||(valmic==C4AF)){
				nleft[valpart]+=1;
			}
			valmicorig=micorig[ix][iy][iz];
			if((valmicorig==C3S)||(valmicorig==C2S)||(valmicorig==C3A)||(valmicorig==C4AF)){
				norig[valpart]+=1;
			}
		}
	}
	}
	}
 
        /* Output results to end of particle hydration file */
	for(ix=100;ix<=partmax;ix++){
                alpart=0.0;
                if(norig[ix]!=0){
			alpart=1.-(float)nleft[ix]/(float)norig[ix];
		}
		fprintf(phydfile,"%d %d %d %.3f\n",ix,norig[ix],nleft[ix],alpart);
	}
	fclose(phydfile);
}
