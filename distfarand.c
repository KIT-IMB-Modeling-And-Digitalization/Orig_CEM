#include <stdio.h>
#include <math.h>

/*###################################################################*/
/* Program distfarand                                                */
/* Purpose: To distribute fly ash phases randomly amongst particles  */
/* Programmer: Dale P. Bentz                                         */
/*                                                                   */
/*###################################################################*/

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

/* Define phase IDs */
#define FLYASH 30
#define C3A 35
#define CACL2 22
#define ASGID 11
#define CAS2ID 12
#define INERTID 9
#define POZZID 8
#define ANHYDRITE 7
#define SYSSIZE 100

int *seed;
#include "ran1.c"

main()
{
        FILE *infile,*outfile;
        char filein[80],fileout[80];
        int ix,iy,iz,valin,valout;
        float probasg,probcacl2,probsio2;
        float prph,probanh,probcas2,probc3a;
        int nseed;

        printf("Enter random number seed value (<0)\n");
        scanf("%d",&nseed);
        printf("%d\n",nseed);
        seed=(&nseed);

        printf("Enter name of file for input \n");
        scanf("%s",filein);
        printf("%s\n",filein);

        printf("Enter name of file for output \n");
        scanf("%s",fileout);
        printf("%s\n",fileout);

        /* Get user input for phase probabilities (volume fractions) */
        printf("Enter probability for fly ash to be aluminosilicate glass \n");
        scanf("%f",&probasg);
        printf("%f\n",probasg);
        printf("Enter probability for fly ash to be calcium aluminodisilicate \n");
        scanf("%f",&probcas2);
        printf("%f\n",probcas2);
        printf("Enter probability for fly ash to be tricalcium aluminate \n");
        scanf("%f",&probc3a);
        printf("%f\n",probc3a);
        printf("Enter probability for fly ash to be calcium chloride \n");
        scanf("%f",&probcacl2);
        printf("%f\n",probcacl2);
        printf("Enter probability for fly ash to be silica \n");
        scanf("%f",&probsio2);
        printf("%f\n",probsio2);
        printf("Enter probability for fly ash to be anhydrite \n");
        scanf("%f",&probanh);
        printf("%f\n",probanh);

        infile=fopen(filein,"r");
        outfile=fopen(fileout,"w");

        /* Convert to cumulative probabilties */
        /* Order must match that used in for loops below */
	probcacl2+=probasg;
	probsio2+=probcacl2;
	probcas2+=probsio2;
	probanh+=probcas2;
	probc3a+=probanh;
        /* Scan each pixel and convert all fly ash pixels to some */
        /* component phase */
        for(ix=0;ix<SYSSIZE;ix++){
        for(iy=0;iy<SYSSIZE;iy++){
        for(iz=0;iz<SYSSIZE;iz++){

                fscanf(infile,"%d",&valin);
                valout=valin;
                if(valin==FLYASH){
                        valout=INERTID;
                        prph=ran1(seed);
                        if(prph<probasg){
                                valout=ASGID;
                        }
                        else if(prph<(probcacl2)){
                                valout=CACL2;
                        }
                        else if(prph<(probsio2)){
                                valout=POZZID;
                        }
                        else if(prph<(probcas2)){
                                valout=CAS2ID;
                        }
                        else if(prph<(probanh)){
                                valout=ANHYDRITE;
                        }
                        else if(prph<(probc3a)){
                                valout=C3A;
                        }
                }

                fprintf(outfile,"%d\n",valout);
        }
        }
        }

        fclose(infile);
        fclose(outfile);
}

