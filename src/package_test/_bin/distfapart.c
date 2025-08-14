#include <stdio.h>
#include <math.h>

/*###################################################################*/
/* Program distfapart to distribute fly ash phases randomly amongst  */
/* monophase particles                                               */
/* Programmer: Dale P. Bentz                                         */
/* Date: May 1997                                                    */
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

/* Define phase identifiers for fly ash components */
#define FLYASH 30
#define C3A 35
#define CACL2 22
#define ASGID 11
#define SLAGID 10
#define INERTID 9
#define POZZID 8
#define ANHYDRITE 7
#define CAS2ID 12
#define SYSSIZE 100
#define NPARTC 12000

int *seed;
#include "ran1.c"

main()
{
        FILE *infile,*partfile,*outfile;
        char filein[80],fileout[80],filepart[80];
        int ix,iy,iz,valin,valout,partin;
        float probasg,probcacl2,probsio2;
        float probc3a,prph,probcas2,probanh;
        static int phase[NPARTC],partid[NPARTC];
        int nseed,count,c1,phnew;
        long int totcnt,ascnt,cacl2cnt,pozzcnt,inertcnt;
        long int anhcnt,cas2cnt,c3acnt;
        long int markc3a,markas,markcacl2,markpozz,markinert,markanh,markcas2;

        ascnt=cacl2cnt=pozzcnt=inertcnt=0;
        cas2cnt=anhcnt=c3acnt=0;
        printf("Enter random number seed value (<0)\n");
        scanf("%d",&nseed);
        printf("%d\n",nseed);
        seed=(&nseed);

        printf("Enter name of file for input \n");
        scanf("%s",filein);
        printf("%s\n",filein);
        printf("Enter name of file for particle IDs \n");
        scanf("%s",filepart);
        printf("%s\n",filepart);

        printf("Enter name of file for output \n");
        scanf("%s",fileout);
        printf("%s\n",fileout);

        printf("Enter total number of fly ash pixels \n");
        scanf("%ld",&totcnt);
        printf("%ld\n",totcnt);

        /* Get user input for phase probabilities (volume fractions */
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

        /* Determine goal counts for each phase */
        markas=(long)(probasg*(float)totcnt);
        markpozz=(long)(probsio2*(float)totcnt);
        markcacl2=(long)(probcacl2*(float)totcnt);
        markanh=(long)(probanh*(float)totcnt);
        markcas2=(long)(probcas2*(float)totcnt);
        markc3a=(long)(probc3a*(float)totcnt);
        markinert=(long)((1.-probasg-probsio2-probcacl2-probanh-
         probcas2-probc3a)*(float)totcnt);
        /* Convert probabilities to cumulative */
        /* Order must be the same as in for loop below */
        probcacl2+=probasg;
        probsio2+=probcacl2;
        probanh+=probsio2;
        probcas2+=probanh;
        probc3a+=probcas2;

        infile=fopen(filein,"r");
        partfile=fopen(filepart,"r");

        for(ix=0;ix<NPARTC;ix++){
                phase[ix]=partid[ix]=0;
        }

        count=0;
        /* First scan-- find each particle and assign phases */
        for(ix=0;ix<SYSSIZE;ix++){
        for(iy=0;iy<SYSSIZE;iy++){
        for(iz=0;iz<SYSSIZE;iz++){

                fscanf(partfile,"%d",&partin);
                fscanf(infile,"%d",&valin);
                if((valin==FLYASH)&&(partid[partin]==0)){
                        count+=1;
                        partid[partin]=count;

                        valout=INERTID;
                        do{
                        prph=ran1(seed);
                        if((prph<probasg)&&(ascnt<markas)){
                                valout=ASGID;
                        }
                        else if((prph<(probcacl2))&&(cacl2cnt<markcacl2)){
                                valout=CACL2;
                        }
                        else if((prph<(probsio2))&&(pozzcnt<markpozz)){
                                valout=POZZID;
                        }
                        else if((prph<(probanh))&&(anhcnt<markanh)){
                                valout=ANHYDRITE;
                        }
                        else if((prph<(probcas2))&&(cas2cnt<markcas2)){
                                valout=CAS2ID;
                        }
                        else if((prph<(probc3a))&&(c3acnt<markc3a)){
                                valout=C3A;
                        }
                        }while((valout==INERTID)&&(inertcnt>markinert));
                        phase[count]=valout;
                }
                if(valin==FLYASH){
                        c1=partid[partin];
                        phnew=phase[c1];
                        if(phnew==ASGID){
                                ascnt+=1;
                        }
                        else if(phnew==CACL2){
                                cacl2cnt+=1;
                        }
                        else if(phnew==POZZID){
                                pozzcnt+=1;
                        }
                        else if(phnew==ANHYDRITE){
                                anhcnt+=1;
                        }
                        else if(phnew==CAS2ID){
                                cas2cnt+=1;
                        }
                        else if(phnew==C3A){
                                c3acnt+=1;
                        }
                        else if(phnew==INERTID){
                                inertcnt+=1;
                        }
               }
        }
        }
        }

        fclose(infile);
        fclose(partfile);

        /* Now distribute phases in second scan */
        infile=fopen(filein,"r");
        partfile=fopen(filepart,"r");
        outfile=fopen(fileout,"w");
        for(ix=0;ix<SYSSIZE;ix++){
        for(iy=0;iy<SYSSIZE;iy++){
        for(iz=0;iz<SYSSIZE;iz++){
                fscanf(partfile,"%d",&partin);
                fscanf(infile,"%d",&valin);
                valout=valin;
                if(valin==FLYASH){
                        count=partid[partin];
                        valout=phase[count];
                }

                fprintf(outfile,"%d\n",valout);

        }
        }
        }
        fclose(infile);
        fclose(partfile);
        fclose(outfile);
}
