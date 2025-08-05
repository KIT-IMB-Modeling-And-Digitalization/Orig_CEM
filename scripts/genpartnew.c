/************************************************************************/
/*                                                                      */
/*      Program genpartnew.c to generate three-dimensional cement       */
/*              and gypsum particles in a 3-D box with periodic         */
/*              boundaries.                                             */
/*      Particles are composed of either cement clinker, gypsum,        */
/*              slag, calcium carbonate, or inert filler,               */
/*              follow a user-specified size distribution, and can      */
/*              be either flocculated, random, or dispersed.            */
/*      Programmer: Dale P. Bentz                                       */
/*                  Building and Fire Research Laboratory               */
/*                  NIST                                                */
/*                  100 Bureau Drive Mail Stop 8615                     */
/*                  Gaithersburg, MD  20899-8615   USA                  */
/*                  (301) 975-5865      FAX: 301-990-6891               */
/*                  E-mail: dale.bentz@nist.gov                         */
/*                                                                      */
/************************************************************************/

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

/* Modified 3/97 to allow placement of pozzolanic, inert and fly ash particles */
/* Modified 9/98 to allow placement of various forms of gypsum */
/* Documented version produced 1/00 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "common.h"


#define SYSSIZE 100     /* system size in pixels per dimension */
#define MAXTRIES 150000 /* maximum number of random tries for sphere placement */

/* phase identifiers */
#define POROSITY 0
/* Note that each particle must have a separate ID to allow for flocculation */
#define CEM 100	        /* and greater */
#define CEMID 1	        /* phase identifier for cement */
#define C2SID 2		/* phase identified for C2S cement */
#define GYPID 5	        /* phase identifier for gypsum */	
#define HEMIHYDRATE 6  /* phase identifier for hemihydrate */
#define ANHYDRITE 7    /* phase identifier for anhydrite */
#define POZZID 8       /* phase identifier for pozzolanic material */
#define INERTID 9      /* phase identifier for inert material */
#define SLAGID 10	/* phase identifier for slag */
#define CACO3 26 	/* phase identifier for calcium carbonate */
#define AGG 28          /* phase identifier for flat aggregate */
#define FLYASH 30       /* phase identifier for all fly ash components */
#define NPARTC 30000    /* maximum number of particles allowed in box*/
#define BURNT 34000     /* this value must be at least 100 > NPARTC */
#define NUMSIZES 200     /* maximum number of different particle sizes */

/* data structure for clusters to be used in flocculation */
struct cluster{
        int partid;	/* index for particle */
        int clustid;  	/* ID for cluster to which this particle belongs */
        int partphase;  /* phase identifier for this particle */
        int x,y,z,r;	/* particle centroid and radius in pixels */
        struct cluster *nextpart;  /* pointer to next particle in cluster */
};

/* 3-D particle structure (each particle has own ID) stored in array cement */
/* 3-D microstructure is stored in 3-D array cemreal */
static unsigned short int cement [SYSSIZE+1] [SYSSIZE+1] [SYSSIZE+1];
static unsigned short int cemreal [SYSSIZE+1] [SYSSIZE+1] [SYSSIZE+1];
int npart,aggsize;     /* global number of particles and size of aggregate */
int *seed;      /* random number seed- global */
int dispdist;   /* dispersion distance in pixels */
int clusleft;   /* number of clusters in system */
/* parameters to aid in obtaining correct sulfate content */
long int n_sulfate=0,target_sulfate=0,n_total=0,target_total=0,volpart[47];
long int n_anhydrite=0,target_anhydrite=0,n_hemi=0,target_hemi=0;
float probgyp,probhem,probanh; /* probability of gypsum particle instead of cement */
			/* and probabilities of anhydrite and hemihydrate */
			/* relative to total sulfate */
struct cluster *clust[NPARTC];	/* limit of NPARTC particles/clusters */

/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */
#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
#define NTAB 32
#define EPS (1.2E-07)
#define MAX(a,b) (a>b)?a:b
#define MIN(a,b) (a<b)?a:b

double ran1(idum)
int *idum;
/* Calls: no routines */
/* Called by: gsphere,makefloc */
{
        int j,k;
        static int iv[NTAB],iy=0;
        void nrerror();
        static double NDIV = 1.0/(1.0+(IM-1.0)/NTAB);
        static double RNMX = (1.0-EPS);
        static double AM = (1.0/IM);

        if ((*idum <= 0) || (iy == 0)) {
                *idum = MAX(-*idum,*idum);
                for(j=NTAB+7;j>=0;j--) {
                        k = *idum/IQ;
                        *idum = IA*(*idum-k*IQ)-IR*k;
                        if(*idum < 0) *idum += IM;
                        if(j < NTAB) iv[j] = *idum;
                }
                iy = iv[0];
        }
        k = *idum/IQ;
        *idum = IA*(*idum-k*IQ)-IR*k;
        if(*idum<0) *idum += IM;
        j = iy*NDIV;
        iy = iv[j];
        iv[j] = *idum;
        return MIN(AM*iy,RNMX);
}
#undef IA 
#undef IM 
#undef IQ
#undef IR
#undef NTAB
#undef EPS 
#undef MAX
#undef MIN

/* routine to add a flat plate aggregate in the microstructure */
void addagg()
/* Calls: no other routines */
/* Called by: main program */
{
        int ix,iy,iz;
        int agglo,agghi;

/* Be sure aggregate size is an even integer */
        do{
             printf("Enter thickness of aggregate to place (an even integer) \n");
             scanf("%d",&aggsize);
             printf("%d\n",aggsize);
        } while (((aggsize%2)!=0)||(aggsize>(SYSSIZE-2)));

        if(aggsize!=0){
                agglo=(SYSSIZE/2)-((aggsize-2)/2);
                agghi=(SYSSIZE/2)+(aggsize/2);

                /* Aggregate is placed in yz plane */
                for(ix=agglo;ix<=agghi;ix++){
                for(iy=1;iy<=SYSSIZE;iy++){
                for(iz=1;iz<=SYSSIZE;iz++){

        /* Mark aggregate in both particle and microstructure images */
                        cement [ix][iy][iz]=AGG;
                        cemreal [ix][iy][iz]=AGG;
                }
                }
                }
        }
}

/* routine to check or perform placement of sphere of ID phasein */
/* centered at location (xin,yin,zin) of radius radd */
/* wflg=1 check for fit of sphere */
/* wflg=2 place the sphere */
/* phasein and phase2 are phases to assign to cement and cemreal images resp. */
int chksph(xin,yin,zin,radd,wflg,phasein,phase2)
        int xin,yin,zin,radd,wflg,phasein,phase2;
/* Calls: no other routines */
/* Called by: gsphere */
{
        int nofits,xp,yp,zp,i,j,k;
        float dist,xdist,ydist,zdist,ftmp;

        nofits=0;       /* Flag indicating if placement is possible */
/* Check all pixels within the digitized sphere volume */
        for(i=xin-radd;((i<=xin+radd)&&(nofits==0));i++){
                xp=i;
                /* use periodic boundary conditions for sphere placement */
                if(xp<1) {xp+=SYSSIZE;}
                else if(xp>SYSSIZE) {xp-=SYSSIZE;}
                ftmp=(float)(i-xin);
                xdist=ftmp*ftmp;
        for(j=yin-radd;((j<=yin+radd)&&(nofits==0));j++){
                yp=j;
                /* use periodic boundary conditions for sphere placement */
                if(yp<1) {yp+=SYSSIZE;}
                else if(yp>SYSSIZE) {yp-=SYSSIZE;}
                ftmp=(float)(j-yin);
                ydist=ftmp*ftmp;
        for(k=zin-radd;((k<=zin+radd)&&(nofits==0));k++){
                zp=k;
                /* use periodic boundary conditions for sphere placement */
                if(zp<1) {zp+=SYSSIZE;}
                else if(zp>SYSSIZE) {zp-=SYSSIZE;}
                ftmp=(float)(k-zin);
                zdist=ftmp*ftmp;

        /* Compute distance from center of sphere to this pixel */
                dist=sqrt(xdist+ydist+zdist);
                if((dist-0.5)<=(float)radd){
                        /* Perform placement */
                        if(wflg==2){
                                cement [xp] [yp] [zp]=phasein;
                                cemreal [xp][yp][zp]=phase2;
                        }
                        /* or check placement */
                        else if((wflg==1)&&(cement [xp] [yp] [zp] !=POROSITY)){
                                nofits=1;
                        }
                }
                /* Check for overlap with aggregate */
      if((wflg==1)&&((fabs(xp-((float)(SYSSIZE+1)/2.0)))<((float)aggsize/2.0))){
                        nofits=1;
                }
        }
        }
        }
	
        /* return flag indicating if sphere will fit */
        return(nofits);
}
					
/* routine to place spheres of various sizes and phases at random */
/* locations in 3-D microstructure */
/* numgen is number of different size spheres to place */
/* numeach holds the number of each size class */
/* sizeeach holds the radius of each size class */
/* pheach holds the phase of each size class */
void gsphere(numgen,numeach,sizeeach,pheach)
        int numgen;
        long int numeach[NUMSIZES];
        int sizeeach[NUMSIZES],pheach[NUMSIZES];
/* Calls: makesph, ran1 */
/* Called by: create */
{
        int count,x,y,z,radius,ig,tries,phnow;
        long int jg;
        float rx,ry,rz,testgyp,typegyp;
        struct cluster *partnew;

/* Generate spheres of each size class in turn (largest first) */
        for(ig=0;ig<numgen;ig++){

		phnow=pheach[ig];        /* phase for this class */
                radius=sizeeach[ig];     /* radius for this class */
                /* loop for each sphere in this size class */
                for(jg=1;jg<=numeach[ig];jg++){

                        tries=0;
                        /* Stop after MAXTRIES random tries */
                        do{
                                tries+=1;
                        /* generate a random center location for the sphere */
                                x=(int)((float)SYSSIZE*ran1(seed))+1;
                                y=(int)((float)SYSSIZE*ran1(seed))+1;
                                z=(int)((float)SYSSIZE*ran1(seed))+1;
                                /* See if the sphere will fit at x,y,z */
                                /* Include dispersion distance when checking */
                            /* to insure requested separation between spheres */
                              count=chksph(x,y,z,radius+dispdist,1,npart+CEM,0);
                                if((tries>MAXTRIES)&&(dispdist==2)){
					tries=0;
					dispdist+=1;
				}
                                if(tries>MAXTRIES){
printf("Could not place sphere %d after %d random attempts \n",npart,MAXTRIES);
                                        printf("Exitting program \n");
                                        exit(1);
                                }
                        } while(count!=0);

                        /* place the sphere at x,y,z */
                        npart+=1;
                        if(npart>=NPARTC){
                                printf("Too many spheres being generated \n");
           printf("User needs to increase value of NPARTC at top of C-code\n");
           printf("Exitting program \n");
                                exit(1);
                        }
                        /* Allocate space for new particle info */
                  clust[npart]=(struct cluster *)malloc(sizeof(struct cluster));
                        clust[npart]->partid=npart;
                        clust[npart]->clustid=npart;
                        /* Default to cement placement */
                        clust[npart]->partphase=CEMID;
                        clust[npart]->x=x;
                        clust[npart]->y=y;
                        clust[npart]->z=z;
                        clust[npart]->r=radius;
                        clusleft+=1;
			if(phnow==1){
                        testgyp=ran1(seed);
                  /* Do not use dispersion distance when placing particle */
                        if(((testgyp>probgyp)&&((target_sulfate-n_sulfate)<(target_total-n_total)))||(n_sulfate>target_sulfate)||(volpart[radius]>(target_sulfate-n_sulfate))||(numeach[ig]<=2)){
                                count=chksph(x,y,z,radius,2,npart+CEM-1,CEMID);
                                n_total+=volpart[radius];
                        }
                        else{
                                /* Place particle as gypsum */
				typegyp=ran1(seed);
                                n_total+=volpart[radius];
				n_sulfate+=volpart[radius];
				if((probanh>=1.0)||((typegyp<probanh)&&(n_anhydrite<target_anhydrite)&&(volpart[radius]<=(target_anhydrite-n_anhydrite)))){
					/* Place particle as anhydrite */
                                        n_anhydrite+=volpart[radius];
					count=chksph(x,y,z,radius,2,npart+CEM-1,ANHYDRITE);
					clust[npart]->partphase=ANHYDRITE;
				}
				else if(((probanh+probhem)>=1.0)||((typegyp<(probanh+probhem))&&(n_hemi<target_hemi)&&(volpart[radius]<=(target_hemi-n_hemi)))){
					/* Place particle as hemihydrate */
                                        n_hemi+=volpart[radius];
					count=chksph(x,y,z,radius,2,npart+CEM-1,HEMIHYDRATE);
					clust[npart]->partphase=HEMIHYDRATE;
				}
				else{
	                                count=chksph(x,y,z,radius,2,npart+CEM-1,GYPID);
       		                         /* Correct phase ID of particle */
       		                         clust[npart]->partphase=GYPID;
				}
                        }
			}
			/* place as inert, CaCO3, C2S, slag, or pozzolanic material */
			else{
				count=chksph(x,y,z,radius,2,npart+CEM-1,phnow);
                                /* Correct phase ID of particle */
                                clust[npart]->partphase=phnow;
			}
                        clust[npart]->nextpart=NULL;
                }
        }
}

/* routine to obtain user input and create a starting microstructure */
void create()
/* Calls: gsphere */
/* Called by: main program */
{
        int numsize,sphrad [NUMSIZES],sphase[NUMSIZES];
        long int sphnum [NUMSIZES],inval1;
        int isph,inval;

        do{
printf("Enter number of different size spheres to use(max. is %d) \n",NUMSIZES);
                scanf("%d",&numsize);
                printf("%d \n",numsize);
        }while ((numsize>NUMSIZES)||(numsize<0));
        do{
                printf("Enter dispersion factor (separation distance in pixels) for spheres (0-2) \n");
                printf("0 corresponds to totally random placement \n");
                scanf("%d",&dispdist);
                printf("%d \n",dispdist);
        }while ((dispdist<0)||(dispdist>2));
        do{
                printf("Enter probability for gypsum particles on a random particle basis (0.0-1.0) \n");
                scanf("%f",&probgyp);
                printf("%f \n",probgyp);
        } while ((probgyp<0.0)||(probgyp>1.0));
	do{
		printf("Enter probabilities for hemihydrate and anhydrite forms of gypsum (0.0-1.0) \n");
		scanf("%f %f",&probhem,&probanh);
		printf("%f %f\n",probhem,probanh);
	} while ((probhem<0.0)||(probhem>1.0)||(probanh<0.0)||(probanh>1.0)||((probanh+probhem)>1.001));

        if((numsize>0)&&(numsize<(NUMSIZES+1))){
       printf("Enter number, radius, and phase ID for each sphere class (largest radius 1st) \n");
	printf("Phases are %d- Cement and (random) calcium sulfate, %d- C2S, %d- Gypsum, %d- hemihydrate %d- anhydrite %d- Pozzolanic, %d- Inert, %d- Slag, %d- CaCO3 %d- Fly Ash \n",CEMID,C2SID,GYPID,HEMIHYDRATE,ANHYDRITE,POZZID,INERTID,SLAGID,CACO3,FLYASH);

        /* Obtain input for each size class of spheres */
                for(isph=0;isph<numsize;isph++){
                        printf("Enter number of spheres of class %d \n",isph+1);
                        scanf("%ld",&inval1);
                        printf("%ld \n",inval1);
                        sphnum[isph]=inval1;
                        do{
                        printf("Enter radius of spheres of class %d \n",isph+1);
                                printf("(Integer <=%d please) \n",SYSSIZE/3);
                                scanf("%d",&inval);
                                printf("%d \n",inval);
                        } while ((inval<1)||(inval>(SYSSIZE/3)));
                        sphrad[isph]=inval;
                        do{
                        printf("Enter phase of spheres of class %d \n",isph+1);
                                scanf("%d",&inval);
                                printf("%d \n",inval);
                        } while ((inval!=CEMID)&&(inval!=C2SID)&&(inval!=GYPID)&&(inval!=HEMIHYDRATE)&&(inval!=ANHYDRITE)&&(inval!=POZZID)&&(inval!=INERTID)&&(inval!=SLAGID)&&(inval!=FLYASH)&&(inval!=CACO3));
                        sphase[isph]=inval;
			if(inval==CEMID){
				target_total+=sphnum[isph]*volpart[sphrad[isph]];
			}
	
                }
		/* Determine target pixel counts for calcium sulfate forms */
                target_sulfate=(int)((float)target_total*probgyp);
		target_anhydrite=(int)((float)target_total*probgyp*probanh);
		target_hemi=(int)((float)target_total*probgyp*probhem);
                gsphere(numsize,sphnum,sphrad,sphase);
        }
}

/* Routine to draw a particle during flocculation routine */
/* See routine chksph for definition of parameters */
void drawfloc(xin,yin,zin,radd,phasein,phase2)
        int xin,yin,zin,radd,phasein,phase2;
/* Calls: no other routines */
/* Called by: makefloc */
{
        int xp,yp,zp,i,j,k;
        float dist,xdist,ydist,zdist,ftmp;

/* Check all pixels within the digitized sphere volume */
        for(i=xin-radd;(i<=xin+radd);i++){
                xp=i;
                /* use periodic boundary conditions for sphere placement */
                if(xp<1) {xp+=SYSSIZE;}
                else if(xp>SYSSIZE) {xp-=SYSSIZE;}
                ftmp=(float)(i-xin);
                xdist=ftmp*ftmp;
        for(j=yin-radd;(j<=yin+radd);j++){
                yp=j;
                /* use periodic boundary conditions for sphere placement */
                if(yp<1) {yp+=SYSSIZE;}
                else if(yp>SYSSIZE) {yp-=SYSSIZE;}
                ftmp=(float)(j-yin);
                ydist=ftmp*ftmp;
        for(k=zin-radd;(k<=zin+radd);k++){
                zp=k;
                /* use periodic boundary conditions for sphere placement */
                if(zp<1) {zp+=SYSSIZE;}
                else if(zp>SYSSIZE) {zp-=SYSSIZE;}
                ftmp=(float)(k-zin);
                zdist=ftmp*ftmp;

        /* Compute distance from center of sphere to this pixel */
                dist=sqrt(xdist+ydist+zdist);
                if((dist-0.5)<=(float)radd){
                        /* Update both cement and cemreal images */
                        cement [xp] [yp] [zp]=phasein;
                        cemreal [xp] [yp] [zp]=phase2;
                }
        }
        }
        }
}

/* Routine to check particle placement during flocculation */
/* for particle of size radd centered at (xin,yin,zin) */
/* Returns flag indicating if placement is possible */
int chkfloc(xin,yin,zin,radd)
        int xin,yin,zin,radd;
/* Calls: no other routines */
/* Called by: makefloc */
{
        int nofits,xp,yp,zp,i,j,k;
        float dist,xdist,ydist,zdist,ftmp;

        nofits=0;       /* Flag indicating if placement is possible */

/* Check all pixels within the digitized sphere volume */
        for(i=xin-radd;((i<=xin+radd)&&(nofits==0));i++){
                xp=i;
                /* use periodic boundary conditions for sphere placement */
                if(xp<1) {xp+=SYSSIZE;}
                else if(xp>SYSSIZE) {xp-=SYSSIZE;}
                ftmp=(float)(i-xin);
                xdist=ftmp*ftmp;
        for(j=yin-radd;((j<=yin+radd)&&(nofits==0));j++){
                yp=j;
                /* use periodic boundary conditions for sphere placement */
                if(yp<1) {yp+=SYSSIZE;}
                else if(yp>SYSSIZE) {yp-=SYSSIZE;}
                ftmp=(float)(j-yin);
                ydist=ftmp*ftmp;
        for(k=zin-radd;((k<=zin+radd)&&(nofits==0));k++){
                zp=k;
                /* use periodic boundary conditions for sphere placement */
                if(zp<1) {zp+=SYSSIZE;}
                else if(zp>SYSSIZE) {zp-=SYSSIZE;}
                ftmp=(float)(k-zin);
                zdist=ftmp*ftmp;

        /* Compute distance from center of sphere to this pixel */
                dist=sqrt(xdist+ydist+zdist);
                if((dist-0.5)<=(float)radd){
                        if((cement [xp] [yp] [zp] !=POROSITY)){
                                /* Record ID of particle hit */
                                nofits=cement [xp] [yp] [zp];
                        }
                }
                /* Check for overlap with aggregate */
                if((fabs(xp-((float)(SYSSIZE+1)/2.0)))<((float)aggsize/2.0)){
                        nofits=AGG;
                }
        }
        }
        }
        /* return flag indicating if sphere will fit */
        return(nofits);
}

/* routine to perform flocculation of particles */
void makefloc()
/* Calls: drawfloc, chkfloc, ran1 */
/* Called by: main program */
{
        int partdo,numfloc;
        int nstart;
        int nleft,ckall;
        int xm,ym,zm,moveran;
        int xp,yp,zp,rp,clushit,valkeep;
        int iclus,cluspart[NPARTC];
        struct cluster *parttmp,*partpoint,*partkeep;
	
        nstart=npart; /* Counter of number of flocs remaining */
        for(iclus=1;iclus<=npart;iclus++){
                  cluspart[iclus]=iclus;
        }
        do{
          printf("Enter number of flocs desired at end of routine (>0) \n");
                scanf("%d",&numfloc);
                printf("%d\n",numfloc);
        } while (numfloc<=0);

        while(nstart>numfloc){
                nleft=0;

                /* Try to move each cluster in turn */
                for(iclus=1;iclus<=npart;iclus++){
                        if(clust[iclus]==NULL){
                                nleft+=1;
                        }
                        else{
                                xm=ym=zm=0;
                /* Generate a random move in one of 6 principal directions */
                moveran=6.*ran1(seed);
                switch(moveran){
                                case 0:
                                        xm=1;
                                        break;
                                case 1:
                                        xm=(-1);
                                        break;
                                case 2:
                                        ym=1;
                                        break;
                                case 3:
                                        ym=(-1);
                                        break;
                                case 4:
                                        zm=1;
                                        break;
                                case 5:
                                        zm=(-1);
                                        break;
                                default:
                                        break;
                }

                /* First erase all particles in cluster */
                partpoint=clust[iclus];
                while(partpoint!=NULL){
                        xp=partpoint->x;
                        yp=partpoint->y;
                        zp=partpoint->z;
                        rp=partpoint->r;
                        drawfloc(xp,yp,zp,rp,0,0);
                        partpoint=partpoint->nextpart;
                }

                ckall=0;
                /* Now try to draw cluster at new location */
                partpoint=clust[iclus];
                while((partpoint!=NULL)&&(ckall==0)){
                        xp=partpoint->x+xm;
                        yp=partpoint->y+ym;
                        zp=partpoint->z+zm;
                        rp=partpoint->r;
                        ckall=chkfloc(xp,yp,zp,rp);
                        partpoint=partpoint->nextpart;
                }

                if(ckall==0){
                /* Place cluster particles at new location */
                        partpoint=clust[iclus];
                        while(partpoint!=NULL){
                                xp=partpoint->x+xm;
                                yp=partpoint->y+ym;
                                zp=partpoint->z+zm;
                                rp=partpoint->r;
                                valkeep=partpoint->partphase;
                                partdo=partpoint->partid;
                                drawfloc(xp,yp,zp,rp,partdo+CEM-1,valkeep);
                                /* Update particle location */
                                partpoint->x=xp;
                                partpoint->y=yp;
                                partpoint->z=zp;
                                partpoint=partpoint->nextpart;
                        }
                }
                else{
                        /* A cluster or aggregate was hit */
                        /* Draw particles at old location */
                        partpoint=clust[iclus];
                        /* partkeep stores pointer to last particle in list */
                        while(partpoint!=NULL){
                                xp=partpoint->x;
                                yp=partpoint->y;
                                zp=partpoint->z;
                                rp=partpoint->r;
                                valkeep=partpoint->partphase;
                                partdo=partpoint->partid;
                                drawfloc(xp,yp,zp,rp,partdo+CEM-1,valkeep);
                                partkeep=partpoint;
                                partpoint=partpoint->nextpart;
                        }
                        /* Determine the cluster hit */
                        if(ckall!=AGG){
                                clushit=cluspart[ckall-CEM+1];
       /* Move all of the particles from cluster clushit to cluster iclus */
                                parttmp=clust[clushit];
                                /* Attach new cluster to old one */
                                partkeep->nextpart=parttmp;
                                while(parttmp!=NULL){
                                        cluspart[parttmp->partid]=iclus;
                          /* Relabel all particles added to this cluster */
                                        parttmp->clustid=iclus;
                                        parttmp=parttmp->nextpart;
                                }
                                /* Disengage the cluster that was hit */
                                clust[clushit]=NULL;
                                nstart-=1;
                        }
                }
                }
        }
  printf("Number left was %d but number of clusters is %d \n",nleft,nstart);
        }      /* end of while loop */
        clusleft=nleft;
}

/* routine to assess global phase fractions present in 3-D system */
void measure()
/* Calls: no other routines */
/* Called by: main program */
{
        long int npor,nc2s,ngyp,ncem,nagg,npozz,ninert,nflyash,nanh,nhem,ncaco3,nslag;
        int i,j,k,valph;

/* counters for the various phase fractions */
        npor=0;
        ngyp=0;
        ncem=0;
        nagg=0;
	ninert=0;
	nslag=0;
	nc2s=0;
	npozz=0;
	nflyash=0;
	nanh=0;
	nhem=0;
	ncaco3=0;

/* Check all pixels in 3-D microstructure */
        for(i=1;i<=SYSSIZE;i++){
        for(j=1;j<=SYSSIZE;j++){
        for(k=1;k<=SYSSIZE;k++){
                valph=cemreal [i] [j] [k];	
                if(valph==POROSITY) {npor+=1;}
                else if(valph==CEMID){ncem+=1;}
                else if(valph==C2SID){nc2s+=1;}
                else if(valph==GYPID){ngyp+=1;}
                else if(valph==ANHYDRITE){nanh+=1;}
                else if(valph==HEMIHYDRATE){nhem+=1;}
                else if(valph==AGG) {nagg+=1;}
                else if(valph==POZZID) {npozz+=1;}
                else if(valph==SLAGID) {nslag+=1;}
                else if(valph==INERTID) {ninert+=1;}
                else if(valph==FLYASH) {nflyash+=1;}
                else if(valph==CACO3) {ncaco3+=1;}
        }
        }
        }

/* Output results */
        printf("\n Phase counts are: \n");
        printf("Porosity= %ld \n",npor);
        printf("Cement= %ld \n",ncem);
        printf("C2S= %ld \n",nc2s);
        printf("Gypsum= %ld \n",ngyp);
	printf("Anhydrite= %ld \n",nanh);
	printf("Hemihydrate= %ld \n",nhem);
	printf("Pozzolan= %ld \n",npozz);
	printf("Inert= %ld \n",ninert);
	printf("Slag= %ld \n",nslag);
	printf("CaCO3= %ld \n",ncaco3);
	printf("Fly Ash= %ld \n",nflyash);
        printf("Aggregate= %ld \n",nagg);
}
	
/* Routine to measure phase fractions as a function of distance from */
/* aggregate surface						     */
void measagg()
/* Calls: no other routines */
/* Called by: main program */
{
        int phase [40],ptot;
        int icnt,ixlo,ixhi,iy,iz,phid,idist;
	FILE *aggfile;

	/* By default, results are sent to output file called agglist.out */
	aggfile=fopen("agglist.out","w");
        printf("Distance  Porosity  Cement C2S  Gypsum  Anhydrite Hemihydrate Pozzolan  Inert   Slag CaCO3  Fly Ash\n");
        fprintf(aggfile,"Distance  Porosity  Cement  C2S Gypsum  Anhydrite Hemihydrate Pozzolan  Inert  Slag  CaCO3   Fly Ash\n");

/* Increase distance from aggregate in increments of one */
        for(idist=1;idist<=(SYSSIZE-aggsize)/2;idist++){
                /* Pixel left of aggregate surface */
                ixlo=((SYSSIZE-aggsize+2)/2)-idist;
                /* Pixel right of aggregate surface */
                ixhi=((SYSSIZE+aggsize)/2)+idist;
	
        /* Initialize phase counts for this distance */
                for(icnt=0;icnt<39;icnt++){
                        phase[icnt]=0;
                }
                ptot=0;

/* Check all pixels which are this distance from aggregate surface */
                for(iy=1;iy<=SYSSIZE;iy++){
                for(iz=1;iz<=SYSSIZE;iz++){
                        phid=cemreal [ixlo] [iy] [iz];
                        ptot+=1;
                        if(phid<=FLYASH){
                                phase[phid]+=1;
                        }
                        phid=cemreal [ixhi] [iy] [iz];
                        ptot+=1;
                        if(phid<=FLYASH){
                                phase[phid]+=1;
                        }
                }
                }

                /* Output results for this distance from surface */
     printf("%d   %d   %d  %d  %d %d %d %d  %d %d %d %d\n",idist,phase[0],phase[CEMID],phase[C2SID],phase[GYPID],phase[ANHYDRITE],phase[HEMIHYDRATE],phase[POZZID],phase[INERTID],phase[SLAGID],phase[CACO3],phase[FLYASH]);
     fprintf(aggfile,"%d   %d   %d  %d  %d %d %d %d  %d %d %d %d\n",idist,phase[0],phase[CEMID],phase[C2SID],phase[GYPID],phase[ANHYDRITE],phase[HEMIHYDRATE],phase[POZZID],phase[INERTID],phase[SLAGID],phase[CACO3],phase[FLYASH]);
	

         }
	fclose(aggfile);
}

/* routine to assess the connectivity (percolation) of a single phase */
/* Two matrices are used here: one for the current burnt locations */
/*                the other to store the newly found burnt locations */
void connect()
/* Calls: no other routines */
/* Called by: main program */
{
        long int inew,ntop,nthrough,ncur,nnew,ntot;
        int i,j,k,nmatx[29000],nmaty[29000],nmatz[29000];
     int xcn,ycn,zcn,npix,x1,y1,z1,igood,nnewx[29000],nnewy[29000],nnewz[29000];
        int jnew,icur;

        do{
                printf("Enter phase to analyze 0) pores 1) Cement  \n");
                scanf("%d",&npix);
                printf("%d \n",npix);
        } while ((npix!=0)&&(npix!=1));

        /* counters for number of pixels of phase accessible from top surface */
        /* and number which are part of a percolated pathway */
        ntop=0;
        nthrough=0;

        /* percolation is assessed from top to bottom only */
        /* and burning algorithm is periodic in x and y directions */
        k=1;
        for(i=1;i<=SYSSIZE;i++){
                for(j=1;j<=SYSSIZE;j++){
                        ncur=0;
                        ntot=0;
                        igood=0;     /* Indicates if bottom has been reached */
if(((cement [i] [j] [k]==npix)&&((cement [i] [j] [SYSSIZE]==npix)||
(cement [i] [j] [SYSSIZE]==(npix+BURNT))))||((cement [i] [j] [SYSSIZE]>=CEM)&&
(cement [i] [j] [k]>=CEM)&&(cement [i] [j] [k]<BURNT)&&(npix==1))){
                                /* Start a burn front */
                                cement [i] [j] [k]+=BURNT;
                                ntot+=1;
                                ncur+=1;
                                /* burn front is stored in matrices nmat* */
                                /* and nnew* */
                                nmatx[ncur]=i;
                                nmaty[ncur]=j;
                                nmatz[ncur]=1;
                            /* Burn as long as new (fuel) pixels are found */
                                do{
                                         nnew=0;
                                         for(inew=1;inew<=ncur;inew++){
                                                xcn=nmatx[inew];
                                                ycn=nmaty[inew];
                                                zcn=nmatz[inew];

                                                /* Check all six neighbors */
                                                for(jnew=1;jnew<=6;jnew++){
                                                        x1=xcn;
                                                        y1=ycn;
                                                        z1=zcn;
                                                        if(jnew==1){
                                                                x1-=1;
                                                                if(x1<1){
                                                                    x1+=SYSSIZE;
                                                                }
                                                        }
                                                        else if(jnew==2){
                                                                x1+=1;
                                                                if(x1>SYSSIZE){
                                                                    x1-=SYSSIZE;
                                                                }
                                                        }
                                                        else if(jnew==3){
                                                                y1-=1;
                                                                if(y1<1){
                                                                    y1+=SYSSIZE;
                                                                }
                                                        }
                                                        else if(jnew==4){
                                                                y1+=1;
                                                                if(y1>SYSSIZE){
                                                                    y1-=SYSSIZE;
                                                                }
                                                        }
                                                        else if(jnew==5){
                                                                z1-=1;
                                                        }
                                                        else if(jnew==6){
                                                                z1+=1;
                                                        }

/* Nonperiodic in z direction so be sure to remain in the 3-D box */
                                                  if((z1>=1)&&(z1<=SYSSIZE)){
if((cement [x1] [y1] [z1]==npix)||((cement [x1] [y1] [z1]>=CEM)&&
(cement [x1] [y1] [z1]<BURNT)&&(npix==1))){
                                                                ntot+=1;
                                                  cement [x1] [y1] [z1]+=BURNT;
                                                                nnew+=1;
                                                          if(nnew>=29000){
                                            printf("error in size of nnew \n");
                                                          }
                                                                nnewx[nnew]=x1;
                                                                nnewy[nnew]=y1;
                                                                nnewz[nnew]=z1;
                                /* See if bottom of system has been reached */
                                                          if(z1==SYSSIZE){
                                                                igood=1;
                                                          }
                                                          }
                                                }
                                                }
                                        }
                                        if(nnew>0){
                                                ncur=nnew;
                                        /* update the burn front matrices */
                                                for(icur=1;icur<=ncur;icur++){
                                                        nmatx[icur]=nnewx[icur];
                                                        nmaty[icur]=nnewy[icur];
                                                        nmatz[icur]=nnewz[icur];
                                                }
                                        }
                                }while (nnew>0);

                                ntop+=ntot;
                                if(igood==1){
                                         nthrough+=ntot;
                                }
                        }
                }
        }

        printf("Phase ID= %d \n",npix);
        printf("Number accessible from top= %ld \n",ntop);
        printf("Number contained in through pathways= %ld \n",nthrough);

        /* return the burnt sites to their original phase values */
        for(i=1;i<=SYSSIZE;i++){
        for(j=1;j<=SYSSIZE;j++){
        for(k=1;k<=SYSSIZE;k++){
                if(cement [i] [j] [k]>=BURNT){
                        cement [i] [j] [k]-=BURNT;
                }
        }
        }
        }
}

/* Routine to output final microstructure to file */
void outmic()
/* Calls: no other routines */
/* Called by: main program */
{
        FILE *outfile,*partfile;
        char filen[80],filepart[80];
        int ix,iy,iz,valout;

        printf("Enter name of file to save microstructure to \n");
        scanf("%s",filen);
        printf("%s\n",filen);

        outfile=fopen(filen,"w");

        printf("Enter name of file to save particle IDs to \n");
        scanf("%s",filepart);
        printf("%s\n",filepart);

        partfile=fopen(filepart,"w");

        for(iz=1;iz<=SYSSIZE;iz++){
        for(iy=1;iy<=SYSSIZE;iy++){
        for(ix=1;ix<=SYSSIZE;ix++){
                valout=cemreal[ix][iy][iz];
                fprintf(outfile,"%1d\n",valout);
                valout=cement[ix][iy][iz];
                if(valout<0){valout=0;}
                fprintf(partfile,"%d\n",valout);
        }
        }
        }
        fclose(outfile);
        fclose(partfile);
}

int main(){
        int userc;      /* User choice from menu */
        int nseed,ig,jg,kg;

	/* Initialize volume array */
        volpart[0]=1;
	volpart[1]=19;
	volpart[2]=81;
	volpart[3]=179;
	volpart[4]=389;
	volpart[5]=739;
	volpart[6]=1189;
	volpart[7]=1791;
	volpart[8]=2553;
	volpart[9]=3695;
	volpart[10]=4945;
	volpart[11]=6403;
	volpart[12]=8217;
	volpart[13]=10395;
	volpart[14]=12893;
	volpart[15]=15515;
	volpart[16]=18853;
	volpart[17]=22575;
	volpart[18]=26745;
	volpart[19]=31103;
	volpart[20]=36137;
	volpart[21]=41851;
	volpart[22]=47833;
	volpart[23]=54435;
	volpart[24]=61565;
	volpart[25]=69599;
	volpart[26]=78205;
	volpart[27]=87271;
	volpart[28]=97233;
	volpart[29]=107783;
        volpart[30]=119009;
	volpart[31]=131155;
	volpart[32]=143761;
	volpart[33]=157563;
	volpart[34]=172317;
	volpart[35]=187511;
	volpart[36]=203965;

        printf("Enter random number seed value (a negative integer) \n");
        scanf("%d",&nseed);
        printf("%d \n",nseed);
        seed=(&nseed);

        /* Initialize counters and system parameters */
        npart=0;
        aggsize=0;
        clusleft=0;

        /* clear the 3-D system to all porosity to start */
        for(ig=1;ig<=SYSSIZE;ig++){
        for(jg=1;jg<=SYSSIZE;jg++){
        for(kg=1;kg<=SYSSIZE;kg++){
                cement [ig] [jg] [kg]=POROSITY;
                cemreal [ig] [jg] [kg]=POROSITY;
        }
        }
        }

        /* present menu and execute user choice */
        do{
printf(" \n Input User Choice \n");
printf("1) Exit \n");
printf("2) Add spherical particles (cement,gypsum, pozzolans, etc.) to microstructure \n");
printf("3) Flocculate system by reducing number of particle clusters \n");
printf("4) Measure global phase fractions \n");
printf("5) Add an aggregate to the microstructure \n");
printf("6) Measure single phase connectivity (pores or solids) \n");
printf("7) Measure phase fractions vs. distance from aggregate surface \n");
printf("8) Output current microstructure to file \n");

                scanf("%d",&userc);
                printf("%d \n",userc);
                fflush(stdout);

                switch (userc) {
                        case 2:
                                create();
                                break;
                        case 3:
                                makefloc();
                                break;
                        case 4:
                                measure();
                                break;
                        case 5:
                                addagg();
                                break;
                        case 6:
                                connect();
                                break;
                        case 7:
                                if(aggsize!=0){
	                                measagg();
		        	}
		        	else{
		        		printf("No aggregate present. \n");
		        	}
                                break;
                        case 8:
                                outmic();
                                break;
                        default:
                                break;
                }
        } while (userc!=1);
}
