#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"

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

/* Random number generator */
#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
#define NTAB 32
#define EPS (1.2E-07)
#define MAX(a,b) (a>b)?a:b
#define MIN(a,b) (a<b)?a:b
#define PI 3.1415926
#define SYSIZE 100
#define SYSSIZE 100
#define MAXCYC 1000    /* maximum sintering cycles to use */
#define MAXSPH 10000   /* maximum number of elements in a spherical template */
#define C3S 1
#define C2S 2
#define C3A 3
#define C4AF 4

int *seed;
static int mask[SYSIZE+1][SYSIZE+1][SYSIZE+1];
static unsigned short int curvature [SYSSIZE+1] [SYSSIZE+1] [SYSSIZE+1];
long int volume[50],surface[50];
int nsph,xsph[MAXSPH],ysph[MAXSPH],zsph[MAXSPH];
long int nsolid[1500],nair[1500];

double ran1(idum)
int *idum;
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

/* routine to create a template for the sphere of interest of radius size */
/* to be used in curvature evaluation */
/* Called by: runsint */
/* Calls no other routines */
int maketemp(size)
        int size;
{
        int icirc,xval,yval,zval;
        float xtmp,ytmp;
        float dist;

/* determine and store the locations of all pixels in the 3-D sphere */
        icirc=0;
        for(xval=(-size);xval<=size;xval++){
                xtmp=(float)(xval*xval);
        for(yval=(-size);yval<=size;yval++){
                ytmp=(float)(yval*yval);
        for(zval=(-size);zval<=size;zval++){
                dist=sqrt(xtmp+ytmp+(float)(zval*zval));
                if(dist<=((float)size+0.5)){
                        icirc+=1;
                        if(icirc>=MAXSPH){
                            printf("Too many elements in sphere \n");
                            printf("Must change value of MAXSPH parameter \n");
                            printf("Currently set at %d \n",MAXSPH);
                            exit(1);
                        }
                        xsph[icirc]=xval;
                        ysph[icirc]=yval;
                        zsph[icirc]=zval;
                }
        }
        }
        }

/* return the number of pixels contained in sphere of radius (size+0.5) */
        return(icirc);
}

/* routine to count phase fractions (porosity and solids) */
/* Called by main routine */
/* Calls no other routines */
void phcount()
{
        long int npore,nsolid [37];
        int ix,iy,iz;

        npore=0;
        for(ix=1;ix<37;ix++){
		nsolid[ix]=0;
	}
        /* check all pixels in the 3-D system */
        for(ix=1;ix<=SYSSIZE;ix++){
        for(iy=1;iy<=SYSSIZE;iy++){
        for(iz=1;iz<=SYSSIZE;iz++){
                if(mask [ix] [iy] [iz]==0){
                        npore+=1;
               	}
                else{
                        nsolid[mask [ix] [iy] [iz]]+=1;
               	}
        }
        }
        }

        printf("Pores are: %ld \n",npore);
        printf("Solids are: %ld %ld %ld %ld %ld %ld\n",nsolid[1],nsolid[2],
             nsolid[3],nsolid[4],nsolid[5],nsolid[6]);
}

/* routine to return number of surface faces exposed to porosity */
/* for pixel located at (xin,yin,zin) */
/* Called by rhcalc */
/* Calls no other routines */
int surfpix(xin,yin,zin)
        int xin,yin,zin;
{
        int npix,ix1,iy1,iz1;

        npix=0;

        /* check each of the six immediate neighbors */
        /* using periodic boundary conditions */
        ix1=xin-1;
        if(ix1<1){ix1+=SYSSIZE;}
        if(mask[ix1][yin][zin]==0){
                npix+=1;
        }
        ix1=xin+1;
        if(ix1>SYSSIZE){ix1-=SYSSIZE;}
        if(mask[ix1][yin][zin]==0){
               	npix+=1;
        }
        iy1=yin-1;
        if(iy1<1){iy1+=SYSSIZE;}
        if(mask[xin][iy1][zin]==0){
                npix+=1;
        }
        iy1=yin+1;
        if(iy1>SYSSIZE){iy1-=SYSSIZE;}
        if(mask[xin][iy1][zin]==0){
                npix+=1;
        }
        iz1=zin-1;
        if(iz1<1){iz1+=SYSSIZE;}
        if(mask[xin][yin][iz1]==0){
                npix+=1;
        }
       	iz1=zin+1;
       	if(iz1>SYSSIZE){iz1-=SYSSIZE;}
       	if(mask[xin][yin][iz1]==0){
                npix+=1;
        }
        return(npix);
}

/* routine to return the current hydraulic radius for phase phin */
/* Calls surfpix */
/* Called by runsint */
float rhcalc(phin)
        int phin;
{
        int ix,iy,iz;
        long int porc,surfc;
        float rhval;

        porc=surfc=0;

        /* Check all pixels in the 3-D volume */
        for(ix=1;ix<=SYSSIZE;ix++){
        for(iy=1;iy<=SYSSIZE;iy++){
        for(iz=1;iz<=SYSSIZE;iz++){
                if(mask [ix] [iy] [iz]==phin){
                        porc+=1;
                        surfc+=surfpix(ix,iy,iz);
                }
        }
        }
        }

        printf("Phase area count is %ld \n",porc);
        printf("Phase surface count is %ld \n",surfc);
        rhval=(float)porc*6./(4.*(float)surfc);
        printf("Hydraulic radius is %f \n",rhval);
        return(rhval);
}

/* routine to return count of pixels in a spherical template which are phase */
/* phin or porosity (phase=0) */
/* Calls no other routines */
/* Called by sysinit */
int countem(xp,yp,zp,phin)
        int xp,yp,zp,phin;
{
        int xc,yc,zc;
        int cumnum,ic;

        cumnum=0;
        for(ic=1;ic<=nsph;ic++){
                xc=xp+xsph[ic];
                yc=yp+ysph[ic];
                zc=zp+zsph[ic];
                /* Use periodic boundaries */
                if(xc<1){xc+=SYSSIZE;}
                else if(xc>SYSSIZE){xc-=SYSSIZE;}
                if(yc<1){yc+=SYSSIZE;}
                else if(yc>SYSSIZE){yc-=SYSSIZE;}
                if(zc<1){zc+=SYSSIZE;}
                else if(zc>SYSSIZE){zc-=SYSSIZE;}
	
                if((xc!=xp)||(yc!=yp)||(zc!=zp)){

                if((mask [xc] [yc] [zc]==phin)||(mask [xc] [yc] [zc]==0)){
                        cumnum+=1;
                }
               	}
        }
        return(cumnum);
}

/* routine to initialize system by determining local curvature */
/* of all phase 1 and phase 2 pixels */
/* Calls countem */
/* Called by runsint */
void sysinit(ph1,ph2)
        int ph1,ph2;
{
        int count,xl,yl,zl;

       	count=0;
       	/* process all pixels in the 3-D box */
       	for(xl=1;xl<=SYSSIZE;xl++){
        for(yl=1;yl<=SYSSIZE;yl++){
        for(zl=1;zl<=SYSSIZE;zl++){

                /* determine local curvature */
                /* For phase 1 want to determine number of porosity pixels */
                /* (phase=0) in immediate neighborhood */
                if(mask [xl] [yl] [zl]==ph1){
                        count=countem(xl,yl,zl,0);
                }
                /* For phase 2 want to determine number of porosity or phase */
                /* 2 pixels in immediate neighborhood */
                if(mask [xl] [yl] [zl]==ph2){
                        count=countem(xl,yl,zl,ph2);
                }
                if((count<0)||(count>=nsph)){
                        printf("Error count is %d \n",count);
                        printf("xl %d  yl %d  zl %d \n",xl,yl,zl);
                }

                /* case where we have a phase 1 surface pixel */
                /* with non-zero local curvature */
                if((count>=0)&&(mask [xl] [yl] [zl]==ph1)){

                        curvature [xl] [yl] [zl]=count;
                       	/* update solid curvature histogram */
                       	nsolid[count]+=1;
                }
			
                /* case where we have a phase 2 surface pixel */
                if((count>=0)&&(mask [xl] [yl] [zl]==ph2)){

                        curvature [xl] [yl] [zl]=count;
                        /* update air curvature histogram */
                       	nair[count]+=1;
                }

        }
        }
       	}	/* end of xl loop */
}

/* routine to scan system and determine nsolid (ph2) and nair (ph1) */
/* histograms based on values in phase and curvature arrays */
/* Calls no other routines */
/* Called by runsint */
void sysscan(ph1,ph2)
        int ph1,ph2;
{
        int xd,yd,zd,curvval;

        /* Scan all pixels in 3-D system */
        for(xd=1;xd<=SYSSIZE;xd++){
        for(yd=1;yd<=SYSSIZE;yd++){
        for(zd=1;zd<=SYSSIZE;zd++){

                curvval=curvature [xd] [yd] [zd];
	
                if(mask [xd] [yd] [zd]==ph2){
                        nair[curvval]+=1;
                }
                else if (mask [xd] [yd] [zd]==ph1){
                        nsolid[curvval]+=1;
                }
        }	
        }
        }
}

/* routine to return how many cells of solid curvature histogram to use */
/* to accomodate nsearch pixels moving */
/* want to use highest values first */
/* Calls no other routines */
/* Called by movepix */
int procsol(nsearch)
        int nsearch;
{
        int valfound,i,stop;
        long int nsofar;

        /* search histogram from top down until cumulative count */
        /* exceeds nsearch */
        valfound=nsph-1;
        nsofar=0;
        stop=0;
        for(i=(nsph-1);((i>=0)&&(stop==0));i--){
                nsofar+=nsolid[i];
                if(nsofar>nsearch){
                        valfound=i;
                       	stop=1;
                }
        }
        return(valfound);
}

/* routine to determine how many cells of air curvature histogram to use */
/* to accomodate nsearch moving pixels */
/* want to use lowest values first */
/* Calls no other routines */
/* Called by movepix */
int procair(nsearch)
        int nsearch;
{
        int valfound,i,stop;
        long int nsofar;

        /* search histogram from bottom up until cumulative count */
        /* exceeds nsearch */
        valfound=0;
        nsofar=0;
        stop=0;
        for(i=0;((i<nsph)&&(stop==0));i++){
                nsofar+=nair[i];
               	if(nsofar>nsearch){
                        valfound=i;
                       	stop=1;
                }
        }
       	return(valfound);
}

/* routine to move requested number of pixels (ntomove) from highest */
/* curvature phase 1 (ph1) sites to lowest curvature phase 2 (ph2) sites */
/* Calls procsol and procair */
/* Called by runsint */
int movepix(ntomove,ph1,ph2)
        int ntomove,ph1,ph2;
{
        int xloc[2100],yloc[2100],zloc[2100];
        int count1,count2,ntot,countc,i,xp,yp,zp;
        int cmin,cmax,cfg;
        int alldone;
        long int nsolc,nairc,nsum,nsolm,nairm,nst1,nst2,next1,next2;
        float pck,plsol,plair;

        alldone=0;
        /* determine critical values for removal and placement */
        count1=procsol(ntomove);
        nsum=0;
        cfg=0;
        cmax=count1;
        for(i=nsph;i>count1;i--){
                if((nsolid[i]>0)&&(cfg==0)){
                        cfg=1;
                       	cmax=i;
                }
               	nsum+=nsolid[i];
        }
        /* Determine movement probability for last cell */
        plsol=(float)(ntomove-nsum)/(float)nsolid[count1];
        next1=ntomove-nsum;
        nst1=nsolid[count1];	

        count2=procair(ntomove);
        nsum=0;
        cmin=count2;
        cfg=0;
        for(i=0;i<count2;i++){
                if((nair[i]>0)&&(cfg==0)){
                        cfg=1;
                       	cmin=i;
                }
               	nsum+=nair[i];
        }
        /* Determine movement probability for last cell */
        plair=(float)(ntomove-nsum)/(float)nair[count2];
        next2=ntomove-nsum;
        nst2=nair[count2];
	
        /* Check to see if equilibrium has been reached --- */
        /* no further increase in hydraulic radius is possible */
        if(cmin>=cmax){
                alldone=1;
                printf("Stopping - at equilibrium \n");
                printf("cmin- %d  cmax- %d \n",cmin,cmax);
                return(alldone); 
        }

        /* initialize counters for performing sintering */
        ntot=0;
        nsolc=0;
        nairc=0;
        nsolm=0;
        nairm=0;

        /* Now process each pixel in turn */
        for(xp=1;xp<=SYSSIZE;xp++){
        for(yp=1;yp<=SYSSIZE;yp++){
        for(zp=1;zp<=SYSSIZE;zp++){
		
                countc=curvature [xp] [yp] [zp];
                /* handle phase 1 case first */
                if(mask [xp] [yp] [zp]==ph1){
                        if(countc>count1){
                                /* convert from phase 1 to phase 2 */
                                mask [xp] [yp] [zp]=ph2;

                                /* update appropriate histogram cells */
                                nsolid[countc]-=1;
                                nair[countc]+=1;
                                /* store the location of the modified pixel */
                                ntot+=1;
                                xloc[ntot]=xp;
                                yloc[ntot]=yp;
                                zloc[ntot]=zp;
                        }
                        if(countc==count1){
                                nsolm+=1;
                             /* generate probability for pixel being removed */
                                pck=ran1(seed);
                                if((pck<0)||(pck>1.0)){pck=1.0;}

                 if(((pck<plsol)&&(nsolc<next1))||((nst1-nsolm)<(next1-nsolc))){
                                        nsolc+=1;
                                        /* convert phase 1 pixel to phase 2 */
                                        mask [xp] [yp] [zp]=ph2;

                                        /* update appropriate histogram cells */
                                        nsolid[count1]-=1;
                                        nair[count1]+=1;
                                  /* store the location of the modified pixel */
                                        ntot+=1;
                                        xloc[ntot]=xp;
                                        yloc[ntot]=yp;
                                        zloc[ntot]=zp;
                                 }
                                 }
                        }

                        /* handle phase 2 case here */
                        else if (mask [xp] [yp] [zp]==ph2){
                                if(countc<count2){
                                        /* convert phase 2 pixel to phase 1 */
                                        mask [xp] [yp] [zp]=ph1;

                                        nsolid[countc]+=1;
                                        nair[countc]-=1;
                                        ntot+=1;
                                        xloc[ntot]=xp;
                                        yloc[ntot]=yp;
                                        zloc[ntot]=zp;
                                }
                                if(countc==count2){
                                        nairm+=1;
                                        pck=ran1(seed);
                                        if((pck<0)||(pck>1.0)){pck=1.0;}

             if(((pck<plair)&&(nairc<next2))||((nst2-nairm)<(next2-nairc))){
                                                nairc+=1;
                                                /* convert phase 2 to phase 1 */
                                                mask [xp] [yp] [zp]=ph1;

                                                nsolid[count2]+=1;
                                                nair[count2]-=1;
                                                ntot+=1;
                                                xloc[ntot]=xp;
                                                yloc[ntot]=yp;
                                                zloc[ntot]=zp;
                                        }
                                }
                        }
		
        } /* end of zp loop */
        } /* end of yp loop */
       	} /* end of xloop */
        printf("ntot is %d \n",ntot);
       	return(alldone);
}

/* routine to execute user input number of cycles of sintering algorithm */
/* Calls maketemp, rhcalc, sysinit, sysscan, and movepix */
/* Called by main routine */
void sinter3d(ph1id,ph2id,rhtarget)
	int ph1id,ph2id;
	float rhtarget;
{
        int natonce,ncyc,i,rade,j,rflag;
        int keepgo;
        long int curvsum1,curvsum2,pixsum1,pixsum2;
        float rhnow,avecurv1,avecurv2;
	
        /* initialize the solid and air count histograms */
        for(i=0;i<=499;i++){
                nsolid[i]=0;
                nair[i]=0;
        }

        /* Obtain needed user input */
	natonce=200;
	rade=3;
        rflag=0;   /* always initialize system */

        nsph=maketemp(rade);
        printf("nsph is %d \n",nsph);
        if(rflag==0){
                sysinit(ph1id,ph2id);
        }
       	else{
                sysscan(ph1id,ph2id);
        }
        i=0;
        rhnow=rhcalc(ph1id);
        while((rhnow<rhtarget)&&(i<MAXCYC)){
                printf("Now: %f  Target: %f \n",rhnow,rhtarget);
                i+=1;
                printf("Cycle: %d \n",i);
                keepgo=movepix(natonce,ph1id,ph2id);
                /* If equilibrium is reached, then return to calling routine */
                if(keepgo==1){
                        return;
                }
                curvsum1=0;
                curvsum2=0;
                pixsum1=0;
                pixsum2=0;
                /* Determine average curvatures for phases 1 and 2 */
                for(j=0;j<=nsph;j++){
                        pixsum1+=nsolid[j];
                       	curvsum1+=(j*nsolid[j]);
                        pixsum2+=nair[j];
                        curvsum2+=(j*nair[j]);
                }
                avecurv1=(float)curvsum1/(float)pixsum1;
                avecurv2=(float)curvsum2/(float)pixsum2;
                printf("Ave. solid curvature: %f \n",avecurv1);
                printf("Ave. air curvature: %f \n",avecurv2);
                rhnow=rhcalc(ph1id);
        } 
}

void stat3d(){
        int valin,ix,iy,iz;
        int ix1,iy1,iz1,k;
        long int voltot,surftot;

        for(ix=0;ix<=42;ix++){
                volume[ix]=surface[ix]=0;
        }

/* Read in image and accumulate volume totals */
        for(iz=1;iz<=SYSIZE;iz++){
        for(iy=1;iy<=SYSIZE;iy++){
        for(ix=1;ix<=SYSIZE;ix++){
                valin=mask[ix][iy][iz];
                volume[valin]+=1;
        }
        }
        }
	

        for(iz=1;iz<=SYSIZE;iz++){
        for(iy=1;iy<=SYSIZE;iy++){
        for(ix=1;ix<=SYSIZE;ix++){
                if(mask [ix] [iy] [iz]!=0){
                        valin=mask [ix] [iy] [iz];
                /* Check six neighboring pixels for porosity */
                        for(k=1;k<=6;k++){

                                 switch (k){
                                        case 1:
                                                ix1=ix-1;
                                                if(ix1<1){ix1+=SYSIZE;}
                                                iy1=iy;
                                                iz1=iz;
                                                break;
                                        case 2:
                                                ix1=ix+1;
                                                if(ix1>SYSIZE){ix1-=SYSIZE;}
                                                iy1=iy;
                                                iz1=iz;
                                                break;
                                        case 3:
                                                iy1=iy-1;
                                                if(iy1<1){iy1+=SYSIZE;}
                                                ix1=ix;
                                                iz1=iz;
                                                break;
                                        case 4:
                                                iy1=iy+1;
                                                if(iy1>SYSIZE){iy1-=SYSIZE;}
                                                ix1=ix;
                                                iz1=iz;
                                                break;
                                        case 5:
                                                iz1=iz-1;
                                                if(iz1<1){iz1+=SYSIZE;}
                                                iy1=iy;
                                                ix1=ix;
                                                break;
                                        case 6:
                                                iz1=iz+1;
                                                if(iz1>SYSIZE){iz1-=SYSIZE;}
                                                iy1=iy;
                                                ix1=ix;
                                                break;
                                        default:
                                                break;
                                }
        if((ix1<1)||(iy1<1)||(iz1<1)||(ix1>SYSIZE)||(iy1>SYSIZE)||(iz1>SYSIZE)){
                                        printf("%d %d %d \n",ix1,iy1,iz1);
                                        exit(1);
                                }
                                if(mask[ix1] [iy1] [iz1]==0){
                                        surface[valin]+=1;
                                }
                        }
                }
        }
        }
        }

        printf("Phase    Volume      Surface     Volume    Surface \n");
        printf(" ID      count        count      fraction  fraction \n");
/* Only include clinker phases in surface area fraction calculation */
        surftot=surface[1]+surface[2]+surface[3]+surface[4];
        voltot=volume[1]+volume[2]+volume[3]+volume[4];
        k=0;
        printf("  %d    %8ld     %8ld  \n",k,volume[0],surface[0]);
        for(k=1;k<=4;k++){
        printf("  %d    %8ld     %8ld     %.5f   %.5f\n",k,volume[k],surface[k],
                    (float)volume[k]/(float)voltot,(float)surface[k]/(float)surftot);
        }
        printf("Total  %8ld     %8ld\n\n\n",voltot,surftot);
        for(k=5;k<=11;k++){
		printf("  %d    %8ld     %8ld\n",k,volume[k],surface[k]);
        }
	printf(" 20    %8ld     %8ld\n",volume[20],surface[20]);
        for(k=24;k<=27;k++){
		printf(" %d    %8ld     %8ld\n",k,volume[k],surface[k]);
        }
	printf(" 28    %8ld     %8ld\n",volume[28],surface[28]);
}

void rand3d(phasein,phaseout,filecorr,xpt)
	int phasein,phaseout;
	char filecorr[80];
	float xpt;
{
        int syssize,ires;
        float snew,s2,ss,sdiff,xtmp,ytmp;
        static float normm[SYSIZE+1][SYSIZE+1][SYSIZE+1];
        static float res[SYSIZE+1][SYSIZE+1][SYSIZE+1];
        static float filter [32][32][32];
        int done,r[61];
        static float s[61],xr[61],sum[502];
        float val2;
        float t1,t2,x1,x2,u1,u2,xrad,resmax,resmin;
        float xtot,filval,radius,sect,sumtot,vcrit;
        int valin,r1,r2,i1,i2,i3,i,j,k,j1,k1;
        int ido,iii,jjj,ix,iy,iz,index;
	FILE *corrfile;

/* Create the Gaussian noise image */
       i1=i2=i3=1;
       for(i=1;i<=((SYSIZE*SYSIZE*SYSIZE)/2);i++){
             u1=ran1(seed);
             u2=ran1(seed);
             t1=2.*PI*u2;
             t2=sqrt(-2.*log(u1));
             x1=cos(t1)*t2;
             x2=sin(t1)*t2;
             normm[i1][i2][i3]=x1;
             i1+=1;
             if(i1>SYSIZE){
                  i1=1;
                  i2+=1;
                  if(i2>SYSIZE){
                        i2=1;
                        i3+=1;
                  }
             }
             normm[i1][i2][i3]=x2;
             i1+=1;
             if(i1>SYSIZE){
                  i1=1;
                  i2+=1;
                  if(i2>SYSIZE){
                        i2=1;
                        i3+=1;
                  }
             }
       }
	
/* Now perform the convolution */
      corrfile=fopen(filecorr,"r");

      fscanf(corrfile,"%d",&ido);
      printf("Number of points in correlation file is %d \n",ido);
      for(i=1;i<=ido;i++){
            fscanf(corrfile,"%d %f",&valin,&val2);
            r[i]=valin;
            s[i]=val2;
            xr[i]=(float)r[i];
       }
       fclose(corrfile);
       ss=s[1];
       s2=ss*ss;
/* Load up the convolution matrix */
       sdiff=ss-s2;
       for(i=0;i<31;i++){
           iii=i*i;
           for(j=0;j<31;j++){
               jjj=j*j;
               for(k=0;k<31;k++){
                  xtmp=(float)(iii+jjj+k*k);
                  radius=sqrt(xtmp);
                  r1=(int)radius+1;
                  r2=r1+1;
                  if(s[r1]<0.0){
                      printf("%d and %d %f and %f with xtmp of %f\n",r1,r2,s[r1],s[r2],xtmp);
                      fflush(stdout); 
                      exit(1);
                  }
                  xrad=radius+1-r1;
                  filval=s[r1]+(s[r2]-s[r1])*xrad;
                  filter[i+1][j+1][k+1]=(filval-s2)/sdiff;
               }
           }
       }
/* Now filter the image maintaining periodic boundaries */
       resmax=0.0;
       resmin=1.0;
       for(i=1;i<=SYSIZE;i++){
       for(j=1;j<=SYSIZE;j++){
       for(k=1;k<=SYSIZE;k++){
            res[i][j][k]=0.0;
            if((float)mask[i][j][k]==phasein){
                  for(ix=1;ix<=31;ix++){
                      i1=i+ix-1;
                      if(i1<1){i1+=SYSIZE;}
                      else if(i1>SYSIZE){i1-=SYSIZE;}
                      for(iy=1;iy<=31;iy++){
                           j1=j+iy-1;
                           if(j1<1){j1+=SYSIZE;}
                           else if(j1>SYSIZE){j1-=SYSIZE;}
                           for(iz=1;iz<=31;iz++){
                                k1=k+iz-1;
                                if(k1<1){k1+=SYSIZE;}
                                else if(k1>SYSIZE){k1-=SYSIZE;}
                           res[i][j][k]+=normm[i1][j1][k1]*filter[ix][iy][iz];
                           }
                      }
                  }
                  if(res[i][j][k]>resmax){resmax=res[i][j][k];}
                  if(res[i][j][k]<resmin){resmin=res[i][j][k];}
            }
       }
       }
       }

/* Now threshold the image */
       sect=(resmax-resmin)/500.;
       for(i=1;i<=500;i++){
           sum[i]=0.0;
       }
       xtot=0.0;
       for(i=1;i<=SYSIZE;i++){
       for(j=1;j<=SYSIZE;j++){
       for(k=1;k<=SYSIZE;k++){
           if((float)mask[i][j][k]==phasein){
                xtot+=1.0;
                index=1+(int)((res[i][j][k]-resmin)/sect);
                if(index>500){index=500;}
                sum[index]+=1.0;
           }
       }
       }
       }
/* Determine which bin to choose for correct thresholding */
       sumtot=vcrit=0.0;
       done=0;
       for(i=1;((i<=500)&&(done==0));i++){
           sumtot+=sum[i]/xtot;
           if(sumtot>xpt){
                ytmp=(float)i;
                vcrit=resmin+(resmax-resmin)*(ytmp-0.5)/500.;
                done=1;
           }
       }
       printf("Critical volume fraction is %f\n",vcrit);
       ires=0;
     
       for(k=1;k<=SYSIZE;k++){
       for(j=1;j<=SYSIZE;j++){
       for(i=1;i<=SYSIZE;i++){
            if((float)mask[i][j][k]==phasein){
                if(res[i][j][k]>vcrit){
                       mask[i][j][k]=phaseout;
                }
            }
       }
       }
       }

}

int main(){
	int i,j,k,iseed,alumflag,alumval,alum2,valin;
	float volin,volf[5],surff[5],rhtest,rdesire;
	char filen[80],fileout[80],filecem[80],filec3s[80],filesil[80],filealum[80];
	FILE *infile,*outfile,*testfile;

	/* Seed the random number generator */
        printf("Enter random number seed (negative integer) \n");
        scanf("%d",&iseed);
        printf("%d\n",iseed);
        seed=(&iseed);

	/* Read in the parameters to use */
        printf("Enter name of cement microstructure image file\n");
        scanf("%s",filen);
        printf("%s\n",filen);

	/* Set up the correlation filenames */
        printf("Enter root name of cement correlation files\n");
        scanf("%s",filecem);
        printf("%s\n",filecem);
	sprintf(filesil,"%s",filecem);
	strcat(filesil,".sil");
	sprintf(filec3s,"%s",filecem);
	strcat(filec3s,".c3s");
	sprintf(filealum,"%s",filecem);
	alumflag=1;
	alumval=4;
	strcat(filealum,".c4f");
	testfile=fopen(filealum,"r");
	if(testfile==NULL){
		alumflag=0;
		sprintf(filealum,"%s",filecem);
		strcat(filealum,".c3a");
		alumval=3;
	}
	else{
		fclose(testfile);
	}

       printf("Enter name of new cement microstructure image file\n");
       scanf("%s",fileout);
       printf("%s\n",fileout);
	for(i=1;i<=4;i++){
		scanf("%f",&volin);
		volf[i]=volin;
		printf("%f\n",volf[i]); 
		scanf("%f",&volin);
		surff[i]=volin; 
		printf("%f\n",surff[i]); 
	}

/* Read in the original microstructure image file */
        infile=fopen(filen,"r");
        if (infile == NULL) {
                perror("Error opening input file");
                return 1;
        }

        for(k=1;k<=SYSIZE;k++){
        for(j=1;j<=SYSIZE;j++){
        for(i=1;i<=SYSIZE;i++){
             fscanf(infile,"%d",&valin);
             mask[i][j][k]=valin;
	     curvature[i][j][k]=0;
        }
        }
        }
        fclose(infile);

	/* First filtering */
	volin=volf[1]+volf[2];
        if(volin<1.0){
		rand3d(1,alumval,filesil,volin);

		/* First sintering */
		stat3d();
		rdesire=(surff[1]+surff[2])*(float)(surface[1]+surface[alumval]);
                if(rdesire!=0.0){
		if((int)rdesire<surface[1]){
			rhtest=(6./4.)*(float)(volume[1])/rdesire;
			sinter3d(1,alumval,rhtest); 
		}
		else{
			rdesire=(surff[3]+surff[4])*(float)(surface[1]+surface[alumval]);
                	if(rdesire!=0.0){
				rhtest=(6./4.)*(float)(volume[alumval])/rdesire;
				sinter3d(alumval,1,rhtest); 
			}
		}
		}
	}

	/* Second filtering */
        if((volf[1]+volf[2])>0.0){
	volin=volf[1]/(volf[1]+volf[2]);
        if(volin<1.0){
		rand3d(1,2,filec3s,volin);

		/* Second sintering */
		stat3d();
		rdesire=(surff[1]/(surff[1]+surff[2]))*(float)(surface[1]+surface[2]);
               	if(rdesire!=0.0){
		if((int)rdesire<surface[1]){
			rhtest=(6./4.)*(float)(volume[1])/rdesire;
			sinter3d(1,2,rhtest); 
		}
		else{
			rdesire=(surff[2]/(surff[1]+surff[2]))*(float)(surface[1]+surface[2]);
                	if(rdesire!=0.0){
				rhtest=(6./4.)*(float)(volume[2])/rdesire;
				sinter3d(2,1,rhtest); 
			}
		}
		}
 	}
        }

	/* Third (final) filtering */
	if(alumval==4){
		volin=volf[4]/(volf[4]+volf[3]);
		alum2=3;
	}
	else{
		volin=volf[3]/(volf[4]+volf[3]);
		alum2=4;
	}
        if(volin<1.0){
	rand3d(alumval,alum2,filealum,volin);

	/* Third (final) sintering */
	stat3d();
	if(alumval==4){
		rdesire=(surff[4]/(surff[3]+surff[4]))*(float)(surface[3]+surface[4]);
        if(rdesire!=0.0){
	if((int)rdesire<surface[4]){
		rhtest=(6./4.)*(float)(volume[4])/rdesire;
		sinter3d(alumval,alum2,rhtest); 
	}
	else{
		rdesire=(surff[3]/(surff[3]+surff[4]))*(float)(surface[3]+surface[4]);
        	if(rdesire!=0.0){
			rhtest=(6./4.)*(float)(volume[3])/rdesire;
			sinter3d(alum2,alumval,rhtest); 
		}
	}
	}
	}
	else{
		rdesire=(surff[3]/(surff[3]+surff[4]))*(float)(surface[3]+surface[4]);
        if(rdesire!=0.0){
		if((int)rdesire<surface[3]){
			rhtest=(6./4.)*(float)(volume[3])/rdesire;
			sinter3d(alumval,alum2,rhtest); 
		}
		else{
			rdesire=(surff[4]/(surff[3]+surff[4]))*(float)(surface[3]+surface[4]);
       	 		if(rdesire!=0.0){
				rhtest=(6./4.)*(float)(volume[4])/rdesire;
				sinter3d(alum2,alumval,rhtest); 
			}
		}
	}
	}
	}

	/* Output final microstructure */
       outfile=fopen(fileout,"w");
     
       for(k=1;k<=SYSIZE;k++){
       for(j=1;j<=SYSIZE;j++){
       for(i=1;i<=SYSIZE;i++){
            fprintf(outfile,"%2d\n",mask[i][j][k]);
       }
       }
       }

}
