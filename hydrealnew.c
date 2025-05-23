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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include "common.h"


#define AGRATE 0.25        /* Probability of gypsum absorption by CSH */

/* routine to select a new neighboring location to (xloc, yloc, zloc) */ 
/* for a diffusing species */
/* Returns a prime number flag indicating direction chosen */
/* Calls ran1 */
/* Called by movecsh, extettr, extfh3, movegyp, extafm, moveettr, */
/* extpozz, movefh3, movech, extc3ah6, movec3a */
/* extfreidel, movecacl2, extstrat, moveas */
int moveone(xloc,yloc,zloc,act,sumold)
        int *xloc,*yloc,*zloc,*act,sumold;
{
        int plok,sumnew,xl1,yl1,zl1,act1;

       	sumnew=1;
        /* store the input values for location */
        xl1=(*xloc);
        yl1=(*yloc);
        zl1=(*zloc);
        act1=(*act);

        /* Choose one of six directions (at random) for the new */
        /* location */
        plok=6.*ran1(seed);
        if((plok>5)||(plok<0)){plok=5;}

        switch (plok){
                case 0: 
                        xl1-=1;
                        act1=1;
                        if(xl1<0){xl1=(SYSIZEM1);}
                        if(sumold%2!=0){sumnew=2;}
                        break;
                case 1:
                        xl1+=1;
                        act1=2;
                        if(xl1>=SYSIZE){xl1=0;}
                        if(sumold%3!=0){sumnew=3;}
                        break;
                case 2:
                        yl1-=1;
                        act1=3;
                        if(yl1<0){yl1=(SYSIZEM1);}
                        if(sumold%5!=0){sumnew=5;}
                        break;
                case 3: 
                        yl1+=1;
                        act1=4;
                        if(yl1>=SYSIZE){yl1=0;}
                        if(sumold%7!=0){sumnew=7;}
                        break;
                case 4:
                        zl1-=1;
                        act1=5;
                        if(zl1<0){zl1=(SYSIZEM1);}
                        if(sumold%11!=0){sumnew=11;}
                        break;
                case 5: 
                        zl1+=1;
                        act1=6;
                        if(zl1>=SYSIZE){zl1=0;}
                        if(sumold%13!=0){sumnew=13;}
                        break;
                default:	
                        break;
        }

        /* Return the new location */
        *xloc=xl1;
        *yloc=yl1;
        *zloc=zl1;
        *act=act1;
        /* sumnew returns a prime number indicating that a specific direction */
        /* has been chosen */
        return(sumnew);
}

/* routine to return count of number of neighboring pixels for pixel */
/* (xck,yck,zck) which are not phase ph1, ph2, or ph3 which are input as */
/* parameters */
/* Calls no other routines */
/* Called by extettr, extfh3, extch, extafm, extpozz, extc3ah6 */
/* extfreidel, extcsh, and extstrat */
int edgecnt(xck,yck,zck,ph1,ph2,ph3)
        int xck,yck,zck,ph1,ph2,ph3;
{
       	int ixe,iye,ize,edgeback,x2,y2,z2,check;

/* counter for number of neighboring pixels which are not ph1, ph2, or ph3 */
        edgeback=0;

/* Examine all pixels in a 3*3*3 box centered at (xck,yck,zck) */
/* except for the central pixel */
       	for(ixe=(-1);ixe<=1;ixe++){
                x2=xck+ixe;
       	for(iye=(-1);iye<=1;iye++){
                y2=yck+iye;
       	for(ize=(-1);ize<=1;ize++){

                 if((ixe!=0)||(iye!=0)||(ize!=0)){
                        z2=zck+ize;
                        /* adjust to maintain periodic boundaries */
                        if(x2<0){x2=(SYSIZEM1);}
                        else if(x2>=SYSIZE){x2=0;}
                        if(y2<0){y2=(SYSIZEM1);}
                        else if(y2>=SYSIZE){y2=0;}
                        if(z2<0){z2=(SYSIZEM1);}
                        else if(z2>=SYSIZE){z2=0;}
                        check=mic[x2][y2][z2];
                        if((check!=ph1)&&(check!=ph2)&&(check!=ph3)){
                                  edgeback+=1;
                        }
                }
        }
       	}
       	}
       	/* return number of neighboring pixels which are not ph1, ph2, or ph3 */
       	return(edgeback);
}

/* routine to add extra CSH when diffusing CSH reacts */
/* Called by movecsh */
/* Calls edgecnt */
void extcsh()
{
        int numnear,sump,xchr,ychr,zchr,fchr,i1,plok,check,msface;
        long int tries;

        fchr=0;
        tries=0;
        /* locate CSH at random location */
        /* in pore space in contact with at least another CSH or C3S or C2S */
        while(fchr==0){
                tries+=1;
                /* generate a random location in the 3-D system */
                xchr=(int)((float)SYSIZE*ran1(seed));
                ychr=(int)((float)SYSIZE*ran1(seed));
                zchr=(int)((float)SYSIZE*ran1(seed));
                if(xchr>=SYSIZE){xchr=0;}
                if(ychr>=SYSIZE){ychr=0;}
                if(zchr>=SYSIZE){zchr=0;}
                check=mic[xchr][ychr][zchr];

                /* if location is porosity, locate the CSH there */
                if(check==POROSITY){
                        numnear=edgecnt(xchr,ychr,zchr,CSH,C3S,C2S);
                        /* be sure that at least one neighboring pixel */
                        /* is C2S, C3S, or diffusing CSH */
                        if((numnear<26)||(tries>5000)){
                                mic[xchr][ychr][zchr]=CSH;
				count[CSH]+=1;
				count[POROSITY]-=1;
				cshage[xchr][ychr][zchr]=cyccnt;
				if(cshgeom==1){
					msface=(int)(3.*ran1(seed)+1.);
					if(msface>3){msface=1;}
					faces[xchr][ychr][zchr]=msface;
					ncshplateinit+=1;
				}
                               	fchr=1;
                        }
                }
        }
}

/* routine to move a diffusing CSH species */
/* Inputs: current location (xcur,ycur,zcur) and flag indicating if final */
/* step in diffusion process */
/* Returns flag indicating action taken (reaction or diffusion/no movement) */
/* Calls moveone,extcsh */
/* Called by hydrate */
int movecsh(xcur,ycur,zcur,finalstep,cycorig)
        int xcur,ycur,zcur,finalstep,cycorig;
{
        int xnew,ynew,znew,plok,action,sumback,sumin,check;
	int msface,mstest,mstest2;
	float prcsh,prcsh1,prcsh2,prtest;

        action=0;
        /* Store current location of species */
        xnew=xcur;
        ynew=ycur;
        znew=zcur;
        sumin=1;
        sumback=moveone(&xnew,&ynew,&znew,&action,sumin);
	if(cshgeom==1){
	/* Determine eligible faces based on direction of move */
		if(xnew!=xcur){
			mstest=1;
			mstest2=2;
		}
		if(ynew!=ycur){
			mstest=2;
			mstest2=3;
		}
		if(znew!=zcur){
			mstest=3;
			mstest2=1;
		}
	}

        if(action==0){printf("Error in value of action \n");}
        check=mic[xnew][ynew][znew];


      /* if new location is solid CSH and plate growth is favorable, */
      /* then convert diffusing CSH species to solid CSH */
       prcsh=ran1(seed);
       if((check==CSH)&&((cshgeom==0)||(faces[xnew][ynew][znew]==0)||(faces[xnew][ynew][znew]==mstest)||(faces[xnew][ynew][znew]==mstest2))){
           /* decrement count of diffusing CSH species */
              count[DIFFCSH]-=1;
           /* and increment count of solid CSH if needed */
  		prtest=molarvcsh[cyccnt]/molarvcsh[cycorig];
                prcsh1=ran1(seed);
		if(prcsh1<=prtest){
                   mic[xcur][ycur][zcur]=CSH;
		   if(cshgeom==1){
			   faces[xcur][ycur][zcur]=faces[xnew][ynew][znew];
		           ncshplategrow+=1;
		   }
         	   cshage[xcur][ycur][zcur]=cyccnt;
                   count[CSH]+=1;
		}
      		else{
			mic[xcur][ycur][zcur]=POROSITY;
			count[POROSITY]+=1;
		}
      /* May need extra solid CSH if temperature goes down with time */
		if(prtest>1.0){
                        prcsh2=ran1(seed);
			if(prcsh2<(prtest-1.0)){
				extcsh();
			}
		}
                action=0;
	}
	/* Changed prcsh limit from 0.1 to 0.01 for CH test 1/27/05 */
       else if((check==SLAGCSH)||(check==POZZCSH)||(finalstep==1)||
       	(((check==C3S)||(check==C2S))&&(prcsh<0.001))||
         (((check==C3A)||(check==C4AF))&&(prcsh<0.2))||
         ((check==CH)&&(prcsh<0.01))||
         (check==CACO3)||(check==INERT)){
           /* decrement count of diffusing CSH species */
              count[DIFFCSH]-=1;
           /* and increment count of solid CSH if needed */
  		prtest=molarvcsh[cyccnt]/molarvcsh[cycorig];
                prcsh1=ran1(seed);
		if(prcsh1<=prtest){
                   mic[xcur][ycur][zcur]=CSH;
         	   cshage[xcur][ycur][zcur]=cyccnt;
		   if(cshgeom==1){
		           msface=(int)(2.*ran1(seed)+1.);
			   if(msface>2){msface=1;}
			   if(msface==1){
			      faces[xcur][ycur][zcur]=mstest;
			   }
			   else{
			      faces[xcur][ycur][zcur]=mstest2;
			   }
		           ncshplateinit+=1;
		   }
                   count[CSH]+=1;
		}
      		else{
			mic[xcur][ycur][zcur]=POROSITY;
			count[POROSITY]+=1;
		}
      /* May need extra solid CSH if temperature goes down with time */
		if(prtest>1.0){
                        prcsh2=ran1(seed);
			if(prcsh2<(prtest-1.0)){
				extcsh();
			}
		}
                action=0;
	}

        if(action!=0){
        /* if diffusion step is possible, perform it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                       	mic[xnew][ynew][znew]=DIFFCSH;
              	 }
                else{
                        /* indicate that diffusing CSH species remained */
                     	/* at original location */
                       	action=7;
                }
        }
        return(action);
}

/* routine to add extra FH3 when gypsum, hemihydrate, anhydrite, CAS2, or */
/* CaCl2 reacts with C4AF at location (xpres,ypres,zpres) */
/* Called by movegyp, moveettr, movecas2, movehem, moveanh, and movecacl2 */
/* Calls moveone and edgecnt */
void extfh3(xpres,ypres,zpres)
        int xpres,ypres,zpres;
{
        int multf,numnear,sump,xchr,ychr,zchr,check,fchr,i1,plok,newact;
       	long int tries;

/* first try 6 neighboring locations until      */
/*	a) successful				*/
/*	b) all 6 sites are tried and full or    */
/*	c) 500 tries are made 			*/
        fchr=0;
       	sump=1;
       	for(i1=1;((i1<=500)&&(fchr==0)&&(sump!=30030));i1++){
		
                /* choose a neighbor at random */
                xchr=xpres;
                ychr=ypres;
                zchr=zpres;
                newact=0;
                multf=moveone(&xchr,&ychr,&zchr,&newact,sump);
                if(newact==0){printf("Error in value of newact in extfh3 \n");}
                check=mic[xchr][ychr][zchr];	 

               	/* if neighbor is porosity   */
                /* then locate the FH3 there */
                if(check==POROSITY){
                        mic[xchr][ychr][zchr]=FH3;
			count[FH3]+=1;
			count[POROSITY]-=1;
                       	fchr=1;
                }
               	else{
                        sump*=multf;
                }
        }

/* if no neighbor available, locate FH3 at random location */
/* in pore space in contact with at least one FH3 */
        tries=0;
        while(fchr==0){
                tries+=1;
                /* generate a random location in the 3-D system */
                xchr=(int)((float)SYSIZE*ran1(seed));
                ychr=(int)((float)SYSIZE*ran1(seed));
                zchr=(int)((float)SYSIZE*ran1(seed));
                if(xchr>=SYSIZE){xchr=0;}
                if(ychr>=SYSIZE){ychr=0;}
                if(zchr>=SYSIZE){zchr=0;}
                check=mic[xchr][ychr][zchr];
                /* if location is porosity, locate the FH3 there */
                if(check==POROSITY){
                        numnear=edgecnt(xchr,ychr,zchr,FH3,FH3,DIFFFH3);
                        /* be sure that at least one neighboring pixel */
                        /* is FH3 or diffusing FH3 */
                        if((numnear<26)||(tries>5000)){
                                mic[xchr][ychr][zchr]=FH3;
				count[FH3]+=1;
				count[POROSITY]-=1;
                               	fchr=1;
                        }
                }
        }
}

/* routine to add extra ettringite when gypsum, anhydrite, or hemihydrate */
/* reacts with aluminates addition adjacent to location (xpres,ypres,zpres) */
/* in a fashion to preserve needle growth */
/* etype=0 indicates primary ettringite */
/* etype=1 indicates iron-rich stable ettringite */
/* Returns flag indicating action taken */
/* Calls moveone and edgecnt */
/* Called by movegyp, movehem, moveanh, and movec3a */
int extettr(xpres,ypres,zpres,etype)
        int xpres,ypres,zpres,etype;
{
        int check,newact,multf,numnear,sump,xchr,ychr,zchr,fchr,i1,plok;
	int numalum,numsil;
        float pneigh,ptest;
        long int tries;

/* first try neighboring locations until        */
/*	a) successful				*/
/*	c) 1000 tries are made 			*/
        fchr=0;
        sump=1;
        /* Note that 30030 = 2*3*5*7*11*13 */
	/* indicating that all six sites have been tried */
        for(i1=1;((i1<=1000)&&(fchr==0));i1++){

/* determine location of neighbor (using periodic boundaries) */
                xchr=xpres;
                ychr=ypres;
                zchr=zpres;
                newact=0;
                multf=moveone(&xchr,&ychr,&zchr,&newact,sump);
                if(newact==0){printf("Error in value of action \n");}

                check=mic[xchr][ychr][zchr];

                /* if neighbor is porosity, and conditions are favorable */
                /* based on number of neighboring ettringite, C3A, or C4AF */
                /* pixels then locate the ettringite there */
                if(check==POROSITY){
                		/* be sure ettringite doesn't touch C3S */
                				numsil=edgecnt(xchr,ychr,zchr,C3S,C2S,C3S);
                           numsil=26-numsil;
			if(etype==0){
	                        numnear=edgecnt(xchr,ychr,zchr,ETTR,ETTR,ETTR);
	                        numalum=edgecnt(xchr,ychr,zchr,C3A,C3A,C3A);
				numalum=26-numalum;
			}
			else{
	                        numnear=edgecnt(xchr,ychr,zchr,ETTRC4AF,ETTRC4AF,ETTRC4AF);
	                        numalum=edgecnt(xchr,ychr,zchr,C4AF,C4AF,C4AF);
				numalum=26-numalum;
			}
                        pneigh=(float)(numnear+1)/26.0;
                        pneigh*=pneigh;
			if(numalum<=1){pneigh=0.0;}
			if(numalum>=2){pneigh+=0.5;}
			if(numalum>=3){pneigh+=0.25;}
			if(numalum>=5){pneigh+=0.25;}
                        ptest=ran1(seed);
                        if(numsil<1){
                        if(pneigh>=ptest){
				if(etype==0){
	                                mic[xchr][ychr][zchr]=ETTR;
					count[ETTR]+=1;
				}
				else{
	                                mic[xchr][ychr][zchr]=ETTRC4AF;
					count[ETTRC4AF]+=1;
				}
                                fchr=1;
				count[POROSITY]-=1;
                        }
                        }
                }
        }

/* if no neighbor available, locate ettringite at random location */
/* in pore space in contact with at least another ettringite */
/* or aluminate surface  */
        tries=0;
        while(fchr==0){
                tries+=1;
                newact=7;
                /* generate a random location in the 3-D system */
                xchr=(int)((float)SYSIZE*ran1(seed));
                ychr=(int)((float)SYSIZE*ran1(seed));
                zchr=(int)((float)SYSIZE*ran1(seed));
                if(xchr>=SYSIZE){xchr=0;}
                if(ychr>=SYSIZE){ychr=0;}
                if(zchr>=SYSIZE){zchr=0;}

                check=mic[xchr][ychr][zchr];
                /* if location is porosity, locate the ettringite there */
                if(check==POROSITY){
                				numsil=edgecnt(xchr,ychr,zchr,C3S,C2S,C3S);
                           numsil=26-numsil;
			if(etype==0){
	                        numnear=edgecnt(xchr,ychr,zchr,ETTR,C3A,C4AF);
			}
			else{
	                        numnear=edgecnt(xchr,ychr,zchr,ETTRC4AF,C3A,C4AF);
			}
                        /* be sure that at least one neighboring pixel */
                        /* is ettringite, or aluminate clinker */
                        if((tries>5000)||((numnear<26)&&(numsil<1))){
				if(etype==0){
	                                mic[xchr][ychr][zchr]=ETTR;
												count[ETTR]+=1;
				}
				else{
	                                mic[xchr][ychr][zchr]=ETTRC4AF;
												count[ETTRC4AF]+=1;
				}
											count[POROSITY]-=1;
                               	fchr=1;
                        }
                }
        }
        return(newact);
}

/* routine to add extra CH when gypsum, hemihydrate, anhydrite, CaCl2, or */
/* diffusing CAS2  reacts with C4AF */
/* Called by movegyp, movehem, moveanh, moveettr, movecas2, and movecacl2 */
/* Calls edgecnt */
void extch()
{
        int numnear,sump,xchr,ychr,zchr,fchr,i1,plok,check;
        long int tries;

        fchr=0;
        tries=0;
        /* locate CH at random location */
        /* in pore space in contact with at least another CH */
        while(fchr==0){
                tries+=1;
                /* generate a random location in the 3-D system */
                xchr=(int)((float)SYSIZE*ran1(seed));
                ychr=(int)((float)SYSIZE*ran1(seed));
                zchr=(int)((float)SYSIZE*ran1(seed));
                if(xchr>=SYSIZE){xchr=0;}
                if(ychr>=SYSIZE){ychr=0;}
                if(zchr>=SYSIZE){zchr=0;}
                check=mic[xchr][ychr][zchr];

                /* if location is porosity, locate the CH there */
                if(check==POROSITY){
                        numnear=edgecnt(xchr,ychr,zchr,CH,DIFFCH,CH);
                        /* be sure that at least one neighboring pixel */
                        /* is CH or diffusing CH */
                        if((numnear<26)||(tries>5000)){
                                mic[xchr][ychr][zchr]=CH;
				count[CH]+=1;
				count[POROSITY]-=1;
                               	fchr=1;
                        }
                }
        }
}

/* routine to add extra gypsum when hemihydrate or anhydrite hydrates */
/* Called by movehem and moveanh */
/* Calls moveone and edgecnt */
void extgyps(xpres,ypres,zpres)
        int xpres,ypres,zpres;
{
        int multf,numnear,sump,xchr,ychr,zchr,check,fchr,i1,plok,newact;
       	long int tries;

/* first try 6 neighboring locations until      */
/*	a) successful				*/
/*	b) all 6 sites are tried and full or    */
/*	c) 500 tries are made 			*/
        fchr=0;
       	sump=1;
       	for(i1=1;((i1<=500)&&(fchr==0)&&(sump!=30030));i1++){
		
                /* choose a neighbor at random */
                xchr=xpres;
                ychr=ypres;
                zchr=zpres;
                newact=0;
                multf=moveone(&xchr,&ychr,&zchr,&newact,sump);
                if(newact==0){printf("Error in value of newact in extfh3 \n");}
                check=mic[xchr][ychr][zchr];	 

               	/* if neighbor is porosity   */
                /* then locate the GYPSUMS there */
                if(check==POROSITY){
                        mic[xchr][ychr][zchr]=GYPSUMS;
			count[GYPSUMS]+=1;
			count[POROSITY]-=1;
                       	fchr=1;
                }
               	else{
                        sump*=multf;
                }
        }

/* if no neighbor available, locate GYPSUMS at random location */
/* in pore space in contact with at least one GYPSUMS */
        tries=0;
        while(fchr==0){
                tries+=1;
                /* generate a random location in the 3-D system */
                xchr=(int)((float)SYSIZE*ran1(seed));
                ychr=(int)((float)SYSIZE*ran1(seed));
                zchr=(int)((float)SYSIZE*ran1(seed));
                if(xchr>=SYSIZE){xchr=0;}
                if(ychr>=SYSIZE){ychr=0;}
                if(zchr>=SYSIZE){zchr=0;}
                check=mic[xchr][ychr][zchr];
                /* if location is porosity, locate the GYPSUMS there */
                if(check==POROSITY){
                        numnear=edgecnt(xchr,ychr,zchr,HEMIHYD,GYPSUMS,ANHYDRITE);
                        /* be sure that at least one neighboring pixel */
                        /* is Gypsum in some form */
                        if((numnear<26)||(tries>5000)){
                                mic[xchr][ychr][zchr]=GYPSUMS;
				count[GYPSUMS]+=1;
				count[POROSITY]-=1;
                               	fchr=1;
                        }
                }
        }
}

/* routine to move a diffusing ANHYDRITE species */
/* Inputs: current location (xcur,ycur,zcur) and flag indicating if final */
/* step in diffusion process */
/* Returns flag indicating action taken (reaction or diffusion/no movement) */
/* Calls moveone */
/* Called by hydrate */
int moveanh(xcur,ycur,zcur,finalstep,nucprgyp)
        int xcur,ycur,zcur,finalstep;
	float nucprgyp;
{
        int xnew,ynew,znew,plok,action,sumback,sumin,check;
	int nexp,iexp,xexp,yexp,zexp,newact,sumold,sumgarb,ettrtype;
	float pgen,pexp,pext,p2diff;

        action=0;
        /* first check for nucleation */
        pgen=ran1(seed);
	p2diff=ran1(seed);
        if((nucprgyp>=pgen)||(finalstep==1)){
                action=0;
                mic[xcur][ycur][zcur]=GYPSUMS;
                count[DIFFANH]-=1;
		count[GYPSUMS]+=1;
               	pexp=ran1(seed);
		if(pexp<0.4){
			extgyps(xcur,ycur,zcur);
		}
        }
	else{
	        /* Store current location of species */
       		 xnew=xcur;
       		 ynew=ycur;
       		 znew=zcur;
       		 sumin=1;
       		 sumback=moveone(&xnew,&ynew,&znew,&action,sumin);
	 
       		 if(action==0){printf("Error in value of action \n");}
       		 check=mic[xnew][ynew][znew];

/* if new location is solid GYPSUM(S) or diffusing GYPSUM, then convert */
/* diffusing ANHYDRITE species to solid GYPSUM */
       	if((check==GYPSUM)||(check==GYPSUMS)||(check==DIFFGYP)){
	                mic[xcur][ycur][zcur]=GYPSUMS;
        	        /* decrement count of diffusing ANHYDRITE species */
               		/* and increment count of solid GYPSUMS */
	                count[DIFFANH]-=1;
       		        count[GYPSUMS]+=1;
	                action=0;
			/* Add extra gypsum as necessary */
               		pexp=ran1(seed);
			if(pexp<0.4){
				extgyps(xnew,ynew,znew);
			}
       		}
        /* if new location is C3A or diffusing C3A, execute conversion */
        /* to ettringite (including necessary volumetric expansion) */
        else if(((check==C3A)&&(p2diff<SOLIDC3AGYP))||((check==DIFFC3A)&&(p2diff<C3AGYP))||((check==DIFFC4A)&&(p2diff<C3AGYP))){
        /* Convert diffusing gypsum to an ettringite pixel */
		ettrtype=0;
                mic[xcur][ycur][zcur]=ETTR;
		if(check==DIFFC4A){
			ettrtype=1;
                	mic[xcur][ycur][zcur]=ETTRC4AF;
		}
                action=0;
                count[DIFFANH]-=1;
		count[check]-=1;

                /* determine if C3A should be converted to ettringite */
                /* 1 unit of hemihydrate requires 0.569 units of C3A */
                /* and should form 4.6935 units of ettringite */
                pexp=ran1(seed);
                nexp=3;
                if(pexp<=0.569){
			if(ettrtype==0){
       		                mic[xnew][ynew][znew]=ETTR;
				count[ETTR]+=1;
			}
			else{
       		                mic[xnew][ynew][znew]=ETTRC4AF;
				count[ETTRC4AF]+=1;
			}
                        nexp=2;
                }
                else{
                        /* maybe someday, use a new FIXEDC3A here */
                        /* so it won't dissolve later */
                        if(check==C3A){
                                mic[xnew][ynew][znew]=C3A;
				count[C3A]+=1;
                        }
                        else{
				if(ettrtype==0){
	                                count[DIFFC3A]+=1;
       		                        mic[xnew][ynew][znew]=DIFFC3A;
				}
				else{
	                                count[DIFFC4A]+=1;
       		                        mic[xnew][ynew][znew]=DIFFC4A;
				}
                        }
                        nexp=3;
                }

/* create extra ettringite pixels to maintain volume stoichiometry */
/* xexp, yexp, and zexp hold coordinates of most recently added ettringite */
/* species as we attempt to grow a needle like structure */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extettr(xexp,yexp,zexp,ettrtype);
                        /* update xexp, yexp and zexp as needed */
                       	switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                        }
                }

                /* probabilistic-based expansion for last ettringite pixel */
                pexp=ran1(seed);
                if(pexp<=0.6935){
                        newact=extettr(xexp,yexp,zexp,ettrtype);
                }
        }

        /* if new location is C4AF execute conversion */
        /* to ettringite (including necessary volumetric expansion) */
        if((check==C4AF)&&(p2diff<SOLIDC4AFGYP)){
                mic[xcur][ycur][zcur]=ETTRC4AF;
		count[ETTRC4AF]+=1;
                count[DIFFANH]-=1;

                /* determine if C4AF should be converted to ettringite */
                /* 1 unit of gypsum requires 0.8174 units of C4AF */
                /* and should form 4.6935 units of ettringite */
                pexp=ran1(seed);
                nexp=3;
                if(pexp<=0.8174){
                        mic[xnew][ynew][znew]=ETTRC4AF;
			count[ETTRC4AF]+=1;
			count[C4AF]-=1;
                        nexp=2;
                        pext=ran1(seed);
                        /* Addition of extra CH */
                        if(pext<0.2584){
                                        extch();
                        }
                        pext=ran1(seed);
                        /* Addition of extra FH3 */
                        if(pext<0.5453){
                                        extfh3(xnew,ynew,znew);
                        }

                }
                else{
                        /* maybe someday, use a new FIXEDC4AF here */
                        /* so it won't dissolve later */
                        mic[xnew][ynew][znew]=C4AF;
                        nexp=3;
                }

/* create extra ettringite pixels to maintain volume stoichiometry */
/* xexp, yexp and zexp hold coordinates of most recently added ettringite */
/* species as we attempt to grow a needle like structure */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extettr(xexp,yexp,zexp,1);
                        /* update xexp, yexp and zexp as needed */
                        switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                       }
                }

                /* probabilistic-based expansion for last ettringite pixel */
                pexp=ran1(seed);
                if(pexp<=0.6935){
                        newact=extettr(xexp,yexp,zexp,1);
                }
                action=0;
        }
	}

        if(action!=0){
        /* if diffusion step is possible, perform it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                       	mic[xnew][ynew][znew]=DIFFANH;
               	}
                else{
                        /* indicate that diffusing ANHYDRITE species remained */
                     	/* at original location */
                       	action=7;
                }
        }
        return(action);
}

/* routine to move a diffusing HEMIHYDRATE species */
/* Inputs: current location (xcur,ycur,zcur) and flag indicating if final */
/* step in diffusion process */
/* Returns flag indicating action taken (reaction or diffusion/no movement) */
/* Calls moveone, extettr, extch, and extfh3 */
/* Called by hydrate */
int movehem(xcur,ycur,zcur,finalstep,nucprgyp)
        int xcur,ycur,zcur,finalstep;
	float nucprgyp;
{
        int xnew,ynew,znew,plok,action,sumback,sumin,check;
	int nexp,iexp,xexp,yexp,zexp,newact,sumold,sumgarb,ettrtype;
	float pgen,pexp,pext,p2diff;

        action=0;
        /* first check for nucleation */
        pgen=ran1(seed);
	p2diff=ran1(seed);
        if((nucprgyp>=pgen)||(finalstep==1)){
                action=0;
                mic[xcur][ycur][zcur]=GYPSUMS;
                count[DIFFHEM]-=1;
		count[GYPSUMS]+=1;
		/* Add extra gypsum as necessary */
               	pexp=ran1(seed);
		if(pexp<0.4){
			extgyps(xcur,ycur,zcur);
		}
        }
	else{
        /* Store current location of species */
	         xnew=xcur;
       		 ynew=ycur;
       		 znew=zcur;
       		 sumin=1;
       		 sumback=moveone(&xnew,&ynew,&znew,&action,sumin);

       		 if(action==0){printf("Error in value of action \n");}
       		 check=mic[xnew][ynew][znew];

/* if new location is solid GYPSUM(S) or diffusing GYPSUM, then convert */
/* diffusing HEMIHYDRATE species to solid GYPSUM */
        	if((check==GYPSUM)||(check==GYPSUMS)||(check==DIFFGYP)){
	                mic[xcur][ycur][zcur]=GYPSUMS;
       		        /* decrement count of diffusing HEMIHYDRATE species */
	                /* and increment count of solid GYPSUMS */
	                count[DIFFHEM]-=1;
       		        count[GYPSUMS]+=1;
       		        action=0;
			/* Add extra gypsum as necessary */
                	pexp=ran1(seed);
			if(pexp<0.4){
				extgyps(xnew,ynew,znew);
			}
	       	}
        /* if new location is C3A or diffusing C3A, execute conversion */
        /* to ettringite (including necessary volumetric expansion) */
        else if(((check==C3A)&&(p2diff<SOLIDC3AGYP))||((check==DIFFC3A)&&(p2diff<C3AGYP))||((check==DIFFC4A)&&(p2diff<C3AGYP))){
        /* Convert diffusing gypsum to an ettringite pixel */
		ettrtype=0;
                mic[xcur][ycur][zcur]=ETTR;
		if(check==DIFFC4A){
			ettrtype=1;
                	mic[xcur][ycur][zcur]=ETTRC4AF;
		}
                action=0;
                count[DIFFHEM]-=1;
		count[check]-=1;

                /* determine if C3A should be converted to ettringite */
                /* 1 unit of hemihydrate requires 0.5583 units of C3A */
                /* and should form 4.6053 units of ettringite */
                pexp=ran1(seed);
                nexp=3;
                if(pexp<=0.5583){
			if(ettrtype==0){
       		                mic[xnew][ynew][znew]=ETTR;
				count[ETTR]+=1;
			}
			else{
       		                mic[xnew][ynew][znew]=ETTRC4AF;
				count[ETTRC4AF]+=1;
			}
                        nexp=2;
                }
                else{
                        /* maybe someday, use a new FIXEDC3A here */
                        /* so it won't dissolve later */
                        if(check==C3A){
                                mic[xnew][ynew][znew]=C3A;
				count[C3A]+=1;
                        }
                        else{
				if(ettrtype==0){
	                                count[DIFFC3A]+=1;
       		                        mic[xnew][ynew][znew]=DIFFC3A;
				}
				else{
	                                count[DIFFC4A]+=1;
       		                        mic[xnew][ynew][znew]=DIFFC4A;
				}
                        }
                        nexp=3;
                }

/* create extra ettringite pixels to maintain volume stoichiometry */
/* xexp, yexp, and zexp hold coordinates of most recently added ettringite */
/* species as we attempt to grow a needle like structure */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extettr(xexp,yexp,zexp,ettrtype);
                        /* update xexp, yexp and zexp as needed */
                       	switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                        }
                }

                /* probabilistic-based expansion for last ettringite pixel */
                pexp=ran1(seed);
                if(pexp<=0.6053){
                        newact=extettr(xexp,yexp,zexp,ettrtype);
                }
        }
	
        /* if new location is C4AF execute conversion */
        /* to ettringite (including necessary volumetric expansion) */
        if((check==C4AF)&&(p2diff<SOLIDC4AFGYP)){
                mic[xcur][ycur][zcur]=ETTRC4AF;
		count[ETTRC4AF]+=1;
                count[DIFFHEM]-=1;

                /* determine if C4AF should be converted to ettringite */
                /* 1 unit of gypsum requires 0.802 units of C4AF */
                /* and should form 4.6053 units of ettringite */
                pexp=ran1(seed);
                nexp=3;
                if(pexp<=0.802){
                        mic[xnew][ynew][znew]=ETTRC4AF;
			count[ETTRC4AF]+=1;
			count[C4AF]-=1;
                        nexp=2;
                        pext=ran1(seed);
                        /* Addition of extra CH */
                        if(pext<0.2584){
                                        extch();
                        }
                        pext=ran1(seed);
                        /* Addition of extra FH3 */
                        if(pext<0.5453){
                                        extfh3(xnew,ynew,znew);
                        }

                }
                else{
                        /* maybe someday, use a new FIXEDC4AF here */
                        /* so it won't dissolve later */
                        mic[xnew][ynew][znew]=C4AF;
                        nexp=3;
                }

/* create extra ettringite pixels to maintain volume stoichiometry */
/* xexp, yexp and zexp hold coordinates of most recently added ettringite */
/* species as we attempt to grow a needle like structure */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extettr(xexp,yexp,zexp,1);
                        /* update xexp, yexp and zexp as needed */
                        switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                       }
                }

                /* probabilistic-based expansion for last ettringite pixel */
                pexp=ran1(seed);
                if(pexp<=0.6053){
                        newact=extettr(xexp,yexp,zexp,1);
                }
                action=0;
        }
	}

        if(action!=0){
        /* if diffusion step is possible, perform it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                       	mic[xnew][ynew][znew]=DIFFHEM;
               	}
                else{
                        /* indicate that diffusing HEMIHYDRATE species */
                     	/* remained at original location */
                       	action=7;
                }
        }
        return(action);
}

/* routine to add extra Freidel's salt when CaCl2 reacts with */
/* C3A or C4AF at location (xpres,ypres,zpres) */
/* Called by movecacl2 and movec3a */
/* Calls moveone and edgecnt */
int extfreidel(xpres,ypres,zpres)
        int xpres,ypres,zpres;
{
        int multf,numnear,sump,xchr,ychr,zchr,check,fchr,i1,plok,newact;
       	long int tries;

/* first try 6 neighboring locations until      */
/*	a) successful				*/
/*	b) all 6 sites are tried and full or    */
/*	c) 500 tries are made 			*/
        fchr=0;
       	sump=1;
       	for(i1=1;((i1<=500)&&(fchr==0)&&(sump!=30030));i1++){
		
                /* choose a neighbor at random */
                xchr=xpres;
                ychr=ypres;
                zchr=zpres;
                newact=0;
                multf=moveone(&xchr,&ychr,&zchr,&newact,sump);
                if(newact==0){printf("Error in value of newact in extfreidel \n");}
                check=mic[xchr][ychr][zchr];	 

               	/* if neighbor is porosity   */
                /* then locate the freidel's salt there */
                if(check==POROSITY){
                        mic[xchr][ychr][zchr]=FREIDEL;
			count[FREIDEL]+=1;
			count[POROSITY]-=1;
                       	fchr=1;
                }
               	else{
                        sump*=multf;
                }
        }

/* if no neighbor available, locate FREIDEL at random location */
/* in pore space in contact with at least one FREIDEL */
        tries=0;
        while(fchr==0){
                tries+=1;
		newact=7;
                /* generate a random location in the 3-D system */
                xchr=(int)((float)SYSIZE*ran1(seed));
                ychr=(int)((float)SYSIZE*ran1(seed));
                zchr=(int)((float)SYSIZE*ran1(seed));
                if(xchr>=SYSIZE){xchr=0;}
                if(ychr>=SYSIZE){ychr=0;}
                if(zchr>=SYSIZE){zchr=0;}
                check=mic[xchr][ychr][zchr];
                /* if location is porosity, locate the FREIDEL there */
                if(check==POROSITY){
                        numnear=edgecnt(xchr,ychr,zchr,FREIDEL,FREIDEL,DIFFCACL2);
                        /* be sure that at least one neighboring pixel */
                        /* is FREIDEL or diffusing CACL2 */
                        if((numnear<26)||(tries>5000)){
                                mic[xchr][ychr][zchr]=FREIDEL;
				count[FREIDEL]+=1;
				count[POROSITY]-=1;
                               	fchr=1;
                        }
                }
        }
	return(newact);
}

/* routine to add extra stratlingite when AS reacts with */
/* CH at location (xpres,ypres,zpres) */
/* or when diffusing CAS2 reacts with aluminates */
/* Called by moveas, movech, and movecas2 */
/* Calls moveone and edgecnt */
int extstrat(xpres,ypres,zpres)
        int xpres,ypres,zpres;
{
        int multf,numnear,sump,xchr,ychr,zchr,check,fchr,i1,plok,newact;
       	long int tries;

/* first try 6 neighboring locations until      */
/*	a) successful				*/
/*	b) all 6 sites are tried and full or    */
/*	c) 500 tries are made 			*/
        fchr=0;
       	sump=1;
       	for(i1=1;((i1<=500)&&(fchr==0)&&(sump!=30030));i1++){
		
                /* choose a neighbor at random */
                xchr=xpres;
                ychr=ypres;
                zchr=zpres;
                newact=0;
                multf=moveone(&xchr,&ychr,&zchr,&newact,sump);
                if(newact==0){printf("Error in value of newact in extstrat \n");}
                check=mic[xchr][ychr][zchr];	 

               	/* if neighbor is porosity   */
                /* then locate the stratlingite there */
                if(check==POROSITY){
                        mic[xchr][ychr][zchr]=STRAT;
			count[STRAT]+=1;
			count[POROSITY]-=1;
                       	fchr=1;
                }
               	else{
                        sump*=multf;
                }
        }

/* if no neighbor available, locate STRAT at random location */
/* in pore space in contact with at least one STRAT */
        tries=0;
        while(fchr==0){
                tries+=1;
		newact=7;
                /* generate a random location in the 3-D system */
                xchr=(int)((float)SYSIZE*ran1(seed));
                ychr=(int)((float)SYSIZE*ran1(seed));
                zchr=(int)((float)SYSIZE*ran1(seed));
                if(xchr>=SYSIZE){xchr=0;}
                if(ychr>=SYSIZE){ychr=0;}
                if(zchr>=SYSIZE){zchr=0;}
                check=mic[xchr][ychr][zchr];
                /* if location is porosity, locate the STRAT there */
                if(check==POROSITY){
                        numnear=edgecnt(xchr,ychr,zchr,STRAT,DIFFCAS2,DIFFAS);
                        /* be sure that at least one neighboring pixel */
                        /* is STRAT, diffusing CAS2, or diffusing AS */
                        if((numnear<26)||(tries>5000)){
                                mic[xchr][ychr][zchr]=STRAT;
				count[STRAT]+=1;
				count[POROSITY]-=1;
                               	fchr=1;
                        }
                }
        }
	return(newact);
}

/* routine to move a diffusing gypsum species */
/* from current location (xcur,ycur,zcur) */
/* Returns action flag indicating response taken */
/* Called by hydrate */
/* Calls moveone, extettr, extch, and extfh3 */
int movegyp(xcur,ycur,zcur,finalstep)
        int xcur,ycur,zcur,finalstep;
{
        int check,xnew,ynew,znew,plok,action,nexp,iexp;
       	int xexp,yexp,zexp,newact,sumold,sumgarb,ettrtype;
       	float pexp,pext,p2diff;

       	sumold=1;

/* First be sure that a diffusing gypsum species is located at xcur,ycur,zcur */
/* if not, return to calling routine */
        if(mic[xcur][ycur][zcur]!=DIFFGYP){
                action=0;
                return(action);
       	}

/* Determine new coordinates (periodic boundaries are used) */
        xnew=xcur;
        ynew=ycur;
        znew=zcur;
        action=0;
        sumgarb=moveone(&xnew,&ynew,&znew,&action,sumold);
        if(action==0){printf("Error in value of action in movegyp \n");}
        check=mic[xnew][ynew][znew];
	p2diff=ran1(seed);
        /* if new location is CSH, check for absorption of gypsum */
        if((check==CSH)&&((float)count[ABSGYP]<(gypabsprob*(float)count[CSH]))){
                pexp=ran1(seed);
                if(pexp<AGRATE){
                        /* update counts for absorbed and diffusing gypsum */
                        count[ABSGYP]+=1;
                        count[DIFFGYP]-=1;
                        mic[xcur][ycur][zcur]=ABSGYP;
                        action=0;
                }
        }

        /* if new location is C3A or diffusing C3A, execute conversion */
        /* to ettringite (including necessary volumetric expansion) */
	/* Use p2diff to try to favor formation of ettringite on */
	/* aluminate surfaces as opposed to in solution */
        else if(((check==C3A)&&(p2diff<SOLIDC3AGYP))||((check==DIFFC3A)&&(p2diff<C3AGYP))||((check==DIFFC4A)&&(p2diff<C3AGYP))){
        /* Convert diffusing gypsum to an ettringite pixel */
		ettrtype=0;
                mic[xcur][ycur][zcur]=ETTR;
		if(check==DIFFC4A){
			ettrtype=1;
                	mic[xcur][ycur][zcur]=ETTRC4AF;
		}
                action=0;
                count[DIFFGYP]-=1;
		count[check]-=1;

                /* determine if C3A should be converted to ettringite */
                /* 1 unit of gypsum requires 0.40 units of C3A */
                /* and should form 3.30 units of ettringite */
                pexp=ran1(seed);
                nexp=2;
                if(pexp<=0.40){
			if(ettrtype==0){
       		                mic[xnew][ynew][znew]=ETTR;
				count[ETTR]+=1;
			}
			else{
       		                mic[xnew][ynew][znew]=ETTRC4AF;
				count[ETTRC4AF]+=1;
			}
                        nexp=1;
                }
                else{
                        /* maybe someday, use a new FIXEDC3A here */
                        /* so it won't dissolve later */
                        if(check==C3A){
                                mic[xnew][ynew][znew]=C3A;
				count[C3A]+=1;
                        }
                        else{
				if(ettrtype==0){
	                                count[DIFFC3A]+=1;
       		                        mic[xnew][ynew][znew]=DIFFC3A;
				}
				else{
	                                count[DIFFC4A]+=1;
       		                        mic[xnew][ynew][znew]=DIFFC4A;
				}
                        }
                        nexp=2;
                }

/* create extra ettringite pixels to maintain volume stoichiometry */
/* xexp, yexp, and zexp hold coordinates of most recently added ettringite */
/* species as we attempt to grow a needle like structure */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extettr(xexp,yexp,zexp,ettrtype);
                        /* update xexp, yexp and zexp as needed */
                       	switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                        }
                }

                /* probabilistic-based expansion for last ettringite pixel */
                pexp=ran1(seed);
                if(pexp<=0.30){
                        newact=extettr(xexp,yexp,zexp,ettrtype);
                }
        }
	
        /* if new location is C4AF execute conversion */
        /* to ettringite (including necessary volumetric expansion) */
        if((check==C4AF)&&(p2diff<SOLIDC4AFGYP)){
                mic[xcur][ycur][zcur]=ETTRC4AF;
		count[ETTRC4AF]+=1;
                count[DIFFGYP]-=1;

                /* determine if C4AF should be converted to ettringite */
                /* 1 unit of gypsum requires 0.575 units of C4AF */
                /* and should form 3.30 units of ettringite */
                pexp=ran1(seed);
                nexp=2;
                if(pexp<=0.575){
                        mic[xnew][ynew][znew]=ETTRC4AF;
			count[ETTRC4AF]+=1;
			count[C4AF]-=1;
                        nexp=1;
                        pext=ran1(seed);
                        /* Addition of extra CH */
                        if(pext<0.2584){
                                        extch();
                        }
                        pext=ran1(seed);
                        /* Addition of extra FH3 */
                        if(pext<0.5453){
                                        extfh3(xnew,ynew,znew);
                        }

                }
                else{
                        /* maybe someday, use a new FIXEDC4AF here */
                        /* so it won't dissolve later */
                        mic[xnew][ynew][znew]=C4AF;
                        nexp=2;
                }

/* create extra ettringite pixels to maintain volume stoichiometry */
/* xexp, yexp and zexp hold coordinates of most recently added ettringite */
/* species as we attempt to grow a needle like structure */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extettr(xexp,yexp,zexp,1);
                        /* update xexp, yexp and zexp as needed */
                        switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                       }
                }

                /* probabilistic-based expansion for last ettringite pixel */
                pexp=ran1(seed);
                if(pexp<=0.30){
                        newact=extettr(xexp,yexp,zexp,1);
                }
                action=0;
        }

        /* if last diffusion step and no reaction, convert back to */
        /* primary solid gypsum */
        if((action!=0)&&(finalstep==1)){
                action=0;
                count[DIFFGYP]-=1;
		count[GYPSUM]+=1;
                mic[xcur][ycur][zcur]=GYPSUM;
        }

        if(action!=0){
                /* if diffusion is possible, execute it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                        mic[xnew][ynew][znew]=DIFFGYP;
                }
                else{
                        /* indicate that diffusing gypsum remained at */
                        /* original location */
                        action=7;
                }
        }
       	return(action);
}

/* routine to move a diffusing CaCl2 species */
/* from current location (xcur,ycur,zcur) */
/* Returns action flag indicating response taken */
/* Called by hydrate */
/* Calls moveone, extfreidel, extch, and extfh3 */
int movecacl2(xcur,ycur,zcur,finalstep)
        int xcur,ycur,zcur,finalstep;
{
        int check,xnew,ynew,znew,plok,action,nexp,iexp;
       	int xexp,yexp,zexp,newact,sumold,sumgarb,keep;
       	float pexp,pext;

       	sumold=1;
	keep=0;

/* First be sure that a diffusing CaCl2 species is located at xcur,ycur,zcur */
/* if not, return to calling routine */
        if(mic[xcur][ycur][zcur]!=DIFFCACL2){
                action=0;
                return(action);
       	}

/* Determine new coordinates (periodic boundaries are used) */
        xnew=xcur;
        ynew=ycur;
        znew=zcur;
        action=0;
        sumgarb=moveone(&xnew,&ynew,&znew,&action,sumold);
        if(action==0){printf("Error in value of action in movecacl2 \n");}
        check=mic[xnew][ynew][znew];

        /* if new location is C3A or diffusing C3A, execute conversion */
        /* to freidel's salt (including necessary volumetric expansion) */
        if((check==C3A)||(check==DIFFC3A)||(check==DIFFC4A)){
        /* Convert diffusing C3A or C3A to a freidel's salt pixel */
                action=0;
                mic[xnew][ynew][znew]=FREIDEL;
		count[FREIDEL]+=1;
                count[check]-=1;

                /* determine if diffusing CaCl2 should be converted to FREIDEL */
                /* 0.5793 unit of CaCl2 requires 1 unit of C3A */
                /* and should form 3.3295 units of FREIDEL */
                pexp=ran1(seed);
                nexp=2;
                if(pexp<=0.5793){
                        mic[xcur][ycur][zcur]=FREIDEL;
			count[FREIDEL]+=1;
			count[DIFFCACL2]-=1;
                        nexp=1;
                }
                else{
			keep=1;
                        nexp=2;
                }

/* create extra Freidel's salt pixels to maintain volume stoichiometry */
/* xexp, yexp, and zexp hold coordinates of most recently added FREIDEL */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extfreidel(xexp,yexp,zexp);
                        /* update xexp, yexp and zexp as needed */
                       	switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                        }
                }

                /* probabilistic-based expansion for last FREIDEL pixel */
                pexp=ran1(seed);
                if(pexp<=0.3295){
                        newact=extfreidel(xexp,yexp,zexp);
                }
        }
	
        /* if new location is C4AF execute conversion */
        /* to freidel's salt (including necessary volumetric expansion) */
        else if(check==C4AF){
                mic[xnew][ynew][znew]=FREIDEL;
		count[FREIDEL]+=1;
                count[C4AF]-=1;

                /* determine if CACL2 should be converted to FREIDEL */
                /* 0.4033 unit of CaCl2 requires 1 unit of C4AF */
                /* and should form 2.3176 units of FREIDEL */
		/* Also 0.6412 units of CH and 1.3522 units of FH3 */
		/* per unit of CACL2 */
                pexp=ran1(seed);
                nexp=1;
                if(pexp<=0.4033){
                        mic[xcur][ycur][zcur]=FREIDEL;
			count[FREIDEL]+=1;
			count[DIFFCACL2]-=1;
                        nexp=0;
                        pext=ran1(seed);
                        /* Addition of extra CH */
                        if(pext<0.6412){
                                        extch();
                        }
                        pext=ran1(seed);
                        /* Addition of extra FH3 */
                        if(pext<0.3522){
                                        extfh3(xnew,ynew,znew);
                        }
			extfh3(xnew,ynew,znew);

                }
                else{
                        nexp=1;
			keep=1;
                }

/* create extra freidel's salt pixels to maintain volume stoichiometry */
/* xexp, yexp and zexp hold coordinates of most recently added FREIDEL */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extfreidel(xexp,yexp,zexp);
                        /* update xexp, yexp and zexp as needed */
                        switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                       }
                }

                /* probabilistic-based expansion for last FREIDEL pixel */
                pexp=ran1(seed);
                if(pexp<=0.3176){
                        newact=extfreidel(xexp,yexp,zexp);
                }
                action=0;
        }

        /* if last diffusion step and no reaction, convert back to */
        /* solid CaCl2 */
        if((action!=0)&&(finalstep==1)){
                action=0;
                count[DIFFCACL2]-=1;
		count[CACL2]+=1;
                mic[xcur][ycur][zcur]=CACL2;
        }

        if(action!=0){
                /* if diffusion is possible, execute it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                        mic[xnew][ynew][znew]=DIFFCACL2;
                }
                else{
                        /* indicate that diffusing CACL2 remained at */
                        /* original location */
                        action=7;
                }
        }
	if(keep==1){action=7;}
       	return(action);
}

/* routine to move a diffusing CAS2 species */
/* from current location (xcur,ycur,zcur) */
/* Returns action flag indicating response taken */
/* Called by hydrate */
/* Calls moveone, extstrat, extch, and extfh3 */
int movecas2(xcur,ycur,zcur,finalstep)
        int xcur,ycur,zcur,finalstep;
{
        int check,xnew,ynew,znew,plok,action,nexp,iexp;
       	int xexp,yexp,zexp,newact,sumold,sumgarb,keep;
       	float pexp,pext;

       	sumold=1;
	keep=0;

/* First be sure that a diffusing CAS2 species is located at xcur,ycur,zcur */
/* if not, return to calling routine */
        if(mic[xcur][ycur][zcur]!=DIFFCAS2){
                action=0;
                return(action);
       	}

/* Determine new coordinates (periodic boundaries are used) */
        xnew=xcur;
        ynew=ycur;
        znew=zcur;
        action=0;
        sumgarb=moveone(&xnew,&ynew,&znew,&action,sumold);
        if(action==0){printf("Error in value of action in movecas2 \n");}
        check=mic[xnew][ynew][znew];

        /* if new location is C3A or diffusing C3A, execute conversion */
        /* to stratlingite (including necessary volumetric expansion) */
        if((check==C3A)||(check==DIFFC3A)||(check==DIFFC4A)){
        /* Convert diffusing CAS2 to a stratlingite pixel */
                action=0;
                mic[xcur][ycur][zcur]=STRAT;
		count[STRAT]+=1;
                count[DIFFCAS2]-=1;

                /* determine if diffusing or solid C3A should be converted to STRAT*/
                /* 1 unit of CAS2 requires 0.886 units of C3A */
                /* and should form 4.286 units of STRAT */
                pexp=ran1(seed);
                nexp=3;
                if(pexp<=0.886){
                        mic[xnew][ynew][znew]=STRAT;
			count[STRAT]+=1;
			count[check]-=1;
                        nexp=2;
                }

/* create extra stratlingite pixels to maintain volume stoichiometry */
/* xexp, yexp, and zexp hold coordinates of most recently added STRAT */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extstrat(xexp,yexp,zexp);
                        /* update xexp, yexp and zexp as needed */
                       	switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                        }
                }

                /* probabilistic-based expansion for last STRAT pixel */
                pexp=ran1(seed);
                if(pexp<=0.286){
                        newact=extstrat(xexp,yexp,zexp);
                }
        }
	
        /* if new location is C4AF execute conversion */
        /* to stratlingite (including necessary volumetric expansion) */
        else if(check==C4AF){
                mic[xnew][ynew][znew]=STRAT;
		count[STRAT]+=1;
                count[C4AF]-=1;

                /* determine if CAS2 should be converted to STRAT */
                /* 0.786 units of CAS2 requires 1 unit of C4AF */
                /* and should form 3.37 units of STRAT */
		/* Also 0.2586 units of CH and 0.5453 units of FH3 */
		/* per unit of C4AF */
                pexp=ran1(seed);
                nexp=2;
                if(pexp<=0.786){
                        mic[xcur][ycur][zcur]=STRAT;
			count[STRAT]+=1;
			count[DIFFCAS2]-=1;
                        nexp=1;
                        pext=ran1(seed);
                        /* Addition of extra CH */
                        /* 0.329= 0.2586/0.786 */
                        if(pext<0.329){
                                        extch();
                        }
                        pext=ran1(seed);
                        /* Addition of extra FH3 */
                        /* 0.6938= 0.5453/0.786 */
                        if(pext<0.6938){
                                        extfh3(xnew,ynew,znew);
                        }

                }
                else{
                        nexp=2;
			keep=1;
                }

/* create extra stratlingite pixels to maintain volume stoichiometry */
/* xexp, yexp and zexp hold coordinates of most recently added STRAT */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extstrat(xexp,yexp,zexp);
                        /* update xexp, yexp and zexp as needed */
                        switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                       }
                }

                /* probabilistic-based expansion for last STRAT pixel */
                pexp=ran1(seed);
                if(pexp<=0.37){
                        newact=extstrat(xexp,yexp,zexp);
                }
                action=0;
        }

        /* if last diffusion step and no reaction, convert back to */
        /* solid CAS2 */
        if((action!=0)&&(finalstep==1)){
                action=0;
                count[DIFFCAS2]-=1;
		count[CAS2]+=1;
                mic[xcur][ycur][zcur]=CAS2;
        }

        if(action!=0){
                /* if diffusion is possible, execute it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                        mic[xnew][ynew][znew]=DIFFCAS2;
                }
                else{
                        /* indicate that diffusing CAS2 remained at */
                        /* original location */
                        action=7;
                }
        }
	if(keep==1){action=7;}
       	return(action);
}

/* routine to move a diffusing AS species */
/* from current location (xcur,ycur,zcur) */
/* Returns action flag indicating response taken */
/* Called by hydrate */
/* Calls moveone, extstrat */
int moveas(xcur,ycur,zcur,finalstep)
        int xcur,ycur,zcur,finalstep;
{
        int check,xnew,ynew,znew,plok,action,nexp,iexp;
       	int xexp,yexp,zexp,newact,sumold,sumgarb,keep;
       	float pexp,pext;

       	sumold=1;
	keep=0;

/* First be sure that a diffusing AS species is located at xcur,ycur,zcur */
/* if not, return to calling routine */
        if(mic[xcur][ycur][zcur]!=DIFFAS){
                action=0;
                return(action);
       	}

/* Determine new coordinates (periodic boundaries are used) */
        xnew=xcur;
        ynew=ycur;
        znew=zcur;
        action=0;
        sumgarb=moveone(&xnew,&ynew,&znew,&action,sumold);
        if(action==0){printf("Error in value of action in moveas \n");}
        check=mic[xnew][ynew][znew];

        /* if new location is CH or diffusing CH, execute conversion */
        /* to stratlingite (including necessary volumetric expansion) */
        if((check==CH)||(check==DIFFCH)){
        /* Convert diffusing CH or CH to a stratlingite pixel */
                action=0;
                mic[xnew][ynew][znew]=STRAT;
		count[STRAT]+=1;
                count[check]-=1;

                /* determine if diffusing AS should be converted to STRAT */
                /* 0.7538 unit of AS requires 1 unit of CH */
                /* and should form 3.26 units of STRAT */
                pexp=ran1(seed);
                nexp=2;
                if(pexp<=0.7538){
                        mic[xcur][ycur][zcur]=STRAT;
			count[STRAT]+=1;
			count[DIFFAS]-=1;
                        nexp=1;
                }
                else{
			keep=1;
                        nexp=2;
                }

/* create extra stratlingite pixels to maintain volume stoichiometry */
/* xexp, yexp, and zexp hold coordinates of most recently added STRAT */
                xexp=xcur;
                yexp=ycur;
                zexp=zcur;
                for(iexp=1;iexp<=nexp;iexp++){
                        newact=extstrat(xexp,yexp,zexp);
                        /* update xexp, yexp and zexp as needed */
                       	switch (newact){
                                case 1:
                                        xexp-=1;
                                        if(xexp<0){xexp=(SYSIZEM1);}
                                       	break;
                                case 2:
                                        xexp+=1;
                                        if(xexp>=SYSIZE){xexp=0;}
                                        break;
                                case 3:
                                        yexp-=1;
                                        if(yexp<0){yexp=(SYSIZEM1);}
                                        break;
                                case 4:
                                        yexp+=1;
                                        if(yexp>=SYSIZE){yexp=0;}
                                        break;
                                case 5:
                                        zexp-=1;
                                        if(zexp<0){zexp=(SYSIZEM1);}
                                        break;
                                case 6:
                                        zexp+=1;
                                        if(zexp>=SYSIZE){zexp=0;}
                                        break;
                                default:
                                        break;
                        }
                }

                /* probabilistic-based expansion for last stratlingite pixel */
                pexp=ran1(seed);
                if(pexp<=0.326){
                        newact=extstrat(xexp,yexp,zexp);
                }
        }
	
        /* if last diffusion step and no reaction, convert back to */
        /* solid ASG */
        if((action!=0)&&(finalstep==1)){
                action=0;
                count[DIFFAS]-=1;
		count[ASG]+=1;
                mic[xcur][ycur][zcur]=ASG;
        }

        if(action!=0){
                /* if diffusion is possible, execute it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                        mic[xnew][ynew][znew]=DIFFAS;
                }
                else{
                        /* indicate that diffusing AS remained at */
                        /* original location */
                        action=7;
                }
        }
	if(keep==1){action=7;}
       	return(action);
}

/* routine to move a diffusing CACO3 species */
/* from current location (xcur,ycur,zcur) */
/* Returns action flag indicating response taken */
/* Called by hydrate */
/* Calls moveone, extettr */
int movecaco3(xcur,ycur,zcur,finalstep)
        int xcur,ycur,zcur,finalstep;
{
        int check,xnew,ynew,znew,plok,action,nexp,iexp;
       	int xexp,yexp,zexp,newact,sumold,sumgarb,keep;
       	float pexp,pext;

       	sumold=1;
	keep=0;

/* First be sure that a diffusing CACO3 species is located at xcur,ycur,zcur */
/* if not, return to calling routine */
        if(mic[xcur][ycur][zcur]!=DIFFCACO3){
                action=0;
                return(action);
       	}

/* Determine new coordinates (periodic boundaries are used) */
        xnew=xcur;
        ynew=ycur;
        znew=zcur;
        action=0;
        sumgarb=moveone(&xnew,&ynew,&znew,&action,sumold);
        if(action==0){printf("Error in value of action in moveas \n");}
        check=mic[xnew][ynew][znew];

        /* if new location is AFM execute conversion */
        /* to carboaluminate and ettringite (including necessary */
        /* volumetric expansion) */
        if(check==AFM){
        /* Convert AFM to a carboaluminate or ettringite pixel */
                action=0;
                pexp=ran1(seed);
                if(pexp<=0.479192){
                      mic[xnew][ynew][znew]=AFMC;
		      count[AFMC]+=1;
                }
                else{
                      mic[xnew][ynew][znew]=ETTR;
		      count[ETTR]+=1;
                }
                count[check]-=1;

                /* determine if diffusing CACO3 should be converted to AFMC */
                /* 0.078658 unit of AS requires 1 unit of AFM */
                /* and should form 0.55785 units of AFMC */
                pexp=ran1(seed);
                if(pexp<=0.078658){
                        mic[xcur][ycur][zcur]=AFMC;
			count[AFMC]+=1;
			count[DIFFCACO3]-=1;
                }
                else{
			keep=1;
                }

/* create extra ettringite pixels to maintain volume stoichiometry */
/* xexp, yexp, and zexp hold coordinates of most recently added ETTR */
                xexp=xnew;
                yexp=ynew;
                zexp=znew;

                /* probabilistic-based expansion for new ettringite pixel */
                pexp=ran1(seed);
                if(pexp<=0.26194){
                        newact=extettr(xexp,yexp,zexp,0);
                }
        }
	
        /* if last diffusion step and no reaction, convert back to */
        /* solid CACO3 */
        if((action!=0)&&(finalstep==1)){
                action=0;
                count[DIFFCACO3]-=1;
		count[CACO3]+=1;
                mic[xcur][ycur][zcur]=CACO3;
        }

        if(action!=0){
                /* if diffusion is possible, execute it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                        mic[xnew][ynew][znew]=DIFFCACO3;
                }
                else{
                        /* indicate that diffusing CACO3 remained at */
                        /* original location */
                        action=7;
                }
        }
	if(keep==1){action=7;}
       	return(action);
}

/* routine to add extra AFm phase when diffusing ettringite reacts */
/* with C3A (diffusing or solid) at location (xpres,ypres,zpres) */
/* Called by moveettr and movec3a */
/* Calls moveone and edgecnt */
void extafm(xpres,ypres,zpres)
        int xpres,ypres,zpres;
{
        int check,sump,xchr,ychr,zchr,fchr,i1,plok,newact,numnear;
       	long int tries;

/* first try 6 neighboring locations until      */
/*	a) successful				*/
/*	b) all 6 sites are tried or             */
/*	c) 100 tries are made 			*/
        fchr=0;
       	sump=1;
        for(i1=1;((i1<=100)&&(fchr==0)&&(sump!=30030));i1++){
		
                /* determine location of neighbor (using periodic boundaries) */
                xchr=xpres;
                ychr=ypres;
                zchr=zpres;
                newact=0;
                sump*=moveone(&xchr,&ychr,&zchr,&newact,sump);
                if(newact==0){printf("Error in value of newact in extafm \n");}
                check=mic[xchr][ychr][zchr];

                /* if neighbor is porosity, locate the AFm phase there */
                if(check==POROSITY){
                        mic[xchr][ychr][zchr]=AFM;
			count[AFM]+=1;
			count[POROSITY]-=1;
                        fchr=1;
                }
         }

         /* if no neighbor available, locate AFm phase at random location */
         /* in pore space */
         tries=0;
         while(fchr==0){
                tries+=1;
                /* generate a random location in the 3-D system */
                xchr=(int)((float)SYSIZE*ran1(seed));
                ychr=(int)((float)SYSIZE*ran1(seed));
                zchr=(int)((float)SYSIZE*ran1(seed));
                if(xchr>=SYSIZE){xchr=0;}
                if(ychr>=SYSIZE){ychr=0;}
                if(zchr>=SYSIZE){zchr=0;}
                check=mic[xchr][ychr][zchr];

                /* if location is porosity, locate the extra AFm there */
                if(check==POROSITY){
                        numnear=edgecnt(xchr,ychr,zchr,AFM,C3A,C4AF);
                        /* Be sure that at least one neighboring pixel is */
                        /* Afm phase, C3A, or C4AF */
                        if((tries>5000)||(numnear<26)){
                                mic[xchr][ychr][zchr]=AFM;
				count[AFM]+=1;
				count[POROSITY]-=1;
                                fchr=1;
                        }
                }
        }
}

/* routine to move a diffusing ettringite species */
/* currently located at (xcur,ycur,zcur) */
/* Called by hydrate */
/* Calls moveone, extch, extfh3, and extafm */
int moveettr(xcur,ycur,zcur,finalstep)
        int xcur,ycur,zcur,finalstep;
{
        int check,xnew,ynew,znew,plok,action,nexp,iexp;
        int xexp,yexp,zexp,newact,sumold,sumgarb;
        float pexp,pafm,pgrow;

/* First be sure a diffusing ettringite species is located at xcur,ycur,zcur */
/* if not, return to calling routine */
        if(mic[xcur][ycur][zcur]!=DIFFETTR){
                action=0;
                return(action);
        }

/* Determine new coordinates (periodic boundaries are used) */
        xnew=xcur;
        ynew=ycur;
        znew=zcur;
        action=0;
        sumold=1;
        sumgarb=moveone(&xnew,&ynew,&znew,&action,sumold);
        if(action==0){printf("Error in value of action in moveettr \n");}
        check=mic[xnew][ynew][znew];

        /* if new location is C4AF, execute conversion */
        /* to AFM phase (including necessary volumetric expansion) */
        if(check==C4AF){
                /* Convert diffusing ettringite to AFM phase */
                mic[xcur][ycur][zcur]=AFM;
		count[AFM]+=1;
                count[DIFFETTR]-=1;

                /* determine if C4AF should be converted to Afm */
                /* or FH3- 1 unit of ettringite requires 0.348 units */
                /* of C4AF to form 1.278 units of Afm, */
                /* 0.0901 units of CH and 0.1899 units of FH3 */
                pexp=ran1(seed);
		
                if(pexp<=0.278){
                        mic[xnew][ynew][znew]=AFM;
			count[AFM]+=1;
			count[C4AF]-=1;
                        pafm=ran1(seed);
                        /* 0.3241= 0.0901/0.278 */
                        if(pafm<=0.3241){
                                extch();
                        }
                        pafm=ran1(seed);
                        /* 0.4313= ((.1899-(.348-.278))/.278)   */
                        if(pafm<=0.4313){
                                extfh3(xnew,ynew,znew);
                        }
                }
                else if (pexp<=0.348){
                        mic[xnew][ynew][znew]=FH3;
			count[FH3]+=1;
			count[C4AF]-=1;
                }
                action=0;
        }

        /* if new location is C3A or diffusing C3A, execute conversion */
        /* to AFM phase (including necessary volumetric expansion) */
        else if((check==C3A)||(check==DIFFC3A)){
                /* Convert diffusing ettringite to AFM phase */
                action=0;
                mic[xcur][ycur][zcur]=AFM;
                count[DIFFETTR]-=1;
		count[AFM]+=1;
		count[check]-=1;	

                /* determine if C3A should be converted to AFm */
                /* 1 unit of ettringite requires 0.2424 units of C3A */
                /* and should form 1.278 units of AFm phase */
                pexp=ran1(seed);
                if(pexp<=0.2424){
                        mic[xnew][ynew][znew]=AFM;
			count[AFM]+=1;
                        pafm=(-0.1);
                }
                else{
                        /* maybe someday, use a new FIXEDC3A here */
                        /* so it won't dissolve later */
                        if(check==C3A){
                                mic[xnew][ynew][znew]=C3A;
				count[C3A]+=1;
                        }
                        else{
                                count[DIFFC3A]+=1;
                                mic[xnew][ynew][znew]=DIFFC3A;
                        }
/*                      pafm=(0.278-0.2424)/(1.0-0.2424);  */
			pafm=0.04699;
                }

                /* probabilistic-based expansion for new AFm phase pixel */
                pexp=ran1(seed);
                if(pexp<=pafm){
                        extafm(xcur,ycur,zcur);
                }
        }

        /* Check for conversion back to solid ettringite */
        else if(check==ETTR){
                pgrow=ran1(seed);
                if(pgrow<=ETTRGROW){
                        mic[xcur] [ycur] [zcur]=ETTR;
			count[ETTR]+=1;
                        action=0;
                        count[DIFFETTR]-=1;
                }
        }

        /* if last diffusion step and no reaction, convert back to */
        /* solid ettringite */
        if((action!=0)&&(finalstep==1)){
                action=0;
                count[DIFFETTR]-=1;
		count[ETTR]+=1;
                mic[xcur][ycur][zcur]=ETTR;
        }

        if(action!=0){
                /* if diffusion is possible, execute it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                        mic[xnew][ynew][znew]=DIFFETTR;
                }
                else{
                        /* indicate that diffusing ettringite remained at */
                        /* original location */
                        action=7;
                }
        }
        return(action);
}

/* routine to add extra pozzolanic CSH when CH reacts at */
/* pozzolanic surface (e.g. silica fume) located at (xpres,ypres,zpres) */
/* Called by movech */
/* Calls moveone and edgecnt */
void extpozz(xpres,ypres,zpres)
        int xpres,ypres,zpres;
{
        int check,sump,xchr,ychr,zchr,fchr,i1,plok,action,numnear;
        long int tries;

/* first try 6 neighboring locations until      */
/*	a) successful				*/
/*	b) all 6 sites are tried or             */
/*	c) 100 tries are made 			*/
        fchr=0;
        sump=1;
        for(i1=1;((i1<=100)&&(fchr==0)&&(sump!=30030));i1++){
		
                /* determine location of neighbor (using periodic boundaries) */
                xchr=xpres;
                ychr=ypres;
                zchr=zpres;
                action=0;
                sump*=moveone(&xchr,&ychr,&zchr,&action,sump);
                if(action==0){printf("Error in value of action in extpozz \n");}
                check=mic[xchr][ychr][zchr];

                /* if neighbor is porosity, locate the pozzolanic CSH there */
                if(check==POROSITY){
                        mic[xchr][ychr][zchr]=POZZCSH;
			count[POZZCSH]+=1;
			count[POROSITY]-=1;
                        fchr=1;
                }
        }

        /* if no neighbor available, locate pozzolanic CSH at random location */
        /* in pore space */
        tries=0;
        while(fchr==0){
                tries+=1;
                /* generate a random location in the 3-D system */
                xchr=(int)((float)SYSIZE*ran1(seed));
                ychr=(int)((float)SYSIZE*ran1(seed));
                zchr=(int)((float)SYSIZE*ran1(seed));
                if(xchr>=SYSIZE){xchr=0;}
                if(ychr>=SYSIZE){ychr=0;}
                if(zchr>=SYSIZE){zchr=0;}
                check=mic[xchr][ychr][zchr];
           /* if location is porosity, locate the extra pozzolanic CSH there */
                if(check==POROSITY){
                        numnear=edgecnt(xchr,ychr,zchr,POZZ,CSH,POZZCSH);
                        /* Be sure that one neighboring species is CSH or */
                        /* pozzolanic material */
                        if((tries>5000)||(numnear<26)){
                                mic[xchr][ychr][zchr]=POZZCSH;
				count[POZZCSH]+=1;
				count[POROSITY]-=1;
                                fchr=1;
                        }
                }
        }
}

/* routine to move a diffusing FH3 species */
/* from location (xcur,ycur,zcur) with nucleation probability nucprob */
/* Called by hydrate */
/* Calls moveone */
int movefh3(xcur,ycur,zcur,finalstep,nucprob)
        int xcur,ycur,zcur,finalstep;
        float nucprob;
{
        int check,xnew,ynew,znew,plok,action,sumold,sumgarb;
        float pgen;

        /* first check for nucleation */
        pgen=ran1(seed);

        if((nucprob>=pgen)||(finalstep==1)){
                action=0;
                mic[xcur][ycur][zcur]=FH3;
		count[FH3]+=1;
                count[DIFFFH3]-=1;
        }
        else{
		
                /* determine new location (using periodic boundaries) */
                xnew=xcur;
                ynew=ycur;
                znew=zcur;
                action=0;
                sumold=1;
                sumgarb=moveone(&xnew,&ynew,&znew,&action,sumold);
                if(action==0){printf("Error in value of action in movefh3 \n");}
                check=mic[xnew][ynew][znew];

               	/* check for growth of FH3 crystal */
                if(check==FH3){
                        mic[xcur][ycur][zcur]=FH3;
			count[FH3]+=1;
                        count[DIFFFH3]-=1;
                        action=0;
                }

                if(action!=0){
                        /* if diffusion is possible, execute it */
                        if(check==POROSITY){
                                mic[xcur][ycur][zcur]=POROSITY;
                                mic[xnew][ynew][znew]=DIFFFH3;
                        }
                        else{
                                /* indicate that diffusing FH3 species */
                                /* remained at original location */
                                action=7;
                        }
                }
        }
        return(action);
}

/* routine to move a diffusing CH species */
/* from location (xcur,ycur,zcur) with nucleation probability nucprob */
/* Called by hydrate */
/* Calls moveone and extpozz */ 
int movech(xcur,ycur,zcur,finalstep,nucprob)
        int xcur,ycur,zcur,finalstep;
        float nucprob;
{
        int check,xnew,ynew,znew,plok,action,sumgarb,sumold;
        float pexp,pgen,pfix;

        /* first check for nucleation */
        pgen=ran1(seed);
        if((nucprob>=pgen)||(finalstep==1)){
                action=0;
                mic[xcur][ycur][zcur]=CH;
                count[DIFFCH]-=1;
		count[CH]+=1;
        }
        else{
		
                /* determine new location (using periodic boundaries) */
                xnew=xcur;
                ynew=ycur;
                znew=zcur;
                action=0;
                sumold=1;
                sumgarb=moveone(&xnew,&ynew,&znew,&action,sumold);
                if(action==0){printf("Error in value of action in movech \n");}
                check=mic[xnew][ynew][znew];

                /* check for growth of CH crystal */
                if((check==CH)&&(pgen<=CHGROW)){
                        mic[xcur][ycur][zcur]=CH;
                        count[DIFFCH]-=1;
			count[CH]+=1;
                        action=0;
                }
              /* check for growth of CH crystal on aggregate or CaCO3 surface */
                /* re suggestion of Sidney Diamond */
                else if(((check==INERTAGG)||(check==CACO3)||(check==INERT))&&(pgen<=CHGROWAGG)&&(chflag==1)){
                        mic[xcur][ycur][zcur]=CH;
                        count[DIFFCH]-=1;
			count[CH]+=1;
                        action=0;
                }

                /* check for pozzolanic reaction */
		/* 36.41 units CH can react with 27 units of S */
                else if((pgen<=ppozz)&&(check==POZZ)&&(npr<=(int)((float)nfill*1.35))){
                        action=0;
                        mic[xcur][ycur][zcur]=POZZCSH;
			count[POZZCSH]+=1;
                        /* update counter of number of diffusing CH */
                        /* which have reacted pozzolanically */
                        npr+=1;
                        count[DIFFCH]-=1;
                        /* Convert pozzolan to pozzolanic CSH as needed */
                        pfix=ran1(seed);
			if(pfix<=(1./1.35)){
				mic[xnew][ynew][znew]=POZZCSH;
				count[POZZ]-=1;
				count[POZZCSH]+=1;
			}
                        /* allow for extra pozzolanic CSH as needed */
                        pexp=ran1(seed);
			/* should form 101.81 units of pozzolanic CSH for */
			/* each 36.41 units of CH and 27 units of S */
			/* 1.05466=(101.81-36.41-27)/36.41 */
			extpozz(xcur,ycur,zcur);
                        if(pexp<=0.05466){
                                extpozz(xcur,ycur,zcur);
                        }
                }
		else if(check==DIFFAS){
			action=0;
			mic[xcur][ycur][zcur]=STRAT;
			count[STRAT]+=1;
			/* update counter of number of diffusing CH */
			/* which have reacted to form stratlingite */
			nasr+=1;
			count[DIFFCH]-=1;
			/* Convert DIFFAS to STRAT as needed */
			pfix=ran1(seed);
			if(pfix<=0.7538){
				mic[xnew][ynew][znew]=STRAT;
				count[STRAT]+=1;
				count[DIFFAS]-=1;
			}
			/* allow for extra stratlingite as needed */
			/* 1.5035=(215.63-66.2-49.9)/66.2 */
			extstrat(xcur,ycur,zcur);
			pexp=ran1(seed);
			if(pexp<=0.5035){
				extstrat(xcur,ycur,zcur);
			}
		}
                
		if(action!=0){
                        /* if diffusion is possible, execute it */
                        if(check==POROSITY){
                                mic[xcur][ycur][zcur]=POROSITY;
                                mic[xnew][ynew][znew]=DIFFCH;
                        }
                        else{
                                /* indicate that diffusing CH species */
                                /* remained at original location */
                                action=7;
                        }
                }
        }
       	return(action);
}

/* routine to add extra C3AH6 when diffusing C3A nucleates or reacts at */
/* C3AH6 surface at location (xpres,ypres,zpres) */
/* Called by movec3a */
/* Calls moveone and edgecnt */
void extc3ah6(xpres,ypres,zpres)
        int xpres,ypres,zpres;
{
        int check,sump,xchr,ychr,zchr,fchr,i1,plok,action,numnear;
        long int tries;

/* First try 6 neighboring locations until      */
/* 	a) successful				*/
/*	b) all 6 sites are tried or             */
/*	c) 100 random attempts are made 	*/
        fchr=0;
        sump=1;
        for(i1=1;((i1<=100)&&(fchr==0)&&(sump!=30030));i1++){
		
                /* determine new coordinates (using periodic boundaries) */
                xchr=xpres;
                ychr=ypres;
                zchr=zpres;
                action=0;
                sump*=moveone(&xchr,&ychr,&zchr,&action,sump);
                if(action==0){printf("Error in action value in extc3ah6 \n");}
                check=mic[xchr][ychr][zchr];

                /* if neighbor is pore space, convert it to C3AH6 */
                if(check==POROSITY){
                        mic[xchr][ychr][zchr]=C3AH6;
			count[C3AH6]+=1;
			count[POROSITY]-=1;
                        fchr=1;
                }
        }

        /* if unsuccessful, add C3AH6 at random location in pore space */
        tries=0;
        while(fchr==0){
                tries+=1;
                xchr=(int)((float)SYSIZE*ran1(seed));
                ychr=(int)((float)SYSIZE*ran1(seed));
                zchr=(int)((float)SYSIZE*ran1(seed));
                if(xchr>=SYSIZE){xchr=0;}
                if(ychr>=SYSIZE){ychr=0;}
                if(zchr>=SYSIZE){zchr=0;}
                check=mic[xchr][ychr][zchr];

                if(check==POROSITY){
                        numnear=edgecnt(xchr,ychr,zchr,C3AH6,C3A,C3AH6);
                        /* Be sure that new C3AH6 is in contact with */
                        /* at least one C3AH6 or C3A */
                        if((tries>5000)||(numnear<26)){
                                mic[xchr][ychr][zchr]=C3AH6;
				count[C3AH6]+=1;
				count[POROSITY]-=1;
                                fchr=1;
                        }
                }
        }
}

/* routine to move a diffusing C3A species */
/* from location (xcur,ycur,zcur) with nucleation probability of nucprob */
/* Called by hydrate */
/* Calls extc3ah6, moveone, extettr, and extafm */
int movec3a(xcur,ycur,zcur,finalstep,nucprob)
        int xcur,ycur,zcur,finalstep;
        float nucprob;
{
        int check,xnew,ynew,znew,plok,action,sumgarb,sumold;
        int xexp,yexp,zexp,nexp,iexp,newact;
        float pgen,pexp,pafm,pgrow,p2diff;

        /* First be sure that a diffusing C3A species is at (xcur,ycur,zcur) */
        if(mic[xcur][ycur][zcur]!=DIFFC3A){
                action=0;
                return(action);
        }

        /* Check for nucleation into solid C3AH6 */
        pgen=ran1(seed);
	p2diff=ran1(seed);

        if((nucprob>=pgen)||(finalstep==1)){
                action=0;
                mic[xcur][ycur][zcur]=C3AH6;
		count[C3AH6]+=1;
                /* decrement count of diffusing C3A species */
                count[DIFFC3A]-=1;
                /* allow for probabilistic-based expansion of C3AH6 */
                /* crystal to account for volume stoichiometry */
                pexp=ran1(seed);
                if(pexp<=0.69){
                        extc3ah6(xcur,ycur,zcur);
                }
        }
        else{
                /* determine new coordinates (using periodic boundaries) */
                xnew=xcur;
                ynew=ycur;
                znew=zcur;
                action=0;
                sumold=1;
                sumgarb=moveone(&xnew,&ynew,&znew,&action,sumold);
                if(action==0){printf("Error in value of action in movec3a \n");}
                check=mic[xnew][ynew][znew];
	
                /* check for growth of C3AH6 crystal */
                if(check==C3AH6){
	                pgrow=ran1(seed);
                        /* Try to slow down growth of C3AH6 crystals to */
                        /* promote ettringite and Afm formation */
                        if(pgrow<=C3AH6GROW){
                                mic[xcur][ycur][zcur]=C3AH6;
				count[C3AH6]+=1;
                                count[DIFFC3A]-=1;
                                action=0;
                         /* allow for probabilistic-based expansion of C3AH6 */
                         /* crystal to account for volume stoichiometry */
                                pexp=ran1(seed);
                                if(pexp<=0.69){
                                        extc3ah6(xcur,ycur,zcur);
                                }
                        }
                }

                /* examine reaction with diffusing gypsum to form ettringite */
                /* Only allow reaction with diffusing gypsum */
                else if((check==DIFFGYP)&&(p2diff<C3AGYP)){
                        /* convert diffusing gypsum to ettringite */
                        mic[xnew][ynew][znew]=ETTR;
			count[ETTR]+=1;
                        /* decrement counts of diffusing gypsum */
                        count[DIFFGYP]-=1;
                        action=0;

/* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
                        pexp=ran1(seed);
                        nexp=2;
                        if(pexp<=0.40){
                                mic[xcur][ycur][zcur]=ETTR;
				count[ETTR]+=1;
				count[DIFFC3A]-=1;
                                nexp=1;
                        }
                        else{
              /* indicate that diffusing species remains in current location */
                                action=7;
                                nexp=2;
                        }

        /* Perform expansion that occurs when ettringite is formed */
        /* xexp, yexp and zexp are the coordinates of the last ettringite */
        /* pixel to be added */
                        xexp=xnew;
                        yexp=ynew;
                        zexp=znew;
                        for(iexp=1;iexp<=nexp;iexp++){
                                newact=extettr(xexp,yexp,zexp,0);
                                /* update xexp, yexp and zexp */
                                switch (newact){
                                        case 1:
                                                xexp-=1;
                                                if(xexp<0){xexp=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xexp+=1;
                                                if(xexp>=SYSIZE){xexp=0;}
                                                break;
                                        case 3:
                                                yexp-=1;
                                                if(yexp<0){yexp=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                yexp+=1;
                                                if(yexp>=SYSIZE){yexp=0;}
                                                break;
                                        case 5:
                                                zexp-=1;
                                                if(zexp<0){zexp=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zexp+=1;
                                                if(zexp>=SYSIZE){zexp=0;}
                                                break;
                                        default:
                                                break;
                                 }
                        }

                  /* probabilistic-based expansion for last ettringite pixel */
                        pexp=ran1(seed);
                        if(pexp<=0.30){
                                newact=extettr(xexp,yexp,zexp,0);
                        }
                }
            /* examine reaction with diffusing hemihydrate to form ettringite */
                /* Only allow reaction with diffusing hemihydrate */
                else if((check==DIFFHEM)&&(p2diff<C3AGYP)){
                        /* convert diffusing hemihydrate to ettringite */
                        mic[xnew][ynew][znew]=ETTR;
			count[ETTR]+=1;
                        /* decrement counts of diffusing hemihydrate */
                        count[DIFFHEM]-=1;
                        action=0;

/* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
                        pexp=ran1(seed);
                        nexp=3;
                        if(pexp<=0.5583){
                                mic[xcur][ycur][zcur]=ETTR;
				count[ETTR]+=1;
				count[DIFFC3A]-=1;
                                nexp=2;
                        }
                        else{
              /* indicate that diffusing species remains in current location */
                                action=7;
                                nexp=3;
                        }

        /* Perform expansion that occurs when ettringite is formed */
        /* xexp, yexp and zexp are the coordinates of the last ettringite */
        /* pixel to be added */
                        xexp=xnew;
                        yexp=ynew;
                        zexp=znew;
                        for(iexp=1;iexp<=nexp;iexp++){
                                newact=extettr(xexp,yexp,zexp,0);
                                /* update xexp, yexp and zexp */
                                switch (newact){
                                        case 1:
                                                xexp-=1;
                                                if(xexp<0){xexp=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xexp+=1;
                                                if(xexp>=SYSIZE){xexp=0;}
                                                break;
                                        case 3:
                                                yexp-=1;
                                                if(yexp<0){yexp=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                yexp+=1;
                                                if(yexp>=SYSIZE){yexp=0;}
                                                break;
                                        case 5:
                                                zexp-=1;
                                                if(zexp<0){zexp=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zexp+=1;
                                                if(zexp>=SYSIZE){zexp=0;}
                                                break;
                                        default:
                                                break;
                                 }
                        }

                  /* probabilistic-based expansion for last ettringite pixel */
                        pexp=ran1(seed);
                        if(pexp<=0.6053){
                                newact=extettr(xexp,yexp,zexp,0);
                        }
                }
             /* examine reaction with diffusing anhydrite to form ettringite */
                /* Only allow reaction with diffusing anhydrite */
                else if((check==DIFFANH)&&(p2diff<C3AGYP)){
                        /* convert diffusing anhydrite to ettringite */
                        mic[xnew][ynew][znew]=ETTR;
			count[ETTR]+=1;
                        /* decrement counts of diffusing anhydrite */
                        count[DIFFANH]-=1;
                        action=0;

/* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
                        pexp=ran1(seed);
                        nexp=3;
                        if(pexp<=0.569){
                                mic[xcur][ycur][zcur]=ETTR;
				count[ETTR]+=1;
				count[DIFFC3A]-=1;
                                nexp=2;
                        }
                        else{
              /* indicate that diffusing species remains in current location */
                                action=7;
                                nexp=3;
                        }

        /* Perform expansion that occurs when ettringite is formed */
        /* xexp, yexp and zexp are the coordinates of the last ettringite */
        /* pixel to be added */
                        xexp=xnew;
                        yexp=ynew;
                        zexp=znew;
                        for(iexp=1;iexp<=nexp;iexp++){
                                newact=extettr(xexp,yexp,zexp,0);
                                /* update xexp, yexp and zexp */
                                switch (newact){
                                        case 1:
                                                xexp-=1;
                                                if(xexp<0){xexp=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xexp+=1;
                                                if(xexp>=SYSIZE){xexp=0;}
                                                break;
                                        case 3:
                                                yexp-=1;
                                                if(yexp<0){yexp=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                yexp+=1;
                                                if(yexp>=SYSIZE){yexp=0;}
                                                break;
                                        case 5:
                                                zexp-=1;
                                                if(zexp<0){zexp=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zexp+=1;
                                                if(zexp>=SYSIZE){zexp=0;}
                                                break;
                                        default:
                                                break;
                                 }
                        }

                  /* probabilistic-based expansion for last ettringite pixel */
                        pexp=ran1(seed);
                        if(pexp<=0.6935){
                                newact=extettr(xexp,yexp,zexp,0);
                        }
                }
                /* examine reaction with diffusing CaCl2 to form FREIDEL */
                /* Only allow reaction with diffusing CaCl2 */
                else if(check==DIFFCACL2){
                        /* convert diffusing C3A to Freidel's salt */
                        mic[xcur][ycur][zcur]=FREIDEL;
			count[FREIDEL]+=1;
                        /* decrement counts of diffusing C3A and CaCl2 */
                        count[DIFFC3A]-=1;
                        action=0;

/* convert diffusing CACL2 to solid FREIDEL or else leave as a diffusing CACL2 */
                        pexp=ran1(seed);
                        nexp=2;
                        if(pexp<=0.5793){
                                mic[xnew][ynew][znew]=FREIDEL;
				count[FREIDEL]+=1;
				count[DIFFCACL2]-=1;
                                nexp=1;
                        }
                        else{
                                nexp=2;
                        }

        /* Perform expansion that occurs when Freidel's salt is formed */
        /* xexp, yexp and zexp are the coordinates of the last FREIDEL */
        /* pixel to be added */
                        xexp=xnew;
                        yexp=ynew;
                        zexp=znew;
                        for(iexp=1;iexp<=nexp;iexp++){
                                newact=extfreidel(xexp,yexp,zexp);
                                /* update xexp, yexp and zexp */
                                switch (newact){
                                        case 1:
                                                xexp-=1;
                                                if(xexp<0){xexp=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xexp+=1;
                                                if(xexp>=SYSIZE){xexp=0;}
                                                break;
                                        case 3:
                                                yexp-=1;
                                                if(yexp<0){yexp=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                yexp+=1;
                                                if(yexp>=SYSIZE){yexp=0;}
                                                break;
                                        case 5:
                                                zexp-=1;
                                                if(zexp<0){zexp=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zexp+=1;
                                                if(zexp>=SYSIZE){zexp=0;}
                                                break;
                                        default:
                                                break;
                                 }
                        }

                  /* probabilistic-based expansion for last FREIDEL pixel */
                        pexp=ran1(seed);
                        if(pexp<=0.3295){
                                newact=extfreidel(xexp,yexp,zexp);
                        }
                }
                /* examine reaction with diffusing CAS2 to form STRAT */
                /* Only allow reaction with diffusing (not solid) CAS2 */
                else if(check==DIFFCAS2){
                        /* convert diffusing CAS2 to stratlingite */
                        mic[xnew][ynew][znew]=STRAT;
			count[STRAT]+=1;
                        /* decrement counts of diffusing C3A and CAS2 */
                        count[DIFFCAS2]-=1;
                        action=0;
	
/* convert diffusing C3A to solid STRAT or else leave as a diffusing C3A */
                        pexp=ran1(seed);
                        nexp=3;
                        if(pexp<=0.886){
                                mic[xcur][ycur][zcur]=STRAT;
				count[STRAT]+=1;
				count[DIFFC3A]-=1;
                                nexp=2;
                        }
                        else{
				action=7;
                                nexp=3;
                        }

        /* Perform expansion that occurs when stratlingite is formed */
        /* xexp, yexp and zexp are the coordinates of the last STRAT */
        /* pixel to be added */
                        xexp=xnew;
                        yexp=ynew;
                        zexp=znew;
                        for(iexp=1;iexp<=nexp;iexp++){
                                newact=extstrat(xexp,yexp,zexp);
                                /* update xexp, yexp and zexp */
                                switch (newact){
                                        case 1:
                                                xexp-=1;
                                                if(xexp<0){xexp=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xexp+=1;
                                                if(xexp>=SYSIZE){xexp=0;}
                                                break;
                                        case 3:
                                                yexp-=1;
                                                if(yexp<0){yexp=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                yexp+=1;
                                                if(yexp>=SYSIZE){yexp=0;}
                                                break;
                                        case 5:
                                                zexp-=1;
                                                if(zexp<0){zexp=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zexp+=1;
                                                if(zexp>=SYSIZE){zexp=0;}
                                                break;
                                        default:
                                                break;
                                 }
                        }

                  /* probabilistic-based expansion for last STRAT pixel */
                        pexp=ran1(seed);
                        if(pexp<=0.286){
                                newact=extstrat(xexp,yexp,zexp);
                        }
                }
/* check for reaction with diffusing or solid ettringite to form AFm */
/* reaction at solid ettringite only possible if ettringite is soluble */
/* and even then on a limited bases to avoid a great formation of AFm */
/* when ettringite first becomes soluble */
                pgrow=ran1(seed);
   if((check==DIFFETTR)||((check==ETTR)&&(soluble[ETTR]==1)&&(pgrow<=C3AETTR))){
                /* convert diffusing or solid ettringite to AFm */
                mic[xnew][ynew][znew]=AFM;
		count[AFM]+=1;
                /* decrement count of ettringite */
		count[check]-=1;
                action=0;
	
                /* convert diffusing C3A to AFm or leave as diffusing C3A */
                pexp=ran1(seed);
                if(pexp<=0.2424){
                        mic[xcur][ycur][zcur]=AFM;
			count[AFM]+=1;
			count[DIFFC3A]-=1;
                        pafm=(-0.1);
                }
                else{
                        action=7;
			pafm=0.04699;
                }

                /* probabilistic-based expansion for new AFm pixel */
                pexp=ran1(seed);
                if(pexp<=pafm){
                        extafm(xnew,ynew,znew);
                }
        }
        if((action!=0)&&(action!=7)){

                /* if diffusion is possible, execute it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                        mic[xnew][ynew][znew]=DIFFC3A;
                }
                else{
                        /* indicate that diffusing C3A remained */
                        /* at original location */
                        action=7;
                }
        }
        }
        return(action);
}

/* routine to move a diffusing C4A species */
/* from location (xcur,ycur,zcur) with nucleation probability of nucprob */
/* Called by hydrate */
/* Calls extc3ah6, moveone, extettr, and extafm */
int movec4a(xcur,ycur,zcur,finalstep,nucprob)
        int xcur,ycur,zcur,finalstep;
        float nucprob;
{
        int check,xnew,ynew,znew,plok,action,sumgarb,sumold;
        int xexp,yexp,zexp,nexp,iexp,newact;
        float pgen,pexp,pafm,pgrow,p2diff;

        /* First be sure that a diffusing C4A species is at (xcur,ycur,zcur) */
        if(mic[xcur][ycur][zcur]!=DIFFC4A){
                action=0;
                return(action);
        }

        /* Check for nucleation into solid C3AH6 */
        pgen=ran1(seed);
	p2diff=ran1(seed);

        if((nucprob>=pgen)||(finalstep==1)){
                action=0;
                mic[xcur][ycur][zcur]=C3AH6;
		count[C3AH6]+=1;
                /* decrement count of diffusing C3A species */
                count[DIFFC4A]-=1;
                /* allow for probabilistic-based expansion of C3AH6 */
                /* crystal to account for volume stoichiometry */
                pexp=ran1(seed);
                if(pexp<=0.69){
                        extc3ah6(xcur,ycur,zcur);
                }
        }
        else{
                /* determine new coordinates (using periodic boundaries) */
                xnew=xcur;
                ynew=ycur;
                znew=zcur;
                action=0;
                sumold=1;
                sumgarb=moveone(&xnew,&ynew,&znew,&action,sumold);
                if(action==0){printf("Error in value of action in movec4a \n");}
                check=mic[xnew][ynew][znew];
	
                /* check for growth of C3AH6 crystal */
                if(check==C3AH6){
	                pgrow=ran1(seed);
                        /* Try to slow down growth of C3AH6 crystals to */
                        /* promote ettringite and Afm formation */
                        if(pgrow<=C3AH6GROW){
                                mic[xcur][ycur][zcur]=C3AH6;
				count[C3AH6]+=1;
                                count[DIFFC4A]-=1;
                                action=0;
                         /* allow for probabilistic-based expansion of C3AH6 */
                         /* crystal to account for volume stoichiometry */
                                pexp=ran1(seed);
                                if(pexp<=0.69){
                                        extc3ah6(xcur,ycur,zcur);
                                }
                        }
                }

                /* examine reaction with diffusing gypsum to form ettringite */
                /* Only allow reaction with diffusing gypsum */
                else if((check==DIFFGYP)&&(p2diff<C3AGYP)){
                        /* convert diffusing gypsum to ettringite */
                        mic[xnew][ynew][znew]=ETTRC4AF;
			count[ETTRC4AF]+=1;
                        /* decrement counts of diffusing gypsum */
                        count[DIFFGYP]-=1;
                        action=0;

/* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
                        pexp=ran1(seed);
                        nexp=2;
                        if(pexp<=0.40){
                                mic[xcur][ycur][zcur]=ETTRC4AF;
				count[ETTRC4AF]+=1;
				count[DIFFC4A]-=1;
                                nexp=1;
                        }
                        else{
              /* indicate that diffusing species remains in current location */
                                action=7;
                                nexp=2;
                        }

        /* Perform expansion that occurs when ettringite is formed */
        /* xexp, yexp and zexp are the coordinates of the last ettringite */
        /* pixel to be added */
                        xexp=xnew;
                        yexp=ynew;
                        zexp=znew;
                        for(iexp=1;iexp<=nexp;iexp++){
                                newact=extettr(xexp,yexp,zexp,1);
                                /* update xexp, yexp and zexp */
                                switch (newact){
                                        case 1:
                                                xexp-=1;
                                                if(xexp<0){xexp=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xexp+=1;
                                                if(xexp>=SYSIZE){xexp=0;}
                                                break;
                                        case 3:
                                                yexp-=1;
                                                if(yexp<0){yexp=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                yexp+=1;
                                                if(yexp>=SYSIZE){yexp=0;}
                                                break;
                                        case 5:
                                                zexp-=1;
                                                if(zexp<0){zexp=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zexp+=1;
                                                if(zexp>=SYSIZE){zexp=0;}
                                                break;
                                        default:
                                                break;
                                 }
                        }

                  /* probabilistic-based expansion for last ettringite pixel */
                        pexp=ran1(seed);
                        if(pexp<=0.30){
                                newact=extettr(xexp,yexp,zexp,1);
                        }
                }
                /* examine reaction with diffusing hemi to form ettringite */
                /* Only allow reaction with diffusing hemihydrate */
                else if((check==DIFFHEM)&&(p2diff<C3AGYP)){
                        /* convert diffusing hemihydrate to ettringite */
                        mic[xnew][ynew][znew]=ETTRC4AF;
			count[ETTRC4AF]+=1;
                        /* decrement counts of diffusing hemihydrate */
                        count[DIFFHEM]-=1;
                        action=0;

/* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
                        pexp=ran1(seed);
                        nexp=3;
                        if(pexp<=0.5583){
                                mic[xcur][ycur][zcur]=ETTRC4AF;
				count[ETTRC4AF]+=1;
				count[DIFFC4A]-=1;
                                nexp=2;
                        }
                        else{
              /* indicate that diffusing species remains in current location */
                                action=7;
                                nexp=3;
                        }

        /* Perform expansion that occurs when ettringite is formed */
        /* xexp, yexp and zexp are the coordinates of the last ettringite */
        /* pixel to be added */
                        xexp=xnew;
                        yexp=ynew;
                        zexp=znew;
                        for(iexp=1;iexp<=nexp;iexp++){
                                newact=extettr(xexp,yexp,zexp,1);
                                /* update xexp, yexp and zexp */
                                switch (newact){
                                        case 1:
                                                xexp-=1;
                                                if(xexp<0){xexp=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xexp+=1;
                                                if(xexp>=SYSIZE){xexp=0;}
                                                break;
                                        case 3:
                                                yexp-=1;
                                                if(yexp<0){yexp=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                yexp+=1;
                                                if(yexp>=SYSIZE){yexp=0;}
                                                break;
                                        case 5:
                                                zexp-=1;
                                                if(zexp<0){zexp=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zexp+=1;
                                                if(zexp>=SYSIZE){zexp=0;}
                                                break;
                                        default:
                                                break;
                                 }
                        }

                  /* probabilistic-based expansion for last ettringite pixel */
                        pexp=ran1(seed);
                        if(pexp<=0.6053){
                                newact=extettr(xexp,yexp,zexp,1);
                        }
                }
             /* examine reaction with diffusing anhydrite to form ettringite */
                /* Only allow reaction with diffusing anhydrite */
                else if((check==DIFFANH)&&(p2diff<C3AGYP)){
                        /* convert diffusing anhydrite to ettringite */
                        mic[xnew][ynew][znew]=ETTRC4AF;
			count[ETTRC4AF]+=1;
                        /* decrement counts of diffusing anhydrite */
                        count[DIFFANH]-=1;
                        action=0;

/* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
                        pexp=ran1(seed);
                        nexp=3;
                        if(pexp<=0.569){
                                mic[xcur][ycur][zcur]=ETTRC4AF;
				count[ETTRC4AF]+=1;
				count[DIFFC4A]-=1;
                                nexp=2;
                        }
                        else{
              /* indicate that diffusing species remains in current location */
                                action=7;
                                nexp=3;
                        }

        /* Perform expansion that occurs when ettringite is formed */
        /* xexp, yexp and zexp are the coordinates of the last ettringite */
        /* pixel to be added */
                        xexp=xnew;
                        yexp=ynew;
                        zexp=znew;
                        for(iexp=1;iexp<=nexp;iexp++){
                                newact=extettr(xexp,yexp,zexp,1);
                                /* update xexp, yexp and zexp */
                                switch (newact){
                                        case 1:
                                                xexp-=1;
                                                if(xexp<0){xexp=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xexp+=1;
                                                if(xexp>=SYSIZE){xexp=0;}
                                                break;
                                        case 3:
                                                yexp-=1;
                                                if(yexp<0){yexp=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                yexp+=1;
                                                if(yexp>=SYSIZE){yexp=0;}
                                                break;
                                        case 5:
                                                zexp-=1;
                                                if(zexp<0){zexp=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zexp+=1;
                                                if(zexp>=SYSIZE){zexp=0;}
                                                break;
                                        default:
                                                break;
                                 }
                        }

                  /* probabilistic-based expansion for last ettringite pixel */
                        pexp=ran1(seed);
                        if(pexp<=0.6935){
                                newact=extettr(xexp,yexp,zexp,1);
                        }
                }
                /* examine reaction with diffusing CaCl2 to form FREIDEL */
                /* Only allow reaction with diffusing CaCl2 */
                else if(check==DIFFCACL2){
                        /* convert diffusing C3A to Freidel's salt */
                        mic[xcur][ycur][zcur]=FREIDEL;
			count[FREIDEL]+=1;
                        /* decrement counts of diffusing C3A and CaCl2 */
                        count[DIFFC4A]-=1;
                        action=0;

/* convert diffusing CACL2 to solid FREIDEL or else leave as a diffusing CACL2 */
                        pexp=ran1(seed);
                        nexp=2;
                        if(pexp<=0.5793){
                                mic[xnew][ynew][znew]=FREIDEL;
				count[FREIDEL]+=1;
				count[DIFFCACL2]-=1;
                                nexp=1;
                        }
                        else{
                                nexp=2;
                        }

        /* Perform expansion that occurs when Freidel's salt is formed */
        /* xexp, yexp and zexp are the coordinates of the last FREIDEL */
        /* pixel to be added */
                        xexp=xnew;
                        yexp=ynew;
                        zexp=znew;
                        for(iexp=1;iexp<=nexp;iexp++){
                                newact=extfreidel(xexp,yexp,zexp);
                                /* update xexp, yexp and zexp */
                                switch (newact){
                                        case 1:
                                                xexp-=1;
                                                if(xexp<0){xexp=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xexp+=1;
                                                if(xexp>=SYSIZE){xexp=0;}
                                                break;
                                        case 3:
                                                yexp-=1;
                                                if(yexp<0){yexp=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                yexp+=1;
                                                if(yexp>=SYSIZE){yexp=0;}
                                                break;
                                        case 5:
                                                zexp-=1;
                                                if(zexp<0){zexp=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zexp+=1;
                                                if(zexp>=SYSIZE){zexp=0;}
                                                break;
                                        default:
                                                break;
                                 }
                        }

                  /* probabilistic-based expansion for last FREIDEL pixel */
                        pexp=ran1(seed);
                        if(pexp<=0.3295){
                                newact=extfreidel(xexp,yexp,zexp);
                        }
                }
                /* examine reaction with diffusing CAS2 to form STRAT */
                /* Only allow reaction with diffusing (not solid) CAS2 */
                else if(check==DIFFCAS2){
                        /* convert diffusing CAS2 to stratlingite */
                        mic[xnew][ynew][znew]=STRAT;
			count[STRAT]+=1;
                        /* decrement counts of diffusing CAS2 */
                        count[DIFFCAS2]-=1;
                        action=0;

/* convert diffusing C3A to solid STRAT or else leave as a diffusing C3A */
                        pexp=ran1(seed);
                        nexp=3;
                        if(pexp<=0.886){
                                mic[xcur][ycur][zcur]=STRAT;
				count[STRAT]+=1;
				count[DIFFC4A]-=1;
                                nexp=2;
                        }
                        else{
				action=7;
                                nexp=3;
                        }

        /* Perform expansion that occurs when stratlingite is formed */
        /* xexp, yexp and zexp are the coordinates of the last STRAT */
        /* pixel to be added */
                        xexp=xnew;
                        yexp=ynew;
                        zexp=znew;
                        for(iexp=1;iexp<=nexp;iexp++){
                                newact=extstrat(xexp,yexp,zexp);
                                /* update xexp, yexp and zexp */
                                switch (newact){
                                        case 1:
                                                xexp-=1;
                                                if(xexp<0){xexp=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xexp+=1;
                                                if(xexp>=SYSIZE){xexp=0;}
                                                break;
                                        case 3:
                                                yexp-=1;
                                                if(yexp<0){yexp=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                yexp+=1;
                                                if(yexp>=SYSIZE){yexp=0;}
                                                break;
                                        case 5:
                                                zexp-=1;
                                                if(zexp<0){zexp=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zexp+=1;
                                                if(zexp>=SYSIZE){zexp=0;}
                                                break;
                                        default:
                                                break;
                                 }
                        }

                  /* probabilistic-based expansion for last STRAT pixel */
                        pexp=ran1(seed);
                        if(pexp<=0.286){
                                newact=extstrat(xexp,yexp,zexp);
                        }
                }
/* check for reaction with diffusing or solid ettringite to form AFm */
/* reaction at solid ettringite only possible if ettringite is soluble */
/* and even then on a limited bases to avoid a great formation of AFm */
/* when ettringite first becomes soluble */
                pgrow=ran1(seed);
   if((check==DIFFETTR)||((check==ETTR)&&(soluble[ETTR]==1)&&(pgrow<=C3AETTR))){
                /* convert diffusing or solid ettringite to AFm */
                mic[xnew][ynew][znew]=AFM;
		count[AFM]+=1;
                /* decrement count of ettringite */
		count[check]-=1;
                action=0;

                /* convert diffusing C4A to AFm or leave as diffusing C4A */
                pexp=ran1(seed);
                if(pexp<=0.2424){
                        mic[xcur][ycur][zcur]=AFM;
			count[AFM]+=1;
			count[DIFFC4A]-=1;
                        pafm=(-0.1);
                }
                else{
                        action=7;
			pafm=0.04699;
                }

                /* probabilistic-based expansion for new AFm pixel */
                pexp=ran1(seed);
                if(pexp<=pafm){
                        extafm(xnew,ynew,znew);
                }
        }
        if((action!=0)&&(action!=7)){

                /* if diffusion is possible, execute it */
                if(check==POROSITY){
                        mic[xcur][ycur][zcur]=POROSITY;
                        mic[xnew][ynew][znew]=DIFFC4A;
                }
                else{
                        /* indicate that diffusing C4A remained */
                        /* at original location */
                        action=7;
                }
        }
        }
        return(action);
}

/* routine to oversee hydration by updating position of all */
/* remaining diffusing species */
/* Calls movech, movec3a, movefh3, moveettr, movecsh, and movegyp */
void hydrate(fincyc,stepmax,chpar1,chpar2,hgpar1,hgpar2,fhpar1,fhpar2,gypar1,gypar2)
        int fincyc,stepmax;
        float chpar1,chpar2,hgpar1,hgpar2,fhpar1,fhpar2,gypar1,gypar2;
{
        int xpl,ypl,zpl,phpl,agepl,xpnew,ypnew,zpnew;
        float chprob,c3ah6prob,fh3prob,gypprob;
        long int icnt,nleft,ntodo,ndale;
        int istep,termflag,reactf;
        float beterm;
        struct ants *curant,*antgone;

        ntodo=nmade;
        nleft=nmade;
        termflag=0;

/* Perform diffusion until all reacted or max. # of diffusion steps reached */
        for(istep=1;((istep<=stepmax)&&(nleft>0));istep++){
                if((fincyc==1)&&(istep==stepmax)){termflag=1;} 

                nleft=0;
                ndale=0;

                /* determine probabilities for CH and C3AH6 nucleation */
                beterm=exp(-(double)(count[DIFFCH])/chpar2);
                chprob=chpar1*(1.-beterm);
                beterm=exp(-(double)(count[DIFFC3A])/hgpar2);
                c3ah6prob=hgpar1*(1.-beterm);
                beterm=exp(-(double)(count[DIFFFH3])/fhpar2);
                fh3prob=fhpar1*(1.-beterm);
                beterm=exp(-(double)(count[DIFFANH]+count[DIFFHEM])/gypar2);
                gypprob=gypar1*(1.-beterm);

                /* Process each diffusing species in turn */
                curant=headant->nextant;
                while(curant!=NULL){
                        ndale+=1;
                        xpl=curant->x;
                        ypl=curant->y;
                        zpl=curant->z;
                        phpl=curant->id;
			agepl=curant->cycbirth;

       /* based on ID, call appropriate routine to process diffusing species */
                        switch (phpl) {
                                case DIFFCSH:
					/* printf("Calling movecsh \n");
					fflush(stdout); */
                                        reactf=movecsh(xpl,ypl,zpl,termflag,agepl);
                                        break;
                                case DIFFANH:
					/* printf("Calling moveanh \n");
					fflush(stdout); */
                                        reactf=moveanh(xpl,ypl,zpl,termflag,gypprob);
                                        break;
                                case DIFFHEM:
					/* printf("Calling movehem \n");
					fflush(stdout); */
                                        reactf=movehem(xpl,ypl,zpl,termflag,gypprob);
                                        break;
                                case DIFFCH:
					/* printf("Calling movech \n");
					fflush(stdout); */
                                  reactf=movech(xpl,ypl,zpl,termflag,chprob);
                                  break;
                                case DIFFFH3:
					/* printf("Calling movefh3 \n");
					fflush(stdout); */
                                  reactf=movefh3(xpl,ypl,zpl,termflag,fh3prob);
                                  break;
                                case DIFFGYP:
					/* printf("Calling movegyp \n");
					fflush(stdout); */
                                        reactf=movegyp(xpl,ypl,zpl,termflag);
                                       	break;
                                case DIFFC3A:
					/* printf("Calling movec3a \n");
					fflush(stdout); */
                                 reactf=movec3a(xpl,ypl,zpl,termflag,c3ah6prob);
                                 break;
                                case DIFFC4A:
					/* printf("Calling movec4a \n");
					fflush(stdout); */
                                 reactf=movec4a(xpl,ypl,zpl,termflag,c3ah6prob);
                                 break;
                                case DIFFETTR:
					/* printf("Calling moveettr \n");
					fflush(stdout); */
                                        reactf=moveettr(xpl,ypl,zpl,termflag);
                                        break;
                                case DIFFCACL2:
					/* printf("Calling movecacl2 \n");
					fflush(stdout); */
                                        reactf=movecacl2(xpl,ypl,zpl,termflag);
                                        break;
                                case DIFFCAS2:
					/* printf("Calling movecas2 \n");
					fflush(stdout); */
                                        reactf=movecas2(xpl,ypl,zpl,termflag);
                                        break;
                                case DIFFAS:
					/* printf("Calling moveas \n");
					fflush(stdout); */
                                        reactf=moveas(xpl,ypl,zpl,termflag);
                                        break;
                                case DIFFCACO3:
					/* printf("Calling movecaco3 \n");
					fflush(stdout); */
                                        reactf=movecaco3(xpl,ypl,zpl,termflag);
                                        break;
                                default:
                                        printf("Error in ID of phase \n");
                                        break;
                        }

                        /* if no reaction */
                        if(reactf!=0){
                                nleft+=1;
                                xpnew=xpl;
                                ypnew=ypl;
                                zpnew=zpl;

                                /* update location of diffusing species */
                                switch (reactf) {
                                        case 1:
                                                xpnew-=1;
                                                if(xpnew<0){xpnew=(SYSIZEM1);}
                                                break;
                                        case 2:
                                                xpnew+=1;
                                                if(xpnew>=SYSIZE){xpnew=0;}
                                                break;
                                        case 3:
                                                ypnew-=1;
                                                if(ypnew<0){ypnew=(SYSIZEM1);}
                                                break;
                                        case 4:
                                                ypnew+=1;
                                                if(ypnew>=SYSIZE){ypnew=0;}
                                                break;
                                        case 5:
                                                zpnew-=1;
                                                if(zpnew<0){zpnew=(SYSIZEM1);}
                                                break;
                                        case 6:
                                                zpnew+=1;
                                                if(zpnew>=SYSIZE){zpnew=0;}
                                                break;
                                        default:
                                                break;
                                }

                                /* store new location of diffusing species */
                                curant->x=xpnew;
                                curant->y=ypnew;
                                curant->z=zpnew;
                                curant->id=phpl;
                                curant=curant->nextant;
                        } /* end of reactf!=0 loop */
                        /* else remove ant from list */
                        else{
                                if(ndale==1){
                                        headant->nextant=curant->nextant;
                                }
                                else{
                                     (curant->prevant)->nextant=curant->nextant;
                                }
                                if(curant->nextant!=NULL){
                                     (curant->nextant)->prevant=curant->prevant;
                                }
                                else{
                                        tailant=curant->prevant;
                                }
                                antgone=curant;
                                curant=curant->nextant;
                                free(antgone);
                                ngoing-=1;
                        }
                } /* end of curant loop */
                ntodo=nleft;
        } /* end of istep loop */
}
