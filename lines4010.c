#include <stdio.h>
#include <signal.h>
#include "grf_obj.h"


typedef struct{
   vec4type *pts;                       /* points array */
   long *lines;                          /* lines array, index into points */
   long nl,np;                           /* # of lines, # of points */
   char name[80];                       /* name of file */
   long white;                           /* white on black ? */
} wiretype;


char nullstr[2];


int onmemfault(sig)
int sig;
{
   printf("Memory fault!\n");
   exit(1);
} 


int onbuserror(sig)
int sig;
{
   printf("Bus error!\n");
   exit(1);
} 


int onfpe(sig)
int sig;
{
   printf("Floating point exception!\n");
   exit(1);
} 



void initstuff(wires)
wiretype *wires;
{
   strcpy(wires->name,"No file loaded");
   wires->pts=NULL;
   wires->white=1;
   signal(SIGSEGV,onmemfault);
   signal(SIGBUS,onbuserror);
   signal(SIGFPE,onfpe);
}


int readscene(fname,wires)
char *fname;
wiretype *wires;
{
   FILE *infile;
   char sceneerr[100];
   scenetype srec;
   int totpts,totlns,objcnt,numpts,numlns,ptcnt,lncnt;
   mat4type view,proj,viewing;

   if((infile=fopen(fname,"r"))==NULL){
      fprintf(stderr, "can't open file %s for input.\n");
      return(0);
   }
   if(!getscene(infile,&srec,sceneerr)){
      sceneerr[strlen(sceneerr)-1]='\0';
      fprintf(stderr, "Error in %s :\n\t'%s'\n",fname, sceneerr);
      fclose(infile);
      return(0);
   } 
   for(totpts=0,totlns=0,objcnt=0;objcnt<srec.numobj;objcnt++){
      predict(srec.objs+objcnt,&numpts,&numlns);
      totpts+=numpts;
      totlns+=numlns;
   }
   if(wires->pts){
      dispose(wires->pts);
      dispose(wires->lines);
   }
   wires->pts=NULL;
   strcpy(wires->name,"No file loaded");
   if((wires->pts=getsome(vec4type,totpts))==NULL)
      printf("out of memory\n");
   else
      if((wires->lines=getsome(int,totlns*2))==NULL)
         printf("out of memory\n");
      else{
         wires->np=totpts;
         wires->nl=totlns;
         for(ptcnt=0,lncnt=0,objcnt=0;objcnt<srec.numobj;objcnt++){
            objtolines(srec.objs+objcnt,wires->pts+ptcnt,wires->lines+lncnt*2,ptcnt);
            predict(srec.objs+objcnt,&numpts,&numlns);
            ptcnt+=numpts;
            lncnt+=numlns;
         }
         /* clip to viewing volume */
         skewtoortho(&view,&srec.eye,&srec.up,srec.uvec,&srec.targ,srec.tvec);
         proj=identmat4;
         proj.a[0][0]=proj.a[1][1]=1/tan(srec.fov*M_PI/360.0);
         proj.a[2][3]=1.0;
         /*multmat4(&view,&proj,&viewing);*/
         multvmblock(wires->pts,&view,wires->pts,wires->np);
         multvmblock(wires->pts,&proj,wires->pts,wires->np);
         /*for(ptcnt=0;ptcnt<wires->np;ptcnt++){*/
            /*printf("xformed: ");*/
            /*writevec(stdout,wires->pts+ptcnt);*/
         /*}*/
         for(ptcnt=0;ptcnt<wires->np;ptcnt++){
            /*printf("before: ");*/
            /*writevec(stdout,wires->pts+ptcnt);*/
            normvec(wires->pts[ptcnt]);
            /*printf("after: %lf %lf\n",wires->pts[ptcnt].a[0],wires->pts[ptcnt].a[1]);*/
         }
         /*homogblock(wires->pts,wires->np);*/
         strcpy(wires->name,fname);
   }
   dispose(srec.objs);
   fclose(infile);
   return(1);
}


void update(wires)
wiretype *wires;
{
   int xsize,ysize;
   register int i,f,t,*iptr,prev=-1;
   register double halfx,halfy,one=1.0;
   register long *sptr;
   long *ipts;

   if(wires->pts){
      xsize=3000;
      ysize=3000;
      if((ipts=getsome(short,wires->np*2))==NULL)
         printf("out of memory\n");
      else{
         halfx=(-xsize)/2.0;
         halfy=(-ysize)/2.0;
         for(sptr=ipts,i=0;i<wires->np;i++,sptr+=2){
            *sptr=(wires->pts[i].a[0]-one)*halfx;
            *(sptr+1)=(wires->pts[i].a[1]-one)*halfy;
         }
         /* should probably take aspect angle into account here. */
         for(iptr=wires->lines,i=0;i<wires->nl;i++,iptr+=2){
            f=(*iptr)<<1; /*f=wires->lines[i*2];*/
            t=(*(iptr+1))<<1; /*t=wires->lines[i*2+1];*/
            if(prev!=f)
               move(ipts[f],3000-ipts[f+1]);
            cont(ipts[t],3000-ipts[t+1]);
            prev=t;
         }
         dispose(ipts);
      }
   }
}


void cleanup()
{
   exit(0);
}


main(argc,argv)
int argc;
char *argv[];
{
   wiretype wires;

   initstuff(&wires);
   if(argc!=2){
      printf("bad arguments :  '%s file.gdf'\n", argv[0]);
      cleanup();
   }
   system("stty 9600");
   openpl();
   readscene(argv[1],&wires);
   update(&wires);
   closepl();
   cleanup();
}
