/*----------------------------------------------------------------------------

   Written by Brad Grantham
   FILE:   ray6.c
   SYSTEM:   MAC II under A/UX using cc and the grf_0.9 library
   PURPOSE:   this is a skeleton for ray-tracing scenes described in the gdf
    format, using the grf_0.9 library.   command line arguments are:

        ray6 scenefile imagefile [xpix ypix [raypix]] 

    scenefile - file to read scene data from
    imagefile - file to store rgb image in
    xpix      - pixels of image in x dimension    (default 400)
    ypix      - pixels of image in y dimension    (default 400)
    raypix    - # of rays to trace per pixel      (default 16)
                ( 0 - one ray, no perturbation
                  1 - one ray, perturbed
                  2 & up - 2 or more, perturbed)

----------------------------------------------------------------------------*/


#include <stdio.h>                      /* standard includes */
#include <math.h>                       /* math includes */
#include <signal.h>                     /* signals */
#include <sys/time.h>                   /* microsecond timing */
#include "grf_math.h"                   /* graphics math */
#include "grf_view.h"                   /* graphics viewing */
#include "grf_obj.h"                    /* graphics objects */


colortype *getbckgnd(ray,scene,color)
vec4type ray[2];
scenetype *scene;
colortype *color;
{
   color->r=0.0;
   color->g=0.0;
   color->b=0.0;
   return(color);
}


int ssec=0,susec=0,
 xsize,ysize,pixrays,
 xcnt,ycnt,
 stop=0;


double tsec=0.0;


FILE *rayimg;


char scenename[80];


/*----------------------------------------------------------------------------
   PURPOSE:   this function performs the tracing of one ray
   PRECONDITIONS:   ray points to the ray to start tracing, scene points to
    the scene data, color points to the final resulting color of the ray 
   POSTCONDITIONS:   color contains the color that the ray finally assumed.
----------------------------------------------------------------------------*/
void raytrace(ray,scene,color)
vec4type ray[2];                        /* ray to trace */
scenetype *scene;                       /* scene to trace inside */
colortype *color;                       /* color of ray when done */
{
   int done=0;                          /* done tracing */
   objtype *objint=NULL;                /* object that was intersected */
   vec4type ipnt[2],                    /* intersection point 0 - local point, */
                                        /*  1 - scene space point */
    myray[2];                           /* my ray (to chuck around ) */
   colortype objcolor;                  /* object color */

   /*printf("ray sent:\n");*/
   /*writevec(stdout,ray);*/
   /*writevec(stdout,ray+1);*/
   myray[0]=ray[0];
   myray[1]=ray[1];
   *color=white_color;
   do{
      if(isblack24(color)) return;
      if((objint=rayintobj(objint,ipnt,myray,scene,ipnt))==NULL){
                                        /* hides ray-object intersection */
                                        /* does ray-to-object-space xform and */
                                        /*  whatever else. */
         getbckgnd(myray,scene,&objcolor);
                                        /* gets background color of scene from ray */
         done=1;
      }else{
         /*printf("well, hit obj '%s' at",objint->name);*/
         /*writevec(stdout,ipnt);*/
         getobjcolor(objint,ipnt,&objcolor);
                                        /* computes surface color, can take into */
                                        /*  account texture mapping or whatever */
         if(objint->radiant)
            done=1;
         else{
            bounceray(myray,objint,ipnt,myray);
                                        /* take care of diffusiveness/reflectiveness */
         }
      }
      multcolors(color,&objcolor,color);
   }while(!done);
   return;
}


void programstop()
{
   FILE *outfile;
   unsigned char ind=255;
   
   fseek(rayimg,0,0);
   fwrite(&ind,1,1,rayimg);
   fseek(rayimg,0,2);
   fprintf(rayimg,"%s",scenename);
   ind=strlen(scenename);
   fwrite(&ind,1,1,rayimg);
   fwrite(&xcnt,4,1,rayimg);
   fwrite(&ycnt,4,1,rayimg);
   fwrite(&pixrays,4,1,rayimg);
   fwrite(&tsec,8,1,rayimg);
   fflush(rayimg);
   stop=1;
}


int intloopflag;


void intset(sig)
int sig;
{
   intloopflag=1;
}


void alarmset(sig)
int sig;
{
   intloopflag=2;
}


void intfunc(sig)
int sig;
{
   char str[80];
   register int sec,t,pix;
   struct timeval tv;

   gettimeofday(&tv,NULL);
   tsec+=(tv.tv_sec-ssec)+.000001*(tv.tv_usec-susec);
   pix=xcnt+ycnt*xsize;
   fprintf(stderr, "*** on pixel (%d, %d), %lf seconds, %lf pixels/sec.\n",xcnt,ycnt,tsec,
    pix/tsec);
   sec=(xsize*ysize-pix)*tsec/pix;
   fprintf(stderr, "*** estimated time to completion: %02d:%02d:%02d\n",sec/3600,
    (sec/60)%60,sec%60);
   fprintf(stderr, "*** press ctrl-c within five seconds to interrupt ray6.\n");
   alarm(5);
   signal(SIGALRM,alarmset);
   signal(SIGINT,intset);
   intloopflag=0;
   while(!intloopflag);
   signal(SIGALRM,SIG_DFL);
   alarm(0);
   if(intloopflag==1){
      programstop();
      return;
   }else
      printf("*** continuing with ray tracing\n");
   signal(SIGINT,intfunc);
   gettimeofday(&tv,NULL);
   ssec=tv.tv_sec;
   susec=tv.tv_usec;
   return;
}


int ray6terminate(sig)
int sig;
{
   printf("*** terminated with a signal 15.\n");
   printf("*** leaving finishable image\n");
   programstop();
   return 0;
}


char myerr[80];


void buserror(sig)
int sig;
{
   struct timeval tv;

   fprintf(stderr,"*** bus error...\n");
   fprintf(stderr,"*** on pixel (%d, %d)\n",xcnt,ycnt);
   fprintf(stderr,"*** (execution must end here)\n");
   gettimeofday(&tv,NULL);
   tsec+=(tv.tv_sec-ssec)+.000001*(tv.tv_usec-susec);
   programstop();
   fclose(rayimg);
   exit(1);
}


void memfault(sig)
int sig;
{
   struct timeval tv;

   fprintf(stderr,"*** memory fault...\n");
   fprintf(stderr,"*** on pixel (%d, %d)\n",xcnt,ycnt);
   fprintf(stderr,"*** (execution will end here)\n");
   gettimeofday(&tv,NULL);
   tsec+=(tv.tv_sec-ssec)+.000001*(tv.tv_usec-susec);
   programstop();
   fclose(rayimg);
   exit(1);
}


void fpefunc(sig)
int sig;
{
   struct timeval tv;

   fprintf(stderr,"*** floating exception...\n");
   fprintf(stderr,"*** indicator: '%s'\n",myerr);
   fprintf(stderr,"*** on pixel (%d, %d)   (execution will end here)\n",xcnt,ycnt);
   gettimeofday(&tv,NULL);
   tsec+=(tv.tv_sec-ssec)+.000001*(tv.tv_usec-susec);
   programstop();
}


#if 0
int matherr(x)
struct exception *x;
{
   fprintf(stderr,"*** math library error, type %d, func '%s'\n",x->type,x->name);
   fpefunc(SIGFPE);
}
#endif


/*----------------------------------------------------------------------------
   PURPOSE:   this function performs the creation of the rays for the ray-
    tracing function, and writing the image file
   PRECONDITIONS:   scene points to the valid scene data to ray-trace, imgfile
    points to the open file to store the image in, xpix and ypix are the x and
    y dimensions of the image, respectively, and rpix is the number of rays per
    pixel.
   POSTCONDITIONS:   the image is ray-traced and written to the file.
----------------------------------------------------------------------------*/
int imagetrace(scene,imgfile,xpix,ypix,rpix,xplace,yplace)
scenetype *scene;                       /* scene to be ray-traced */
FILE *imgfile;                          /* file to hold image data */
int xpix,ypix,                          /* x pixels, y pixels */ 
 rpix,                                  /* rays per pixel */
 xplace,yplace;                         /* x and y starting pixels */
{
   mat4type iviewing;                   /* inverse viewing matrix */
   double field;                        /* field of view size (or 1/2) */
   register double xinc,yinc,           /* x and y increments */
    xloop,yloop,                        /* x and y loop counters per coord */
    colfrac;                            /* color fraction per ray */
   register int raycnt;                 /* x and y and rays counters */
   vec4type pixray[2],                  /* pixel ray */
    spaceray[2];                        /* space ray, xformed out of image */
                                        /*  space */
   colortype pixelcol,raycol;           /* pixel color and ray color;  ray */
                                        /*  colors get averaged into pixel */
                                        /*  color */
   unsigned char rgbc[3];               /* pixel r,g,b data for image file */
   double drand48();
   struct timeval tv;

   gettimeofday(&tv,NULL);
   ssec=tv.tv_sec;
   susec=tv.tv_usec;
   ysize=ypix;
   xsize=xpix;
   pixrays=rpix;
   rayimg=imgfile;
   signal(SIGINT,intfunc);
   sprintf(myerr,"none");
   signal(SIGFPE,fpefunc);
   signal(SIGSEGV,memfault);
   signal(SIGBUS,buserror);
   signal(SIGTERM,ray6terminate);
   pixray[0].a[0]=pixray[0].a[1]=pixray[0].a[2]=0.0;
   pixray[0].a[0]=pixray[0].a[1]=pixray[0].a[2]=0.0;
   pixray[0].a[0]=pixray[0].a[1]=pixray[0].a[2]=0.0;
   pixray[0].a[3]=1.0;
   orthotoskew(&iviewing,&scene->eye,&scene->up,scene->uvec,&scene->targ,
    scene->tvec);
   /*printf("iviewing:\n");*/
   /*writemat(stdout,&iviewing);*/
   field=tan(scene->fov/2.0/180.0*M_PI);
   yinc=field*2.0/ypix;
   xinc=field*2.0/xpix;
   /*printf("field: %lf\n",field);*/
#if PICFORMAT
   if((xplace==0)&&(yplace==0)){
      if(fwrite(&xpix,4,1,imgfile)!=1){
         fprintf(stderr,"*** ERROR: cannot write x dimension to image file\n");
         return(0);
      }
      if(fwrite(&ypix,4,1,imgfile)!=1){
         fprintf(stderr,"*** ERROR: cannot write y dimension to image file\n");
         return(0);
      }
   }
#else /* PPM6 */
	fprintf(imgfile, "P6 %d %d 255 ", xpix, ypix);
#endif

   for(ycnt=yplace,yloop=field;ycnt<ypix;ycnt++,yloop-=yinc){
      for(xcnt=xplace,xloop=field;xcnt<xpix;xcnt++,xloop-=xinc){
         pixray[1].a[0]=xloop;
         pixray[1].a[1]=yloop;
         pixray[1].a[2]=pixray[1].a[3]=1.0;
         normalize(pixray+1,pixray+1);
         if(!rpix){
            multraymat(pixray,&iviewing,spaceray);
            raytrace(spaceray,scene,&pixelcol);
         }else{
            pixelcol=black_color;
            colfrac=1.0/rpix;
            for(raycnt=0;raycnt<rpix;raycnt++){
               spaceray[0]=pixray[0];
               spaceray[1]=pixray[1];
               spaceray[1].a[0]+=(drand48()-0.5)*xinc;
               spaceray[1].a[1]+=(drand48()-0.5)*yinc;
               /*translvec((drand48()-0.5)*xinc,(drand48()-0.5)*yinc,0.0,pixray+1,*/
                /*spaceray+1);*/
               multraymat(spaceray,&iviewing,spaceray);
               raytrace(spaceray,scene,&raycol);
               addcolors(&pixelcol,scalecolor(&raycol,colfrac,&raycol),
                &pixelcol);
            }
         }
         if(stop)
            return(0);
         rgbc[0]=(unsigned char)(pixelcol.r*255.0);
         rgbc[1]=(unsigned char)(pixelcol.g*255.0);
         rgbc[2]=(unsigned char)(pixelcol.b*255.0);
         if(fwrite(rgbc,3,1,imgfile)!=1){
            fprintf(stderr,"*** ERROR: cannot write pixel (%d,%d) "
		"to image file\n", xcnt,ycnt);
            return(0);
         }
      }
      xplace=0;
   }
   gettimeofday(&tv,NULL);
   tsec+=(tv.tv_sec-ssec)+.000001*(tv.tv_usec-susec);
   return(1);
}


/*----------*/
/*  main()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function (which is also main()) takes care of initial
    comm-line arg checking, reads the scene from its file, pre-processes the
    data, and calls imagetrace() to do the ray-tracing.
   PRECONDITIONS:   argc contains # of comm-line args, argv points to pointers
    to comm-line args
   POSTCONDITIONS:   the ray-traced image is in a file named after the second
    argument, with x and y dimensions after arg 3 and 4, and was traced with
    arg 5 rays per pixel.
----------------------------------------------------------------------------*/
main(argc,argv)
int argc;                               /* # of command line args */
char *argv[];                           /* command line arg strings */
{
   int xpix=400,ypix=400,pixrays=16,    /* image parameters and defaults */
    xcnt=0,ycnt=0,                      /* xcnt and ycnt from unfinished .pic */
    strcnt,
    midstart=0;                         /* starting in the middle? */
   scenetype scene;                     /* scene data structure */
   FILE *inf,                           /* input file for scene */
    *outf;                              /* output file for images */
   unsigned char inbyte;
   
   if(argc<3){
      fprintf(stderr,"bad args: '%s scenefile imagefile [xpix [ypix [pixrays]]]'\n",
       argv[0]);
      exit(0);
   }
   if(!strcmp(argv[1],"-f")){
      midstart=1;
      if((outf=fopen(argv[2],"r+"))==NULL){
         fprintf(stderr,"*** ERROR: could not open image file for completion\n");
         exit(0);
      }
      if(fread(&inbyte,1,1,outf)!=1){
         fprintf(stderr,"*** ERROR: could not read marker byte from uncompleted\n");
         fprintf(stderr,"***  image file\n");
         exit(0);
      }
      if(inbyte!=255){
         fprintf(stderr,"*** ERROR: marker byte wrong value for unfinished image\n");
         fprintf("%d\n",(int)(inbyte-255));
         exit(0);
      }
      fseek(outf,0,0);
      inbyte=0;
      fwrite(&inbyte,1,1,outf);         /* fix that marker byte */
      fseek(outf,0,0);
      fread(&xpix,4,1,outf);            /* get image size */
      fread(&ypix,4,1,outf);
                                        /* someday I'll actually put error */
                                        /*  checking in here... */
      fseek(outf,-21,2);                /* go to the end of the file to read */
      fread(&inbyte,1,1,outf);          /* length of scenename */
      fread(&xcnt,4,1,outf);            /* x and y pixel counters */
      fread(&ycnt,4,1,outf);
      fread(&pixrays,4,1,outf);         /* rays per pixel */
      fread(&tsec,8,1,outf);            /* total seconds of rendering */
      fseek(outf,-21-inbyte,2);         /* reset to get scenename */
      fread(scenename,(int)inbyte,1,outf);
      scenename[(int)inbyte]='\0';
                                        /* get scene name */
      fseek(outf,-21-inbyte,2);         /* reset for imagetrace */
   }else
      strcpy(scenename,argv[1]);
   if((inf=fopen(scenename,"r"))==NULL){
      fprintf(stderr,"*** ERROR: could not open scene file for input.\n");
      exit(0);
   }
   if(!getscene(inf,&scene)){
      fprintf(stderr,"*** ERROR: could not read scene from file.\n");
      exit(0);
   }
   fclose(inf);

   /*pre-process objects (bboxes, subdivision, etc...)*/

   if(!midstart){
      if(argc>3)
         if(sscanf(argv[3],"%d",&xpix)!=1){
            fprintf(stderr,"bad args: '%d scenefile xpix ypix pixrays'\n",argv[0]);
            /*free(scene->objs);*/
            /*free(scene->lights);*/
            exit(0);
         }
      if(argc>4)
         if(sscanf(argv[4],"%d",&ypix)!=1){
            fprintf(stderr,"bad args: '%d scenefile xpix ypix pixrays'\n",argv[0]);
            /*free(scene->objs);*/
            /*free(scene->lights);*/
            exit(0);
         }
      if(argc>5)
         if(sscanf(argv[5],"%d",&pixrays)!=1){
            fprintf(stderr,"bad args: '%d scenefile xpix ypix pixrays'\n",argv[0]);
            /*free(scene->objs);*/
            /*free(scene->lights);*/
            exit(0);
         }
      if((outf=fopen(argv[2],"w"))==NULL){
         fprintf(stderr,"*** ERROR: cannot open image file '%s' for output.\n",argv[1]);
         exit(0);
      }
   }
   if(!imagetrace(&scene,outf,xpix,ypix,pixrays,xcnt,ycnt))
      fprintf(stderr,"*** WARNING: image not completed\n");
   else{
      printf("*** Rendering Statistics: %d rays computed, %d pixels shaded\n",
       pixrays*xpix*ypix,ypix*xpix);
      printf("***  rendering time: %02d:%02d:%02lf, rays/sec: %lf, pixels/sec: %lf\n",
       ((int)tsec)/3600,(((int)tsec)/60)%60,tsec-((int)tsec/60)*60,pixrays*xpix*ypix/
       tsec,xpix*ypix/tsec);
   }
   fclose(outf);
}
