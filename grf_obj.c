/*----------------------------------------------------------------------------
   Written by Brad Grantham
   FILE:   grf_obj.h 
   PURPOSE:   This file contains the object manipulation functions.  These 
    routines handle intersections, object processing, and manipulation of 
    object lists.

----------------------------------------------------------------------------*/


#include <stdio.h>
#include "grf_obj.h"


int getscene(infile,srec,errstr)
FILE *infile;
scenetype *srec;
char *errstr;
{
   int objcnt,
    plcnt;
   double red,green,blue;
   objtype *tobj;
   char code[5];
   int (* printfun)(); /*  sprintf(),fprintf(); */
   char *outobj;

   if(errstr){
      printfun=sprintf;
      outobj=(char *)errstr;
   }else{
      printfun=fprintf;
      outobj=(char *)stderr;
   }
   if(fgets(code,5,infile)==NULL){
      printfun(outobj,"error reading 4-character id code\n");
      return(0);
   }
   if(strcmp(code,"gdf1")){
      printfun(outobj,"don't recognize id code ('%s')\n",code);
      return(0);
   }
   if(!readvec3(infile,&srec->eye)){
      printfun(outobj,"error reading eye vector\n");
      return(0);
   }
   if(!readvec3(infile,&srec->targ)){
      printfun(outobj,"error reading target vector\n");
      return(0);
   }
   if(fscanf(infile,"%d",&srec->tvec)!=1){
      printfun(outobj,"error reading target vec/pt flag\n");
      return(0);
   }
   if(!readvec3(infile,&srec->up)){
      printfun(outobj,"error reading up vector\n");
      return(0);
   }
   if(fscanf(infile,"%d",&srec->uvec)!=1){
      printfun(outobj,"error reading up vec/pt flag\n");
      return(0);
   }
   if(fscanf(infile,"%lf %lf %lf",&srec->fov,&srec->near,&srec->far)!=3){
      printfun(outobj,"error reading fov, near, or far\n");
      return(0);
   }
   if(fscanf(infile,"%d",&srec->numobj)!=1){
      printfun(outobj,"error reading number of objects\n");
      return(0);
   }
   if((srec->objs=getsome(objtype,srec->numobj))==NULL){
      printfun(outobj,"not enough memory for objects\n");
      return(0);
   }
   for(objcnt=0;objcnt<srec->numobj;objcnt++){
      tobj=srec->objs+objcnt;
      if(fscanf(infile,"%s",tobj->name)!=1){
         printfun(outobj,"object #%d : can't read name string\n",objcnt);
         dispose(srec->objs);
         return(0);
      }
      if(fscanf(infile,"%d",&tobj->objtag)!=1){
         printfun(outobj,"object #%d ('%s'): can't read type\n",
          objcnt,tobj->name);
         dispose(srec->objs);
         return(0);
      }
      switch(tobj->objtag){
        case SPHERE:
        case CYLINDER:
        case DISC:
        case CONE:
         break;
        case TORUS:
         if(fscanf(infile,"%lf",tobj->lesser)!=1){
            printfun(outobj,"object #%d ('%s'): can't read lesser radius\n",
             objcnt,tobj->name);
            dispose(srec->objs);
            return(0);
         }
         break;
        case RING:
         if(fscanf(infile,"%lf",tobj->inner)!=1){
            printfun(outobj,"object #%d ('%s'): can't read inner radius\n",
             objcnt,tobj->name);
            dispose(srec->objs);
            return(0);
         }
         break;
        case TRIANGLE:
         if(!readvec3(infile,&tobj->pts[0])){
            printfun(outobj,"object #%d ('%s'): can't read first vertex\n",
             objcnt,tobj->name);
            dispose(srec->objs);
            return(0);
         }
         if(!readvec3(infile,&tobj->pts[1])){
            printfun(outobj,"object #%d ('%s'): can't read second vertex\n",
             objcnt,tobj->name);
            dispose(srec->objs);
            return(0);
         }
         if(!readvec3(infile,&tobj->pts[2])){
            printfun(outobj,"object #%d ('%s'): can't read third vertex\n",
             objcnt,tobj->name);
            dispose(srec->objs);
            return(0);
         }
         break;
      }
      if(!readmat(infile,&tobj->xform)){
         printfun(outobj,"object #%d ('%s'): can't read xform matrix\n",
          objcnt,tobj->name);
         dispose(srec->objs);
         return(0);
      }
      if(!invertmat(&tobj->xform,&tobj->ixform)){
         printfun(outobj,"object $%d ('%s'): xform has no inverse\n",objcnt,
          tobj->name);
         tobj->ixform=identmat4;
      }
      if(fscanf(infile,"%lf",&red)!=1){
         printfun(outobj,"object #%d ('%s'): can't read red/texmap-flag\n",
          objcnt,tobj->name);
         dispose(srec->objs);
         return(0);
      }
      if(red!=(-1.0)){
         if(fscanf(infile,"%lf %lf",&green,&blue)!=2){
            printfun(outobj,"object #%d ('%s'): can't read green or blue\n",
             objcnt,tobj->name);
            dispose(srec->objs);
            return(0);
         }
         /*if((red>1.0)||(green<0.0)||(green>1.0)||(blue<0.0)||(blue>1.0)){*/
            /*printfun(outobj,"object #%d ('%s'): r,g,b values not in [0.0,1.0]\n",*/
             /*objcnt,tobj->name);*/
            /*dispose(srec->objs);*/
            /*return(0);*/
         /*}*/
         tobj->color.r=red;
         tobj->color.g=green;
         tobj->color.b=blue;
         tobj->texmap[0]=0;
      }else{
         if(fscanf(infile,"%s",tobj->texmap)!=1){
            printfun(outobj,"object #%d ('%s'): can't read texmap filename\n",
             objcnt,tobj->name);
            dispose(srec->objs);
            return(0);
         }
         /*----------------------------------------*/
         /* read texture mapping parameters here...*/
         /*----------------------------------------*/
      }
      if(fscanf(infile,"%lf",&tobj->trans)!=1){
         printfun(outobj,"object #%d ('%s'): can't read transparency\n",
          objcnt,tobj->name);
         dispose(srec->objs);
         return(0);
      }
      if(fscanf(infile,"%lf",&tobj->refl)!=1){
         printfun(outobj,"object #%d ('%s'): can't read reflectivity\n",
          objcnt,tobj->name);
         dispose(srec->objs);
         return(0);
      }
      if(fscanf(infile,"%d",&tobj->radiant)!=1){
         printfun(outobj,"object #%d ('%s'): can't read radiance flag\n",
          objcnt,tobj->name);
         dispose(srec->objs);
         return(0);
      }
   }
   return(1);
}


void predict(obj,npts,nlns)
objtype *obj;
int *npts,
 *nlns;
{
   switch(obj->objtag){
     case SPHERE: 
      *npts=42;
      *nlns=51;
      break;
     case DISC: 
      *npts=16;
      *nlns=16;
      break;
     case CYLINDER: 
      *npts=32;
      *nlns=36;
      break;
     case CONE: 
      *npts=17;
      *nlns=20;
      break;
     case RING: 
      *npts=32;
      *nlns=36;
      break;
     case TORUS: 
      *npts=80;
      *nlns=80;
      break;
     case TRIANGLE: 
      *npts=3;
      *nlns=3;
      break;
   }
}


#include "sphere.dat"
#include "disc.dat"
#include "cylinder.dat"
#include "cone.dat"
#include "ring.dat"
#include "torus.dat"


void objtolines(obj,points,lines,ptcnt)
objtype *obj;
vec4type points[];
int lines[][2],
 ptcnt;
{
   register int lncnt;
   mat4type xform;

   switch(obj->objtag){
     case SPHERE:
      multvmblock(sphere_pts,&obj->xform,points,42);
      for(lncnt=0;lncnt<51;lncnt++){
         lines[lncnt][0]=ptcnt+sphere_lns[lncnt][0];
         lines[lncnt][1]=ptcnt+sphere_lns[lncnt][1];
      }
      break;
     case DISC:
      multvmblock(disc_pts,&obj->xform,points,16);
      for(lncnt=0;lncnt<18;lncnt++){
         lines[lncnt][0]=ptcnt+disc_lns[lncnt][0];
         lines[lncnt][1]=ptcnt+disc_lns[lncnt][1];
      }
      break;
     case CYLINDER:
      multvmblock(cylinder_pts,&obj->xform,points,32);
      for(lncnt=0;lncnt<36;lncnt++){
         lines[lncnt][0]=ptcnt+cylinder_lns[lncnt][0];
         lines[lncnt][1]=ptcnt+cylinder_lns[lncnt][1];
      }
      break;
     case CONE:
      multvmblock(cone_pts,&obj->xform,points,17);
      for(lncnt=0;lncnt<20;lncnt++){
         lines[lncnt][0]=ptcnt+cone_lns[lncnt][0];
         lines[lncnt][1]=ptcnt+cone_lns[lncnt][1];
      }
      break;
     case RING:
      scalemat(obj->inner,1.0,obj->inner,&identmat4,&xform);
      multmat4(&xform,&obj->xform,&xform);
      multvmblock(ring_pts,&obj->xform,points,16);
      multvmblock(ring_pts,&xform,points+16,16);
      for(lncnt=0;lncnt<36;lncnt++){
         lines[lncnt][0]=ptcnt+ring_lns[lncnt][0];
         lines[lncnt][1]=ptcnt+ring_lns[lncnt][1];
      }
      break;
     case TORUS:
      multvmblock(torus_pts,&obj->xform,points,16);

      scalemat(obj->lesser,1.0,obj->lesser,&identmat4,&xform);
      rotmat(M_PI/2.0,'X',&xform,&xform);
      translmat(1.0,0.0,0.0,&xform,&xform);
      multmat4(&xform,&obj->xform,&xform);
      multvmblock(torus_pts,&xform,points+16,16);

      scalemat(obj->lesser,1.0,obj->lesser,&identmat4,&xform);
      rotmat(M_PI/2.0,'X',&xform,&xform);
      translmat(-1.0,0.0,0.0,&xform,&xform);
      multmat4(&xform,&obj->xform,&xform);
      multvmblock(torus_pts,&xform,points+32,16);

      scalemat(obj->lesser,1.0,obj->lesser,&identmat4,&xform);
      rotmat(M_PI/2.0,'Z',&xform,&xform);
      translmat(0.0,0.0,1.0,&xform,&xform);
      multmat4(&xform,&obj->xform,&xform);
      multvmblock(torus_pts,&xform,points+48,16);

      scalemat(obj->lesser,1.0,obj->lesser,&identmat4,&xform);
      rotmat(M_PI/2.0,'Z',&xform,&xform);
      translmat(0.0,0.0,1.0,&xform,&xform);
      multmat4(&xform,&obj->xform,&xform);
      multvmblock(torus_pts,&xform,points+64,16);

      for(lncnt=0;lncnt<80;lncnt++){
         lines[lncnt][0]=ptcnt+torus_lns[lncnt][0];
         lines[lncnt][1]=ptcnt+torus_lns[lncnt][1];
      }
      break;
     case TRIANGLE:
      multvecmat(obj->pts,&obj->xform,points);
      multvecmat(obj->pts+1,&obj->xform,points+1);
      multvecmat(obj->pts+2,&obj->xform,points+2);
      lines[0][0]=ptcnt+0;   lines[0][1]=ptcnt+1;
      lines[1][0]=ptcnt+1;   lines[1][1]=ptcnt+2;
      lines[2][0]=ptcnt+2;   lines[2][1]=ptcnt+0;
      break;
   }
   return;
}


double raytoxz(ray)
vec4type ray[2];
{
   register double t;

   if(fabs(ray[1].a[1])<0.00000001)
      return(-1.0);
   t=ray[0].a[1]/ray[1].a[1];
   if(t<0.0)
      return(-1.0);
   return(t);
}


objtype *rayintobj(lastobj,lastpnt,ray,scene,ipnt)
objtype *lastobj;
vec4type lastpnt[2],
 ray[2];
scenetype *scene;
vec4type ipnt[2];
{
   register objtype *intrsct=NULL,*tmp;
   vec4type newray[2],intpnt,tmpvec;
   register int objcnt;
   register double tval=-1.0,t,mag;

   for(tmp=scene->objs,objcnt=0;objcnt<scene->numobj;tmp++,objcnt++){
      multraymat(ray,&tmp->ixform,newray);
      mag=magnitude(newray+1);
      normalize(newray+1,newray+1);
      t=-1.0;
      if(tmp==lastobj)
         addvec(newray,sizevec(newray+1,.00001,&tmpvec),newray);
      switch(tmp->objtag){
        case SPHERE:{
         register double m,n,o,b,root,t1,t2;

         m=newray[0].a[0];
         n=newray[0].a[1];
         o=newray[0].a[2];
         b=2.0*(newray[1].a[0]*m+newray[1].a[1]*n+newray[1].a[2]*o);
         root=b*b-4.0*(m*m+n*n+o*o-1.0);
         if(root>0.00000001){
            root=sqrt(root);
            t1=(-b-root)/2.0;
            t2=(-b+root)/2.0;
            if((t1>0.0)||(t2>0.0))
               t=(t1<0.0)?t2:t1;

            /*if(-root>b)
               t=(-b-root)/2.0;
            else
               if(root>b)
                  t=(-b+root)/2.0;*/
         }
         break;
        }
        case DISC:
         if((t=raytoxz(newray))>0.0){
            addvec(newray,sizevec(newray+1,t,&intpnt),&intpnt);
            if(intpnt.a[0]*intpnt.a[0]+intpnt.a[1]*intpnt.a[1]+intpnt.a[2]*
             intpnt.a[2]>1.0)
               t=-1.0;
         }else
            t=-1.0;
         break;
        case CYLINDER:{
         register double a,b,c,root, t1, t2;

         a=newray[1].a[0]*newray[1].a[0]+newray[1].a[2]*newray[1].a[2];
         b=2.0*(newray[0].a[0]*newray[1].a[0]+newray[0].a[2]*newray[1].a[2]);
         c=newray[0].a[0]*newray[0].a[0]+newray[0].a[2]*newray[0].a[2]-1.0;
         if((root=b*b-4.0*a*c)>0.0){
            root=sqrt(root);
            t1=(-b-root)/(2.0*a);
            t2=(-b+root)/(2.0*a);
            t = t1;
            if((t1>0.0)&&(fabs(newray[0].a[1]+t1*newray[1].a[1]-0.5)<=0.5))
               t=t1;
            else
               if((t2>0.0)&&(fabs(newray[0].a[1]+t2*newray[1].a[1]-0.5)<=0.5))
                  t = t2;
               else 
                  t = -1.0;
         }
         break;
        }
        case CONE: /*?*/
         t = -1.0;
         break;
        case RING:{
         register double d;

         if((t=raytoxz(newray))>0.0){
            addvec(newray,sizevec(newray+1,t,&intpnt),&intpnt);
            d=intpnt.a[0]*intpnt.a[0]+intpnt.a[1]*intpnt.a[1]+intpnt.a[2]*
             intpnt.a[2];
            if((d>1.0)||(d<tmp->inner*tmp->inner))
               t=-1.0;
         }else
            t=-1.0;
         break;
        }
        case TORUS: /*?*/
         t = -1.0;
         break;
    
        case TRIANGLE:{
         register double da,db,dc;
         vec4type ma[2],mb[2],mc[2];

         /* don't forget; triangle has three arbitrary viewpoints, and */
         /*  THEN a xform matrix. */

	 /* this code is wrong */
         if((t=raytoxz(newray))>0.0){
            addvec(newray,sizevec(newray+1,t,&intpnt),&intpnt);
            sizevec(addvec(tmp->pts,tmp->pts+1,ma),0.5,ma);
            sizevec(addvec(tmp->pts+1,tmp->pts+2,mb),0.5,mb);
            sizevec(addvec(tmp->pts+2,tmp->pts,mc),0.5,mc);
            normalize(subvec(mb,ma,ma+1),ma+1);
            normalize(subvec(mc,mb,mb+1),mb+1);
            normalize(subvec(ma,mc,mc+1),mc+1);
            da=disttoray(tmp->pts+2,ma);
            db=disttoray(tmp->pts,mb);
            dc=disttoray(tmp->pts+1,mc);
            if((disttoray(&intpnt,ma)>da)||(disttoray(&intpnt,mb)>db)||(
             disttoray(&intpnt,mc)>dc))
               t=-1.0;
         }else
            t=-1.0;
         break;
        }
      }
      t/=mag;
      if((t>0.0)&&((tval==-1.0)||(t<tval))){
         tval=t;
         intrsct=tmp;
      }
   }
   if(intrsct){
      addvec(ray,sizevec(ray+1,tval,&intpnt),ipnt);
      multvecmat(ipnt,&intrsct->ixform,ipnt+1);
   }
   return(intrsct);
}


colortype *getobjcolor(obj,ipnt,color)
objtype *obj;
vec4type ipnt[2];
colortype *color;
{
   (*color)=obj->color;
   return(color);
}


vec4type upnorm={{0.0,1.0,0.0,1.0}};


vec4type *bounceray(ray,obj,ipnt,refl)
vec4type ray[2];
objtype *obj;
vec4type ipnt[2],
 refl[2];
{
   vec4type inc[2],
    norm,
    ref[2];

   multraymat(ray,&obj->ixform,inc);
   normalize(inc+1,inc+1);
   switch(obj->objtag){
     case SPHERE:
      norm=ipnt[1];
      break;
     case DISC:
      norm = upnorm;
      break;
     case RING:
      norm = upnorm;
      break;
     case TRIANGLE:
      norm = upnorm;
      break;
     case CYLINDER:
      norm.a[0] = ipnt[1].a[0];
      norm.a[1] = 0.0;
      norm.a[2] = ipnt[1].a[2];
      norm.a[3] = 1.0;
     case CONE: /*?*/
      break;
     case TORUS: /*?*/
      break;
   }
   ref[0]=ipnt[1];
   reflect(inc+1,&norm,ref+1);
   multraymat(ref,&obj->xform,refl);
   normalize(refl+1,refl+1);
   return(refl);
}
