/*----------------------------------------------------------------------------
  
   r6obj.h 
   Created by Brad Grantham

   This file contains the object manipulation functions.  These routines
    handle intersections, object processing, and manipulation of object lists.

----------------------------------------------------------------------------*/


#ifndef GRF_OBJ_H
#define GRF_OBJ_H


#include "grf_math.h"
#include "grf_color.h"


/*-------------------*/
/*  typedef objtype  */
/*-------------------*/
typedef struct {                        /* object structure */
   char name[40];
   int objtag;                          /* type of object */
      vec4type pts[3];                  /* triangle points */
      double lesser;                    /* size of lesser radius in torus */
                                        /*  with greater radius of one unit */
      double inner;                     /* inner radius of ring with outer */
                                        /*  radius of one unit */
   mat4type xform,                      /* transformation matrix */
    ixform;                             /* inverse of xform */
   colortype color;                     /* object color */
   char texmap[80];                     /* if len>0, get color from tex map */
         /* put mapping parameters here */
   double trans,                        /* ratio total transparency */
    refl;                               /* ratio total reflectivity */
   int radiant;                         /* does object radiate light? */
         /* optimization and computed parameters here */
   vec4type C;                          /* center of bounding sphere */
   double R;                            /* radius of bounding sphere */
} objtype;


/*----------------------*/
/*  define object types */
/*----------------------*/
#define SPHERE 0
#define DISC 1
#define CYLINDER 2
#define CONE 3
#define RING 4
#define TORUS 5
#define TRIANGLE 6
#define PATCH 7                         /* what should this be? */
#define LATHING 8                       /* this is a collection of cone sections */
#define EXTENSION 9                     /* this is a line extension into space */


typedef struct {
   vec4type eye,targ,up;
   int tvec,uvec;
   double fov,near,far;
   int numobj,numpl;
   objtype *objs;
} scenetype;


int getscene(/*infile,srec*/);
void predict(/*obj,npts,nlns*/);
void objtolines(/*obj,points,lines*/);
objtype *rayintobj(/*lastobj,lastpnt,ray,scene,ipnt*/);
colortype *getobjcolor(/*obj,ipnt,color*/);
vec4type *bounceray(/*ray,obj,ipnt,refl*/);
scenetype *forraytracing(/*scene*/);


#endif
