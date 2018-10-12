/*----------------------------------------------------------------------------

   Written by Brad Grantham

   opened:   12/26/89
   last modification:   01/09/90

   FILE:       grf_math.h
   LANGUAGE:   C
   SYSTEM:     Mac II, under A/UX
   PURPOSE:    contains definitions and some declarations for use with the
    grf_math.o object and grf_math.c code

----------------------------------------------------------------------------*/


#ifndef GRF_MATH_H
#define GRF_MATH_H


#include <stdio.h>
#include <math.h>
#include "brads.h"

   
/*----------------*/
/*  #define sq()  */
/*----------------*/
#define sq(a) ((a)*(a))                 /* square macro: warning, ARGUMENT */
                                        /*  EVALUATED TWICE */


/*-----------------*/
/*  #define DINKY  */
/*-----------------*/
#define DINKY (0.0000000001)            /* really small number */


/*----------------*/
/*  typedef vec4  */
/*----------------*/
typedef struct{
   double a[4];                          /* 4 element vector */
} vec4type;


/*-----------------*/
/*  typedef mat4  */
/*-----------------*/
typedef struct{
   double a[4][4];                      /* 4 by 4 matrix (transformations) */
                                        /* my convention is to use form */
                                        /*  mat[row][col] */
} mat4type;


/*--------------------*/
/*  typedef vec2type  */
/*--------------------*/
typedef struct {                        /* 2 element vector */
   double a[2];                         /* perhaps screen coordinates before */
                                        /*  being rounded to integers */
} vec2type;


/*-----------------------------------------*/
/*  struct raytypestruct; typedef raytype  */
/*-----------------------------------------*/
typedef struct raytypestruct{           /* ray type */ 
   vec4type o,                          /* origin of ray */
    d;                                  /* direction of ray */
} raytype;


/*------------------------*/
/*  mat4type identmat4  */
/*------------------------*/
extern mat4type identmat4;            /* identity matrix */


#define normvec(v) {register double n_z=(v).a[3];(v).a[0]/=n_z;(v).a[1]/=n_z;\
 (v).a[2]/=n_z;v.a[3]=1.0;}

#define homogblock(v,n) {register int h_c;for(h_c=0;h_c<(n);h_c++)normvec((v)[h_c]);};


#define CLIP_NONE 0                     /* defines for clipsegment routine */
#define CLIP_WHOLE 1
#define CLIP_LOSTP2 2
#define CLIP_LOSTP1 3


vec4type *addvec(/*vec1,vec2,sumvec*/);
vec4type *subvec(/*vec1,vec2,vec3*/);
vec4type *negate(/*vec,negvec*/);
vec4type *sizevec(/*vec,factor,scaled*/);
vec4type *multvecmat(/*vec,mat,prod*/);
vec4type *multvmblock(/*vec,mat,prod,total*/);
mat4type *multmat4(/*mat1,mat2,prod*/);
mat4type *invertmat(/*mat,inv*/);
mat4type *transposemat(/*mat,trans*/);
vec4type *normalize(/*vec1,vec2*/);
vec4type *reflect(/*inc,norm,ref*/);
double dot(/*vec1,vec2*/);
vec4type *cross(/*vec1,vec2,ortho*/);
double magnitude(/*vector*/);
double disttoray(/*point,ray*/);
double disttoplane(/*point,ray*/);
vec4type *angtovec(/*alpha,beta,vector*/);
vec4type *pointonray(/*point,ray,proj*/);
vec4type *pointonplane(/*point,ray,proj*/);
vec4type *multraymat(/*ray,mat,prod*/);
vec4type *multrmblock(/*rays,mat,prods,total*/);
vec4type *rotvec(/*angle,axis,vec,result*/);
mat4type *rotmat(/*angle,axis,mat,result*/);
vec4type rotvecvec(/*angle,vec,src,res*/);
mat4type rotmatvec(/*angle,vec,src,res*/);
vec4type *translvec(/*x,y,z,vec,res*/);
mat4type *translmat(/*x,y,z,mat,res*/);
vec4type *scalevec(/*x,y,z,vec,res*/);
mat4type *scalemat(/*x,y,z,mat,res*/);
mat4type *skewtoortho(/*mat,src,ypt,yvec,zpt,zvec*/);
mat4type *orthotoskew(/*mat,src,ypt,yvec,zpt,zvec*/);
int writemat(/*outf,matrix*/);
int readmat(/*outf,matrix*/);
int writevec(/*outf,vec*/);
int readvec(/*outf,vec*/);
int writevec3(/*outf,vec*/);
int readvec3(/*outf,vec*/);
double angletween(/*vec1,vec2*/);


#endif
