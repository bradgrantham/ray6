/*----------------------------------------------------------------------------
   Written By Brad Grantham
   Copyright 1990
   FILE:   grf_math.c
   SYSTEM/COMPILER:   MAC II A/UX using gcc
   DESCRIPTION:   this is the source for my graphics vector and matrix math
    routines.   The rotation rule is that rotation around an axis by a positive
    angle is counter-clockwise if looking from the origin down the positive axis.
    These functions will work if the same vector or matrix pointer is used for
    more than one parameter; i.e. addvec(vec,vec,vec) will add vec to itself
    properly.  These functions are designed to return the result parameter of
    the operation (or NULL if unsuccessful).
----------------------------------------------------------------------------*/


#include "grf_math.h"                /* graphics math header */


static vec4type tempvec;             /* temporary vector for this library */


static mat4type tempmat;             /* temporary matrix for this library */


mat4type identmat4={{{1.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0},{0.0,0.0,1.0,0.0},
    {0.0,0.0,0.0,1.0}}};              /* identity matrix */


/*------------*/
/*  addvec()  */
/*----------------------------------------------------------------------------
   DESCRIPTION:   this function adds two vectors together and stores the sum
    in another vector
   PRECONDITIONS:   vec1 and vec2 must point to the vectors to add, sumvec
    must point to storage to place the sum in.
   POSTCONDITIONS:   the sum of vec1 and vec2 is stored in sumvec
----------------------------------------------------------------------------*/
vec4type *addvec(vec1,vec2,sumvec)
vec4type *vec1,*vec2,                /* vectors to sum */
 *sumvec;                            /* storage for the sum */
{
   sumvec->a[0]=vec1->a[0]+vec2->a[0];
   sumvec->a[1]=vec1->a[1]+vec2->a[1];
   sumvec->a[2]=vec1->a[2]+vec2->a[2];
   sumvec->a[3]=1.0;
   return(sumvec);
}


/*------------*/
/*  subvec()  */
/*----------------------------------------------------------------------------
   DESCRIPTION:   this function subtracts the second vector passed from the
    first vector and returne the result in the third vector.
   PRECONTITIONS:   vec2 points to the vector to subtract from the vector
    pointed to by vec1
   POSTCONDITIONS:   vec3=vec1-vec2
----------------------------------------------------------------------------*/
vec4type *subvec(vec1,vec2,diffvec)
vec4type *vec1,*vec2,                /* vectors to subtract */
 *diffvec;                           /* storage for result */
{
   diffvec->a[0]=vec1->a[0]-vec2->a[0];
   diffvec->a[1]=vec1->a[1]-vec2->a[1];
   diffvec->a[2]=vec1->a[2]-vec2->a[2];
   diffvec->a[3]=1.0;
   return(diffvec);
}


/*---------------*/
/*  negatevec()  */
/*----------------------------------------------------------------------------
   DESCRIPTION:   this function accepts a function, and places its negation in
    another vector.
   PRECONDITIONS:   vec points to the vector to negate, and negvec has storage
    for the result.
   POSTCONDITIONS:   the vector passed comes back with all elements negated
----------------------------------------------------------------------------*/
vec4type *negatevec(vec,negvec)
vec4type *vec,                       /* vector negate */
 *negvec;                            /* result of negation */
{
   negvec->a[0]=(-vec->a[0]);
   negvec->a[1]=(-vec->a[1]);
   negvec->a[2]=(-vec->a[2]);
   negvec->a[3]=vec->a[3];
   return(negvec);
}


/*----------------*/
/*  enlargevec()  */
/*----------------------------------------------------------------------------
   DESCRIPTION:   this function scales a given vector by a factor (all
    components of the vector are scaled)
   PRECONDITIONS:   vec points to the vector to scale, factor is the value to
    scale the vector by, and scaled points to storage for the scaled vector
   POSTCONDITIONS:   scaled equals vec scaled by factor
----------------------------------------------------------------------------*/
vec4type *sizevec(vec,factor,scaled)
vec4type *vec;                          /* vector to scale */
double factor;                          /* factor to scale by  */
vec4type *scaled;                       /* scaled version of vector */
{
   scaled->a[0]=vec->a[0]*factor;
   scaled->a[1]=vec->a[1]*factor;
   scaled->a[2]=vec->a[2]*factor;
   scaled->a[3]=vec->a[3];
   return(scaled);
}


/*----------------*/
/*  multvecmat()  */
/*----------------------------------------------------------------------------
   DESCRIPTION:   this function multiplies a vector by a matrix and returns the
    vector product. (i.e. the vector is transformed by the matrix)
   PRECONDITIONS:   vec points to the vector, mat points to the matrix, prod
    points to storage to hold the result.
   POSTCONDITIONS:   prod contains the vector after being transformed by the
    matrix
----------------------------------------------------------------------------*/
vec4type *multvecmat(vec,mat,prod)
vec4type *vec;                          /* vector to multiply */
mat4type *mat;                          /* matrix to multiply by */
vec4type *prod;                         /* product vector */
{
   register double a,b,c,d;             /* temps to keep array access down */

   /*printf("multvecmat\n");*/
   /*writevec(stdout,vec);*/
   /*writemat(stdout,mat);*/
   a=vec->a[0];
   b=vec->a[1];
   c=vec->a[2];
   d=vec->a[3];
   /*printf("registers : %lf %lf %lf %lf\n",a,b,c,d);*/
   prod->a[0]=a*mat->a[0][0]+b*mat->a[1][0]+c*mat->a[2][0]+d*mat->a[3][0];
   prod->a[1]=a*mat->a[0][1]+b*mat->a[1][1]+c*mat->a[2][1]+d*mat->a[3][1];
   prod->a[2]=a*mat->a[0][2]+b*mat->a[1][2]+c*mat->a[2][2]+d*mat->a[3][2];
   prod->a[3]=a*mat->a[0][3]+b*mat->a[1][3]+c*mat->a[2][3]+d*mat->a[3][3];
   return(prod);
}


/*-----------------*/
/*  multvmblock()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function performs exactly as multvecmat (above), except it
    transforms an array of vectors by a matrix
   PRECONDITIONS:   vec is the array of vectors, mat points to the matrix,
    prod points to enough storage to hold the transformed vectors, and total
    is the number of vectors to transform
   POSTCONDITIONS:   prod holds the transformed versions of the vectors
----------------------------------------------------------------------------*/
vec4type *multvmblock(vec,mat,prod,total)
vec4type vec[];                         /* vector array */
mat4type *mat;                          /* matrix to transform by */
vec4type prod[];                        /* product vector array */
int total;                              /* total vectors to transform */
{
   register double a,b,c,d;             /* temps to keep down array access */
   register int veccnt;                 /* vector counter */

   for(veccnt=0;veccnt<total;veccnt++){
      a=vec[veccnt].a[0];
      b=vec[veccnt].a[1];
      c=vec[veccnt].a[2];
      d=vec[veccnt].a[3];
      prod[veccnt].a[0]=a*mat->a[0][0]+b*mat->a[1][0]+c*mat->a[2][0]+d*mat->a[3][0];
      prod[veccnt].a[1]=a*mat->a[0][1]+b*mat->a[1][1]+c*mat->a[2][1]+d*mat->a[3][1];
      prod[veccnt].a[2]=a*mat->a[0][2]+b*mat->a[1][2]+c*mat->a[2][2]+d*mat->a[3][2];
      prod[veccnt].a[3]=a*mat->a[0][3]+b*mat->a[1][3]+c*mat->a[2][3]+d*mat->a[3][3];
   }
   return(prod);
}


/*---------------*/
/*  multmat4()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function multiplies two 4x4 matrices together
   PRECONDITIONS:   mat1 points to the matrix to multiply by the matrix pointed
    to by mat2 and prod points to storage for the product matrix
   POSTCONDITIONS:   prod contains the product matrix
----------------------------------------------------------------------------*/
mat4type *multmat4(mat1,mat2,prod)
mat4type *mat1,                         /* a         */
 *mat2,                                 /*   * b     */
 *prod;                                 /*       = c */
{
   register int i,j,k;                  /* indices for array traversal */
   mat4type tempmat;

   for(i=0;i<4;i++)
      for(j=0;j<4;j++){
         tempmat.a[i][j]=0;
         for(k=0;k<4;k++){
            tempmat.a[i][j]+=mat1->a[i][k]*mat2->a[k][j];
         }
   }
   *prod=tempmat;
   return(prod);
}


double matdet(mat)
mat4type *mat;
{

   double A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P;

   A=mat->a[0][0];
   B=mat->a[0][1];
   C=mat->a[0][2];
   D=mat->a[0][3];
   E=mat->a[1][0];
   F=mat->a[1][1];
   G=mat->a[1][2];
   H=mat->a[1][3];
   I=mat->a[2][0];
   J=mat->a[2][1];
   K=mat->a[2][2];
   L=mat->a[2][3];
   M=mat->a[3][0];
   N=mat->a[3][1];
   O=mat->a[3][2];
   P=mat->a[3][3];
   return((A*F-B*E)*(K*P-L*O)+
    (C*E-A*G)*(J*P-L*N)+
    (A*H-D*E)*(J*O-K*N)+
    (B*G-C*F)*(I*P-L*M)+
    (D*F-B*H)*(I*O-K*M)+
    (C*H-D*G)*(I*N-J*M));
}


/*----------------------------------------------------------------------------
   PURPOSE:   this function calculates the inverse of a 4x4 matrix
   PRECONDITIONS:   mat points to the matrix to invert, and inv points to 
    storage for the result of the inversion
   POSTCONDITIONS:   inv contains the inversion of mat, or NULL is returned
    to indicate that the matrix has no inverse
----------------------------------------------------------------------------*/
mat4type *invertmat(mat,inv)
mat4type *mat,                          /* matrix to invert */
 *inv;                                  /* inversion */
{
   register int i,j,rswap;
   register double det,div,swap;
   mat4type hold,imat;

   hold=*mat;
   imat=identmat4;
   det=matdet(mat);
   if(fabs(det)<DINKY)
      return(NULL);

   if(fabs(hold.a[0][0])<DINKY){
      if(fabs(hold.a[1][0])>DINKY) rswap=1;
      else if(fabs(hold.a[2][0])>DINKY) rswap=2;
           else if(fabs(hold.a[3][0])>DINKY) rswap=3;
      for(i=0;i<4;i++){
         swap=hold.a[0][i];
         hold.a[0][i]=hold.a[rswap][i];
         hold.a[rswap][i]=swap;
         swap=imat.a[0][i];
         imat.a[0][i]=imat.a[rswap][i];
         imat.a[rswap][i]=swap;
      }
   }
      
   div=hold.a[0][0];
   for(i=0;i<4;i++){
      hold.a[0][i]/=div;
      imat.a[0][i]/=div;
   }

   div=hold.a[1][0];
   for(i=0;i<4;i++){
      hold.a[1][i]-=div*hold.a[0][i];
      imat.a[1][i]-=div*imat.a[0][i];
   }
   div=hold.a[2][0];
   for(i=0;i<4;i++){
      hold.a[2][i]-=div*hold.a[0][i];
      imat.a[2][i]-=div*imat.a[0][i];
   }
   div=hold.a[3][0];
   for(i=0;i<4;i++){
      hold.a[3][i]-=div*hold.a[0][i];
      imat.a[3][i]-=div*imat.a[0][i];
   }

   if(fabs(hold.a[1][1])<DINKY){
      if(fabs(hold.a[2][1])>DINKY) rswap=2;
      else if(fabs(hold.a[3][1])>DINKY) rswap=3;
      for(i=0;i<4;i++){
         swap=hold.a[1][i];
         hold.a[1][i]=hold.a[rswap][i];
         hold.a[rswap][i]=swap;
         swap=imat.a[1][i];
         imat.a[1][i]=imat.a[rswap][i];
         imat.a[rswap][i]=swap;
      }
   }

   div=hold.a[1][1];
   for(i=0;i<4;i++){
      hold.a[1][i]/=div;
      imat.a[1][i]/=div;
   }

   div=hold.a[0][1];
   for(i=0;i<4;i++){
      hold.a[0][i]-=div*hold.a[1][i];
      imat.a[0][i]-=div*imat.a[1][i];
   }
   div=hold.a[2][1];
   for(i=0;i<4;i++){
      hold.a[2][i]-=div*hold.a[1][i];
      imat.a[2][i]-=div*imat.a[1][i];
   }
   div=hold.a[3][1];
   for(i=0;i<4;i++){
      hold.a[3][i]-=div*hold.a[1][i];
      imat.a[3][i]-=div*imat.a[1][i];
   }

   if(fabs(hold.a[2][2])<DINKY){
      for(i=0;i<4;i++){
         swap=hold.a[2][i];
         hold.a[2][i]=hold.a[3][i];
         hold.a[3][i]=swap;
         swap=imat.a[2][i];
         imat.a[2][i]=imat.a[3][i];
         imat.a[3][i]=swap;
      }
   }

   div=hold.a[2][2];
   for(i=0;i<4;i++){
      hold.a[2][i]/=div;
      imat.a[2][i]/=div;
   }

   div=hold.a[0][2];
   for(i=0;i<4;i++){
      hold.a[0][i]-=div*hold.a[2][i];
      imat.a[0][i]-=div*imat.a[2][i];
   }
   div=hold.a[1][2];
   for(i=0;i<4;i++){
      hold.a[1][i]-=div*hold.a[2][i];
      imat.a[1][i]-=div*imat.a[2][i];
   }
   div=hold.a[3][2];
   for(i=0;i<4;i++){
      hold.a[3][i]-=div*hold.a[2][i];
      imat.a[3][i]-=div*imat.a[2][i];
   }

   div=hold.a[3][3];
   for(i=0;i<4;i++){
      hold.a[3][i]/=div;
      imat.a[3][i]/=div;
   }

   div=hold.a[0][3];
   for(i=0;i<4;i++){
      hold.a[0][i]-=div*hold.a[3][i];
      imat.a[0][i]-=div*imat.a[3][i];
   }
   div=hold.a[1][3];
   for(i=0;i<4;i++){
      hold.a[1][i]-=div*hold.a[3][i];
      imat.a[1][i]-=div*imat.a[3][i];
   }
   div=hold.a[2][3];
   for(i=0;i<4;i++){
      hold.a[2][i]-=div*hold.a[3][i];
      imat.a[2][i]-=div*imat.a[3][i];
   }

   *inv=imat;
   return(inv);
}


/*------------------*/
/*  transposemat()  */
/*----------------------------------------------------------------------------
   PURPOSE:   takes a matrix and calculates the transpose of the matrix
   PRECONDITIONS:   mat points to the matrix to transpose and trans points to
    storage for the result
   POSTCONDITIONS:   the transpose of mat is in trans, or NULL if the function
    failed
----------------------------------------------------------------------------*/
mat4type *transposemat(mat,trans)
mat4type *mat,                          /* matrix to transpose */
 *trans;                                /* storage for transposition */
{
   register double t;
   
   t=mat->a[0][1];
   trans->a[0][1]=mat->a[1][0];
   trans->a[1][0]=t;
   
   t=mat->a[0][2];
   trans->a[0][2]=mat->a[2][0];
   trans->a[2][0]=t;
   
   t=mat->a[0][3];
   trans->a[0][3]=mat->a[3][0];
   trans->a[3][0]=t;
   
   t=mat->a[1][2];
   trans->a[1][2]=mat->a[2][1];
   trans->a[2][1]=t;
   
   t=mat->a[1][3];
   trans->a[1][3]=mat->a[3][1];
   trans->a[3][0]=t;
   
   t=mat->a[2][3];
   trans->a[2][3]=mat->a[3][2];
   trans->a[3][2]=t;

   return(trans);
}


/*---------------*/
/*  normalize()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function normalizes a vector
   PRECONDITIONS:   vec1 points to the vector, vec2 points to storage for the
    result
   POSTCONDITIONS:   vec2 holds the results of vec1 now with a length of 1,
    or NULL is returned if the vector originally had no length
----------------------------------------------------------------------------*/
vec4type *normalize(vec1,vec2)
vec4type *vec1,                         /* vec to normalize */
 *vec2;                                 /* normalized */
{
   register double len;                 /* the length of the vector */

   len=sqrt(vec1->a[0]*vec1->a[0]+vec1->a[1]*vec1->a[1]+vec1->a[2]*vec1->a[2]);
   if(len!=0.0){
      vec2->a[0]=vec1->a[0]/len;
      vec2->a[1]=vec1->a[1]/len;
      vec2->a[2]=vec1->a[2]/len;
      vec2->a[3]=1.0;
      return(vec2);
   }else
      return(NULL);
}


/*----------------------------------------------------------------------------

note that the two-vector representation of a plane can be created from the
 three-point representation by these operations:

	origin=point1;
        normal=normalize(cross(subvec(point1,point2),subvec(point1,point3)))

also note that this is pseudo-code, not actual functions from this library.

----------------------------------------------------------------------------*/


/*-------------*/
/*  reflect()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function calculates the reflection of a vector off of a
    plane with a given normal
   PRECONDITIONS:   inc points to the incident vector, norm points to the
    plane normal, and ref points to storage for the reflection.  both inc and
    norm must be normalized
   POSTCONDITIONS:   ref holds the normalized reflection of inc off of norm
----------------------------------------------------------------------------*/
vec4type *reflect(inc,norm,ref)
vec4type *inc,                          /* incident vector */
 *norm,                                 /* normal vector */
 *ref;                                  /* reflected vector */
{
   register double dxsq,dysq,dzsq,      /* these variables hold numbers to */
    dxrn,dyrn,dzrn;                     /*  keep things from being */
                                        /*  calcualted twice */

   dxsq=2*sq(norm->a[0]);
   dysq=2*sq(norm->a[1]);
   dzsq=2*sq(norm->a[2]);
   dxrn=norm->a[0]*inc->a[0];
   dyrn=norm->a[1]*inc->a[1];
   dzrn=norm->a[2]*inc->a[2];
   tempvec.a[0]=inc->a[0]*(dysq+dzsq-1)-2*norm->a[0]*(dyrn+dzrn);
   tempvec.a[1]=inc->a[1]*(dxsq+dzsq-1)-2*norm->a[1]*(dxrn+dzrn);
   tempvec.a[2]=inc->a[2]*(dysq+dxsq-1)-2*norm->a[2]*(dyrn+dxrn);
   tempvec.a[3]=1.0;
   return(normalize(&tempvec,ref));
}


/*---------*/
/*  dot()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function returns the dot product of two vectors 
   PRECONDITIONS:   vec1 and vec2 point to the vectors to dot
   POSTCONDITIONS:   the dot product is returned
----------------------------------------------------------------------------*/
double dot(vec1,vec2)
vec4type *vec1,                         /* vectors to dot */
 *vec2;
{
    return(vec1->a[0]*vec2->a[0]+vec1->a[1]*vec2->a[1]+vec1->a[2]*vec2->a[2]);
}


/*-----------*/
/*  cross()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function calculates the cross product of two vectors
   PRECONDITIONS:   vec1 and vec2 point to the vectors to cross, and ortho
    points to storage for the cross product
   POSTCONDITIONS:   ortho contains the cross product, or NULL is returned if
    the vectors are parallel
----------------------------------------------------------------------------*/
vec4type *cross(vec1,vec2,ortho)
vec4type *vec1,                         /* 1st vector to cross */
 *vec2,                                 /* 2nd vector to cross */
 *ortho;                                /* cross product */
{
   tempvec.a[0]=vec1->a[1]*vec2->a[2]-vec1->a[2]*vec2->a[1];
   tempvec.a[1]=vec1->a[2]*vec2->a[0]-vec1->a[0]*vec2->a[2];
   tempvec.a[2]=vec1->a[0]*vec2->a[1]-vec1->a[1]*vec2->a[0];
   tempvec.a[3]=1.0;
   return(normalize(&tempvec,ortho));
}


/*---------------*/
/*  magnitude()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function returns the length (magnitude) of the vector passed.
   PRECONDITIONS:   vector points to the vector to measure
   POSTCONDITIONS:   returns the length of the passed vector
----------------------------------------------------------------------------*/
double magnitude(vector)
vec4type *vector;                       /* vector to measure */
{
   return(sqrt(sq(vector->a[0])+sq(vector->a[1])+sq(vector->a[2])));
}


/*-------------*/
/*  disttoray  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function returns the shortest distance from a point to a 
    ray (represented by an origin and direction)
   PRECONDITIONS:   point points to the point, ray points to the origin and
    then direction vector of the ray
   POSTCONDITIONS:   returns the distance fro the point to the ray
----------------------------------------------------------------------------*/
double disttoray(point,ray)
vec4type *point,                        /* point to measure from */
 ray[2];                                /* ray to measure to */
{
   register double dx,dy,dz,            /* direction from point to ray */
                                        /*  origin */
    t;                                  /* direction down ray to closest */
                                        /*  approach to point */

   dx=point->a[0]-ray[0].a[0];
   dy=point->a[1]-ray[0].a[1];
   dz=point->a[2]-ray[0].a[2];
   t=ray[1].a[0]*dx+ray[1].a[1]*dy+ray[1].a[2]*dz;
   return(sqrt(dx*dx+dy*dy+dz*dz-t*t));
}


/*-----------------*/
/*  disttoplane()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function finds the distance from a point to a plane
    (represented by a ray as the normal and point on the plane)
   PRECONDITIONS:   point points to the point, and ray points to the normal,
    which must be normalized.
   POSTCONDITIONS:   returns the distance from the point to the plane.
    a negative result implies that the point lies on the opposite side of the
    plane from the plane normal.
----------------------------------------------------------------------------*/
double disttoplane(point,ray)
vec4type *point,                        /* point */
 ray[2];                                /* normal of and point on the plane */
{
   vec4type toplane,                    /* a vector from the point to the */
                                        /*  plane */
    tonorm;                             /* vector from the point to the */
                                        /*  origin of the plane normal */

   negatevec(ray+1,&toplane);
   subvec(ray,point,&tonorm);
   return(dot(&toplane,&tonorm));
}


/*--------------*/
/*  angtovec()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function produces a vector equal to the unit z vector 
    rotated alpha radians around the x-axis, and then rotated beta radians 
    around the y-axis
   PRECONDITIONS:   alpha and beta are in radians, vector points to storage
    for the new vector 
   POSTCONDITIONS:   vector contains a vector through alpha and beta angles
----------------------------------------------------------------------------*/
vec4type *angtovec(alpha,beta,vector)
double alpha,                           /* angle to rotate around x */
 beta;                                  /* angle to rotate around beta */
vec4type *vector;                       /* vector created */
{
   vector->a[0]=cos(beta)*cos(alpha);
   vector->a[1]=sin(alpha);
   vector->a[2]=sin(beta)*cos(alpha);
   return(vector);
}


/*--------------*/
/*  pointonray  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function projects a point onto a ray at its closest
   PRECONDITIONS:   point points to the point, ray points to the ray, and proj
    points to storage for the projection
   POSTCONDITIONS:   proj holds the projection of the point onto the ray
----------------------------------------------------------------------------*/
vec4type *pointonray(point,ray,proj)
vec4type *point,                        /* point to project */
 ray[2],                                /* ray to project onto */
 *proj;                                 /* projection */
{
   vec4type normal,                     /* normalized direction of ray */
    toorigin;                           /* vector from point to ray origin */

   normalize(&ray[1],&normal);
   subvec(point,&ray[0],&toorigin);
   sizevec(&normal,dot(&normal,&toorigin),&normal);
   return(addvec(&ray[0],&normal,proj));
}


/*------------------*/
/*  pointonplane()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function projects a point onto a plane at its closest
   PRECONDITIONS:   point points to the point, ray points to a two-vector
    representation of the plane (origin, normal), and proj points to storage
    for the projection
   POSTCONDITIONS:   proj contains the projected point 
----------------------------------------------------------------------------*/
vec4type *pointonplane(point,ray,proj)
vec4type *point,                        /* point to project */
 ray[2],                                /* plane to project int onto */
 *proj;                                 /* projection */
{
   register double dist;                /* distance from point to plane */
   vec4type toplane;                    /* vector from the point to the */
                                        /*  plane */

   /*writevec(stdout,point);*/
   /*writevec(stdout,ray);*/
   /*writevec(stdout,ray+1);*/
   dist=disttoplane(point,ray);
   printf("dist is %lf\n",dist);
   sizevec(ray+1,(-dist),&toplane);
   addvec(point,&toplane,proj);
   return(proj);
}


/*----------------------------------------------------------------------------
   PURPOSE:   this function calculates the intersection of a ray with a plane
    and returns a point.
   PRECONDITIONS:   ray is a 2-element vector consisting of the origin of the
    ray and direction numbers, plane is a 2-el vector consisting of a point
    on the plane and the normal vector, and proj points to storage for the
    intersection.
   POSTCONDITIONS:   on success, the intersection point is stored in proj,
    and the address of proj is returned, or NULL returned on failure (the ray
    does not intersect the plane.)
----------------------------------------------------------------------------*/
vec4type *raytoplane(ray,plane,proj)
vec4type ray[2],                      /* ray intersecting... */
 plane[2],                            /* ... plane */
 *proj;                               /* intersection point */
{
   vec4type t1,t2,t3,                 /* temporary calculation vectors */
    i;                                /* projection of ray origin onto plane */
   double prop;                       /* dot product of ray direction and */
                                      /*  directions from ray origin to */
                                      /*  plane */

   pointonplane(ray,plane,&i);
   printf("ray origin projected onto plane: ");
   writevec(stdout,&i);
   normalize(ray+1,&t1);
   subvec(&i,ray,&t2);
   normalize(&t2,&t3);
   prop=dot(&t1,&t3);
   if(prop==0.0)
      return(NULL);
   sizevec(&t1,magnitude(&t2)/prop,&t3);
   addvec(ray,&t3,proj);
   return(proj);
}


/*----------------------------------------------------------------------------
   PURPOSE:   this function clips a line segment to a plane
   PRECONDITIONS:   plane is a two-element array consisting of a point on the
    plane and the normal vector, pt1 and pt2 are the endpoints of the segment,
    and new points to storage for the new endpoint (if there is one.)
   POSTCONDITIONS:   the value returned can be one of the following; note that
    for the last two, new contains the new endpoint for the segment.
              0 - no segment remains after clip
              1 - whole segment after clip
              2 - segment from pt1 to new remains
              3 - segment from pt2 to new remains
   NOTE:   this would be a hell of a lot faster with some forethought and no
    function calls.
----------------------------------------------------------------------------*/
int clipsegment(plane,pt1,pt2,new)
vec4type plane[2],                    /* plane array */
 *pt1,*pt2,                           /* endpoints of the segment to clip */ 
 *new;                                /* possible new endpoint of segment */
{
   vec4type r[2],                     /* ray representing the line segment */
    t1,t2;                            /* temporary calculation vectors */
   double a,b;                        /* location indicators for endpoints */

   a=dot(plane+1,subvec(pt1,plane,&t1));
   b=dot(plane+1,subvec(pt2,plane,&t2));
   if((a<0.0)&&(b<0.0))
      return(CLIP_NONE);
   if((a>0.0)&&(b>0.0))
      return(CLIP_WHOLE);
   r[0]=*pt1;
   subvec(pt2,pt1,r+1);
   normalize(r+1);
   raytoplane(r,plane,new);
   if(a>0.0)
      return(CLIP_LOSTP2);
   return(CLIP_LOSTP1);
}


/*----------------------------------------------------------------------------
   PURPOSE:  this function transforms a ray through a matrix (this operation
    is really the multiplication of two vectors by the matrix; one being the
    ray origin, and the other being the origin plus the direction vector.)
   PRECONDITIONS:   ray points to the ray to transform, mat points to the ,
    matrix, and prod points to storage for the transformed ray 
   POSTCONDITIONS:   prod contains the transformation of the passed ray 
    though the passed matrix
----------------------------------------------------------------------------*/
vec4type *multraymat(ray,mat,prod)
vec4type ray[2];                        /* the ray to transform */
mat4type *mat;                          /* matrix to transform by */
vec4type prod[2];                       /* product of transformation */
{
   vec4type rayvec;                     /* vector from origin + direction */

   addvec(&ray[0],&ray[1],&rayvec);
   /*translvec(ray[0].a[0],ray[0].a[1],ray[0].a[2],ray+1,&rayvec);*/
   multvecmat(&ray[0],mat,&prod[0]);
   multvecmat(&rayvec,mat,&rayvec);
   /*translvec(-ray[0].a[0],-ray[0].a[1],-ray[0].a[2],&rayvec,prod+1);*/
   subvec(&rayvec,&prod[0],&prod[1]);
   return(prod);
}


/*-----------------*/
/*  multrmblock()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function performs exactly like multraymat() (above), except
    that it operates on an array of rays
   PRECONDITIONS:   rays points to the ray array, mat points to the matrix,
    prods points to storage for the results, and total is the number of rays
    to transform
   POSTCONDITIONS:   prods holds the transformed array
----------------------------------------------------------------------------*/
vec4type *multrmblock(rays,mat,prods,total)
vec4type rays[][2];                     /* rays array */
mat4type *mat;                          /* matrix to transform by */
vec4type prods[][2];                    /* product of transformation */
double total;                           /* total rays to transform */
{
   register int rayloop;                /* ray counter */
   vec4type rayvec;                     /* vector from origin + direction */

   for(rayloop=0;rayloop<total;rayloop++){
      addvec(&rays[rayloop][0],&rays[rayloop][1],&rayvec);
      multvecmat(&rays[rayloop][0],mat,&prods[rayloop][0]);
      multvecmat(&rayvec,mat,&rayvec);
      subvec(&rayvec,&rays[rayloop][0],&prods[rayloop][1]);
   }
   return(&prods[0][0]);
}


/*------------*/
/*  rotvec()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function rotates a vector around a given axis
   PRECONDITIONS:   angle is the rotation angle in radians, axis is a
    character representing the axis ('X','x','Y','y','Z','z'), vec points to
    the vector to rotate, and result points to storage for the rotated vector
   POSTCONDITIONS:   result contains the result of the rotation, or NULL if
    the character is invalid
----------------------------------------------------------------------------*/
vec4type *rotvec(angle,axis,vec,result)
double angle;                           /* angle to rotate */
char axis;                              /* axis of rotation */
vec4type *vec,                          /* vector to rotate */
 *result;                               /* result of rotation */
{
   register double newx,newy,newz,      /* temp storage for new values */
    angcos,angsin;                      /* angle sine and cosine */

   angcos=cos(angle);
   angsin=sin(angle);
   switch(axis){
     case 'x':  case 'X':
      newx=vec->a[0];
      newy=vec->a[1]*angcos+vec->a[2]*angsin;
      newz=vec->a[2]*angcos-vec->a[1]*angsin;
      break;
     case 'y':  case 'Y':
      newx=vec->a[0]*angcos-vec->a[2]*angsin;
      newy=vec->a[1];
      newz=vec->a[2]*angcos+vec->a[0]*angsin;
      break;
     case 'z':  case 'Z':
      newx=vec->a[0]*angcos+vec->a[1]*angsin;
      newy=vec->a[1]*angcos-vec->a[0]*angsin;
      newz=vec->a[2];
      break;
     default:
      return(NULL);
   }
   result->a[0]=newx;
   result->a[1]=newy;
   result->a[2]=newz;
   result->a[3]=vec->a[3];
   return(result);
}


/*------------*/
/*  rotmat()  */
/*----------------------------------------------------------------------------
   PURPOSE:   multiply a matrix M by a rotation matrix R so that a vector
    transformed by the result is equivalent to v*M*R
   PRECONDITIONS:   angle is the angle to rotate in radians, axis is a
    character represting the axis (see rotvec, above), mat points to the
    matrix to multiply, and result points to storage for the product
   POSTCONDITIONS:   result holds the matrix rotated, or NULL is returned if
    the character passed is invalid 
----------------------------------------------------------------------------*/
mat4type *rotmat(angle,axis,mat,result)
double angle;                           /* angle to rotate by */
char axis;                              /* axis of rotation */
mat4type *mat,                          /* matrix to rotate */
 *result;                               /* rotation */
{
   register double angcos,angsin,       /* angle sine and cosine */
    val1,val2;                          /* temp storage for numbers */
   register int rowcnt;                 /* row counter */

   angcos=cos(angle);
   angsin=sin(angle);
   *result=(*mat);
   switch(axis){
     case 'x':  case 'X':
      for(rowcnt=0;rowcnt<4;rowcnt++){
         val1=mat->a[rowcnt][1]*angcos+mat->a[rowcnt][2]*angsin;
         val2=mat->a[rowcnt][2]*angcos-mat->a[rowcnt][1]*angsin;
         result->a[rowcnt][1]=val1;
         result->a[rowcnt][2]=val2;
      }
      break;
     case 'y':  case 'Y':
      for(rowcnt=0;rowcnt<4;rowcnt++){
         val1=mat->a[rowcnt][0]*angcos-mat->a[rowcnt][2]*angsin;
         val2=mat->a[rowcnt][2]*angcos+mat->a[rowcnt][0]*angsin;
         result->a[rowcnt][0]=val1;
         result->a[rowcnt][2]=val2;
      }
      break;
     case 'z':  case 'Z':
      for(rowcnt=0;rowcnt<4;rowcnt++){
         val1=mat->a[rowcnt][0]*angcos+mat->a[rowcnt][1]*angsin;
         val2=mat->a[rowcnt][1]*angcos-mat->a[rowcnt][0]*angsin;
         result->a[rowcnt][0]=val1;
         result->a[rowcnt][1]=val2;
      }
      break;
     default:
      return(NULL);
   }
   return(result);
}
         


/*---------------*/
/*  rotvecvec()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function rotates a vector around a ray
   PRECONDITIONS:   angle is the rotation in radians, vec is the axis ray
    of rotation (direction must be normalized), src is the vector to rotate, 
    and res is the result of the rotation 
   POSTCONDITIONS:   res holds the result of the rotation
----------------------------------------------------------------------------*/
vec4type rotvecvec(angle,vec,src,res)
double angle;                           /* angle of rotation */
vec4type vec[2],                        /* axis of rotation vector */
 *src,                                  /* vector to rotate */
 *res;                                  /* result of rotation */
{
   register double angcos,angsin,foocos,foosin,barcos,barsin;
}


/*---------------*/
/*  rotmatvec()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function is the same as rotmat() (above), except that the
    matrix is rotated around an arbitrary ray
   PRECONDITIONS:   angle is the rotation in radians, vec points to the axis
    of rotation ray, mat points to the matrix to multiply, and res points to
    storage for the results
   POSTCONDITIONS:   res holds the result of the rotation
----------------------------------------------------------------------------*/
mat4type rotmatvec(angle,vec,src,res)
double angle;                           /* angle of rotation */
vec4type vec[2];                        /* axis of rotation vector */
mat4type *src,                          /* matrix to rotate */
 *res;                                  /* result */
{
}


/*---------------*/
/*  translvec()  */
/*----------------------------------------------------------------------------
   PURPOSE:   translate a vector by arbitrary displacements
   PRECONDITIONS:   x,y,z are doubles representing the translation, vec points
    to the vector to translate, res points to storage for the result
   POSTCONDITIONS:   res holds the result of the translation
----------------------------------------------------------------------------*/
vec4type *translvec(x,y,z,vec,res)
double x,y,z;                           /* translation values */
vec4type *vec,                          /* vector to translate */
 *res;                                  /* result of translation */
{
   res->a[0]=vec->a[0]+x;
   res->a[1]=vec->a[1]+y;
   res->a[2]=vec->a[2]+z; 
   res->a[3]=vec->a[3];
   return(res);
}


/*---------------*/
/*  translmat()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function multiplies a matrix such that a vector multiplied
    by the product is equivalent to the vector transformed by the original
    matrix and then translated by the given values.
   PRECONDITIONS:   x,y,z are doubles representing the displacement, mat
    points to the matrix to translate, and res points to storage for the result
   POSTCONDITIONS:   res holds the result of the multiplication/translation
----------------------------------------------------------------------------*/
mat4type *translmat(x,y,z,mat,res)
double x,y,z;                           /* translation values */
mat4type *mat,                          /* matrix to translate */
 *res;                                  /* result matrix */
{
   register int rowcnt;                 /* row counter */

   *res=(*mat);
   for(rowcnt=0;rowcnt<4;rowcnt++){
      res->a[rowcnt][0]+=res->a[rowcnt][3]*x;
      res->a[rowcnt][1]+=res->a[rowcnt][3]*y;
      res->a[rowcnt][2]+=res->a[rowcnt][3]*z;
   }
   return(res);
}


/*--------------*/
/*  scalevec()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function scales a vector by given x, y, and z components
   PRECONDITIONS:   x,y,z are doubles representing the scaling, vec points to
    the vector to scale, and res points to storage for the result
   POSTCONDITIONS:   res holds the result of the scaling
----------------------------------------------------------------------------*/
vec4type *scalevec(x,y,z,vec,res)
double x,y,z;                           /* scale values */
vec4type *vec,                          /* vector to scale */
 *res;                                  /* result of scaling */
{
   res->a[0]=vec->a[0]*x;
   res->a[1]=vec->a[1]*y;
   res->a[2]=vec->a[2]*z;
   res->a[3]=vec->a[3];
   return(res);
}


/*--------------*/
/*  scalemat()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function multiplies a matrix so that a vector transformed
    by the product is equal to the vector transformed by the original matrix
    and then scaled by the given values
   PRECONDITIONS:   x,y,z are doubles representing the scale values, mat
    points to the matrix to multiply, and res points to storage for the result
   POSTCONDITIONS:   res holds the result of the transformation
----------------------------------------------------------------------------*/
mat4type *scalemat(x,y,z,mat,res)
double x,y,z;                           /* scale values */
mat4type *mat,                          /* matrix to multiply */
 *res;                                  /* result of multiplication */
{
   register int rowcnt;                 /* row counter */

   *res=(*mat);
   for(rowcnt=0;rowcnt<4;rowcnt++){
      res->a[rowcnt][0]*=x;
      res->a[rowcnt][1]*=y;
      res->a[rowcnt][2]*=z;
   }
   return(res);
}


/*-----------------*/
/*  skewtoortho()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function creates a matrix that will transform a pair of
    vectors or points or both onto the y and z axes, through a series of
    translations and rotations
   PRECONDITIONS:   mat points to storage for the result, src points to the
    point that will become the origin, ypt points to either a vector
    indicating the direction of the y axis or a point that represents the axis
    as a line going through src and ypt, and yvec is true if it's a vector,
    false if it's a point, and finally, likewise for zpt.
   POSTCONDITIONS:   mat contains the resulting matrix
----------------------------------------------------------------------------*/
mat4type *skewtoortho(mat,src,ypt,yvec,zpt,zvec)
mat4type *mat;                          /* storage for the result */
vec4type *src,                          /* origin point */
 *ypt;                                  /* y axis representation */
int yvec;                               /* y a vector or point */
vec4type *zpt;                          /* z axis representation */
int zvec;                               /* z a vector or point */
{
   vec4type ztemp,ytemp;                /* temp vectors during rotations */
   mat4type tempmat;

   translmat(-src->a[0],-src->a[1],-src->a[2],&identmat4,mat);

   if(!yvec){
      translvec(-src->a[0],-src->a[1],-src->a[2],ypt,&ytemp);
      normalize(&ytemp,&ytemp);
   }else
      ytemp=(*ypt);

   if(!zvec){
      translvec(-src->a[0],-src->a[1],-src->a[2],zpt,&ztemp);
      normalize(&ztemp,&ztemp);
   }else
      ztemp=(*zpt);

   if((fabs(ztemp.a[0])>DINKY)||(fabs(ztemp.a[2])>DINKY)){
      rotmat(atan2(ztemp.a[0],ztemp.a[2]),'y',&identmat4,&tempmat);
      multmat4(mat,&tempmat,mat);
      multvecmat(&ztemp,&tempmat,&ztemp);
      multvecmat(&ytemp,&tempmat,&ytemp);
   }

   rotmat(-atan2(ztemp.a[1],ztemp.a[2]),'x',&identmat4,&tempmat);
   multmat4(mat,&tempmat,mat);
   multvecmat(&ytemp,&tempmat,&ytemp);

   rotmat(-atan2(ytemp.a[0],ytemp.a[1]),'z',&identmat4,&tempmat);
   multmat4(mat,&tempmat,mat);
   return(mat);
}


/*-----------------*/
/*  orthotoskew()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function performs the same function as skewtoortho()
    (above), except in reverse; the matrix produced will transform the origin
    coordinate system to the given representation
   PRECONDITIONS:   mat points to storage for the result, src points to the
    point that will become the origin, ypt points to either a vector
    indicating the direction of the y axis or a point that represents the axis
    as a line going through src and ypt, and yvec is true if it's a vector,
    false if it's a point, and finally, likewise for zpt.
   POSTCONDITIONS:   mat contains the resulting matrix
----------------------------------------------------------------------------*/
mat4type *orthotoskew(mat,src,ypt,yvec,zpt,zvec)
mat4type *mat;                          /* storage for the result */
vec4type *src,                          /* origin point */
 *ypt;                                  /* y axis representation */
int yvec;                               /* y a vector or point */
vec4type *zpt;                          /* z axis representation */
int zvec;                               /* z a vector or point */
{
   vec4type ztemp,ytemp;                /* temp vectors during rotations */
   mat4type tempmat,tmat2;
   double yang,xang,zang;

   skewtoortho(&tempmat,src,ypt,yvec,zpt,zvec);
   return(invertmat(&tempmat,mat));
   if(!yvec)
      translvec(-src->a[0],-src->a[1],-src->a[2],ypt,&ytemp);
   else
      ytemp=(*ypt);

   if(!zvec)
      translvec(-src->a[0],-src->a[1],-src->a[2],zpt,&ztemp);
   else
      ztemp=(*zpt);

   if((fabs(ztemp.a[0])>DINKY)||(fabs(ztemp.a[2])>DINKY)){
      rotmat((yang=atan2(ztemp.a[0],ztemp.a[2])),'y',&identmat4,&tempmat);
      multvecmat(&ztemp,&tempmat,&ztemp);
      multvecmat(&ytemp,&tempmat,&ytemp);
   }

   rotmat((xang=-atan2(ztemp.a[1],ztemp.a[2])),'x',&identmat4,&tempmat);
   multvecmat(&ytemp,&tempmat,&ytemp);

   zang=-atan2(ytemp.a[0],ytemp.a[1]);
   rotmat(-zang,'z',&identmat4,mat);
   rotmat(-xang,'x',mat,mat);
   rotmat(-yang,'y',mat,mat);
   translmat(src->a[0],src->a[1],src->a[2],mat,mat);

   return(mat);
}


/*--------------*/
/*  writemat()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function writes a 4x4 matrix to a given file
   PRECONDITIONS:   outf points to the output file, and matrix points to the
    matrix to write
   POSTCONDITIONS:   returns 1 on success, 0 on write failure
----------------------------------------------------------------------------*/
int writemat(outf,matrix)
FILE *outf;                             /* output file */
mat4type *matrix;                       /* matrix to write */
{
   if(fprintf(outf,"%lf %lf %lf %lf\n",matrix->a[0][0],matrix->a[0][1],matrix->a[0][2],
    matrix->a[0][3])<0)
      return(0);
   if(fprintf(outf,"%lf %lf %lf %lf\n",matrix->a[1][0],matrix->a[1][1],matrix->a[1][2],
    matrix->a[1][3])<0)
      return(0);
   if(fprintf(outf,"%lf %lf %lf %lf\n",matrix->a[2][0],matrix->a[2][1],matrix->a[2][2],
    matrix->a[2][3])<0)
      return(0);
   if(fprintf(outf,"%lf %lf %lf %lf\n",matrix->a[3][0],matrix->a[3][1],matrix->a[3][2],
    matrix->a[3][3])<0)
      return(0);
   return(1);
}


/*-------------*/
/*  readmat()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function reads a 4x4 matrix from a given file
   PRECONDITIONS:   outf points to the input file and matrix points to the
    matrix to read into
   POSTCONDITIONS:   matrix contains the input and function returns 1, or
    returns 0 on read failure
----------------------------------------------------------------------------*/
int readmat(outf,matrix)
FILE *outf;                             /* input file */
mat4type *matrix;                       /* storage for input */
{
   if(fscanf(outf,"%lf %lf %lf %lf\n",&matrix->a[0][0],&matrix->a[0][1],&matrix->a[0][2],
    &matrix->a[0][3])!=4)
      return(0);
   if(fscanf(outf,"%lf %lf %lf %lf\n",&matrix->a[1][0],&matrix->a[1][1],&matrix->a[1][2],
    &matrix->a[1][3])!=4)
      return(0);
   if(fscanf(outf,"%lf %lf %lf %lf\n",&matrix->a[2][0],&matrix->a[2][1],&matrix->a[2][2],
    &matrix->a[2][3])!=4)
      return(0);
   if(fscanf(outf,"%lf %lf %lf %lf\n",&matrix->a[3][0],&matrix->a[3][1],&matrix->a[3][2],
    &matrix->a[3][3])!=4)
      return(0);
   return(1);
}


/*--------------*/
/*  writevec()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function writes a vec4type to a given output file
   PRECONDITIONS:   outf points to the output file, and vec points to the
    vector to write
   POSTCONDITIONS:   returns 1 on success, 0 on write failure
----------------------------------------------------------------------------*/
int writevec(outf,vec)
FILE *outf;                             /* output file */
vec4type *vec;                          /* vector to output */
{
   if(fprintf(outf,"%lf %lf %lf %lf\n",vec->a[0],vec->a[1],vec->a[2],
    vec->a[3])<0)
      return(0);
   return(1);
}


/*-------------*/
/*  readvec()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function reads a vec4type from a file
   PRECONDITIONS:   outf points to the input file, vec points to the vector
    to be read in
   POSTCONDITIONS:   vector contains the input and function returns 1, or
    returns 0 on read failure
----------------------------------------------------------------------------*/
int readvec(outf,vec)
FILE *outf;                             /* input file */
vec4type *vec;                          /* storage for input */
{
   if(fscanf(outf,"%lf %lf %lf %lf\n",&vec->a[0],&vec->a[1],&vec->a[2],
    &vec->a[3])!=4)
      return(0);
   return(1);
}


/*---------------*/
/*  writevec3()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function writes a three element vector from a vec4type
    (by ignoring element 4) to a given output file
   PRECONDITIONS:   outf points to the output file, and vec points to the
    vector to write
   POSTCONDITIONS:   returns 1 on success, 0 on write failure
----------------------------------------------------------------------------*/
int writevec3(outf,vec)
FILE *outf;                             /* output file */
vec4type *vec;                          /* vector to output */
{
   if(fprintf(outf,"%lf %lf %lf\n",vec->a[0],vec->a[1],vec->a[2])<0)
      return(0);
   return(1);
}


/*--------------*/
/*  readvec3()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function reads a 3-element vector from a file and converts
    it into a 4-element vector with 1.0 as the fourth element
   PRECONDITIONS:   outf points to the input file, vec points to the vector
    to be read in
   POSTCONDITIONS:   vector contains the input and function returns 1, or
    returns 0 on read failure
----------------------------------------------------------------------------*/
int readvec3(outf,vec)
FILE *outf;                             /* input file */
vec4type *vec;                          /* storage for input */
{
   if(fscanf(outf,"%lf %lf %lf\n",&vec->a[0],&vec->a[1],&vec->a[2])!=3)
      return(0);
   vec->a[3]=1.0;
   return(1);
}


/*----------------*
/*  angletween()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this routine calculates the angle between two vectors in three
    dimensions
   PRECONDITIONS:   vec1 and vec2 point to the vectors to calculate the vector
    between
   POSTCONDITIONS:   returns the angle between vec1 and vec2
----------------------------------------------------------------------------*/
double angletween(vec4type *vec1, vec4type *vec2)
{
   return(acos(dot(vec1,vec2)/magnitude(vec1)/magnitude(vec2)));
}


vec4type *xformnorm(norm,mat,res)
vec4type *norm;
mat4type *mat;
vec4type *res;
{
   
}
