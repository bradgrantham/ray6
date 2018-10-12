/*----------------------------------------------------------------------------
   Written by Brad Grantham
   FILE:   grf_color.c
   ENVIRONMENT:   MAC II under A/UX using gcc
   PURPOSE:   this file has the includes for the graphics colors unit.  these
    definitions and functions allow for color modularity and ease of 
    programming.
----------------------------------------------------------------------------*/


#include "grf_color.h"                  /* color includes */


colortype black_color={0.0,0.0,0.0},    /* no color */
 white_color={1.0,1.0,1.0};             /* all color */


#define TINY_COLOR 0.0000000596         /* really dinky */


#define HUGE_COLOR 0.9999999404         /* really big */


/*---------------*/
/*  isblack24()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function returns true if the color passed to it is
    satisfactorily dark;  i.e. less than 1/(2^24th) of white color
   PRECONDITIONS:   color points to the color to check
   POSTCONDITIONS:   returns true if color is pretty much black
----------------------------------------------------------------------------*/
int isblack24(color)
colortype *color;                       /* color to check */
{
   return(((color->r<TINY_COLOR)&&(color->g<TINY_COLOR)&&(color->b<
    TINY_COLOR))?1:0);
}


/*---------------*/
/*  iswhite24()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function returns true if the color passed to it is
    satisfactorily bright;  i.e. more than 1-1/(2^24th) of white color
   PRECONDITIONS:   color points to the color to check
   POSTCONDITIONS:   returns true if color is pretty much white
----------------------------------------------------------------------------*/
int iswhite24(color)
colortype *color;                       /* color to check */
{
   return(((color->r>HUGE_COLOR)&&(color->g>HUGE_COLOR)&&(color->b>
    HUGE_COLOR))?1:0);
}


/*----------------*/
/*  multcolors()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function basically multiplies two colors together, which I
    define as color1*color2/(maximum value of a color)
   PRECONDITIONS:   color1 and color2 point to the colors to multiply, product
    points to storage for the result
   POSTCONDITIONS:   result in product, and the function returns the pointer
    to the product
----------------------------------------------------------------------------*/
colortype *multcolors(color1,color2,product)
colortype *color1,*color2,              /* colors to multiply */
 *product;                              /* product of multiplication */
{
   product->r=color1->r*color2->r;
   product->g=color1->g*color2->g;
   product->b=color1->b*color2->b;
   return(product);
}


/*----------------*/
/*  scalecolor()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function basically scales a color by a certain factor
   PRECONDITIONS:   color1 points to the color to scale, result points to 
    storage for the result
   POSTCONDITIONS:   result in result, and the function returns the pointer
    to the result
----------------------------------------------------------------------------*/
colortype *scalecolor(color,factor,result)
colortype *color;                       /* color to scale */
double factor;                          /* factor value */
colortype *result;                      /* result of scaling */
{
   result->r=color->r*factor;
   result->g=color->g*factor;
   result->b=color->b*factor;
   return(result);
}


/*---------------*/
/*  addcolors()  */
/*----------------------------------------------------------------------------
   PURPOSE:   this function adds two colors together
   PRECONDITIONS:   color1 and color2 point to the colors to add, sum
    points to storage for the result
   POSTCONDITIONS:   result in sum, and the function returns the pointer
    to the sum
----------------------------------------------------------------------------*/
colortype *addcolors(color1,color2,sum)
colortype *color1,*color2,              /* colors to add */
 *sum;                                  /* sum of colors */
{
   sum->r=color1->r+color2->r;
   sum->g=color1->g+color2->g;
   sum->b=color1->b+color2->b;
   return(sum);
}
