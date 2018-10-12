/*----------------------------------------------------------------------------
   Written by Brad Grantham
   FILE:   grf_color.h
   ENVIRONMENT:   MAC II under A/UX using gcc
   PURPOSE:   this file has the includes for the graphics colors unit.  these
    definitions and functions allow for color modularity and ease of 
    programming.
----------------------------------------------------------------------------*/


/*---------------------*/
/*  typedef colortype  */
/*---------------------*/
typedef struct {
   double r,g,b;                        /* red, green, and blue */
} colortype;


extern colortype black_color,           /* no color */
 white_color;                           /* all color */


int isblack24(/*colortype *col*/);
int iswhite24(/*colortype *col*/);
colortype *multcolors(/*color1,color2,product*/);
colortype *scalecolor(/*color,factor,result*/);
colortype *addcolors(/*color1,color2,sum*/);
