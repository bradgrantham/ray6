/*----------------------------------------------------------------------------

   Written by Brad Grantham

   opened:   01/31/90
   last modification:   02/15/90

   FILE:       brads.h
   LANGUAGE:   C
   SYSTEM:     Mac II, under A/UX
   PURPOSE:    this file contains miscellaneous macros, functions, data
    structures that I use on a regular basis.

----------------------------------------------------------------------------*/

#ifndef BRADS_H                         /* this ifndef keeps the file from */
                                        /*  being compiled if it has already */
                                        /*  been read */
#define BRADS_H


/*--------------------*/
/*  #define getone()  */
/*--------------------*/
#define getone(s) ((s *)malloc(sizeof(s)))
                                        /* this macro allocates one thing of */
                                        /*  a given type */


/*---------------------*/
/*  #define getsome()  */
/*---------------------*/
#define getsome(s,n) ((s *)malloc(sizeof(s)*(n)))
                                        /* this macro allocates n of a given */
                                        /*  type */

/*---------------------*/
/*  #define dispose()  */
/*---------------------*/
#define dispose(s) free((char *)(s))    /* this macro takes care of the */
                                        /*  (char *) casting to use free(). */


/*--------------------*/
/*  #define highbyte  */
/*--------------------*/
#define highbyte(i) ((i)>>24)           /* macro for chopping the most  */
                                        /*  significant byte off of a 4-byte */
                                        /*  word. */


/*-------------------*/
/*  #define lowbyte  */
/*-------------------*/
#define lowbyte(i) ((i)&0xff)           /* macro for chopping the least  */
                                        /*  significant byte off of a 4-byte */
                                        /*  word. */


/*----------------------*/
/*  #define secondbyte  */
/*----------------------*/
#define secondbyte(i) (((i)>>16)&0xff)  /* removes the 2nd most significant */
                                        /*  byte from a 4-byte word */


/*---------------------*/
/*  #define thirdbyte  */
/*---------------------*/
#define thirdbyte(i) (((i)>>8)&0xff)    /* removes the 2nd least significant */
                                        /*  byte from a 4-byte word */


/*-------------------*/
/*  typedef funcptr  */
/*-------------------*/
typedef int (* funcptr)();              /* this is the type of a function */
                                        /*  that returns an integer.  it */
                                        /*  allows you to pass functions as */
                                        /*  parameters. */

typedef unsigned char byte;            /* why isn't this a standard C type? */


int strcasecmp(/*str1,str2*/);          /* case-insensitive string compare. */


#endif                                  /* end of BRADS_H stuff */
