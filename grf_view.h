/*----------------------------------------------------------------------------

   Written by Brad Grantham

   opened:   01/09/90
   last modification:   01/09/90

   FILE:       grf_view.h
   LANGUAGE:   C
   SYSTEM:     Mac II, under A/UX
   PURPOSE:    this file contains functions and declarations for use in
    graphics, specifically for viewplane definition, viewing transformations,
    and viewplane-ray computations

----------------------------------------------------------------------------*/


#ifndef GRF_VIEW_H
#define GRF_VIEW_H


#include "grf_math.h"


void polarview();
void planeview();
void directview();
void perspproj();
void orthoproj();


#endif
