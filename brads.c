#include "brads.h"


int strcasecmp(str1,str2)
char *str1,*str2;
{
   register int check,
    chcnt,
    minlen,
    len1,len2;

   len1=strlen(str1);
   len2=strlen(str2);
   minlen=(len1<len2)?len1:len2;
   for(chcnt=0;chcnt<minlen;chcnt++)
      if((check=toupper(str1[chcnt])-toupper(str2[chcnt]))!=0)
         return(check);
   if(len1!=len2)
      return((len1<len2)?-1:1);
   else
      return(0);
}
