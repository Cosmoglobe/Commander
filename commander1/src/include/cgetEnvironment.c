#include <stdlib.h>
#include <string.h>

#if defined(RS6000) || defined(__xlC__)
void cgetenvironment
#else
void cgetenvironment_
#endif
  (char *e, char *v, int ne, int nv) {
  int i,ns;
  char *s,*p,*q;
  s=getenv(e);
  if (s==NULL) s="";
  ns=strlen(s);
  if (ns>nv) ns=nv;
  p=&s[0];
  q=&v[0];
  for (i=0; i<ns; i++) {*q++=*p++;}
  *q='\0';
}
