#include "icoord.h"
#include "utils.h"
using namespace std;

#if 0
//not working
double randomf(double a, double b){

  timeval t1;
  gettimeofday(&t1, NULL);
//  srand(t1.tv_usec*t1.tv_sec);
  srand( time(NULL) );

  printf(" RAND: %i \n",rand());

  double range = b-a;
  double randn = a + double(range*rand()/(RAND_MAX+1));
  return randn;
}
#endif

int sign(double x){

  if (x>0) return 1;
  else if (x<=0) return -1;

}

void cross(double* m, double* r1, double* r2){

  m[0] = r1[1]*r2[2]-r2[1]*r1[2];
  m[1] = -r1[0]*r2[2]+r2[0]*r1[2];
  m[2] = r1[0]*r2[1]-r2[0]*r1[1];

  return;
}

int close_value(double a, double b, double c) {

  int close = false;

  if (b<=a && b+c>a)
    close = true;
  else if (a<=b && a+c>b)
    close = true;

  return close;
}
