#include <sys/time.h>
#include <sys/resource.h>

void getrusage_sec(double& utime, double& stime)
{
   struct rusage t;
   struct timeval tvu, tvs;
   getrusage(RUSAGE_SELF, &t);
   tvu = t.ru_utime;
   tvs = t.ru_stime;

   utime = (double)(tvu.tv_sec) + (double)(tvu.tv_usec)*1e-6;
   stime = (double)(tvs.tv_sec) + (double)(tvs.tv_usec)*1e-6;
   return;

}

double getrusage_usec()
{
   struct rusage t;
   struct timeval tvu;
   getrusage(RUSAGE_SELF, &t);
   tvu = t.ru_utime;

   return (double)(tvu.tv_sec) + (double)(tvu.tv_usec)*1e-6;

}

double getrusage_ssec()
{
   struct rusage t;
   struct timeval tvs;
   getrusage(RUSAGE_SELF, &t);
   tvs = t.ru_stime;

   return (double)(tvs.tv_sec) + (double)(tvs.tv_usec)*1e-6;

}

