//
// Tomas Oliveira e Silva
//
// rand_city_coords()
// * generates "random" city coodinates for a maximum of max_n_cities cities
// * the coordinates are stored in the arrays city_x[] and city_y[]
// * use your student number as its argument
//
// print_tour()
// * create a postscript file with a drawing of a given tour
// * its first argument is the number of cities
// * its second argument is the a[] array holding the order the cities are visited
// * its third argument is the file name (use the extension .ps) where is figure will be stored
//     on a linux terminal, that file can be converted to pdf using the command epspdf
//
// see also the commented out main() function at the end of this file for one possible
// way of using these functions; to compile this as a standalone program use the command
//   cc -Wall -O2 -DINCLUDE_MAIN tsp_aux.c -lm
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(_WIN32)
#include <Windows.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

#define max_n_cities 15

double city_x[max_n_cities];
double city_y[max_n_cities];

void rand_city_coords(int seed)
{
  int i,j;

  srand(seed);

  for(i = 0;i < max_n_cities;i++)
    for(;;)
    {
      city_x[i] = 0.01 * floor(1000.0 * (double)rand() / (double)RAND_MAX); // range 0.00:0.01:10.00
      city_y[i] = 0.01 * floor(1000.0 * (double)rand() / (double)RAND_MAX); // range 0.00:0.01:10.00
      for(j = 0;j < i;j++)
        if(fabs(city_x[i] - city_x[j]) + fabs(city_y[i] - city_y[j]) < 2.0)
          break; // cities too close together (Manhattan distance)
      if(j == i)
        break; // accept city coordinates
    }
}

void print_tour(int n_cities,int *a,char *file_name)
{
  FILE *fp;
  int i;

  if(n_cities < 2 || n_cities > max_n_cities)
  {
    fprintf(stderr,"bad number of cities\n");
    exit(1);
  }
  fp = fopen(file_name,"w");
  if(fp == NULL)
  {
    fprintf(stderr,"unable to create file %s\n",file_name);
    exit(1);
  }
  fprintf(fp,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(fp,"%%%%Creator: TOS\n");
  fprintf(fp,"%%%%BoundingBox: 0 0 312 312\n");
  fprintf(fp,"72 25.4 div dup scale 1 setlinejoin 1 setlinecap\n");
  fprintf(fp,"/Helvetica-Bold findfont 4 scalefont setfont\n");
  fprintf(fp,"0.5 setlinewidth 0.7 setgray\n");
  fprintf(fp,"5 10 105 { 5 moveto 0 100 rlineto } for stroke\n");
  fprintf(fp,"5 10 105 { 5 exch moveto 100 0 rlineto } for stroke\n");
  fprintf(fp,"1 setlinewidth 1 0 0 setrgbcolor\n");
  for(i = 0;i < n_cities;i++)
    fprintf(fp,"%.2f %.2f %s\n",5.0 + 10.0 * city_x[a[i]],5.0 + 10.0 * city_y[a[i]],(i == 0) ? "moveto" : "lineto");
  fprintf(fp,"closepath stroke 0.5 setlinewidth\n");
  for(i = 0;i < n_cities;i++)
  { // display city numbers
    fprintf(fp,"1 setgray %.2f %.2f 3 0 360 arc fill\n",5.0 + 10.0 * city_x[i],5.0 + 10.0 * city_y[i]);
    fprintf(fp,"0 setgray %.2f %.2f 3 0 360 arc closepath stroke\n",5.0 + 10.0 * city_x[i],5.0 + 10.0 * city_y[i]);
    fprintf(fp,"(%d) dup stringwidth pop 2 div neg %.2f add %.2f 1.4 sub moveto show\n",i,5.0 + 10.0 * city_x[i],5.0 + 10.0 * city_y[i]);
  }
  fprintf(fp,"showpage\n");
  fprintf(fp,"%%%%EOF\n");
  fclose(fp);
//fprintf(stderr,"run the command\n  epspdf %s\nto convert the encapsulated postscript file to pdf\n",file_name);
}

double distMatrix[max_n_cities][max_n_cities];

#if 0
double dist(int i,int j)
{
  double dx = city_x[i] - city_x[j];
  double dy = city_y[i] - city_y[j];
  return sqrt(dx * dx + dy * dy);
}
#else
# define dist(i,j) distMatrix[i][j]
#endif

double min_d;
int min_a[max_n_cities];

double max_d;
int max_a[max_n_cities];

void calculate_all_distances(unsigned int n)
{
    int i,j;
    double dx, dy;
    
    for(i = 0;i < n;i++){
		for(j = 0;j < n;j++){
		    dx = city_x[i] - city_x[j];
		    dy = city_y[i] - city_y[j];
			distMatrix[i][j] = sqrt(dx * dx + dy * dy);
		}
	}
}
void generate_all_permutations(unsigned int n,unsigned int m,int a[])  // int a[] is a synonym of int *a
{
  unsigned int i;
  if(m < n - 1)
  { // not yet at the end
    for(i = m;i < n;i++)
    {
      #define swap(i,j)  do { int t = a[i]; a[i] = a[j]; a[j] = t; } while(0)
      swap(i,m);                            // exchange a[i] with a[m]
      generate_all_permutations(n,m + 1,a); // recurse
      swap(i,m);                            // undo the exchange of a[i] with a[m]
      #undef swap
    }
  }
  else
  { // visit permutation
    double d = dist(a[n-1],a[0]);
    for(i = 0;i < n - 1;i++)
      d += dist(a[i],a[i+1]);
    if(d < min_d) { min_d = d; for(i = 0;i < n;i++) min_a[i] = a[i]; }
    if(d > max_d) { max_d = d; for(i = 0;i < n;i++) max_a[i] = a[i]; }
  }
}

double cpu_time(void)
{
	#if defined(_WIN32)
	FILETIME createTime;
	FILETIME exitTime;
	FILETIME kernelTime;
	FILETIME userTime;
	if ( GetProcessTimes( GetCurrentProcess( ),
		&createTime, &exitTime, &kernelTime, &userTime ) != -1 )
	{
		SYSTEMTIME userSystemTime;
		if ( FileTimeToSystemTime( &userTime, &userSystemTime ) != -1 )
			return (double)userSystemTime.wHour * 3600.0 +
				(double)userSystemTime.wMinute * 60.0 +
				(double)userSystemTime.wSecond +
				(double)userSystemTime.wMilliseconds / 1000.0;
	}
	#else
	  struct rusage r;

	  getrusage(RUSAGE_SELF,&r);
	  return (double)r.ru_utime.tv_sec + 0.000001 * (double)r.ru_utime.tv_usec;
	#endif
	return -1;
}

int main(void)
{
  int n_cities,i,a[max_n_cities];
  char file_name[32];
  double dt;
  
  rand_city_coords(83087);
  
  calculate_all_distances(max_n_cities);
  
  for(n_cities = 3;n_cities <= max_n_cities;n_cities++)
  {
    for(i = 0;i < n_cities;i++)
      a[i] = i;
    min_d = 1.0e+100;
    max_d = -1.0;
    
    dt = cpu_time();
    generate_all_permutations(n_cities,1,a);
    dt = cpu_time() - dt;
    printf("%d [%.3f]\n",n_cities,dt);
    printf("  min %.6f",min_d); for(i = 0;i < n_cities;i++) printf(" %d",min_a[i]); printf("\n");
    printf("  max %.6f",max_d); for(i = 0;i < n_cities;i++) printf(" %d",max_a[i]); printf("\n");
    sprintf(file_name,"tsp_%02d_min.ps",n_cities); print_tour(n_cities,min_a,file_name);
    sprintf(file_name,"tsp_%02d_max.ps",n_cities); print_tour(n_cities,max_a,file_name);
  }
  return 0;
}
