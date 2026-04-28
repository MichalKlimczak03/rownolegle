/**
@file=sito_omp.c
Wersja OpenMP bazująca na sicie z materiałów kursowych.
Grafy są przetwarzane partiami, równolegle na CPU.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NMAX 20
#define BUFSIZE 1024
#define MAX_GRAPHS 1024

typedef int ta[NMAX][NMAX];

int eigensymmatrix(int n, long double* a, int k1, int k2, long double* x)
{
  int i,j,k,k3,k4,L,L1,z;
  long double eps,g,h,ma,mn,norm,s,t,u,w;
  int cond;
  long double d[NMAX], e[NMAX], e2[NMAX], Lb[NMAX];

  if ((1<=k1) && (k1<=k2) && (k2<=n))
   {
    i = 0;
    for (L=1;L<=n;L++) { i += L; d[L] = a[i]; }
    for (L=n;L>=2;L--)
     {
      i--; j = i; h = a[j]; s = 0;
      for (k=L-2;k>=1;k--) { i--; g = a[i]; s += g*g; }
      i--;
      if (s == 0) { e[L] = h; e2[L] = h*h; a[j] = 0.0; }
       else
        {
          s += h*h; e2[L] = s; g = sqrt(s); if (h>=0.0) g=-g;
          e[L] = g;
          s = 1.0 / (s-h*g);
          a[j] = h - g; h = 0.0; L1 = L - 1; k3 = 1;
          for (j=1;j<=L1;j++)
           {
             k4 = k3; g = 0;
             for (k=1;k<=L1;k++) { g +=a[k4]*a[i+k];
                                   if (k<j)  z = 1; else z = k;
                                   k4 += z; }
             k3 += j; g *= s; e[j] = g; h += a[i+j]*g;
           }
          h *= 0.5*s; k3 = 1;
          for (j=1;j<=L1;j++)
           {
             s = a[i+j]; g = e[j]-h*s; e[j] = g;
             for (k=1;k<=j;k++) { a[k3] += -s*e[k]-a[i+k]*g; k3++; }
           }
        }
      h = d[L]; d[L] = a[i+L]; a[i+L] = h;
     }
    h = d[1]; d[1] = a[1]; a[1] = h; e[1] = 0.0; e2[1] = 0.0; s = d[n];
    t = fabs(e[n]); mn = s - t; ma = s + t;
    for (i=n-1;i>=1;i--)
     {
      u = fabs(e[i]); h = t + u; t = u; s = d[i]; u = s - h;
      if (u < mn) mn = u;
      u = s + h;
      if (u > ma) ma = u;
     }
    for (i=1;i<=n;i++) { Lb[i] = mn; x[i] = ma; }
    norm = fabs(mn); s = fabs(ma);
    if (s>norm) norm = s;
    w = ma; eps = 7.28e-17*norm;
    for (k=k2;k>=k1;k--)
     {
      s = mn; i = k;
      do {cond = 0; g = Lb[i];
         if (s < g) s = g; else { i--; if (i>=k1) cond = 1; }
      } while (cond);
      g = x[k];
      if (w>g) w = g;
      while (w-s>2.91e-16*(fabs(s)+fabs(w))+eps)
       {
         if (floor(w+10e-5)<s-10e-5) return 0;
         L1 = 0; g = 1.0; t = 0.5*(s+w);
         for (i=1;i<=n;i++)
          {
            if (g!=0)  g = e2[i] / g; else g = fabs(6.87e15*e[i]);
            g = d[i]-t-g;
            if (g<0) L1++;
          }
         if (L1<k1) { s = t; Lb[k1] = s; }
          else
           { if (L1<k)
               {
                 s = t; Lb[L1+1] = s;
                 if (x[L1]>t) x[L1] = t;
               }
              else w = t;
           }
      }
      u = 0.5*(s+w); x[k] = u;
      if  (!(( ceil(u) - u  < 10e-5 ) || ( u - floor(u) < 10e-5 ))) { return 0; };
    }
  }
  return 1;
}

void AToa(int N, ta A, long double * a) {
 int poz,i,j;
 for (i = 0 ; i < N ; i++) A[i][i] = 0;
 a[0] = 0.0;
 poz = 1;
 for ( i = 0 ; i < N; i++)
 {
  for ( j = 0; j <= i ; j++)
   { if (A[i][j]) a[poz++]=1.0; else a[poz++]=0.0;}
 }
}

void BMKdecode(char * BUFFOR, int *N, ta A)
{
  int bit, poz, i, j;
  bit = 32;
  poz = 1;

  *N = BUFFOR[0] - 63;

  for (i = 1; i < *N; i++)
   for (j = 0; j< i; j++)
    {
      if (bit == 0) { bit = 32;  poz++; }
      if ((BUFFOR[poz] - 63) & bit)
                { A[i][j] = A[j][i] = 1; }
               else
                { A[i][j] = A[j][i] = 0; }
      bit = bit >> 1;
    }
}

int main(int argc, char *argv[])
{
  char BUFFOR[MAX_GRAPHS][BUFSIZE];
  int valid[MAX_GRAPHS];

  unsigned long long checked = 0;
  unsigned long long found = 0;

  int threads = omp_get_max_threads();

  if (argc > 1)
    threads = atoi(argv[1]);

  omp_set_num_threads(threads);

  double start = omp_get_wtime();

  while (1)
  {
    int n_graphs = 0;

    while (n_graphs < MAX_GRAPHS && fgets(BUFFOR[n_graphs], BUFSIZE, stdin))
    {
      valid[n_graphs] = 0;
      n_graphs++;
    }

    if (n_graphs == 0)
      break;

    checked += n_graphs;

    #pragma omp parallel for schedule(dynamic) reduction(+:found)
    for (int i = 0; i < n_graphs; i++)
    {
      ta A;
      int N;
      long double a[NMAX*NMAX];
      long double x[NMAX+1];

      BMKdecode(BUFFOR[i], &N, A);
      AToa(N, A, a);

      if (eigensymmatrix(N, a, 1, N, x))
      {
        valid[i] = 1;
        found++;
      }
    }

    for (int i = 0; i < n_graphs; i++)
    {
      if (valid[i])
        printf("%s", BUFFOR[i]);
    }
  }

  double end = omp_get_wtime();
  double time = end - start;

  fprintf(stderr, "\n--- STATYSTYKI OPENMP ---\n");
  fprintf(stderr, "Watki: %d\n", threads);
  fprintf(stderr, "Sprawdzone grafy: %llu\n", checked);
  fprintf(stderr, "Znalezione grafy calkowite: %llu\n", found);
  fprintf(stderr, "Czas wykonania: %.6f s\n", time);

  if (checked > 0)
    fprintf(stderr, "Sredni czas na graf: %.12f s\n", time / checked);

  return EXIT_SUCCESS;
}