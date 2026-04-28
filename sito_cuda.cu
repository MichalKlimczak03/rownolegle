/**
@file=sito_cuda.cu
Wersja CUDA bazująca na sicie8.cu z materiałów kursowych.
Jeden wątek CUDA sprawdza jeden graf.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <cuda_runtime.h>

#define BUFSIZE 1024
#define NMAX 20
#define MAX_GRAPHS 1024
#define THREADS 256

#define CUDA_CHECK(call) do { \
  cudaError_t err = call; \
  if (err != cudaSuccess) { \
    fprintf(stderr, "CUDA error: %s (%s:%d)\n", cudaGetErrorString(err), __FILE__, __LINE__); \
    exit(EXIT_FAILURE); \
  } \
} while (0)

double get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec / 1000000.0;
}

__global__ void test(char * BUFFOR, int * valid, int n_graphs) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n_graphs) return;

  valid[idx] = 0;

  char *SBUFFER = BUFFOR + idx * BUFSIZE;

  int i,j,k,k3,k4,L,L1,z;
  double eps,g,h,ma,mn,norm,s,t,u,w;
  int cond;
  double d[NMAX+1], e[NMAX+1], e2[NMAX+1], Lb[NMAX+1];
  double x[NMAX+1];
  double a[NMAX*(NMAX-1)/2 + NMAX + 1];
  int n;
  int bit, poz, poz2;

  bit = 32;
  poz = 1;
  poz2 = 1;
  n = SBUFFER[0] - 63;

  a[0] = 0.0;

  for (i = 0; i < n; i++)
   for (j = 0; j<=i; j++)
    {
      if (i==j) {a[poz2++] = 0; }
      else {
        if (bit == 0) { bit = 32;  poz++; }
        if ((SBUFFER[poz] - 63) & bit)
                { a[poz2++] = 1; }
               else
                { a[poz2++] = 0; }
        bit = bit >> 1;
      }
    }

  int k1 = 1;
  int k2 = n;

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

    w = ma;
    eps = 7.28e-17*norm;

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
         if (floor(w+10e-5)<s-10e-5) return;

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

      if  (!(( ceil(u) - u  < 10e-5 ) || ( u - floor(u) < 10e-5 ))) { return; };
    }
  }

  valid[idx] = 1;
}

int main(int argc, char *argv[])
{
  char (*BUFFOR)[BUFSIZE];
  int *valid;

  char *cuda_bufor;
  int *cuda_valid;

  unsigned long long checked = 0;
  unsigned long long found = 0;

  BUFFOR = (char (*)[BUFSIZE])malloc(MAX_GRAPHS * BUFSIZE);
  valid = (int*)malloc(MAX_GRAPHS * sizeof(int));

  CUDA_CHECK(cudaMalloc((void**)&cuda_bufor, MAX_GRAPHS * BUFSIZE));
  CUDA_CHECK(cudaMalloc((void**)&cuda_valid, MAX_GRAPHS * sizeof(int)));

  double start = get_time();

  while (1)
  {
    int n_graphs = 0;

    while (n_graphs < MAX_GRAPHS && fgets(BUFFOR[n_graphs], BUFSIZE-1, stdin))
    {
      n_graphs++;
    }

    if (n_graphs == 0)
      break;

    checked += n_graphs;

    CUDA_CHECK(cudaMemcpy(cuda_bufor, BUFFOR, n_graphs * BUFSIZE, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(cuda_valid, 0, n_graphs * sizeof(int)));

    int blocks = (n_graphs + THREADS - 1) / THREADS;

    test<<<blocks, THREADS>>>(cuda_bufor, cuda_valid, n_graphs);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(valid, cuda_valid, n_graphs * sizeof(int), cudaMemcpyDeviceToHost));

    for (int i = 0; i < n_graphs; i++)
    {
      if (valid[i])
      {
        found++;
        printf("%s", BUFFOR[i]);
      }
    }
  }

  double end = get_time();
  double time = end - start;

  fprintf(stderr, "\n--- STATYSTYKI CUDA ---\n");
  fprintf(stderr, "Watki na blok: %d\n", THREADS);
  fprintf(stderr, "Rozmiar batcha: %d\n", MAX_GRAPHS);
  fprintf(stderr, "Sprawdzone grafy: %llu\n", checked);
  fprintf(stderr, "Znalezione grafy calkowite: %llu\n", found);
  fprintf(stderr, "Czas wykonania: %.6f s\n", time);

  if (checked > 0)
    fprintf(stderr, "Sredni czas na graf: %.12f s\n", time / checked);

  CUDA_CHECK(cudaFree(cuda_bufor));
  CUDA_CHECK(cudaFree(cuda_valid));

  free(BUFFOR);
  free(valid);

  return EXIT_SUCCESS;
}