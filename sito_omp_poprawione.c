/**
@file sito_omp_poprawione.c

Wersja OpenMP sita grafow calkowitych.
Bazuje bezposrednio na kodzie prowadzacego, ale dodaje:
- jawne ustawienie liczby watkow,
- jawne ustawienie rozmiaru paczki,
- tryb -q do uczciwszych pomiarow bez kosztu wypisywania grafow.

Kompilacja:
gcc -O3 sito_omp_poprawione.c -o sito_omp -fopenmp -lm

Uruchomienie standardowe:
./geng -c 9 2>/dev/null | ./sito_omp 4 4096 | wc -l

Uruchomienie do pomiarow:
./geng -c 9 2>/dev/null > grafy9.g6
./sito_omp 4 4096 -q < grafy9.g6

Argumenty:
1. liczba watkow, np. 4
2. rozmiar paczki grafow, np. 4096
-q oznacza brak wypisywania grafow na stdout
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define NMAX 20
#define BUFSIZE 1024
#define DEFAULT_BATCH_SIZE 4096

/**
@brief Testowanie czy graf jest calkowity
@param BUFOR - lancuch tekstowy kodujacy graf w formacie graph6
@return 1 jesli graf jest calkowity, 0 w przeciwnym razie
*/
int eigensymmatrix(char *BUFOR)
{
  int i,j,k,k3,k4,L,L1,z;
  long double eps,g,h,ma,mn,norm,s,t,u,w;
  int cond;
  long double d[NMAX+1], e[NMAX+1], e2[NMAX+1], Lb[NMAX+1];
  long double x[NMAX+1];
  long double a[NMAX*(NMAX-1)/2 + NMAX + 1];
  int n;
  int bit, poz, poz2;

  bit = 32;
  poz = 1;
  poz2 = 1;
  n = BUFOR[0] - 63;

  if (n < 1 || n > NMAX) return 0;

  a[0] = 0.0;

  for (i = 0; i < n; i++)
   for (j = 0; j <= i; j++)
    {
        if (i == j) {
            a[poz2++] = 0;
        }
        else {
            if (bit == 0) {
                bit = 32;
                poz++;
            }

            if ((BUFOR[poz] - 63) & bit)
                a[poz2++] = 1;
            else
                a[poz2++] = 0;

            bit = bit >> 1;
        }
    }

  int k1 = 1;
  int k2 = n;

  i = 0;
  for (L=1; L<=n; L++) {
      i += L;
      d[L] = a[i];
  }

  for (L=n; L>=2; L--)
   {
    i--;
    j = i;
    h = a[j];
    s = 0;

    for (k=L-2; k>=1; k--) {
        i--;
        g = a[i];
        s += g*g;
    }

    i--;

    if (s == 0) {
        e[L] = h;
        e2[L] = h*h;
        a[j] = 0.0;
    }
    else
    {
      s += h*h;
      e2[L] = s;
      g = sqrt(s);
      if (h >= 0.0) g = -g;

      e[L] = g;
      s = 1.0 / (s - h*g);
      a[j] = h - g;
      h = 0.0;
      L1 = L - 1;
      k3 = 1;

      for (j=1; j<=L1; j++)
       {
         k4 = k3;
         g = 0;

         for (k=1; k<=L1; k++) {
             g += a[k4] * a[i+k];
             if (k < j) z = 1;
             else z = k;
             k4 += z;
         }

         k3 += j;
         g *= s;
         e[j] = g;
         h += a[i+j] * g;
       }

      h *= 0.5 * s;
      k3 = 1;

      for (j=1; j<=L1; j++)
       {
         s = a[i+j];
         g = e[j] - h*s;
         e[j] = g;

         for (k=1; k<=j; k++) {
             a[k3] += -s*e[k] - a[i+k]*g;
             k3++;
         }
       }
    }

    h = d[L];
    d[L] = a[i+L];
    a[i+L] = h;
   }

  h = d[1];
  d[1] = a[1];
  a[1] = h;
  e[1] = 0.0;
  e2[1] = 0.0;
  s = d[n];

  t = fabs(e[n]);
  mn = s - t;
  ma = s + t;

  for (i=n-1; i>=1; i--)
   {
    u = fabs(e[i]);
    h = t + u;
    t = u;
    s = d[i];

    u = s - h;
    if (u < mn) mn = u;

    u = s + h;
    if (u > ma) ma = u;
   }

  for (i=1; i<=n; i++) {
      Lb[i] = mn;
      x[i] = ma;
  }

  norm = fabs(mn);
  s = fabs(ma);
  if (s > norm) norm = s;

  w = ma;
  eps = 7.28e-17 * norm;

  for (k=k2; k>=k1; k--)
   {
    s = mn;
    i = k;

    do {
       cond = 0;
       g = Lb[i];

       if (s < g)
           s = g;
       else {
           i--;
           if (i >= k1) cond = 1;
       }
    } while (cond);

    g = x[k];
    if (w > g) w = g;

    while (w - s > 2.91e-16 * (fabs(s) + fabs(w)) + eps)
     {
       if (floor(w + 10e-5) < s - 10e-5)
           return 0;

       L1 = 0;
       g = 1.0;
       t = 0.5 * (s + w);

       for (i=1; i<=n; i++)
        {
          if (g != 0)
              g = e2[i] / g;
          else
              g = fabs(6.87e15 * e[i]);

          g = d[i] - t - g;

          if (g < 0)
              L1++;
        }

       if (L1 < k1) {
           s = t;
           Lb[k1] = s;
       }
       else {
         if (L1 < k) {
             s = t;
             Lb[L1+1] = s;

             if (x[L1] > t)
                 x[L1] = t;
         }
         else {
             w = t;
         }
       }
    }

    u = 0.5 * (s + w);
    x[k] = u;

    if (!((ceil(u) - u < 10e-5) || (u - floor(u) < 10e-5)))
        return 0;
   }

  return 1;
}

int main(int argc, char *argv[])
{
  int threads = omp_get_max_threads();
  int batch_size = DEFAULT_BATCH_SIZE;
  int quiet = 0;

  if (argc > 1 && argv[1][0] != '-')
      threads = atoi(argv[1]);

  if (argc > 2 && argv[2][0] != '-')
      batch_size = atoi(argv[2]);

  for (int i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-q") == 0)
          quiet = 1;
  }

  if (threads < 1)
      threads = 1;

  if (batch_size < 1)
      batch_size = DEFAULT_BATCH_SIZE;

  omp_set_num_threads(threads);

  char **bufory = (char **)malloc((size_t)batch_size * sizeof(char *));
  int *valid = (int *)malloc((size_t)batch_size * sizeof(int));

  if (bufory == NULL || valid == NULL) {
      fprintf(stderr, "Blad alokacji pamieci.\n");
      free(bufory);
      free(valid);
      return EXIT_FAILURE;
  }

  for (int i = 0; i < batch_size; i++) {
      bufory[i] = (char *)malloc(BUFSIZE * sizeof(char));

      if (bufory[i] == NULL) {
          fprintf(stderr, "Blad alokacji pamieci.\n");

          for (int j = 0; j < i; j++)
              free(bufory[j]);

          free(bufory);
          free(valid);
          return EXIT_FAILURE;
      }
  }

  unsigned long long checked = 0;
  unsigned long long found = 0;

  double start = omp_get_wtime();

  while (1)
  {
      int graphs = 0;

      for (graphs = 0; graphs < batch_size; graphs++) {
          if (fgets(bufory[graphs], BUFSIZE, stdin) == NULL)
              break;

          valid[graphs] = 0;
      }

      if (graphs == 0)
          break;

      checked += (unsigned long long)graphs;

      #pragma omp parallel for schedule(static) reduction(+:found)
      for (int gid = 0; gid < graphs; gid++) {
          if (eigensymmatrix(bufory[gid])) {
              valid[gid] = 1;
              found++;
          }
      }

      if (!quiet) {
          for (int gid = 0; gid < graphs; gid++) {
              if (valid[gid])
                  printf("%s", bufory[gid]);
          }
      }

      if (graphs < batch_size)
          break;
  }

  double end = omp_get_wtime();
  double elapsed = end - start;

  fprintf(stderr, "tryb: openmp\n");
  fprintf(stderr, "watki: %d\n", threads);
  fprintf(stderr, "rozmiar_paczki: %d\n", batch_size);
  fprintf(stderr, "sprawdzone_grafy: %llu\n", checked);
  fprintf(stderr, "znalezione_grafy_calkowite: %llu\n", found);
  fprintf(stderr, "czas_s: %.6f\n", elapsed);

  if (checked > 0)
      fprintf(stderr, "sredni_czas_na_graf_s: %.12f\n", elapsed / checked);

  for (int i = 0; i < batch_size; i++)
      free(bufory[i]);

  free(bufory);
  free(valid);

  return EXIT_SUCCESS;
}
