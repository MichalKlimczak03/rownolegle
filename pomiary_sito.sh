#!/bin/bash

# Skrypt do pomiaru czasu sita sekwencyjnego i OpenMP.
# Zaklada, ze w katalogu sa pliki:
# - sito_seq_poprawione.c
# - sito_omp_poprawione.c
# oraz program geng z pakietu nauty.
#
# Przyklad uruchomienia:
# chmod +x pomiary_sito.sh
# ./pomiary_sito.sh
#
# Wyniki zostana zapisane do:
# wyniki_sito.csv

set -e

SEQ_SRC="sito_seq_poprawione.c"
OMP_SRC="sito_omp_poprawione.c"

SEQ_EXE="sito_seq"
OMP_EXE="sito_omp"

OUT="wyniki_sito.csv"

# Liczby wierzcholkow grafow do testowania.
# Dla n=10 czas moze byc juz znacznie dluzszy.
N_VALUES="5 6 7 8 9"

# Liczby watkow dla OpenMP.
THREADS_VALUES="1 2 4 8"

# Rozmiar paczki dla wersji OpenMP.
BATCH_SIZE=4096

# Ile razy powtorzyc kazdy pomiar.
REPEATS=5

echo "Kompilacja programow..."

gcc -O3 "$SEQ_SRC" -o "$SEQ_EXE" -lm
gcc -O3 "$OMP_SRC" -o "$OMP_EXE" -fopenmp -lm

echo "algorytm,n,watki,rozmiar_paczki,powtorzenie,sprawdzone_grafy,znalezione_grafy_calkowite,czas_s,sredni_czas_na_graf_s" > "$OUT"

for N in $N_VALUES
do
    GRAPH_FILE="grafy_${N}.g6"

    echo "Generowanie grafow dla n=$N..."
    ./geng -c "$N" 2>/dev/null > "$GRAPH_FILE"

    echo "Pomiary sekwencyjne dla n=$N..."

    for REP in $(seq 1 $REPEATS)
    do
        RESULT=$(./"$SEQ_EXE" -q < "$GRAPH_FILE" 2>&1 >/dev/null)

        CHECKED=$(echo "$RESULT" | awk -F': ' '/sprawdzone_grafy/ {print $2}')
        FOUND=$(echo "$RESULT" | awk -F': ' '/znalezione_grafy_calkowite/ {print $2}')
        TIME=$(echo "$RESULT" | awk -F': ' '/czas_s/ {print $2}')
        AVG=$(echo "$RESULT" | awk -F': ' '/sredni_czas_na_graf_s/ {print $2}')

        echo "seq,$N,1,0,$REP,$CHECKED,$FOUND,$TIME,$AVG" >> "$OUT"
    done

    for THREADS in $THREADS_VALUES
    do
        echo "Pomiary OpenMP dla n=$N, watki=$THREADS..."

        for REP in $(seq 1 $REPEATS)
        do
            RESULT=$(./"$OMP_EXE" "$THREADS" "$BATCH_SIZE" -q < "$GRAPH_FILE" 2>&1 >/dev/null)

            CHECKED=$(echo "$RESULT" | awk -F': ' '/sprawdzone_grafy/ {print $2}')
            FOUND=$(echo "$RESULT" | awk -F': ' '/znalezione_grafy_calkowite/ {print $2}')
            TIME=$(echo "$RESULT" | awk -F': ' '/czas_s/ {print $2}')
            AVG=$(echo "$RESULT" | awk -F': ' '/sredni_czas_na_graf_s/ {print $2}')

            echo "openmp,$N,$THREADS,$BATCH_SIZE,$REP,$CHECKED,$FOUND,$TIME,$AVG" >> "$OUT"
        done
    done
done

echo ""
echo "Gotowe."
echo "Wyniki zapisano do pliku: $OUT"
echo ""
echo "Podglad:"
head "$OUT"
