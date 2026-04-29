#!/bin/bash

GENG="../nauty2_8_9/geng"
N=16
K=88

SIZES=(100 200 400 800 1600 3200 6400 12800 25600 51200 102400)

SEQ="./sito_seq"
OMP="./sito_omp"

OMP_THREADS=4
BATCH_SIZE=4096

OUT="wyniki_pomiarow.csv"

echo "algorytm,grafy,czas,znalezione" > $OUT

for SIZE in "${SIZES[@]}"
do
    echo "Test dla $SIZE grafów"

    # --- SEQ ---
    $GENG -c $N $K:$K | head -n $SIZE | $SEQ -q > /dev/null 2> tmp_seq.txt

    TIME=$(grep "czas_s" tmp_seq.txt | awk -F': ' '{print $2}')
    FOUND=$(grep "znalezione_grafy_calkowite" tmp_seq.txt | awk -F': ' '{print $2}')

    echo "seq,$SIZE,$TIME,$FOUND" >> $OUT

    # --- OMP ---
    $GENG -c $N $K:$K | head -n $SIZE | $OMP $OMP_THREADS $BATCH_SIZE -q > /dev/null 2> tmp_omp.txt

    TIME=$(grep "czas_s" tmp_omp.txt | awk -F': ' '{print $2}')
    FOUND=$(grep "znalezione_grafy_calkowite" tmp_omp.txt | awk -F': ' '{print $2}')

    echo "openmp,$SIZE,$TIME,$FOUND" >> $OUT
done

rm -f tmp_seq.txt tmp_omp.txt

echo "Gotowe. Wyniki zapisane w $OUT"
