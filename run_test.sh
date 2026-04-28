#!/bin/bash

GENG="../nauty2_8_9/geng"
N=16
K=88

SIZES=(1000 5000 10000 50000 100000)

SEQ="./sito_seq"
OMP="./sito_omp"
CUDA="./sito_cuda"

OMP_THREADS=4

OUT="wyniki_pomiarow.csv"

echo "algorytm,grafy,czas,znalezione" > $OUT

for SIZE in "${SIZES[@]}"
do
    echo "Test dla $SIZE grafów"

    $GENG -c $N $K:$K | head -n $SIZE | $SEQ > /dev/null 2> tmp_seq.txt
    TIME=$(grep "Czas wykonania" tmp_seq.txt | awk '{print $(NF-1)}')
    FOUND=$(grep "Znalezione" tmp_seq.txt | awk '{print $NF}')
    echo "seq,$SIZE,$TIME,$FOUND" >> $OUT

    $GENG -c $N $K:$K | head -n $SIZE | $OMP $OMP_THREADS > /dev/null 2> tmp_omp.txt
    TIME=$(grep "Czas wykonania" tmp_omp.txt | awk '{print $(NF-1)}')
    FOUND=$(grep "Znalezione" tmp_omp.txt | awk '{print $NF}')
    echo "openmp,$SIZE,$TIME,$FOUND" >> $OUT

    $GENG -c $N $K:$K | head -n $SIZE | $CUDA > /dev/null 2> tmp_cuda.txt
    TIME=$(grep "Czas wykonania" tmp_cuda.txt | awk '{print $(NF-1)}')
    FOUND=$(grep "Znalezione" tmp_cuda.txt | awk '{print $NF}')
    echo "cuda,$SIZE,$TIME,$FOUND" >> $OUT
done

rm -f tmp_seq.txt tmp_omp.txt tmp_cuda.txt

echo "Gotowe. Wyniki zapisane w $OUT"