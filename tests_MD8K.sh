#!/bin/bash
export LC_NUMERIC=C

SRC="mD8K_p2.c"
EXE="mD8K_p2"
N_SIZE=8000

mpicc -O3 $SRC -o $EXE

if [ $? -ne 0 ]; then
    echo "ERROR: No s'ha pogut compilar $SRC."
    exit 1
fi

# Temps sequencial marcat per l'enunciat per a N=8000
SEQ_JUTGES=19.5
ROW_FMT="%-20s | %-12s | %-16s | %-8s | %-14s |\n"

echo "================================================================================================"
echo " PROVES MPI EN MÀQUINES JUTGES (Temps seq base: 19.5s | N=$N_SIZE)"
echo "================================================================================================"

# ==============================================================================
echo "Generant Referencia Mestra (1 proceso) per validació estricta..."
salloc -p jutjat -N 1 --exclusive srun -N 1 -n 1 --ntasks-per-node=1 --distribution=block:block ./$EXE $N_SIZE > full_ref.out 2>/dev/null
grep 'Suma dels elements' full_ref.out > ref.out
rm -f full_ref.out

if [ ! -s ref.out ]; then
    echo "ERROR: No s'ha pogut generar la referencia. Revisa el codi C."
    exit 1
fi
echo "Referència generada correctament."
echo "------------------------------------------------------------------------------------------------"

sum_inv_S_time=0
count=0

printf "$ROW_FMT" "Nodes / PPN (Procs)" "Time (s)" "Speedup (Time)" "Valid" "Suma C"
echo "------------------------------------------------------------------------------------------------"

CONFIGS=(
    "2 1" "4 1" "8 1"
    "16 1" "16 2" "16 4"
    "16 8" "16 16" "16 32"
)

for config in "${CONFIGS[@]}"; do
    read -r nodes ppn <<< "$config"
    procs=$((nodes * ppn))
    
    if [ "$ppn" -gt "$nodes" ]; then
        OVERCOMMIT_FLAG="--overcommit"
    else
        OVERCOMMIT_FLAG=""
    fi
    
    # Execucio controlada enviant la sortida a run.out
    salloc -p jutjat -N $nodes --exclusive srun -N $nodes -n $procs --ntasks-per-node=$ppn --distribution=block:block $OVERCOMMIT_FLAG /usr/bin/time -f "TEMPS_REAL:%e" ./$EXE $N_SIZE > run.out 2>&1
    OUTPUT=$(cat run.out)
    rm -f run.out
    
    # Extraiem el temps real i la suma final
    TIME_REAL=$(echo "$OUTPUT" | grep "TEMPS_REAL:" | cut -d':' -f2 | grep -oE '[0-9.]+')
    SUMA_TEXT=$(echo "$OUTPUT" | grep "Suma dels elements de C" | grep -oE '[0-9]+')
    
    # Verificacio estricta contra la referencia mestre
    echo "$OUTPUT" | grep "Suma dels elements" > test.out
    if diff -q ref.out test.out > /dev/null; then
        VALID="OK"
    else
        VALID="MAL"
    fi     
    
    if [ -z "$TIME_REAL" ] || [ "$TIME_REAL" == "" ]; then
        TIME_REAL_FMT="ERR"
        SPEEDUP_TIME="0.00000"
        SUMA_TEXT="ERROR"
        VALID="ERR"
    else
        SPEEDUP_TIME=$(awk -v seq="$SEQ_JUTGES" -v par="$TIME_REAL" 'BEGIN { printf "%.5f", seq/par }')
        sum_inv_S_time=$(awk -v sum="$sum_inv_S_time" -v seq="$SEQ_JUTGES" -v par="$TIME_REAL" 'BEGIN { print sum + (par/seq) }')
        TIME_REAL_FMT=$(awk -v real="$TIME_REAL" 'BEGIN { printf "%.3f", real }')
    fi
    
    printf "$ROW_FMT" "$nodes / $ppn ($procs)" "${TIME_REAL_FMT}" "${SPEEDUP_TIME}" "$VALID" "${SUMA_TEXT}"
    ((count++))
done

HMEAN_TIME=$(awk -v n="$count" -v sum_inv="$sum_inv_S_time" 'BEGIN { printf "%.5f", n/sum_inv }')

echo "------------------------------------------------------------------------------------------------"
echo "-> MITJANA HARMÒNICA SPEEDUP (Time): $HMEAN_TIME"
echo "================================================================================================"

# Neteja de fitxers temporals finals
rm -f $EXE ref.out test.out