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

# ==============================================================================
# TRUCO HPC: Autoreserva de SLURM
# ==============================================================================
if [ -z "$SLURM_JOB_ID" ]; then
    echo "========================================================================================================================="
    echo " Solicitant reserva exclusiva de 16 nodes de cop. Esperant a SLURM..."
    echo "========================================================================================================================="
    salloc -p jutjat -N 16 --exclusive bash "$0" "$@"
    exit $?
fi

# ==============================================================================
# A PARTIR D'AQUÍ: Ja som amos de 16 nodes exclusius. Tot serà instantani.
# ==============================================================================

# Temps seqüencial marcat per l'enunciat per a N=8000
SEQ_JUTGES=19.5
ROW_FMT="%-20s | %-12s | %-16s | %-8s | %-14s | %-12s | %-14s |\n"

echo "========================================================================================================================="
echo " PROVES MPI EN MÀQUINES JUTGES (Temps seq base: 19.5s | N=$N_SIZE)"
echo "========================================================================================================================="

echo "Generant Referencia Mestra (1 proceso) per validació estricta..."
srun -N 1 -n 1 --ntasks-per-node=1 --distribution=block:block ./$EXE $N_SIZE > full_ref.out 2>/dev/null
grep 'Suma dels elements' full_ref.out > ref.out
rm -f full_ref.out

if [ ! -s ref.out ]; then
    echo "ERROR: No s'ha pogut generar la referencia. Revisa el codi C."
    exit 1
fi
echo "Referència generada correctament."
echo "-------------------------------------------------------------------------------------------------------------------------"

sum_inv_S_time=0
count=0

printf "$ROW_FMT" "Nodes / PPN (Procs)" "Time (s)" "Speedup (Time)" "Valid" "Suma C" "Elements CD" "Operacions"
echo "-------------------------------------------------------------------------------------------------------------------------"

CONFIGS=(
    "2 1" "4 1" "8 1"
    "16 1" "16 2" "16 4"
    "16 8" "16 16" "16 32"
)

for config in "${CONFIGS[@]}"; do
    read -r nodes ppn <<< "$config"
    procs=$((nodes * ppn))
    
    # Overcommit
    if [ "$ppn" -gt 16 ]; then
        OVERCOMMIT_FLAG="--overcommit"
    else
        OVERCOMMIT_FLAG=""
    fi
    
    # Tornem a ficar el temps DINS de l'srun perquè no ens mesuri l'overhead de SLURM
    OUTPUT=$(srun -N $nodes -n $procs --ntasks-per-node=$ppn --distribution=block:block $OVERCOMMIT_FLAG /usr/bin/time -f "TEMPS_REAL:%e" ./$EXE $N_SIZE 2>&1)
    
    # Extraiem les dades forçant a agafar només la primera línia (seguretat anti-errors)
    TIME_REAL=$(echo "$OUTPUT" | grep "TEMPS_REAL:" | cut -d':' -f2 | grep -oE '[0-9.]+' | head -n 1)
    SUMA_TEXT=$(echo "$OUTPUT" | grep "Suma dels elements de C" | grep -oE '[0-9]+' | head -n 1)
    NELEC_TEXT=$(echo "$OUTPUT" | grep "Numero elements de la matriu dispersa C" | grep -oE '[0-9]+' | head -n 1) 
    OPS_TEXT=$(echo "$OUTPUT" | grep "Total operacions multiplicacio" | grep -oE '[0-9]+' | head -n 1)
    
    if [ -z "$NELEC_TEXT" ]; then
        NELEC_TEXT="ERR"
    fi
    
    if [ -z "$OPS_TEXT" ]; then
        OPS_TEXT="ERR"
    fi

    # Verificacio estricta contra la referencia mestre
    echo "$OUTPUT" | grep "Suma dels elements de C" | grep -oE '[0-9]+' | tail -n 1 > test.out
    if [ -s ref.out ] && [ -s test.out ] && diff -q <(grep -oE '[0-9]+' ref.out) test.out > /dev/null; then 
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
     
    printf "$ROW_FMT" "$nodes / $ppn ($procs)" "${TIME_REAL_FMT}" "${SPEEDUP_TIME}" "$VALID" "${SUMA_TEXT}" "${NELEC_TEXT}" "${OPS_TEXT}"
    ((count++))
done

# Hem canviat el 'bc' per 'awk' per evitar l'error de comanda no trobada
if [ "$count" -gt 0 ] && [ "$(awk -v sum="$sum_inv_S_time" 'BEGIN { print (sum > 0 ? 1 : 0) }')" -eq 1 ]; then
    HMEAN_TIME=$(awk -v n="$count" -v sum_inv="$sum_inv_S_time" 'BEGIN { printf "%.5f", n/sum_inv }')
else
    HMEAN_TIME="ERR"
fi

echo "-------------------------------------------------------------------------------------------------------------------------"
echo "-> MITJANA HARMÒNICA SPEEDUP (Time): $HMEAN_TIME"
echo "========================================================================================================================="

# Neteja de fitxers temporals finals
rm -f $EXE ref.out test.out