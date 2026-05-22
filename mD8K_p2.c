#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <strings.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>

typedef struct {
    int i,j,v;
} tmd;

int *A, *B; 
int *C_local, *C1_local, *C2_local;
int *jBD, *VCcol, *VBcol;
tmd *AD, *BD, *CD;

long long Suma_local = 0, Suma_global = 0;
int neleC_local = 0, neleC_global = 0;
long long ops_local = 0, ops_global = 0;

int cmp_fil(const void *pa, const void *pb)
{
tmd * a = (tmd*)pa;
tmd * b = (tmd*)pb;

  if (a->i > b->i) return(1);
  else if (a->i < b->i) return (-1);
  else return (a->j - b->j);
}

int cmp_col(const void *pa, const void *pb)
{
tmd * a = (tmd*)pa;
tmd * b = (tmd*)pb;

  if (a->j > b->j) return(1);
  else if (a->j < b->j) return (-1);
  else return (a->i - b->i);
}

int main(int argc, char *argv[])
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) printf("ERROR: Falta la mida de la matriu N.\n");
        MPI_Finalize();
        return 1;
    }
    long N = atol(argv[1]);
    long ND = N * N / 100;

    int i,j,k;

    // REPARTIMENT DE LA FEINA (COLUMNES DE LES MATRIUS RESULTANTS)
    int base_count = N / size;
    int remainder = N % size;
    int N_local = base_count + (rank < remainder ? 1 : 0);

    int offset = 0;
    for (int p = 0; p < rank; p++) {
        offset += base_count + (p < remainder ? 1 : 0);
    }
    
    // RESERVA DE MEMÒRIA DINÀMICA AL HEAP (per evitar stack overflow amb N=8000)
    A = (int *)calloc(N * N, sizeof(int));
    B = (int *)calloc(N * N, sizeof(int));

    C_local  = (int *)calloc(N * N_local, sizeof(int));
    C1_local = (int *)calloc(N * N_local, sizeof(int));
    C2_local = (int *)calloc(N * N_local, sizeof(int));

    jBD = (int *)malloc((N + 1) * sizeof(int));
    VCcol = (int *)calloc(N, sizeof(int));
    VBcol = (int *)calloc(N, sizeof(int));

    AD = (tmd *)malloc(ND * sizeof(tmd));
    BD = (tmd *)malloc(ND * sizeof(tmd));
    CD = (tmd *)malloc(N * N * sizeof(tmd));
    
    // Generació a l'ombra amb llavor fixa per garantir coherència    
    srand(1);
    for(k=0;k<ND;k++)
    {
        AD[k].i=rand()%(N-1);
        AD[k].j=rand()%(N-1);
        AD[k].v=rand()%100+1;
        while (A[AD[k].i * N + AD[k].j]) {
            if(AD[k].i < AD[k].j)
                AD[k].i = (AD[k].i + 1)%N;
            else 
                AD[k].j = (AD[k].j + 1)%N;
        }
        A[AD[k].i * N + AD[k].j] = AD[k].v;
    }
    qsort(AD,ND,sizeof(tmd),cmp_fil);

    for(k=0;k<ND;k++)
    {
        BD[k].i=rand()%(N-1);
        BD[k].j=rand()%(N-1);
        BD[k].v=rand()%100+1;
        while (B[BD[k].i * N + BD[k].j]) {
            if(BD[k].i < BD[k].j)
                BD[k].i = (BD[k].i + 1)%N;
            else 
                BD[k].j = (BD[k].j + 1)%N;
        }
        B[BD[k].i * N + BD[k].j] = BD[k].v;
    }
    qsort(BD,ND,sizeof(tmd),cmp_col);
    
    k=0;
    for (j=0; j<N+1; j++)
     {
      while (k < ND && j>BD[k].j) k++;
      jBD[j] = k;
     }
 
    // Matriu dispersa per matriu densa (cache-friendly amb inversió de bucles)
    // Guardem en variables locals registrades els límits de C1
    // per evitar recalcular l'offset contínuament dins dels registres L1 de la CPU.
    for (k = 0; k < ND; k++) {
        int row_A = AD[k].i * N_local;
        int row_B = AD[k].j * N;
        int val_A = AD[k].v;
        int end_i = offset + N_local;
        for (i = offset; i < end_i; i++) {
            C1_local[row_A + (i - offset)] += val_A * B[row_B + i];
            ops_local++;
        }
    }

    // Matriu dispersa per matriu dispersa amb optimitzacions de localitat i neteja esparsa
    bzero(VBcol, sizeof(int) * N);

    for(i = offset; i < offset + N_local; i++) {
        // Reduïm operacions traient la resta d'índex local fora del bucle intern
        int col_idx = i - offset;
        
        // Expandir columna de B[*][i] a VBcol
        for (k = jBD[i]; k < jBD[i+1]; k++)
            VBcol[BD[k].i] = BD[k].v;
        
        for (k = 0; k < ND; k++) {
            int target_j = AD[k].j;
            if (VBcol[target_j] != 0) {
                C2_local[AD[k].i * N_local + col_idx] += AD[k].v * VBcol[target_j];
                ops_local++;
            }
        }
        
        // NETEJA ESPARSA: Només posem a zero el que hem modificat
        for (k = jBD[i]; k < jBD[i+1]; k++) {
            VBcol[BD[k].i] = 0;
        }
    }
                
    // COMPRESSIÓ: Generació de la llista CD amb neteja esparsa integrada
    neleC_local = 0;
    bzero(VBcol, sizeof(int) * N);
    bzero(VCcol, sizeof(int) * N);

    for(i = offset; i < offset + N_local; i++) {
        // Expandir columna de B[*][i] a VBcol
        for (k = jBD[i]; k < jBD[i+1]; k++)
            VBcol[BD[k].i] = BD[k].v;
        
        // Càlcul de la columna dispersa de C i neteja esparsa de VBcol al mateix temps
        for (k = 0; k < ND; k++) {
            int target_j = AD[k].j;
            if (VBcol[target_j] != 0) {
                VCcol[AD[k].i] += AD[k].v * VBcol[target_j];
            }
        }
        
        // Neteja esparsa de VBcol
        for (k = jBD[i]; k < jBD[i+1]; k++) {
            VBcol[BD[k].i] = 0;
        }
        
        // Compressió i neteja esparsa de VCcol (només si la posició té dades)
        for (k = 0; k < ND; k++) {
            int row_A = AD[k].i;
            if (VCcol[row_A] != 0) {
                CD[neleC_local].i = row_A;
                CD[neleC_local].j = i;
                CD[neleC_local].v = VCcol[row_A];
                VCcol[row_A] = 0;
                neleC_local++;
            }
        }
    }

    // COMPROVACIONS I REDUCCIONS 
    for (i = offset; i < offset + N_local; i++) {
        int col_idx = i - offset;
        for (j = 0; j < N; j++) {
            if (C2_local[j * N_local + col_idx] != C1_local[j * N_local + col_idx]) {
                printf("Diferencies C1 i C2 pos %d,%d: %d != %d (Rank %d)\n", j, i, C1_local[j * N_local + col_idx], C2_local[j * N_local + col_idx], rank);
            }
        }
    }

    // Comprovacio MD X MD -> M i MD x MD -> MD
    Suma_local = 0;
    for(k = 0; k < neleC_local; k++) {
        Suma_local += CD[k].v;
        
        // Calculem la posició dins del buffer local per a la columna global guardada a CD[k].j
        int col_idx = CD[k].j - offset;
        // Comprovem només si col_idx cau dins del rang local assignat
        if (col_idx >= 0 && col_idx < N_local) {
            if (CD[k].v != C1_local[CD[k].i * N_local + col_idx]) {
                printf("Diferencies C1 i CD a i:%d,j:%d,v%d, k:%d, vd:%d (Rank %d)\n", CD[k].i, CD[k].j, C1_local[CD[k].i * N_local + col_idx], k, CD[k].v, rank);
            }
        }
    }

    MPI_Reduce(&Suma_local, &Suma_global, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&neleC_local, &neleC_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ops_local, &ops_global, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("\nNumero elements de la matriu dispersa C %d\n", neleC_global);   
        printf("Suma dels elements de C %lld\n", Suma_global);
        printf("Total operacions multiplicacio %lld\n", ops_global);
        fflush(stdout);
    }

    free(A); free(B); free(C_local); free(C1_local); free(C2_local);
    free(jBD); free(VCcol); free(VBcol);
    free(AD); free(BD); free(CD);

    MPI_Finalize();
    return 0;
}