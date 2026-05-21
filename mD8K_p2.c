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

    // CÀLCUL DEL REPARTIMENT DE COLUMNES
    int base_count = N / size;
    int remainder = N % size;
    int N_local = base_count + (rank < remainder ? 1 : 0);

    int offset = 0;
    for (int p = 0; p < rank; p++) {
        offset += base_count + (p < remainder ? 1 : 0);
    }
    
    // RESERVA DE MEMÒRIA DINÀMICA
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
 
    // Matriu dispersa per matriu densa
    for(i = offset; i < offset + N_local; i++) {
        int col_idx = i - offset;
        for (k = 0; k < ND; k++) {
            C1_local[AD[k].i * N_local + col_idx] += AD[k].v * B[AD[k].j * N + i];
        }
    }

    // Matriu dispersa per matriu dispersa
    bzero(VBcol, sizeof(int) * N);
    for(i = offset; i < offset + N_local; i++) {
        int col_idx = i - offset;
        for (k = jBD[i]; k < jBD[i+1]; k++)
            VBcol[BD[k].i] = BD[k].v;
        
        for (k = 0; k < ND; k++) {
            C2_local[AD[k].i * N_local + col_idx] += AD[k].v * VBcol[AD[k].j];
        }
        for (j = 0; j < N; j++)
            VBcol[j] = 0;
    }
                
    // Compressió i càlcul de la Suma de CD
    neleC_local = 0;
    for (j=0; j<N; j++) VBcol[j] = VCcol[j] = 0;

    for(i = offset; i < offset + N_local; i++) {
        for (k = jBD[i]; k < jBD[i+1]; k++)
            VBcol[BD[k].i] = BD[k].v;
        
        for (k = 0; k < ND; k++)
            VCcol[AD[k].i] += AD[k].v * VBcol[AD[k].j];
            
        for (j = 0; j < N; j++) {
            VBcol[j] = 0;
            if (VCcol[j]) {
                CD[neleC_local].i = j;
                CD[neleC_local].j = i;
                CD[neleC_local].v = VCcol[j];
                VCcol[j] = 0;
                neleC_local++;
            }
        }
    }

    // Comprovacio MD x M -> M i MD x MD -> M (Cada un comprova les seves columnes locals)
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

    // Reduïm les sumes i els elements totals cap al procés 0
    MPI_Reduce(&Suma_local, &Suma_global, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&neleC_local, &neleC_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("\nNumero elements de la matriu dispersa C %d\n", neleC_global);   
        printf("Suma dels elements de C: %lld\n", Suma_global);
        fflush(stdout);
    }

    free(A); free(B); free(C_local); free(C1_local); free(C2_local);
    free(jBD); free(VCcol); free(VBcol);
    free(AD); free(BD); free(CD);

    MPI_Finalize();
    return 0;
}