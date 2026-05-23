#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <strings.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>

typedef struct {
    int i,j,v;
} __attribute__((aligned(16))) tmd;

int *A, *B; 
int *C_local, *C1_local, *C2_local;
int *jBD, *VCcol, *VBcol;
tmd *AD, *BD, *CD;

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
    long long Suma_local = 0, Suma_global = 0;
    int neleC_local = 0, neleC_global = 0;
    long long ops_local = 0, ops_global = 0;

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

    // Càlcul de distribució de columnes per procés
    int base_count = N / size;
    int remainder = N % size;
    int N_local = base_count + (rank < remainder ? 1 : 0);

    int offset = 0;
    for (int p = 0; p < rank; p++) {
        offset += base_count + (p < remainder ? 1 : 0);
    }
    
    // Assignació de memòria per a estructures de dades
    // Reserva de memòria avançada amb posix_memalign.
    // Forcem que l'adreça d'inici de cada matriu sigui múltiple de 64 bytes.
    posix_memalign((void**)&A, 64, N * N * sizeof(int)); 
    bzero(A, N * N * sizeof(int));
    
    posix_memalign((void**)&B, 64, N * N * sizeof(int)); 
    bzero(B, N * N * sizeof(int));

    posix_memalign((void**)&C_local, 64, N * N_local * sizeof(int)); 
    bzero(C_local, N * N_local * sizeof(int));
    
    posix_memalign((void**)&C1_local, 64, N * N_local * sizeof(int)); 
    bzero(C1_local, N * N_local * sizeof(int));
    
    posix_memalign((void**)&C2_local, 64, N * N_local * sizeof(int)); 
    bzero(C2_local, N * N_local * sizeof(int));

    posix_memalign((void**)&jBD, 64, (N + 1) * sizeof(int));
    
    posix_memalign((void**)&VCcol, 64, N * sizeof(int)); 
    bzero(VCcol, N * sizeof(int));
    
    posix_memalign((void**)&VBcol, 64, N * sizeof(int)); 
    bzero(VBcol, N * sizeof(int));

    posix_memalign((void**)&AD, 64, ND * sizeof(tmd));
    posix_memalign((void**)&BD, 64, ND * sizeof(tmd));
    posix_memalign((void**)&CD, 64, N * N * sizeof(tmd));
    
    // Inicialització pseudoaleatòria de les matrius disperses
    srand(1);
    for(k=0; k<ND; k++)
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

    for(k=0; k<ND; k++)
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
 
    // 1. Matriu Dispersa per Matriu Densa (MD x M -> M)
    long long ops_phase1 = 0;
    for (k = 0; k < ND; k++) {
        int row_A = AD[k].i * N_local;
        int row_B = AD[k].j * N;
        int val_A = AD[k].v;
        int end_i = offset + N_local;
        
        for (i = offset; i < end_i; i++) {
            C1_local[row_A + (i - offset)] += val_A * B[row_B + i];
        }
        ops_phase1 += (end_i - offset);
    }
    ops_local += ops_phase1;

    // 2. Matriu Dispersa per Matriu Dispersa (MD x MD -> M)
    bzero(VBcol, sizeof(int) * N);
    long long ops_phase2 = 0;

    for(i = offset; i < offset + N_local; i++) {
        int col_idx = i - offset;
        
        for (k = jBD[i]; k < jBD[i+1]; k++)
            VBcol[BD[k].i] = BD[k].v;
        
        for (k = 0; k < ND; k++) {
            int target_j = AD[k].j;
            int val_B = VBcol[target_j]; 
            if (val_B != 0) {
                C2_local[AD[k].i * N_local + col_idx] += AD[k].v * val_B;
                ops_phase2++;
            }
        }
        
        for (k = jBD[i]; k < jBD[i+1]; k++) {
            VBcol[BD[k].i] = 0;
        }
    }
    ops_local += ops_phase2;
                
    // 3. Matriu Dispersa per Matriu Dispersa amb compressió (MD x MD -> MD)
    neleC_local = 0;
    bzero(VBcol, sizeof(int) * N);
    bzero(VCcol, sizeof(int) * N);

    for(i = offset; i < offset + N_local; i++) {
        for (k = jBD[i]; k < jBD[i+1]; k++)
            VBcol[BD[k].i] = BD[k].v;
        
        for (k = 0; k < ND; k++) {
            int target_j = AD[k].j;
            int val_B = VBcol[target_j];
            if (val_B != 0) {
                VCcol[AD[k].i] += AD[k].v * val_B;
            }
        }
        
        for (k = jBD[i]; k < jBD[i+1]; k++) {
            VBcol[BD[k].i] = 0;
        }
        
        for (j = 0; j < N; j++) {
            int val_C = VCcol[j];
            if (val_C != 0) {
                CD[neleC_local].i = j;
                CD[neleC_local].j = i;
                CD[neleC_local].v = val_C;
                VCcol[j] = 0; 
                neleC_local++;
            }
        }
    }

    // Comprovacions d'integritat de les matrius resultants
    for (i = offset; i < offset + N_local; i++) {
        int col_idx = i - offset;
        for (j = 0; j < N; j++) {
            if (C2_local[j * N_local + col_idx] != C1_local[j * N_local + col_idx]) {
                printf("Diferencies C1 i C2 pos %d,%d: %d != %d (Rank %d)\n", j, i, C1_local[j * N_local + col_idx], C2_local[j * N_local + col_idx], rank);
            }
        }
    }

    Suma_local = 0;
    for(k = 0; k < neleC_local; k++) {
        Suma_local += CD[k].v;
        
        int col_idx = CD[k].j - offset;
        if (col_idx >= 0 && col_idx < N_local) {
            if (CD[k].v != C1_local[CD[k].i * N_local + col_idx]) {
                printf("Diferencies C1 i CD a i:%d,j:%d,v%d, k:%d, vd:%d (Rank %d)\n", CD[k].i, CD[k].j, C1_local[CD[k].i * N_local + col_idx], k, CD[k].v, rank);
            }
        }
    }

    // Reducció MPI de mètriques d'execució
    long long stats_local[3] __attribute__((aligned(64))) = {Suma_local, (long long)neleC_local, ops_local};
    long long stats_global[3] __attribute__((aligned(64)));

    MPI_Reduce(stats_local, stats_global, 3, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        Suma_global = stats_global[0];
        neleC_global = (int)stats_global[1];
        ops_global = stats_global[2];

        printf("\nNumero elements de la matriu dispersa C %d\n", neleC_global);   
        printf("Suma dels elements de C %lld\n", Suma_global);
        printf("Total operacions multiplicacio %lld\n", ops_global);
        fflush(stdout);
    }

    // Alliberament de memòria
    free(A); free(B); free(C_local); free(C1_local); free(C2_local);
    free(jBD); free(VCcol); free(VBcol);
    free(AD); free(BD); free(CD);

    MPI_Finalize();
    return 0;
}