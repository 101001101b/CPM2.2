#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <strings.h>
#include <assert.h>

#define N 16000L
#define ND N*N/100

typedef struct {
    int i,j,v;
} tmd;

int *A[N],*B[N],*C[N],*C1[N],*C2[N];
int jBD[N+1],VCcol[N],VBcol[N];
tmd AD[ND],BD[ND],*CD;

long long Suma;

void reserva_matriu(int nn,int *pi[nn])
{
	pi[0]=malloc(nn*nn*sizeof(int));
    assert(pi[0]);
    for (int i=1;i<nn;i++) pi[i]=pi[i-1]+nn;
}

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

int main(int na, char* arg[] )
{
	int nn,nnd;
    int i,j,k,neleC;

	assert(na == 2);
	nn = atoi(arg[1]);
	assert((nn <= N) && (nn >= 500));
	nnd = nn * nn /100;

	reserva_matriu(nn,A);
	reserva_matriu(nn,B);
	reserva_matriu(nn,C);
	reserva_matriu(nn,C1);
	reserva_matriu(nn,C2);
	CD=malloc(nn*nn*sizeof(tmd));
	assert(CD);
    
    bzero(C[0],sizeof(int)*(nn*nn));
    bzero(C1[0],sizeof(int)*(nn*nn));
    bzero(C2[0],sizeof(int)*(nn*nn));
     
    for(k=0;k<nnd;k++)
    {
        AD[k].i=rand()%(nn-1);
        AD[k].j=rand()%(nn-1);
        AD[k].v=rand()%100+1;
        while (A[AD[k].i][AD[k].j]) {
            if(AD[k].i < AD[k].j)
                AD[k].i = (AD[k].i + 1)%nn;
            else 
                AD[k].j = (AD[k].j + 1)%nn;
        }
        A[AD[k].i][AD[k].j] = AD[k].v;
    }
    qsort(AD,nnd,sizeof(tmd),cmp_fil); // ordenat per files

    for(k=0;k<nnd;k++)
    {
        BD[k].i=rand()%(nn-1);
        BD[k].j=rand()%(nn-1);
        BD[k].v=rand()%100+1;
        while (B[BD[k].i][BD[k].j]) {
            if(BD[k].i < BD[k].j)
                BD[k].i = (BD[k].i + 1)%nn;
            else 
                BD[k].j = (BD[k].j + 1)%nn;
        }
        B[BD[k].i][BD[k].j] = BD[k].v;
    }

    qsort(BD,nnd,sizeof(tmd),cmp_col); // ordenat per columnes
    
    // calcul dels index de les columnes
    k=0;
    for (j=0; j<nn+1; j++)
     {
      while (k < nnd && j>BD[k].j) k++;
      jBD[j] = k;
     }

    ////Matriu x matriu original (recorregut de C per columnes)
    //for (i=0;i<nn;i++)
    //    for (j=0;j<nn;j++)
    //        for (k=0;k<nn;k++)
    //            C[j][i] += A[j][k] * B[k][i];
 
    //Matriu dispersa per matriu
    for(i=0;i<nn;i++)
        for (k=0;k<nnd;k++)
            C1[AD[k].i][i] += AD[k].v * B[AD[k].j][i];
            
    //Matriu dispersa per matriu dispersa
    for (j=0;j<nn;j++)
        VBcol[j] = 0;

    for(i=0;i<nn;i++)
      {
        // expandir Columna de B[*][i]
        for (k=jBD[i];k<jBD[i+1];k++)
                VBcol[BD[k].i] = BD[k].v;
        // Calcul de tota una columna de C
        for (k=0;k<nnd;k++)
            C2[AD[k].i][i] += AD[k].v * VBcol[AD[k].j];
        // neteja vector de B[*][i]
        for (j=0;j<nn;j++)
            VBcol[j] = 0;
      }
                
    //Matriu dispersa per matriu dispersa -> dona matriu Dispersa
    neleC=0;
    for (j=0;j<nn;j++)
        VBcol[j] = VCcol[j] = 0;

    for(i=0;i<nn;i++)
      {
        // expandir Columna de B[*][i]
        for (k=jBD[i];k<jBD[i+1];k++)
                VBcol[BD[k].i] = BD[k].v;
        // Calcul de tota una columna de C
        for (k=0;k<nnd;k++)
            VCcol[AD[k].i] += AD[k].v * VBcol[AD[k].j];
        for (j=0;j<nn;j++)
         {
            // neteja vector de B[*][i]
            VBcol[j] = 0;
            // Compressio de C
            if (VCcol[j])
             {
                CD[neleC].i = j;
                CD[neleC].j = i;
                CD[neleC].v = VCcol[j];
                VCcol[j] = 0;
                neleC++;
             }
         }
      }

    // Comprovacio MD x M -> M i MD x MD -> M
    for (i=0;i<nn;i++)
        for(j=0;j<nn;j++)
            if (C2[i][j] != C1[i][j])
                printf("Diferencies C1 i C2 pos %d,%d: %d != %d\n",i,j,C1[i][j],C2[i][j]);

    // Comprovacio MD X MD -> M i MD x MD -> MD
    Suma = 0;
    for(k=0;k<neleC;k++)
     {
        Suma += CD[k].v;
        if (CD[k].v != C1[CD[k].i][CD[k].j])
            printf("Diferencies C1 i CD a i:%d,j:%d,v%d, k:%d, vd:%d\n",CD[k].i,CD[k].j,C1[CD[k].i][CD[k].j],k,CD[k].v);
     }
     
    printf ("\nNumero elements de la matriu dispersa C %d\n",neleC);   
    printf("Suma dels elements de C %lld \n",Suma);
    exit(0);
}
