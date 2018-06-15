#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define ABS(x) (((x) > (0)) ? (x) : -(x))

typedef struct organismo
{
    int *dna;
    unsigned long int avaliacao;
} Organismo;

void AlgoritmoGenetico(int tam_dna, int qtd_ind, int geracoes, int prob_mutacao, int individuos_salvos);
int avaliacao (Organismo org[], int qtd_ind, int tam_dna, int individuos_salvos);
unsigned long int definicao_problema (int *dna, int tam_dna);
Organismo *orgmalloc (int tam, int tam_dna);
void orgfree (Organismo org[], int tam);
void gera_populacao (Organismo org[], int qtd_ind, int tam_dna);
void mutacao (Organismo org[], int individuo, int tam_dna);
void cruzamento (Organismo Origem[], Organismo Destino[], int pai, int mae, int filho, int tam_dna);
void Merge (Organismo vet[], int inicio, int meio, int fim);
void MergeSort (Organismo vet[], int inicio, int fim);
void ini_ExponentialRanking (double c, int *S_selection, int qtd_ind);
void ini_LinearRanking(double tx_rep_pior, int *S_selection, int qtd_ind);
int RankingSelection  (int *S_selection, int qtd_ind);
int buscabinaria(int x, int n, int v[]);

int BIT[31] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824};

int main (void)
{
    int tam_dna = 1000;
    int qtd_ind = 100;
    int geracoes = 1000;
    int prob_mutacao = 60;
    int individuos_salvos = 20;

    AlgoritmoGenetico(tam_dna, qtd_ind, geracoes, prob_mutacao, individuos_salvos);

    return 0;
}

void AlgoritmoGenetico(int tam_dna, int qtd_ind, int geracoes, int prob_mutacao, int individuos_salvos)
{
    srand((unsigned int) time(NULL));
    int *S_selection = (int *) malloc (qtd_ind * sizeof(int));
    //ini_LinearRanking(0.1, S_selection, qtd_ind);
    ini_ExponentialRanking(0.95, S_selection, qtd_ind);
    int geracao, pai, mae, filho, termino, idna;

    Organismo *org = orgmalloc(qtd_ind, tam_dna);
    Organismo *org2 = orgmalloc(qtd_ind, tam_dna);

    gera_populacao(org, qtd_ind, tam_dna);
    termino = avaliacao(org, qtd_ind, tam_dna, 0); // Analise inicial de todos (cuidado com o zero)
    MergeSort(org, 0, qtd_ind - 1);

    for (geracao = 0; geracao < geracoes && termino == 0; geracao++)
    {
        for (filho = 0; filho < qtd_ind; filho++)
            for (idna = 0; idna < tam_dna; idna++)
                org2[filho].dna[idna] = org[filho].dna[idna];

        for (filho = 0; filho < qtd_ind - individuos_salvos; filho++)
        {
            pai = RankingSelection(S_selection, qtd_ind);
            mae = RankingSelection(S_selection, qtd_ind);
            cruzamento(org2, org, pai, mae, filho, tam_dna);
            if (rand() % 100 < prob_mutacao)
                mutacao(org, filho, tam_dna);
        }
        termino = avaliacao(org, qtd_ind, tam_dna, individuos_salvos);
        MergeSort(org, 0, qtd_ind - 1);
        printf("G: %d\tMelhor: %d\n", geracao, org[qtd_ind - 1].avaliacao);
    }

    orgfree(org, qtd_ind);
    orgfree(org2, qtd_ind);
    free(S_selection);
}

int avaliacao (Organismo org[], int qtd_ind, int tam_dna, int individuos_salvos)
{
    int i, j;
    int retorno = 0;

    for (i = 0; i < (qtd_ind - individuos_salvos) && !retorno; i++)
    {
        org[i].avaliacao = definicao_problema(org[i].dna, tam_dna);
        // para o problema do one-max fitnes maximo = tam_dna (criterio de parada)
        if (org[i].avaliacao == tam_dna)
            retorno = 1;
    }

    return retorno;
}

unsigned long int definicao_problema (int *dna, int tam_dna)
{
    int j;
    unsigned long int total = 0;

    // problema do one-max (maximizar quantidades de uns no dna)
    for (j = 0; j < tam_dna; j++)
    {
        total += (unsigned long int) dna[j];
    }

    return total;
}

Organismo *orgmalloc (int tam, int tam_dna)
{
    Organismo *org = (Organismo *) malloc(tam * sizeof(Organismo));
    int i;
    for (i = 0; i < tam; i++)
        org[i].dna = (int *) malloc(tam_dna * sizeof(int));
    return org;
}

void orgfree (Organismo org[], int tam)
{
    int i;
    for (i = 0; i < tam; i++)
        free(org[i].dna);
    free(org);
}

void gera_populacao (Organismo org[], int qtd_ind, int tam_dna)
{
    int i, j;

    for (i = 0; i < qtd_ind; i++)
    {
        for (j = 0; j < tam_dna; j++)
        {
            org[i].dna[j] = rand() % 2;
        }
    }
}

void mutacao (Organismo org[], int individuo, int tam_dna)
{
    int gene = rand() % tam_dna;
    org[individuo].dna[gene] = !org[individuo].dna[gene];
}

void cruzamento (Organismo Origem[], Organismo Destino[], int pai, int mae, int filho, int tam_dna)
{
    int i;

    for (i = 0; i < tam_dna; i++)
    {
        if(Origem[mae].dna[i] == Origem[pai].dna[i])
            Destino[filho].dna[i] = Origem[mae].dna[i];
        else
            Destino[filho].dna[i] = rand() % 2;
    }
}

void Merge (Organismo vet[], int inicio, int meio, int fim)
{
    Organismo *w = (Organismo*) malloc((fim - inicio + 1) * sizeof(Organismo));
    int i = inicio, j = meio + 1, k = 0;
    while(i <= meio && j <= fim)
    {
        if(vet[i].avaliacao <= vet[j].avaliacao) w[k++] = vet[i++];
        else w[k++] = vet[j++];
    }
    while(i <= meio) w[k++] = vet[i++];
    while(j <= fim) w[k++] = vet[j++];
    for(i = 0; i < k; i++) vet[i + inicio] = w[i];
    free(w);
}

void MergeSort (Organismo vet[], int inicio, int fim)
{
    if(inicio < fim)
    {
        int meio = (inicio + fim) / 2;
        MergeSort(vet, inicio, meio);
        MergeSort(vet, meio + 1, fim);
        Merge(vet, inicio, meio, fim);
    }
}

void ini_ExponentialRanking (double c, int *S_selection, int qtd_ind)
{
    int i;
    double k = (c - 1.0) / (pow(c, qtd_ind) - 1.0);
    double *S_tmp = (double *) malloc (qtd_ind * sizeof(double));

    i = 1;
    S_tmp[i - 1] = k * pow(c, qtd_ind - i);

    for (i = 2; i <= qtd_ind; i++)
        S_tmp[i - 1] = S_tmp[i - 2] + k * pow(c, qtd_ind - i);

    double diff1 = (double) (S_tmp[1] - S_tmp[0]) / (double) 2.0;
    double diff2 = (double) (S_tmp[qtd_ind - 1] - S_tmp[qtd_ind - 2]) / (double) 2.0;
    double minimo = S_tmp[0];
    double maximo = S_tmp[qtd_ind - 1];
    double a = (double)(((double) 1.0 - (double) diff1 - (double) diff2) / ((double) maximo - (double) minimo));
    double b = ((double) diff1 - (double) a * (double) minimo);

    a = (double) a * (double) 32767.0;
    b = (double) b * (double) 32767.0;

    for (i = 0; i < qtd_ind; i++)
        S_selection[i] = round((double) a * (double) S_tmp[i] + (double) b);

    free(S_tmp);
}

void ini_LinearRanking(double tx_rep_pior, int *S_selection, int qtd_ind)
{
    int i;
    double tx_rep_melhor = 2 - tx_rep_pior;
    double *S_tmp = (double *) malloc (qtd_ind * sizeof(double));

    i = 1;
    S_tmp[i - 1] = ((double) 1.0 / (double) qtd_ind) * (tx_rep_pior + (tx_rep_melhor - tx_rep_pior) *  (double) (i - 1) / (double)(qtd_ind - 1));

    for (i = 2; i <= qtd_ind; i++)
        S_tmp[i - 1] = S_tmp[i - 2] + ((double) 1.0 / (double) qtd_ind) * (tx_rep_pior + (tx_rep_melhor - tx_rep_pior) *  (double) (i - 1) / (double)(qtd_ind - 1));

    double diff1 = (double) (S_tmp[1] - S_tmp[0]) / (double) 2.0;
    double diff2 = (double) (S_tmp[qtd_ind - 1] - S_tmp[qtd_ind - 2]) / (double) 2.0;
    double minimo = S_tmp[0];
    double maximo = S_tmp[qtd_ind - 1];
    double a = (double)(((double) 1.0 - (double) diff1 - (double) diff2) / ((double) maximo - (double) minimo));
    double b = ((double) diff1 - (double) a * (double) minimo);

    a = (double) a * (double) 32767.0;
    b = (double) b * (double) 32767.0;

    for (i = 0; i < qtd_ind; i++)
        S_selection[i] = round((double) a * (double) S_tmp[i] + (double) b);

    free(S_tmp);
}

int RankingSelection  (int *S_selection, int qtd_ind)
{
    return buscabinaria(rand(), qtd_ind, S_selection);
}

int buscabinaria(int x, int n, int v[])
{
    int inicio, meio, fim;
    inicio = 0;
    fim = n - 1;
    while (fim - inicio != 1)
    {
        meio = (inicio + fim)/2;
        if (v[meio] < x)
            inicio = meio;
        else
            fim = meio;
    }

    if (abs(v[inicio] - x) < abs(v[fim] - x))
        meio = inicio;
    else
        meio = fim;

    return meio;
}
