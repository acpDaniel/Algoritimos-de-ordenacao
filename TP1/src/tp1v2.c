#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

typedef struct
{
    int chave;
    char valor[99];
} Item;

int findMax(Item *v, int n);
void swap(Item *a, Item *b);
void reversev(Item *v, int size);
void generateRandomv(Item *v, int size, int sorted, int unique, int ascending);
void bubbleSort(Item *v, int n);
void insertionSort(Item *v, int n);
void selectionSort(Item *v, int n);
void mergeSort(Item *v, int l, int r);
void merge(Item *v, int l, int m, int r);
void quickSort(Item *v, int l, int r);
void partition(Item *v, int l, int r, int *i, int *j);
void shellSort(Item *v, int n);
void countingSort(Item *v, int n);
void bucketSort(Item *v, int size);
void radixsort(Item *v, int n);
unsigned bits(unsigned x, int k, int j);
int calculateMaxBitPosition(int max);
void radixexchange(Item *v, int l, int r, int b);
int compare(const void *a, const void *b);

int compare(const void *a, const void *b)
{
    Item *itemA = (Item *)a;
    Item *itemB = (Item *)b;
    return (itemA->chave - itemB->chave);
}

int findMax(Item *v, int n)
{
    int max = v[0].chave;
    for (int i = 1; i < n; i++)
    {
        if (v[i].chave > max)
        {
            max = v[i].chave;
        }
    }
    return max;
}

void swap(Item *a, Item *b)
{
    Item temp = *a;
    *a = *b;
    *b = temp;
}

void reversev(Item *v, int size)
{
    int i;
    Item temp;
    for (i = 0; i < size / 2; i++)
    {
        temp = v[i];
        v[i] = v[size - i - 1];
        v[size - i - 1] = temp;
    }
}

void generateRandomv(Item *v, int size, int sorted, int unique, int ascending)
{
    int i, j, num;
    srand(37);

    for (i = 0; i < size; i++)
    {
        if (unique)
        {
            do
            {
                num = rand() % (size * 10);
                for (j = 0; j < i; j++)
                {
                    if (v[j].chave == num)
                        break;
                }
            } while (j != i);
        }
        else
        {
            num = rand() % (size * 10);
        }
        v[i].chave = num;
        memset(v[i].valor, 'a', sizeof(v[i].valor) - 1);
        v[i].valor[sizeof(v[i].valor) - 1] = '\0';
    }

    if (sorted)
    {
        // Implementar a ordenação se necessário
    }
}

void bubbleSort(Item *v, int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (v[j].chave > v[j + 1].chave)
            {
                swap(&v[j], &v[j + 1]);
            }
        }
    }
}

void insertionSort(Item *v, int n)
{
    int i, j;
    Item key;
    for (i = 1; i < n; i++)
    {
        key = v[i];
        j = i - 1;

        while (j >= 0 && v[j].chave > key.chave)
        {
            v[j + 1] = v[j];
            j = j - 1;
        }
        v[j + 1] = key;
    }
}

void selectionSort(Item *v, int n)
{
    int i, j, min;
    for (i = 0; i < n - 1; i++)
    {
        min = i;
        for (j = i + 1; j < n; j++)
        {
            if (v[j].chave < v[min].chave)
            {
                min = j;
            }
        }

        // Trocar o menor elemento encontrado com o primeiro elemento do subvetor não ordenado
        if (min != i)
        {
            swap(&v[i], &v[min]);
        }
    }
}

void merge(Item *v, int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    Item *L = (Item *)malloc(n1 * sizeof(Item));
    Item *R = (Item *)malloc(n2 * sizeof(Item));

    for (i = 0; i < n1; i++)
        L[i] = v[l + i];
    for (j = 0; j < n2; j++)
        R[j] = v[m + 1 + j];

    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2)
    {
        if (L[i].chave <= R[j].chave)
        {
            v[k] = L[i];
            i++;
        }
        else
        {
            v[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1)
    {
        v[k] = L[i];
        i++;
        k++;
    }

    while (j < n2)
    {
        v[k] = R[j];
        j++;
        k++;
    }

    free(L);
    free(R);
}

void mergeSort(Item *v, int l, int r)
{
    if (l < r)
    {
        int m = l + (r - l) / 2;
        mergeSort(v, l, m);
        mergeSort(v, m + 1, r);
        merge(v, l, m, r);
    }
}

void partition(Item *v, int l, int r, int *i, int *j)
{
    *i = l;
    *j = r;
    int pivo = v[(*i + *j) / 2].chave;

    do
    {

        while (*i <= r && v[*i].chave < pivo)
        {
            (*i)++;
        }

        while (*j >= l && v[*j].chave > pivo)
        {
            (*j)--;
        }

        if (*i <= *j)
        {
            swap(&v[*i], &v[*j]);
            (*i)++;
            (*j)--;
        }

    } while (*i <= *j);
}

void quickSort(Item *v, int l, int r)
{
    if (l < r)
    {
        int i = 0;
        int j = 0;
        partition(v, l, r, &i, &j);
        if (l < j)
            quickSort(v, l, j);

        if (i < r)
            quickSort(v, i, r);
    }
}

void shellSort(Item *v, int n)
{
    for (int h = n / 2; h > 0; h /= 2)
    {
        for (int i = h; i < n; i++)
        {
            Item aux = v[i];
            int j = i;

            for (j = i; j >= h && v[j - h].chave > aux.chave; j -= h)
            {
                v[j] = v[j - h];
            }

            v[j] = aux;
        }
    }
}

// verificar
void countingSort(Item *v, int n)
{
    int i;
    int max = v[0].chave;
    for (i = 1; i < n; i++)
    {
        if (v[i].chave > max)
        {
            max = v[i].chave;
        }
    }

    // Alocar espaço para max + 1 elementos
    int *counts = (int *)calloc(max + 1, sizeof(int));
    Item *output = (Item *)malloc(n * sizeof(Item));

    for (i = 0; i < n; i++)
    {
        counts[v[i].chave]++;
    }

    for (i = 1; i <= max; i++)
    {
        counts[i] += counts[i - 1];
    }

    for (i = n - 1; i >= 0; i--)
    {
        output[counts[v[i].chave] - 1] = v[i];
        counts[v[i].chave]--;
    }

    for (i = 0; i < n; i++)
    {
        v[i] = output[i];
    }

    free(counts);
    free(output);
}

void bucketSort(Item *v, int size)
{
    int numBuckets = 10;
    int i, j;

    // Encontra o valor máximo no array para determinar a escala dos buckets
    int maxValue = v[0].chave;
    for (i = 1; i < size; i++)
    {
        if (v[i].chave > maxValue)
        {
            maxValue = v[i].chave;
        }
    }

    // Cria os buckets
    Item **buckets = (Item **)malloc(numBuckets * sizeof(Item *));
    int *bucketSizes = (int *)malloc(numBuckets * sizeof(int));
    int *bucketCounts = (int *)malloc(numBuckets * sizeof(int));

    for (i = 0; i < numBuckets; i++)
    {
        bucketSizes[i] = size / numBuckets;
        buckets[i] = (Item *)malloc(bucketSizes[i] * sizeof(Item));
        bucketCounts[i] = 0;
    }

    // Distribui os elementos do array nos buckets
    for (i = 0; i < size; i++)
    {
        int bucketIndex = (v[i].chave * numBuckets) / (maxValue + 1);
        if (bucketCounts[bucketIndex] == bucketSizes[bucketIndex])
        {
            // Redimensiona o bucket se necessário
            bucketSizes[bucketIndex] *= 2;
            buckets[bucketIndex] = (Item *)realloc(buckets[bucketIndex], bucketSizes[bucketIndex] * sizeof(Item));
        }
        buckets[bucketIndex][bucketCounts[bucketIndex]++] = v[i];
    }

    // Ordena os elementos em cada bucket usando qsort
    for (i = 0; i < numBuckets; i++)
    {
        if (bucketCounts[i] > 0)
        {
            qsort(buckets[i], bucketCounts[i], sizeof(Item), compare);
        }
    }

    // Concatena todos os buckets ordenados no array original
    int index = 0;
    for (i = 0; i < numBuckets; i++)
    {
        for (j = 0; j < bucketCounts[i]; j++)
        {
            v[index++] = buckets[i][j];
        }
        free(buckets[i]);
    }

    free(buckets);
    free(bucketSizes);
    free(bucketCounts);
}

void radixsort(Item *v, int n)
{
    int i;
    Item *b;
    int exp = 1;

    b = (Item *)calloc(n, sizeof(Item));

    int max = v[0].chave;
    for (i = 1; i < n; i++)
    {
        if (v[i].chave > max)
        {
            max = v[i].chave;
        }
    }

    while (max / exp > 0)
    {
        int bucket[10] = {0};
        for (i = 0; i < n; i++)
            bucket[(v[i].chave / exp) % 10]++;
        for (i = 1; i < 10; i++)
            bucket[i] += bucket[i - 1];
        for (i = n - 1; i >= 0; i--)
            b[--bucket[(v[i].chave / exp) % 10]] = v[i];
        for (i = 0; i < n; i++)
            v[i] = b[i];
        exp *= 10;
    }

    free(b);
}

unsigned bits(unsigned x, int k, int j)
{
    return (x >> k) & ~(~0 << j);
}

int calculateMaxBitPosition(int max)
{
    return (int)floor(log2(max)) + 1;
}

void radixexchange(Item *v, int l, int r, int b)
{
    if (r <= l || b < 0)
    {
        return;
    }

    int i = l, j = r;
    while (i < j)
    {
        while (i < j && bits(v[i].chave, b, 1) == 0)
            i++;
        while (i < j && bits(v[j].chave, b, 1) == 1)
            j--;
        if (i < j)
        {
            Item t = v[i];
            v[i] = v[j];
            v[j] = t;
        }
    }

    if (bits(v[r].chave, b, 1) == 0)
    {
        j++;
    }

    radixexchange(v, l, j - 1, b - 1);
    radixexchange(v, j, r, b - 1);
}

int main(int argc, char *argv[])
{
    if (argc < 6)
    {
        printf("Uso: %s <tamanho do vetor> <numeros unicos 0|1> <pre-ordenado 0|1> <ordem ascendente 0|1> <algoritmo 1-9>\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    int unique = atoi(argv[2]);
    int sorted = atoi(argv[3]);
    int ascending = atoi(argv[4]);
    int choice = atoi(argv[5]);

    Item *v = (Item *)malloc(n * sizeof(Item));

    // Geração do vetor
    generateRandomv(v, n, sorted, unique, ascending);

    // printf("Vetor desordenado:\n");
    // for (int i = 0; i < n; i++)
    // {
    //     printf("Chave: %d, Valor: %s\n", v[i].chave, v[i].valor);
    // }
    // printf("\n");

    struct timeval start, end;
    long seconds, useconds;
    double cpu_time_used;

    int max = max = findMax(v, n);
    int maxBitPosition = maxBitPosition = calculateMaxBitPosition(max);

    gettimeofday(&start, NULL);

    // Execução do algoritmo escolhido
    switch (choice)
    {
    case 1:
        bubbleSort(v, n);
        break;
    case 2:
        insertionSort(v, n);
        break;
    case 3:
        selectionSort(v, n);
        break;
    case 4:
        mergeSort(v, 0, n - 1);
        break;
    case 5:
        quickSort(v, 0, n - 1);
        break;
    case 6:
        shellSort(v, n);
        break;
    case 7:
        countingSort(v, n);
        break;
    case 8:
        bucketSort(v, n);
        break;
    case 9:
        radixexchange(v, 0, n - 1, maxBitPosition - 1);
        break;
    default:
        printf("Escolha inválida de algoritmo.\n");
        free(v);
        return 1;
    }

    gettimeofday(&end, NULL);

    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    cpu_time_used = seconds + useconds / 1e6;

    printf("%f,", cpu_time_used);

    // printf("Vetor ordenado:\n");
    // for (int i = 0; i < n; i++)
    // {
    //     printf("Chave: %d, Valor: %s\n", v[i].chave, v[i].valor);
    // }
    // printf("\n");

    free(v);
    return 0;
}