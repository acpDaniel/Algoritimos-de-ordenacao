#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

int findMax(int *v, int n);
void swap(int *a, int *b);
void reversev(int *v, int size);
void generateRandomv(int *v, int size, int sorted, int unique, int ascending);
void bubbleSort(int *v, int n);
void insertionSort(int *v, int n);
void selectionSort(int *v, int n);
void mergeSort(int *v, int l, int r);
void merge(int *v, int indiceEsq, int indiceM, int indiceDir);
void quickSort(int *v, int l, int r);
void partition(int *v, int l, int r, int *i, int *j);
void shellSort(int *v, int n);
void countingSort(int *v, int n);
void bucketSort(int *v, int n);
void radixsort(int *v, int n);
unsigned bits(unsigned x, int k, int j);
int calculateMaxBitPosition(int max);
void radixexchange(int a[], int l, int r, int b);
int compare(const void *a, const void *b);

int compare(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

int findMax(int *v, int n)
{
    int max = v[0];
    for (int i = 1; i < n; i++)
    {
        if (v[i] > max)
        {
            max = v[i];
        }
    }
    return max;
}

void swap(int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

void reversev(int *v, int size)
{
    int i = 0;
    int j = size - 1;
    while (i < j)
    {
        swap(&v[i], &v[j]);
        i++;
        j--;
    }
}

void generateRandomv(int *v, int size, int sorted, int unique, int ascending)
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
                    if (v[j] == num)
                        break;
                }
            } while (j != i);
        }
        else
        {
            num = rand() % (size * 10);
        }
        v[i] = num;
    }

    if (sorted)
    {
        selectionSort(v, size);
        if (!ascending)
        {
            reversev(v, size);
        }
    }
}

void bubbleSort(int *v, int n)
{
    int i = 0;
    int j = 0;

    for (i = 0; i < n - 1; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (v[j] > v[j + 1])
            {
                swap(&v[j], &v[j + 1]);
            }
        }
    }
}

void insertionSort(int *v, int n)
{
    int i = 0;
    int j = 0;
    int aux;
    for (i = 1; i < n; i++)
    {
        aux = v[i];
        j = i - 1;

        while (j >= 0 && v[j] > aux)
        {
            v[j + 1] = v[j];
            j--;
        }
        v[j + 1] = aux;
    }
}

void selectionSort(int *v, int n)
{
    int i = 0;
    int j = 0;
    int min = 0;

    for (i = 0; i < n - 1; i++)
    {
        min = i;
        for (j = i + 1; j < n; j++)
        {
            if (v[j] < v[min])
            {
                min = j;
            }
        }

        // trocar o menor elemento encontrado com o primeiro elemento do subv nao ordenado
        if (min != i)
        {
            swap(&v[min], &v[i]);
        }
    }
}

void merge(int arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    int L[n1], R[n2];

    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    i = 0;

    j = 0;

    k = l;
    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void mergeSort(int *v, int l, int r)
{
    if (l < r)
    {
        int m = l + (r - l) / 2;
        mergeSort(v, l, m);
        mergeSort(v, m + 1, r);
        merge(v, l, m, r);
    }
}

void partition(int *v, int l, int r, int *i, int *j)
{
    *i = l;
    *j = r;
    int pivo = v[(*i + *j) / 2];

    do
    {

        while (*i <= r && v[*i] < pivo)
        {
            (*i)++;
        }

        while (*j >= l && v[*j] > pivo)
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

void quickSort(int *v, int l, int r)
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

void shellSort(int *v, int n)
{
    for (int h = n / 2; h > 0; h /= 2)
    {
        for (int i = h; i < n; i++)
        {
            int aux = v[i];
            int j = i;

            for (j = i; j >= h && v[j - h] > aux; j = j - h)
            {
                v[j] = v[j - h];
            }

            v[j] = aux;
        }
    }
}

// verificar
void countingSort(int *v, int n)
{
    int max = v[0];
    for (int i = 1; i < n; i++)
    {
        if (v[i] > max)
        {
            max = v[i];
        }
    }

    // Alocar espaço para max + 1 elementos
    int *counts = (int *)calloc(max + 1, sizeof(int));

    for (int i = 0; i < n; i++)
    {
        counts[v[i]]++;
    }

    int index = 0;
    for (int j = 0; j <= max; j++)
    {
        while (counts[j] > 0)
        {
            v[index++] = j;
            counts[j]--;
        }
    }

    free(counts);
}

void bucketSort(int *v, int size)
{
    int numBuckets = 10;
    int i, j;

    // Encontra o valor máximo no array para determinar a escala dos buckets
    int maxValue = v[0];
    for (i = 1; i < size; i++)
    {
        if (v[i] > maxValue)
        {
            maxValue = v[i];
        }
    }

    // Cria os buckets
    int **buckets = (int **)malloc(numBuckets * sizeof(int *));
    int *bucketSizes = (int *)malloc(numBuckets * sizeof(int));
    int *bucketCounts = (int *)malloc(numBuckets * sizeof(int));

    for (i = 0; i < numBuckets; i++)
    {
        bucketSizes[i] = size / numBuckets;
        buckets[i] = (int *)malloc(bucketSizes[i] * sizeof(int));
        bucketCounts[i] = 0;
    }

    // Distribui os elementos do array nos buckets
    for (i = 0; i < size; i++)
    {
        int bucketIndex = (v[i] * numBuckets) / (maxValue + 1);
        if (bucketCounts[bucketIndex] == bucketSizes[bucketIndex])
        {
            // Redimensiona o bucket se necessário
            bucketSizes[bucketIndex] *= 2;
            buckets[bucketIndex] = (int *)realloc(buckets[bucketIndex], bucketSizes[bucketIndex] * sizeof(int));
        }
        buckets[bucketIndex][bucketCounts[bucketIndex]++] = v[i];
    }

    // Ordena os elementos em cada bucket usando quicksort
    for (i = 0; i < numBuckets; i++)
    {
        if (bucketCounts[i] > 0)
        {
            qsort(buckets[i], bucketCounts[i], sizeof(int), compare);
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

void radixsort(int *v, int n)
{
    int i;
    int *b;
    int exp = 1;

    b = (int *)calloc(n, sizeof(int));

    int max = v[0];
    for (int i = 1; i < n; i++)
    {
        if (v[i] > max)
        {
            max = v[i];
        }
    }

    while (max / exp > 0)
    {
        int bucket[10] = {0};
        for (i = 0; i < n; i++)
            bucket[(v[i] / exp) % 10]++;
        for (i = 1; i < 10; i++)
            bucket[i] += bucket[i - 1];
        for (i = n - 1; i >= 0; i--)
            b[--bucket[(v[i] / exp) % 10]] = v[i];
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

void radixexchange(int *v, int l, int r, int b)
{
    if (r <= l || b < 0)
    {
        return;
    }

    int i = l, j = r;
    while (i < j)
    {
        while (i < j && bits(v[i], b, 1) == 0)
            i++;
        while (i < j && bits(v[j], b, 1) == 1)
            j--;
        if (i < j)
        {
            int t = v[i];
            v[i] = v[j];
            v[j] = t;
        }
    }

    if (bits(v[r], b, 1) == 0)
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

    int *v = (int *)malloc(n * sizeof(int));

    // Geração do vetor
    generateRandomv(v, n, sorted, unique, ascending);

    // printf("Vetor desordenado:\n");
    // for (int i = 0; i < n; i++)
    // {
    //     printf("%d \n", v[i]);
    // }
    // printf("\n");

    struct timeval start, end;
    long seconds, useconds;
    double cpu_time_used;

    // int max = max = findMax(v, n);
    // int maxBitPosition = maxBitPosition = calculateMaxBitPosition(max);

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
        // radixexchange(v, 0, n - 1, maxBitPosition - 1);
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
    //     printf("%d \n", v[i]);
    // }
    // printf("\n");

    free(v);
    return 0;
}