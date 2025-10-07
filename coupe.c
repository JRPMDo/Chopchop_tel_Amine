#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct int_list {
    int* line;
    struct int_list* next;
    
}int_list;

typedef struct matrice {
    int* values;
    int size;
    
}matrice;


void printmat(matrice* mat){
    
    for (int i = 0; i < mat->size; i++) {
        for (int j = 0; j < mat->size; j++) {
            printf("%3d ", mat->values[i * mat->size + j]);
        }
        printf("\n");
    }
}

matrice* copymat(matrice* mat){
    matrice* res = malloc(sizeof(matrice));
    res->size = mat->size;
    res->values = malloc(sizeof(int)*(res->size)*(res->size));
    for(int i = 0; i <(mat->size)*(mat->size);i++){
        (res->values)[i] = (mat->values)[i];
    }
}

matrice* multmat(const matrice* a, const matrice* b) {
    int n = a->size;
    matrice* res = malloc(sizeof(matrice));
    res->size = n;
    res->values = calloc(n * n, sizeof(int));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int sum = 0;
            for (int k = 0; k < n; k++) {
                sum += a->values[i * n + k] * b->values[k * n + j];
            }
            res->values[i * n + j] = sum;
        }
    }
    return res;
}
 

void freemat(matrice* m) {
    if (!m) return;
    free(m->values);
    free(m);
}


matrice* identite(int n) {
    matrice* I = malloc(sizeof(matrice));
    I->size = n;
    I->values = calloc(n * n, sizeof(int));
    for (int i = 0; i < n; i++) {
        I->values[i * n + i] = 1;
    }
    return I;
}


matrice* powmat(matrice* mat, int pow) {
    if (pow < 0) {
        fprintf(stderr, "Erreur : puissance négative non gérée.\n");
        return NULL;
    }

    matrice* res = identite(mat->size); 
    matrice* base = malloc(sizeof(matrice));
    base->size = mat->size;
    base->values = malloc(sizeof(int) * mat->size * mat->size);
    for (int i = 0; i < mat->size * mat->size; i++)
        base->values[i] = mat->values[i];

    while (pow > 0) {
        if (pow % 2 == 1) { 
            matrice* tmp = multmat(res, base);
            freemat(res);
            res = tmp;
        }
        pow /= 2;
        if (pow > 0) { 
            matrice* tmp2 = multmat(base, base);
            freemat(base);
            base = tmp2;
        }
    }

    freemat(base);
    return res;
}


int checkcycle(matrice* mat){
    for(int i =0; i< mat->size ; i++){
        if((mat->values)[i+(mat->size)*i]) return i ;
    }
    return -1;
}


matrice* coupe(matrice* mat, int nbr) {
    int n = mat->size;

    matrice* res = malloc(sizeof(matrice));
    res->size = n - 1;
    res->values = malloc(sizeof(int) * (res->size) * (res->size));

    int ii = 0;

    for (int i = 0; i < n; i++) {
        if (i == nbr) continue;

        int jj = 0;
        for (int j = 0; j < n; j++) {
            if (j == nbr) continue; 

            res->values[ii * (n - 1) + jj] = mat->values[i * n + j];

            jj++;
        }

        ii++;
    }

    return res;
}

matrice* matcoupee(matrice* mat) {
    matrice* res = copymat(mat);  

    while (1) {
        matrice* tampon = powmat(res, res->size);
        int idx = checkcycle(tampon);
        freemat(tampon);

        if (idx == -1) break;  

        matrice* old = res;
        res = coupe(res, idx);  
        freemat(old);           
    }

    return res;
}

int is_acyclic(int* subset, int subset_size, matrice* mat) {
    matrice* test = copymat(mat);

    for (int i = subset_size - 1; i >= 0; i--) {
        matrice* old = test;
        test = coupe(test, subset[i]);
        freemat(old);
    }
    matrice* powmat_test = powmat(test, test->size);
    int c = checkcycle(powmat_test);
    freemat(powmat_test);
    freemat(test);
    return c == -1;
}


void combinations(int n, int k, int start, int* combo, int depth, matrice* mat, int* best_subset, int* best_size) {
    if (depth == k) {
        if (is_acyclic(combo, k, mat)) {
            if (k < *best_size) {
                *best_size = k;
                for (int i = 0; i < k; i++)
                    best_subset[i] = combo[i];
            }
        }
        return;
    }

    for (int i = start; i < n; i++) {
        combo[depth] = i;
        combinations(n, k, i + 1, combo, depth + 1, mat, best_subset, best_size);
    }
}


void find_min_cut(matrice* mat) {
    int n = mat->size;
    int best_size = n + 1;
    int* best_subset = malloc(sizeof(int) * n);
    int* combo = malloc(sizeof(int) * n);

    for (int k = 0; k <= n; k++) { 
        combinations(n, k, 0, combo, 0, mat, best_subset, &best_size);
        if (best_size != n + 1) break; 
    }

    printf("Nombre minimal de sommets à supprimer : %d\n", best_size);
    printf("Sommets à supprimer : ");
    for (int i = 0; i < best_size; i++)
        printf("%d ", best_subset[i]);
    printf("\n");

    free(best_subset);
    free(combo);
}
