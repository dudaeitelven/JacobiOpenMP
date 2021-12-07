#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#pragma pack(1)

#define ITERACOES 1000
#define ERRO 1e-10
#define N 5
#define N_MAX 10

/*
Para compilar:
1 - Abrir o local do fonte
2 - Digitar para compilar: gcc -o <programa> <programa.c> -fopenmp -lm
3 - Digitar para rodar: ./<programa> <numero_threads>
*/

void gerarMatrizAleatoria(float *matrizInicial) {
	int n;
	int i, j;
	int posMatriz;

	for(i=0; i<N; i++) {
		for(j=0; j<N; j++) {
			n = rand() % N_MAX + 1;

			posMatriz = (i * N) + j;
			matrizInicial[posMatriz] = n;
		}
	}
}

void gerarVetorAleatorio(float *vetor) {
	int n;
	int i;
	
	for(i=0; i<N; i++) {
		n = rand() % N_MAX + 1;
		vetor[i] = n;
	}
}

void gerarVetorZerado(float *vetor) {
	int i;
	
	for(i=0; i<N; i++) {
		vetor[i] = 0;
	}
}

void mostrarMatriz(float *matriz) {
	int i, j;
	int posMatriz;

	for(i=0; i<N; i++) {
		for(j=0; j<N; j++) {
			posMatriz = (i * N) + j;
			printf("%.2f ", matriz[posMatriz]);
		}
		printf("\n");
	}
	printf("\n");
}

void mostrarVetor(float *vetor) {
	int i;
	
	for(i=0; i<N; i++) {
		printf("%.2f ", vetor[i]);
	}
	printf("\n\n");
}


void calcularJacobi(int ini, int contador, int nThreads, float *matrizInicial, float *vetorInicial, float *vetorCalculado, float *vetorCalculadoAnt) {
	int 	posMatriz;
	int 	i, j;
	float 	soma = 0;
	float 	dp = 0;

	contador = 0;
	
	while (contador < ITERACOES) {
		for(i = ini; i < N; i+=nThreads) {
			soma = 0;

			for(j=0; j<N; j++) {
				posMatriz = (i * N) + j;

				if (i != j) {
					soma += matrizInicial[posMatriz] * vetorCalculado[i];	        		
				}
				else {
					dp =  matrizInicial[posMatriz];
				}
			}
			vetorCalculado[i] = (vetorInicial[i] - soma) / dp;
		}	
		
		// if ((calcularNorma(vetorCalculadoAnt) - calcularNorma(vetorCalculado)) < ERRO) {
		// 	contador = ITERACOES;
		// }
		// else {
		// 	vetorCalculadoAnt = vetorCalculado;
		// }
			
		contador++;
	}
}


float calcularNorma(float *vetor) {
	int i;
	double soma = 0;
	
	for(i=0; i<N; i++) {
		soma += vetor[i] * vetor[i];
	}
	
	soma = sqrt(soma);
	
	return soma;
}	

int main(int argc, char **argv ){
	int 	nThreads;
	int		contador;
	float 	*matrizInicial;
	float 	*vetorInicial;
	float 	*vetorCalculado;
	float 	*vetorCalculadoAnt;
	
	if ( argc != 2){
		printf("%s <num_threads>\n", argv[0]);
		exit(0);
	}

	nThreads = atoi(argv[1]);

	//Alocar memoria dinamicamente
	matrizInicial 		= (float *) malloc((N * N) * sizeof(float));
	vetorInicial 		= (float *) malloc(N * sizeof(float));
	vetorCalculado 		= (float *) malloc(N * sizeof(float));
	vetorCalculadoAnt 	= (float *) malloc(N * sizeof(float));

	//Iniciar as variaveis
	contador = 0;

	gerarMatrizAleatoria(matrizInicial);
	gerarVetorAleatorio(vetorInicial);
	gerarVetorZerado(vetorCalculado);
	gerarVetorZerado(vetorCalculadoAnt);

	//Mostrar valores iniciais
	printf("Matriz:\n");
	mostrarMatriz(matrizInicial);
	
	printf("Vetor inicial:\n");
	mostrarVetor(vetorInicial);

	//Calcular Jacobi
	int ini = 0; //Remover
	calcularJacobi(ini, contador, nThreads, matrizInicial,vetorInicial,vetorCalculado,vetorCalculadoAnt);

	//Mostrar valores calculados
	printf("Vetor calculado:\n");
	mostrarVetor(vetorCalculado);

	//Limpar memoria
	free(matrizInicial);
	free(vetorInicial);
	free(vetorCalculado);
	free(vetorCalculadoAnt);
}