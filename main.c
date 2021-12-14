#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#pragma pack(1)

#define ITERACOES 1000
#define ERRO 1.0e-10
#define N 4

/*
Para compilar:
1 - Abrir o local do fonte
2 - Digitar para compilar: gcc -o <programa> <programa.c> -fopenmp -lm
3 - Digitar para rodar: ./<programa> <numero_threads>
*/

void gerarMatriz(float *matrizInicial) {
	matrizInicial[0] = 10;
	matrizInicial[1] = -1;
	matrizInicial[2] =  2;
	matrizInicial[3] =  0;

	matrizInicial[4] = -1;
	matrizInicial[5] = 11;
	matrizInicial[6] = -1;
	matrizInicial[7] =  3;

	matrizInicial[8]  =  2;
	matrizInicial[9]  = -1;
	matrizInicial[10] = 10;
	matrizInicial[11] = -1;

	matrizInicial[12] =  0;
	matrizInicial[13] =  3;
	matrizInicial[14] = -1;
	matrizInicial[15] =  8;
}

void gerarVetor(float *vetor) {
	vetor[0] =   6;
	vetor[1] =  25;
	vetor[2] = -11;
	vetor[3] =  15;
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

void copiarVetor(float *copiado, float *original)
{
    int i;

    for(i=0 ; i < N ; i++)
        copiado[i] = original[i];
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

void calcularJacobi(int contador, int nThreads, float *matrizInicial, float *vetorInicial, float *vetorX, float *vetorXNovo) {
	int 	posMatriz;
	int 	i, j;
	float 	soma = 0;
	float 	dp = 0;
	int 	id;
	int 	iter = 0;

	contador = 0;

	omp_set_num_threads(nThreads);

	#pragma omp parallel private (i, j, id, soma, posMatriz)
	{
		while (contador < ITERACOES) {
			id = omp_get_thread_num();

			for(i = id+1; i < N; i+=nThreads) {
				soma = 0;

				for(j = 0; j<N; j++) {
					posMatriz = (i * N) + j;

					if (i != j) {
						soma += matrizInicial[posMatriz] * vetorX[j];
					}
					else {
						dp =  matrizInicial[posMatriz];
					}
				
				}
				vetorXNovo[i] = (vetorInicial[i] - soma) / dp;
			}

			#pragma omp critical
			{
				if (fabs(calcularNorma(vetorX)- calcularNorma(vetorXNovo)) < ERRO) {
					iter = contador;
					contador = ITERACOES;
				}
				else {
					copiarVetor(vetorX, vetorXNovo);
				}
			}

			contador++;
		}
	}

	printf("Iteracoes: %d \n", iter);
}

double tempoCorrente(void){
     struct timeval tval;

     gettimeofday(&tval, NULL);

     return (tval.tv_sec + tval.tv_usec/1000000.0);
}

int main(int argc, char **argv ){
	int 	nThreads;
	int		contador;
	float 	*matrizInicial;
	float 	*vetorInicial;
	float 	*vetorCalculado;
	float 	*vetorCalculadoAnt;
	double 	ti,tf;

	if ( argc != 2){
		printf("%s <num_threads>\n", argv[0]);
		exit(0);
	}

	nThreads = atoi(argv[1]);

	ti = tempoCorrente();

	//Alocar memoria dinamicamente
	matrizInicial 		= (float *) malloc((N * N) * sizeof(float));
	vetorInicial 		= (float *) malloc(N * sizeof(float));
	vetorCalculado 		= (float *) malloc(N * sizeof(float));
	vetorCalculadoAnt 	= (float *) malloc(N * sizeof(float));

	//Iniciar as variaveis
	contador = 0;

	gerarMatriz(matrizInicial);
	gerarVetor(vetorInicial);
	gerarVetorZerado(vetorCalculado);
	gerarVetorZerado(vetorCalculadoAnt);

	//Mostrar valores iniciais
	printf("Matriz:\n");
	mostrarMatriz(matrizInicial);

	printf("Vetor inicial:\n");
	mostrarVetor(vetorInicial);

	//Calcular Jacobi
	calcularJacobi(contador, nThreads, matrizInicial,vetorInicial,vetorCalculado,vetorCalculadoAnt);

	//Mostrar valores calculados
	printf("Vetor calculado:\n");
	mostrarVetor(vetorCalculadoAnt);

	tf = tempoCorrente();
	printf("Tempo = %f\n", tf - ti );

	//Limpar memoria
	free(matrizInicial);
	free(vetorInicial);
	free(vetorCalculado);
	free(vetorCalculadoAnt);
}