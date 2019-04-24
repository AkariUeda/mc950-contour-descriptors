#include "MO445.h"
#include <stdio.h>
#include<string.h>

typedef struct distance{
	double salience_dist;
	double fractals_dist;
	char *name;
}distance;

int cmpfractals (const void * a, const void * b) {
	const struct distance *da = a, *db = b;

   return ( db->fractals_dist - da->fractals_dist );
}

int cmpsaliences (const void * a, const void * b) {
	const struct distance *da = a, *db = b;

   return ( db->salience_dist - da->salience_dist );
}


int main(int argc,char **argv){
	// Variaveis e variaveis
	char *filename = NULL,*dir_name = NULL, aux_file[100], name[1400][30];
	unsigned int len;
	int i = 0,j, n=0;

	FILE * fractals_file,  *saliences_file;
	Image *img1 = NULL;

	FeatureVector1D	*fvMS1 = NULL, *fvMS2 = NULL;
	FeatureVector1D *fvSS1 = NULL, *fvSS2 = NULL;
	FeatureVector1D **fractals, **saliences;
	fractals = calloc(1400,sizeof(FeatureVector1D*));
	saliences = calloc(1400,sizeof(FeatureVector1D*));

	fractals_file = fopen ("fractals.txt","w");
	saliences_file = fopen ("saliences.txt","w");

	// Diretorio com o dataset
	dir_name = argv[1];

	// Gera os fractais e as saliencias
	while(getline(&filename, &len, stdin) != EOF){
		// Bruxaria pra pegar o nome da imagem
		filename[strlen(filename)-1] = '\0';
		strcpy(name[i], filename);
		strcpy(aux_file,dir_name);
		strcat(aux_file, "/");
		strcat(aux_file, filename);
		printf("%s\n", aux_file);

		// Le a imagem
		img1 = ReadImage(aux_file);

		// Calcula os fractais
		fvMS1 = MS_ExtractionAlgorithm(img1);
		fractals[i] = fvMS1;

		//Calcula as saliencias
		fvSS1 = SS_ExtractionAlgorithm(img1);
		saliences[i] = fvSS1;

		DestroyImage(&img1);
		i++;
	}

	n = i;
	for(i=0;i<n;i++){
		// Cria os vetor de distancias
		distance *distances;
		distances = calloc(n,sizeof(distance));

		// Calcula a distancia par a par
		for(j=0;j<n;j++){
			distances[j].name = name[j];
			distances[j].fractals_dist = MS_DistanceAlgorithm(fractals[i], fractals[j]);
			distances[j].salience_dist = SS_DistanceAlgorithm(saliences[i], saliences[j]);
		}

		// Ordena os vetores em ordem decrescente e escreve as distancias
		qsort(distances, n, sizeof(double),cmpfractals);
		for(j=0;j<n;j++){
			fprintf(fractals_file, "%s Q0 %s %d %lf STANDARD\n", name[i], distances[j].name, i*n + j, distances[j].fractals_dist);
		}
		qsort(distances, n, sizeof(double),cmpsaliences);
		for(j=0;j<n;j++){
			fprintf(saliences_file, "%s Q0 %s %d %lf STANDARD\n", name[i], distances[j].name, i*n + j, distances[j].salience_dist);	
		}
	}
	
	fclose(fractals_file);
	fclose(saliences_file);

	return 0;
}
