#include "MO445.h"
#include <stdio.h>
#include<string.h>
#include <pthread.h> 

typedef struct distance{
	double salience_dist;
	double fractals_dist;
	char *name;
}distance;

distance distmat[1500][1500];
bool been[1500][1500];

// Cria os vetor de distancias
distance distances[1400];

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
	char *filename = NULL, *dir_name = NULL;
	char aux_file[100], name[1400][30];
	unsigned int len;
	int i = 0,j, n=0;
	char output_fractals[1400][70];
	char output_saliences[1400][70];

	FILE *fractals_file,  *saliences_file;
	Image *img1 = NULL;

	FeatureVector1D	*fvMS1 = NULL,*fvSS1 = NULL;
	FeatureVector1D **fractals, **saliences;
	fractals = calloc(1400,sizeof(FeatureVector1D*));
	saliences = calloc(1400,sizeof(FeatureVector1D*));

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
		printf("i: %d\n", i);
	}

	n = i;
	fractals_file = fopen("fractals.txt","w");
	saliences_file = fopen("saliences.txt","w");




	for(i=0;i<n;i++){
		// Calcula a distancia par a par
		double max_frac = 0;
		double max_sal = 0;
		for(j=0;j<n;j++){
			distances[j].name = name[j];
			if(been[i][j]){
				distances[j] = distmat[i][j];
			}else{
				been[i][j] = true;
				distances[j].fractals_dist = MS_DistanceAlgorithm(fractals[i], fractals[j]);
				distances[j].salience_dist = SS_DistanceAlgorithm(saliences[i], saliences[j]);
				distmat[i][j] = distmat[j][i] = distances[j];
			}
			
			if(distances[j].fractals_dist > max_frac)
				max_frac = distances[j].fractals_dist;
			if(distances[j].salience_dist > max_sal)
				max_sal = distances[j].salience_dist;
		}



		// Ordena os vetores em ordem decrescente e escreve as distancias
		// qsort(distances, n, sizeof(distance),cmpfractals);
		printf("%s\n", name[i]);
		for(j=0;j<n;j++){
			sprintf(output_fractals[j], "%s Q0 %s %d %lf STANDARD\n", name[i], distances[j].name, i*n + j, 1.0-(distances[j].fractals_dist/max_frac));
		}
		// qsort(distances, n, sizeof(distance),cmpsaliences);
		for(j=0;j<n;j++){
			sprintf(output_saliences[j], "%s Q0 %s %d %lf STANDARD\n", name[i], distances[j].name, i*n + j, 1.0-(distances[j].salience_dist/max_sal));	
		}
	}

	for(j=0;j<n;j++){
		fprintf(fractals_file, "%s", output_fractals[j]);	
	}
	for(j=0;j<n;j++){
		fprintf(saliences_file, "%s", output_saliences[j]);	
	}
	fclose(fractals_file);
	fclose(saliences_file);

	return 0;
}
