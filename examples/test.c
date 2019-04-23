#include "MO445.h"
#include <stdio.h>
#include<string.h>

int cmpfunc (const void * a, const void * b) {
   return ( *(double*)b - *(double*)a );
}


int main(int argc,char **argv){
	fprintf(stderr,"\nInstructions: ");
	fprintf(stderr,"\nFor the shape descriptors (Moments Invariant, Fourier Descriptor, BAS Descriptor abd Tensor Scale Descriptor) the image \nneeds to be in binary format (0 - background 1 - object)\n");
	char *filename = NULL,*dir_name = NULL, aux_file[100];
	unsigned int len;
	FILE * fractals_file,  *saliences_file;

	Image *img1 = NULL;
	FeatureVector1D	*fvMS1 = NULL, *fvMS2 = NULL;
	FeatureVector1D *fvSS1 = NULL, *fvSS2 = NULL;
	FeatureVector1D **fractals, **saliences;
	fractals = calloc(1400,sizeof(FeatureVector1D*));
	saliences = calloc(1400,sizeof(FeatureVector1D*));
	char name[1400][30];
	int i = 0,j, n=0;
	fprintf(stderr,"Extracting Multiscale Fractal Dimension and Segment Saliences...\n");
	fractals_file = fopen ("fractals.txt","w");
	saliences_file = fopen ("saliences.txt","w");

	dir_name = argv[1];
	printf("%s\n", dir_name);
	while(getline(&filename, &len, stdin) != EOF){
		filename[strlen(filename)-1] = '\0';
		strcpy(name[i], filename);
		strcpy(aux_file,dir_name);
		strcat(aux_file, "/");
		strcat(aux_file, filename);

		printf("%s\n", aux_file);
		img1 = ReadImage(aux_file);

		fvMS1 = MS_ExtractionAlgorithm(img1);
		fractals[i] = fvMS1;

		fvSS1 = SS_ExtractionAlgorithm(img1);
		saliences[i] = fvSS1;

		DestroyImage(&img1);
		i++;
	}
	n = i;
	for(i=0;i<n;i++){
		double fractals_dist[1400];
		double salience_dist[1400];

		for(j=0;j<n;j++){
			fractals_dist[j] = MS_DistanceAlgorithm(fractals[i], fractals[j]);
			salience_dist[j]= SS_DistanceAlgorithm(saliences[i], saliences[j]);
		}
		qsort(fractals_dist, n, sizeof(double),cmpfunc );
		qsort(salience_dist, n, sizeof(double),cmpfunc );
		for(j=0;j<n;j++){
			fprintf(fractals_file, "%s Q0 %s %d %lf STANDARD\n", name[i], name[j], j+1, fractals_dist[j]);

			fprintf(saliences_file, "%s Q0 %s %d %lf STANDARD\n", name[i], name[j], j+1, salience_dist[j]);	
		}
	}
	

	// fprintf(stderr,"\nExtracting Segment Saliences ... ");
	// WriteFeatureVector1D(fvSS1, "results/segmentsaliences_mpeg7-bat.txt");

	// fprintf(stderr,"\nExtracting Multiscale Fractal Dimension ... ");
	// fvMS2 = MS_ExtractionAlgorithm(img2);
	// WriteFeatureVector1D(fvMS2, "results/multiscale_mpeg7-bird.txt");

	// fprintf(stderr,"\nExtracting Segment Saliences ... ");
	// fvSS2 = SS_ExtractionAlgorithm(img2);
	// WriteFeatureVector1D(fvSS2, "results/segmentsaliences_mpeg7-bird.txt");

	// // CS and Euclidean Distance
	// fprintf(stderr,"%lf",SS_DistanceAlgorithm(fvSS1, fvSS2));
	// fprintf(stderr,"%lf",Fourier_DistanceAlgorithm(fvMS1, fvMS2));


	// DestroyCImage(&cimg1);	
	// DestroyCImage(&cimg2);	DestroyImage(&img2);
	// DestroyFeatureVector1D(&fvBIC1); DestroyFeatureVector1D(&fvBIC2);
	// DestroyFeatureVector1D(&fvMoments1); DestroyFeatureVector1D(&fvMoments2);
	// DestroyFeatureVector1D(&fvFourier1); DestroyFeatureVector1D(&fvFourier2);
	// DestroyFeatureVector1D(&fvBAS1); DestroyFeatureVector1D(&fvBAS2);
	// DestroyFeatureVector1D(&fvTensor1); DestroyFeatureVector1D(&fvTensor2);
	// DestroyFeatureVector1D(&fvMS1); DestroyFeatureVector1D(&fvMS2);
	// DestroyFeatureVector2D(&fvCS1); DestroyFeatureVector2D(&fvCS2);
	// DestroyFeatureVector1D(&fvSS1); DestroyFeatureVector1D(&fvSS2);


	return 0;
}
