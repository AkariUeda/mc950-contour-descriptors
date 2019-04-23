#include "MO445.h"
#include <stdio.h>

int main(int argc,char **argv){
	fprintf(stderr,"\nInstructions: ");
	fprintf(stderr,"\nFor the shape descriptors (Moments Invariant, Fourier Descriptor, BAS Descriptor abd Tensor Scale Descriptor) the image \nneeds to be in binary format (0 - background 1 - object)\n");
	char *filename;
	Image *img1 = NULL;
	FeatureVector1D	*fvMS1 = NULL, *fvMS2 = NULL;
	FeatureVector1D *fvSS1 = NULL, *fvSS2 = NULL;
	FeatureVector1D **features;
	features = calloc(1400,sizeof(FeatureVector1D*);
	int i;

	while(scanf("%s\n", &filename) != EOF){
		printf("%s\n", filename);
	}
	// img1 = ReadImage(argv[1]);
	// img2 = ReadImage("figs/mpeg7-bird.pgm");
		
	// fprintf(stderr,"\nExtracting Multiscale Fractal Dimension ... ");
	// fvMS1 = MS_ExtractionAlgorithm(img1);
	// WriteFeatureVector1D(fvMS1, "results/multiscale_mpeg7-bat.txt");

	// for (i=0; i < desc->n; i++)
    // printf(fp,"%f\n",desc->X[i]);

	// fprintf(stderr,"\nExtracting Segment Saliences ... ");
	// fvSS1 = SS_ExtractionAlgorithm(img1);
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


	// DestroyCImage(&cimg1);	DestroyImage(&img1);
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
