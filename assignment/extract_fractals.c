#include <stdio.h>
#include "MO445.h"

int main(){
    char *filename = "./pgm_imgs/beetle-4.pgm";
    Image *img;
    img = ReadImage(filename);
    FeatureVector1D *fv;
    
    fv = MS_ExtractionAlgorithm(img);

    for(int i = 0;i< (*fv).n;i++) printf("%lf ",(*fv).X[i]);

    return 0;
}