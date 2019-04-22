echo $1;
cd $1;
for file in *; do
    echo $file
    convert $file -compress None '../pgm_imgs/'${file%.*}'.pgm'
done;