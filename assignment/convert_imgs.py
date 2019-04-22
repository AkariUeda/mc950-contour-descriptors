import cv2
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib import cm
import os
import sys

directory = sys.argv[1]
for file in os.listdir(directory):
    filename = directory + '/' + file
    im_in =Image.open(filename)
    im_in = im_in.convert('L')
    im_in = np.array(im_in)
    
    # Threshold.
    # Set values equal to or above 220 to 0.
    # Set values below 220 to 255.
    last_x = im_in.shape[0]
    print(im_in.shape)
    print(last_x)

    # Create a border to prevent images that touch the border
    # to become full white
    bordersize=10
    mean= 0
    im_in = cv2.copyMakeBorder(im_in, top=bordersize, bottom=bordersize, left=bordersize, right=bordersize, borderType= cv2.BORDER_CONSTANT, value=[mean,mean,mean] )
    im_in = cv2.GaussianBlur(im_in,(5,5),0)

    th, im_in = cv2.threshold(im_in, 220, 255, cv2.THRESH_BINARY);

    # Copy the thresholded image.
    im_floodfill = im_in.copy()
    # Mask used to flood filling.
    # Notice the size needs to be 2 pixels than the image.
    h, w = im_floodfill.shape[:2]
    mask = np.zeros((h+2, w+2), np.uint8)
    # Floodfill from point (0, 0)
    cv2.floodFill(im_floodfill, mask, (0,0), 255);
    
    # Invert floodfilled image
    im_floodfill_inv = cv2.bitwise_not(im_floodfill)
    
    # Combine the two images to get the foreground.
    im_out = im_in | im_floodfill_inv

    im_out = Image.fromarray(np.uint8(cm.gist_earth(im_out)*255))
    # im_out = Image.fromarray((im_out * 255).astype(np.uint8))
    im_out = im_out.convert('L')
    # plt.imshow( im_out, cmap='gray')
    im_out.save('jpeg_imgs/' + file[:-4] + '.jpeg')
