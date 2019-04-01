

## Assignment

1. Convert .gif to .jpeg files and for each image:
2. Extract Multi-scale fractals and calculate the Euclidean distance
3. Calculate the segment saliences and calculate the contour salience(?)
4. For each image, create a list ranking the other 1399 images based on the euclidean distance and on contour salience
5. Use trec_eval to calculate precision and recall of the ranked lists (the ground truth is an input for trec_eval)

## Some notes

* We are working with contour based descriptors.
* Some images that have contours inside the object may present issues when running the code
* Due date: last week of April
* Submission: Explain what we have done, and plot the precision x recall curves (one page max.)
* Work in pairs
* Bull's Eye Score: Recall at the first 20 retrieved images.

