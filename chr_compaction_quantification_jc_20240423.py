# this script:
# finds bounding box of the chromatin reagion
# rotate the image according to the rotation angle of the bounding box
# finds chromatin area' centroid in the rotated image
# in the rotated image, calculates inner compartment intensity proportion. Inner compartment is a circle whose center is the chromatin centroid (see line 4) and whose area is 11.1% of the total ROI area (ROI != chromatin != image)
# in the rotated image, profiles intensity along a vertical line that goes through the chromatin centroid

# the intensity profiling part is modified from https://stackoverflow.com/questions/72227233/extract-intensity-profile-along-a-line-from-image
# bounding box find is modified from https://colab.research.google.com/drive/1tBnU6Xdc58dLHzBZPwKz_Y_mCb1Pp40e?usp=sharing#scrollTo=BiD_w9dAVoUo

# directory structure requirement for this script: main folder -> subfolders whose names have - Copy - Copy in the end -> subfolders whose names contain pgc (pgc1 or pgc2)

# input for this script: in the pgc folder, you need to two tiffs whose name are *_ori.tif and *_fill.tif
# *_ori.tif will be used for intensity calculation, rotation, etc.
# *_fill.tif indicates the ROI

# output of this script:
# in the pgc folder: a bounding box tif, a thresholded rotated tif, a intensity profiling tiff on rotated image, a CSV with line intensity profiling result
# in the main folder: a CSV with all quantified images' metrics (RawIntDen, ROI_area, CV, etc.)

###################################################################################################################################################
# IF YOU RUN INTO ERRORS, YOU MIGHT WANT TO CHECK IF YOUR WHOLE DIRECTORY NAME IS TOO LONG#
###################################################################################################################################################

# load modules
# If you exited python shell to install some modules via python3 -m pip --upgrade <module-name>,
# to run this script, you will need to reenter python shell by entering python3 in the terminal
from PIL import Image
import numpy as np
import os
import re
import pandas
import math
import matplotlib.pyplot as plt
import cv2
import skimage.measure as ski
from scipy import ndimage, datasets

# change dir
# os.chdir('C:\\Users\\chen3\\OneDrive\\Downloads\\test\\test - Copy - Copy')

os.chdir('C:\\Users\\chen3\\Dropbox (Duke Bio_Ea)\\Baugh Lab\\2_baugh_lab_users\\Jingxian\\confocal_Zeiss880')
os.chdir('.\\20240502_LRB663_IAA_36h_post_bleach_airyscan_zstack\\processed_autostrength6_3D') # change this accordingly
os.chdir('.\\Image 3_Out TIFFs - Copy - Copy\\pgc2')

# set up metrics you want to get
cols = ['directory', 'image', 'RawIntDen', 'ROI_area', 'MeanGrayValue', 'std', 'CV', 'inner_comp_IntDen', 'prop_Int_inner_comp']
lst = []

file = 'Image-0006.tif_1_ori.tif'

file_index = re.sub("_ori.tif", '', file) # pattern, replaced with what, input string
roi_file = ('%s_fill.tif' % file_index)
                                        
im = Image.open(file)
imarray = np.array(im)

ROI = Image.open(roi_file)
ROIarray = np.array(ROI)

ROI_indices = np.where(ROIarray > 0)
ROI_values = imarray[ROI_indices]

mean = np.mean(ROI_values)
standard_deviation = np.std(ROI_values)
area = len(ROI_values) # same result as np.count_nonzero(ROIarray)
radius = np.sqrt(0.111 * area / math.pi)

image = cv2.imread(file) # cv2.imread reads BGR (blue green red) even if image is single-channel. three channels will have the same values
image_red = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

# plt.rcParams["figure.figsize"] = (10, 10) # Increase size of the visualization. rc stands for runtime configuration
# plt.imshow(image_red)
# plt.title('Preview of the original image') # add title to the canvas
# plt.show()
# plt.clf() # IMPORTANT!! without clearing the canvas, all future plots will be on top of each other on the same canvas

# convert the picture to single channel
image_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

# thresholding result will be used for contouring
# contouring works best on binary pictures so here let's convert the pic to grayscale first
# and then convert everything within the ROI to no color (zero)
# and everything outside ROI full color (255)
# read more about cv2.THRESH_OTSU adaptive thresholding here https://docs.opencv.org/4.x/d7/d1b/group__imgproc__misc.html#ggaa9e58d2860d4afa658ef70a9b1115576a95251923e8e22f368ffa86ba8bce87ff
# currently adaptive thresholding cv2.THRESH_OTSU only works on 8-bit
# cv2.THRESH_OTSU outputs an optimal thresholding value different from the specified value (here, 0)
thresh_result = cv2.threshold(image_gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
thresh_img = thresh_result[1] # second element of the tuple is the thresholded image

# find contours
# read more about what hierarchy here https://opencvpython.blogspot.com/2013/01/contours-5-hierarchy.html
contours, hierarchy = cv2.findContours(thresh_img, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE) # contour extraction mode is RETR_TREE, returning a tree/hierarchy of contours; contour approximation method is CHAIN_APPROX_SIMPLE

# find bounding box and rotation angle
len(contours) > 1


contours
rect = cv2.minAreaRect(contours[1]) # second element of the tuple with index 1 is the chromatin (index 0 is everything except for chromatin)
x_y = rect[0] # x_y is centroid of rectangle in terms of x, y instead of nrow, ncol
angle_of_rotation = rect[2]
box = np.int0(cv2.boxPoints(rect))
cv2.drawContours(image_red, [box], 0, (36, 255, 12), 1) # second argument needs a list and [] is used to construct a list. third argument is index of the list, fourth and fifth arguments are contours color (here is RGB not BGR...) and thickness
plt.imshow(image_red)
plt.scatter(x_y[0], x_y[1]) # input is x, y instead of nrow, ncol
plt.show()
box_file = ('%s_bounding_box.tif' % file_index)
plt.savefig(box_file)
plt.clf()
                            
# rotate image
image_red = cv2.cvtColor(image, cv2.COLOR_BGR2RGB) # important! need a clean reload of image_red
image_rotated = ndimage.rotate(image_red, angle_of_rotation, reshape = True, order = 1) # reshape = True contains all the input array values. order = 1 does linear spline interpolation
imarray_rotated = ndimage.rotate(imarray, angle_of_rotation, reshape = True, order = 1)

plt.imshow(image_rotated)
plt.show()
plt.clf()

# after the interpolation, IntDen varies a little bit
np.sum(imarray)
np.sum(imarray_rotated) 

# use the integrated intensity in the rotated image
IntDen = np.sum(imarray_rotated)

# binarilize and threshold to paint chromatin region in the rotated image
image_gray_rotated = cv2.cvtColor(image_rotated, cv2.COLOR_RGB2GRAY)
thresh_result_rotated = cv2.threshold(image_gray_rotated, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
thresh_img_rotated = thresh_result_rotated[1]

# get chromatin centroid
chr_indices = np.where(thresh_img_rotated == 0)
chr_centroid_nrow = (np.min(chr_indices[0]) + np.max(chr_indices[0]))/2
chr_centroid_ncol = (np.min(chr_indices[1]) + np.max(chr_indices[1]))/2

plt.imshow(thresh_img_rotated)
plt.scatter(chr_centroid_ncol, chr_centroid_nrow) # input is x, y instead of nrow, ncol
rotated_thresholded_file = ('%s_rot_thres.tif' % file_index)
plt.savefig(rotated_thresholded_file)
plt.clf()

# calculate proportion of IntDen in the inner compartment
total_row, total_col = np.shape(image_rotated)[:2]
row, col = np.ogrid[:total_row, :total_col]
distance_from_centroid = np.sqrt((row - chr_centroid_nrow)**2 + (col - chr_centroid_ncol)**2)

outer_region_indices = distance_from_centroid > radius
circlearray = imarray_rotated.copy()
circlearray[outer_region_indices] = 0

circle_IntDen = np.sum(circlearray)

# profile intensity along a vertical line that goes through the chromatin centroid in the rotated image
start = (0, chr_centroid_ncol) # order is nrow, ncol instead of x, y 
end = (image_rotated.shape[0], chr_centroid_ncol) # order is nrow, ncol instead of x, y 

scan_width = 2
profile = ski.profile_line(imarray_rotated, start, end, linewidth = scan_width) # linewidth is scan width; start and end are specified by nrow, ncol instead of x, y

# plot image, profile line position, and intensity histogram at the same time
# plt.grid(False)
plt.imshow(image_rotated)
plt.plot([start[1], end[1]], [start[0], end[0]], 'w--', lw = scan_width) # white dashed line with linewidth being 2; order of arguments: [x1, x2], [y1, y2] instead of [nrow1, nrow2], [ncol1, ncol2]
plt.plot(profile, list(range(total_row + 1)), "w-") # white solid line
# plt.plot(profile, list(range(total_row + 1)), "k-") # black solid line
plt.xlim(0, imarray_rotated.shape[1])
plt.ylim(imarray_rotated.shape[0], 0)
plt.show()
profile_img_file = ('%s_intensity_profile.tif' % file_index)
plt.savefig(profile_img_file)
plt.clf()

# save intensity profiling result as a dataframe
profile_df = pandas.DataFrame({'raw_int' : profile, 'nrow' : list(range(total_row + 1))})
profile_df_file = ('%s_intensity_profile.csv' % file_index)
profile_df.to_csv(profile_df_file)

# append metrics to the metric list
lst.append([final_dir, file_index, IntDen, area, mean, standard_deviation, standard_deviation/mean, circle_IntDen, circle_IntDen/IntDen])

                            
df = pandas.DataFrame(lst, columns = cols) 
df

os.chdir(rootdir)
df.to_csv('chr_compaction_metrics.csv')



