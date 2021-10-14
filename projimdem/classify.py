'''
Set of functions to classify pixels in images in 2 or more categories. Original intent is the classify snow vs. no snow area.

Technique 0, Manual thresholding

Technique 1, MultuOsu (scikit-image)

Technique 2, Salzano et al (2019)

Technique 3, PCA

Technique 4, Clustering

'''

import matplotlib.pyplot as plt

from skimage.color import rgb2hsv
from skimage.filters import threshold_multiotsu

#======================================
# Technique 0, Manual Thresholding

def classify_threshold(img_rgb, threshold, mask=None):
    img_value = rgb2hsv(img_rgb)[:, :, 2]
    
    regions = img_value < threshold
    if mask is not None:
        regions = regions * mask
    return regions

def classify_blue_threshold(img_rgb, threshold, mask=None):
        
    regions = img_rgb[:,:,1] < threshold
    if mask is not None:
        regions = regions * mask
    return regions
    
#======================================
# Technique 1, MultiOtsu Thresholding
#
# reference: https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.threshold_multiotsu

def classify_multiotsu(img_rgb, nb_classes=2, mask=None):
    thresholds = threshold_multiotsu(img_rgb, classes=nb_classes)
    regions = np.digitize(img_rgb, bins=thresholds)
    
    if mask is not None:
        regions = regions * mask
    
    return regions


#======================================
# Technique 2, Salzano et al (2019)
#
# reference: https://www.mdpi.com/2076-3263/9/2/97/htm#B32-geosciences-09-00097

def classify_snow_salzano(img_rgb):
    # See with Kris to understand method detailed in Salzano et al. 2019 (https://www.mdpi.com/2076-3263/9/2/97/htm#B32-geosciences-09-00097)
    


