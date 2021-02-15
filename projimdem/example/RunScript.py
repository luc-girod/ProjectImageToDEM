# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 14:36:49 2021

@author: lucg
"""


# Input Finse
#

from projimdem import resection as rs
from projimdem import ProjectImageOnDEM as pi
import numpy as np

cam_file = './FinseFromPhoto4D/CamFinseInit.json'
GCP_file = './FinseFromPhoto4D/GCPs.csv'
dem_file = './FinseFromPhoto4D/DEM_2m.tif'
viewshed_file = './FinseFromPhoto4D/viewshed_test.tif'
image_file = './FinseFromPhoto4D/DSC03111.JPG'
output_file = './FinseFromPhoto4D/finse_proj_2m.tif'

#================================
# 1. Compute transformation
finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,4000], [3.15,3.15,3.15,np.inf,np.inf,np.inf,6000]))
#finse.GCPs
finse.estimate_cam()
finse.print_residuals()

#================================
# 2. Visualize transform
finse.proj_GCPs2img()

#================================
# 3. Projection
finse_proj = pi.ProjIm2dem(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=finse.new_cam.proj_param,
                           output_file=output_file)
finse_proj.ProjectImage2DEM(return_raster=True, epsg=32632)



'''
cam_file = './FinseDemoData/CamFinseInit.json'
GCP_file = './FinseDemoData/GCPs.csv'
dem_file = './FinseDemoData/time_lapse_finse_DSM_mid.tif'
viewshed_file = './FinseDemoData/viewshed_test.tif'
image_file = './FinseDemoData/2019-05-24_12-00.jpg'
output_file = './FinseDemoData/finse_proj_1m.tif'
