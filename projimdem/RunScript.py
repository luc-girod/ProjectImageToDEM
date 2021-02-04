# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 14:36:49 2021

@author: lucg
"""


# Input Finse
#

from projimdem import resection as rs
from projimdem import ProjectImageOnDEM as pi

cam_file = './example/FinseDemoData/CamFinseInit.json'
GCP_file = './example/FinseDemoData/GCPs.csv'
dem_file = './example/FinseDemoData/time_lapse_finse_DSM_mid.tif'
viewshed_file = './example/FinseDemoData/viewshed_test.tif'
image_file = './example/FinseDemoData/2019-05-24_12-00.jpg'
output_file = './example/FinseDemoData/finse_proj_1m.tif'

finse = rs.resection(cam_file, GCP_file, image_file)

#finse.GCPs
finse.estimate_cam()
finse.print_residuals()
finse.proj_GCPs2img()


finse_proj = pi.ProjIm2dem(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=finse.estimate.new_cam,
                           output_file=output_file
                          )
finse_proj.ProjectImage2DEM(return_raster=True, epsg=32632)