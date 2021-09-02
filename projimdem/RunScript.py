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

cam_file = './example/FinseDemoData/CamFinseInit.json'
GCP_file = './example/FinseDemoData/GCPs_pointagev4.csv'
dem_file = './example/FinseDemoData/time_lapse_finse_DSM_midfilled.tif'
viewshed_file = './example/FinseDemoData/viewshed_mid.tif'
image_file = './example/FinseDemoData/2019-05-24_12-00_ori_marked.jpg'
output_file = './example/FinseDemoData/finse_proj_4m.tif'

finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=(
                        [-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,1000], [3.15,3.15,3.15,np.inf,np.inf,np.inf,2000]))
finse.estimate_cam(method='trf', loss='soft_l1')
#finse.proj_GCPs2img()

#finse.ChangeFreeParams(free_param=[ 'Foc','DCx', 'DCy', 'R1','R3', 'R5'],
#                    param_bounds=([-1200,-50,-50,-10,-10,-10], [1500,50,50,10,10,10]))


finse.ChangeFreeParams(free_param=['DCx', 'DCy', 'K1', 'K2', 'K3', 'K4','K5', 'P1', 'P2','P3','P4','P5','P6','P7'],
                    param_bounds=([-50,-50,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1], [50,50,1,1,1,1,1,1,1,1,1,1,1,1]))


finse.estimate_cam('trf',xtol=1e-16, loss='soft_l1', ftol=1e-12)
#%matplotlib
#finse.proj_GCPs2img()

finse_proj = pi.ProjIm2dem(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=finse.new_cam.proj_param,
                           output_file=output_file
                          )
finse_proj.ProjectImage2DEM(return_raster=True, epsg=32632)












finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,1000], [3.15,3.15,3.15,np.inf,np.inf,np.inf,2000]))
#finse.GCPs
finse.estimate_cam()
finse.print_residuals()
finse.proj_GCPs2img()


# With Distortion
finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc', 'DCx', 'DCy', 'R1', 'R3', 'R5'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,1000,930,480,-0.00000001,-0.000000000000001,-0.0000000000000000000001], [3.15,3.15,3.15,np.inf,np.inf,np.inf,2000,990,540,0.00000001,0.000000000000001,0.0000000000000000000001]))

# With Distortion large numbers
finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc', 'DCx', 'DCy', 'R1', 'R3', 'R5'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,1000,-50,-50,-10,-10,-10], [3.15,3.15,3.15,np.inf,np.inf,np.inf,2000,50,50,10,10,10]))

# With Distortion mid numbers
finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc', 'DCx', 'DCy', 'R1', 'R3', 'R5'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,1000,930,480,-1,-1,-1], [3.15,3.15,3.15,np.inf,np.inf,np.inf,2000,990,540,1,1,1]))

# With Distortion mid numbers
finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc', 'DCx', 'DCy', 'R1'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,1000,-50,-50,-1], [3.15,3.15,3.15,np.inf,np.inf,np.inf,2000,50,50,1]))

# With Distortion mid numbers
finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc', 'R1'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,1000,-1000], [3.15,3.15,3.15,np.inf,np.inf,np.inf,2000,1000]))


# With Distortion Model 1
finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=(
                        [-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,1000], [3.15,3.15,3.15,np.inf,np.inf,np.inf,2000]))
finse.estimate_cam(method='trf', loss='soft_l1')


finse.ChangeFreeParams(free_param=['FOC', 'DCx', 'DCy', 'K1', 'K2', 'P1','P2','P3','P4'],
                    param_bounds=([1000,-50,-50,-10,-10,-10,-10,-10,-10], [2000,50,50,10,10,10,10,10,10]))


finse.estimate_cam(method='trf', loss='soft_l1')
finse.print_residuals()
finse.proj_GCPs2img()




cam_file = './example/FinseFromPhoto4D/CamFinseInit.json'
GCP_file = './example/FinseFromPhoto4D/GCPs.csv'
dem_file = './example/FinseFromPhoto4D/DEM_2m.tif'
viewshed_file = './example/FinseFromPhoto4D/viewshed_test.tif'
image_file = './example/FinseFromPhoto4D/DSC03111.JPG'
output_file = './example/FinseFromPhoto4D/finse_proj_2m.tif'

finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,4000], [3.15,3.15,3.15,np.inf,np.inf,np.inf,6000]))

# With Distortion
cam_file = './example/FinseFromPhoto4D/CamFinseAfterNoDist.json'
finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc', 'DCx', 'DCy', 'R1', 'R3', 'R5'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,4000,2728-30,1816-30,-0.00000001,-0.000000000000001,-0.0000000000000000000001], [3.15,3.15,3.15,np.inf,np.inf,np.inf,6000,2728+30,1816+30,0.00000001,0.000000000000001,0.0000000000000000000001]))
# With Distortion mid numbers
finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc', 'DCx', 'DCy', 'R1'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,4000,2728-30,1816-30,-1], [3.15,3.15,3.15,np.inf,np.inf,np.inf,6000,2728+30,1816+30,1]))




#finse.GCPs
finse.estimate_cam(method='trf')
finse.print_residuals()
finse.proj_GCPs2img()


finse_proj = pi.ProjIm2dem(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=finse.new_cam.proj_param,
                           output_file=output_file
                          )
finse_proj.ProjectImage2DEM(return_raster=True, epsg=32632)


