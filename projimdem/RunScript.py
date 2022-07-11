# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 14:36:49 2021

@author: lucg
"""
#########################################
####        Input Cucza Drone        ####
#########################################

from projimdem import resection as rs
from projimdem import projection as pi
import numpy as np

cam_file       = './example/CuczaDemoData/CamCucza.json'
GCP_file       = './example/CuczaDemoData/GCPs.csv'
dem_file       = './example/CuczaDemoData/DEM.tif'
viewshed_file  = './example/CuczaDemoData/viewshed.tif'
image_file     = './example/CuczaDemoData/Abbey-IMG_0209_RGB.jpg'
output_file    = './example/CuczaDemoData/Ortho.tif'
GCP_to_img_file= './example/CuczaDemoData/GCPs_to_img.png'
                   
cucza = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=   ['omega', 'phi', 'kappa', 'X_ini', 'Y_ini', 'Z_ini', 'Foc', 'DCx', 'DCy', 'K1', 'K2', 'K3'],
                    param_bounds=([-6.30  , -6.30, -6.30 ,  -np.inf, -np.inf, -np.inf,  1900,   680,   430, -100, -100, -100],
                                  [ 6.30  ,  6.30,  6.30 ,   np.inf,  np.inf,  np.inf,  2100,   720,   510,  100,  100,  100])
                    )
cucza.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)

# Plot resection results
cucza.project_GCPs_to_img()
from matplotlib import pyplot as plt
plt.savefig(GCP_to_img_file, dpi=300)

proj = pi.Projection(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=cucza.new_cam.proj_param,
                          output_file=output_file
                          )
proj.project_img_to_DEM(return_raster=True, epsg=32632)

plt.figure(figsize=(10,10))
[DCx, DCy] = proj.cam_param[3]
deltaX = proj.X_undistort.min()
deltaY = proj.Y_undistort.min()
plt.imshow(proj.image_undistort)
plt.scatter(proj.pt_proj.X_img - DCx - deltaX, proj.pt_proj.Y_img - DCy - deltaY, c=proj.pt_proj.Z_cam, alpha=0.05)
plt.savefig('./example/CuczaDemoData/DEM_to_image.png', dpi=300)




#########################################
####       Input Finse Photo4D       ####
#########################################

from projimdem import resection as rs
from projimdem import projection as pi
import numpy as np

cam_file = './example/FinseFromPhoto4D/CamFinseInit.json'
GCP_file = './example/FinseFromPhoto4D/GCPs.csv'
dem_file = './example/FinseFromPhoto4D/DEM_2m.tif'
viewshed_file = './example/FinseFromPhoto4D/viewshed.tif'
image_file = './example/FinseFromPhoto4D/DSC03111.JPG'
output_file = './example/FinseFromPhoto4D/finseP4D_proj_2m.tif'
GCP_to_img_file= './example/FinseFromPhoto4D/GCPs_to_img.png'
                   
finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=   ['omega', 'phi', 'kappa', 'X_ini', 'Y_ini', 'Z_ini', 'Foc', 'DCx', 'DCy', 'K1', 'K2', 'K3', 'P1', 'P2'],
                    param_bounds=([-6.30  , -6.30, -6.30 ,  -np.inf, -np.inf, -np.inf,  4800,   2700,   1780, -100, -100, -100, -100, -100],
                                  [ 6.30  ,  6.30,  6.30 ,   np.inf,  np.inf,  np.inf,  4900,   2750,   1840,  100,  100,  100,  100,  100])
                    )
finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)

# Plot resection results
finse.project_GCPs_to_img()
from matplotlib import pyplot as plt
plt.savefig(GCP_to_img_file, dpi=300)

finse_proj = pi.Projection(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=finse.new_cam.proj_param,
                           output_file=output_file
                          )
finse_proj.project_img_to_DEM(return_raster=True, epsg=32632)



#########################################
####        Input Finse Roof         ####
#########################################

from projimdem import resection as rs
from projimdem import projection as pi
import numpy as np

cam_file = './example/FinseDemoData/CamFinseInit.json'
GCP_file = './example/FinseDemoData/GCPs_pointagev4.csv'
dem_file = './example/FinseDemoData/time_lapse_finse_DSM_midfilled.tif'
viewshed_file = './example/FinseDemoData/viewshed.tif'
image_file = './example/FinseDemoData/2019-05-24_12-00_ori_marked.jpg'
#image_file = './example/FinseDemoData/2022-07-08.jpg'
output_file = './example/FinseDemoData/FinseRoof_2019-05-24_12-00.tif'


GCP_to_img_file= './example/FinseDemoData/GCPs_to_img_FinseRoof_2019-05-24_12-00.png'
finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa', 'X_ini','Y_ini','Z_ini', 'Foc','DCx', 'DCy', 'K1', 'K2', 'K3', 'P1', 'P2'],
                    param_bounds=([-6.30, -6.30, -6.30,-np.inf,-np.inf,-np.inf, 1480,900,  500,-100,-100,-100, -100, -100],
                                  [ 6.30,  6.30,  6.30, np.inf, np.inf, np.inf, 1490,1100, 600, 100, 100,  100,  100,  100])
                    )

finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-15)

# Plot resection results
finse.project_GCPs_to_img()
from matplotlib import pyplot as plt
plt.savefig(GCP_to_img_file, dpi=300)

finse_proj = pi.Projection(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=finse.new_cam.proj_param,
                           output_file=output_file
                          )
finse_proj.project_img_to_DEM(return_raster=True, epsg=32632)


plt.figure(figsize=(10,10))

[DCx, DCy] = finse_proj.cam_param[3]
deltaX = finse_proj.X_undistort.min()
deltaY = finse_proj.Y_undistort.min()
plt.imshow(finse_proj.image_undistort)
plt.scatter(finse_proj.pt_proj.X_img - DCx - deltaX, finse_proj.pt_proj.Y_img - DCy - deltaY, c= finse_proj.pt_proj.Z_cam,alpha=0.05)
plt.savefig('./example/FinseDemoData/DEM_to_image.png', dpi=300)

