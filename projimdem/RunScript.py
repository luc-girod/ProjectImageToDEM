# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 14:36:49 2021

@author: lucg
"""


# Input Finse
#

from projimdem import resection as rs
from projimdem import projection as pi
import numpy as np

cam_file = './example/FinseDemoData/CamFinseInit.json'
GCP_file = './example/FinseDemoData/GCPs_pointagev4.csv'
dem_file = './example/FinseDemoData/time_lapse_finse_DSM_midfilled.tif'
viewshed_file = './example/FinseDemoData/viewshed_mid.tif'
image_file = './example/FinseDemoData/2019-05-24_12-00_ori_marked.jpg'
output_file = './example/FinseDemoData/finse_proj_4m.tif'

finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=(
                        [-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,1000], [3.15,3.15,3.15,np.inf,np.inf,np.inf,2000]))
finse.estimate_cam(method='trf', loss='soft_l1')
#finse.proj_GCPs2img()



# trying an iterative appraoach...not sure if wights are any good
for i in range(5):
    finse.change_free_params(free_param=['DCx', 'DCy', 'K1', 'K2', 'K3', 'P1', 'P2','P3'],
                        param_bounds=([-50,-50,-1e-5,-1e-10,-1e-15,-1e-2,-1e-2,-1e-2], 
                                      [ 50, 50, 1e-5, 1e-10, 1e-15, 1e-2, 1e-2, 1e-2]))
    
    finse.estimate_cam('trf',xtol=1e-16, loss='soft_l1', ftol=1e-12)
    
    finse.change_free_params(free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                        param_bounds=([finse.new_cam.omega-0.1, finse.new_cam.phi-0.1, finse.new_cam.kappa-0.1,finse.new_cam.X_ini-2,finse.new_cam.Y_ini-2,finse.new_cam.Z_ini-2, finse.new_cam.Foc*0.95],
                                      [finse.new_cam.omega+0.1, finse.new_cam.phi+0.1, finse.new_cam.kappa+0.1,finse.new_cam.X_ini+2,finse.new_cam.Y_ini+2,finse.new_cam.Z_ini+2, finse.new_cam.Foc*1.05]))
    
    finse.estimate_cam('trf',xtol=1e-16, loss='soft_l1', ftol=1e-12)


%matplotlib
finse.project_GCPs_to_img()

from projimdem import projection as pi
finse_proj = pi.Projection(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=finse.new_cam.proj_param,
                           output_file=output_file
                          )
finse_proj.project_img_to_DEM(return_raster=False, epsg=32632)

kwargs={'alpha':.7,'cmap':pyplot.cm.terrain}
finse_proj.plot_DEM_on_img(**kwargs)



'''







finse.change_free_params(free_param=['DCx', 'DCy', 'K1', 'K2', 'K3', 'K4','K5', 'P1', 'P2','P3','P4','P5','P6','P7'],
                    param_bounds=([-50,-50,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1], [50,50,1,1,1,1,1,1,1,1,1,1,1,1]))


finse.change_free_params(free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc','DCx', 'DCy', 'K1', 'K2', 'K3', 'P1', 'P2','P3'],
                    param_bounds=([finse.new_cam.omega-0.1, finse.new_cam.phi-0.1, finse.new_cam.kappa-0.1,finse.new_cam.X_ini-2,finse.new_cam.Y_ini-2,finse.new_cam.Z_ini-2, finse.new_cam.Foc*0.95, -50,-50,-1e-5,-1e-10,-1e-15,-1e-5,-1e-5,-1e-5],
                                  [finse.new_cam.omega+0.1, finse.new_cam.phi+0.1, finse.new_cam.kappa+0.1,finse.new_cam.X_ini+2,finse.new_cam.Y_ini+2,finse.new_cam.Z_ini+2, finse.new_cam.Foc*1.05, 50,50,1e-5,1e-10,1e-15,1e-5,1e-5,1e-5]))

finse.change_free_params(free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=([finse.new_cam.omega-0.1, finse.new_cam.phi-0.1, finse.new_cam.kappa-0.1,finse.new_cam.X_ini-2,finse.new_cam.Y_ini-2,finse.new_cam.Z_ini-2, finse.new_cam.Foc*0.95],
                                  [finse.new_cam.omega+0.1, finse.new_cam.phi+0.1, finse.new_cam.kappa+0.1,finse.new_cam.X_ini+2,finse.new_cam.Y_ini+2,finse.new_cam.Z_ini+2, finse.new_cam.Foc*1.05]))









finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,1000], [3.15,3.15,3.15,np.inf,np.inf,np.inf,2000]))
#finse.GCPs
finse.estimate_cam()
finse.print_residuals()
finse.project_GCPs_to_img()





cam_file = './example/FinseFromPhoto4D/CamFinseInit.json'
GCP_file = './example/FinseFromPhoto4D/GCPs.csv'
dem_file = './example/FinseFromPhoto4D/DEM_2m.tif'
viewshed_file = './example/FinseFromPhoto4D/viewshed_test.tif'
image_file = './example/FinseFromPhoto4D/DSC03111.JPG'
output_file = './example/FinseFromPhoto4D/finse_proj_2m.tif'

finse = rs.resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=([-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf,4000], [3.15,3.15,3.15,np.inf,np.inf,np.inf,6000]))




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


