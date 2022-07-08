# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 14:36:49 2021

@author: lucg
"""


# Input Finse Photo4D
#

from projimdem import resection as rs
from projimdem import projection as pi
import numpy as np

cam_file = './example/FinseFromPhoto4D/CamFinseInit.json'
GCP_file = './example/FinseFromPhoto4D/GCPs.csv'
dem_file = './example/FinseFromPhoto4D/DEM_2m.tif'
viewshed_file = './example/FinseFromPhoto4D/viewshed_test.tif'
image_file = './example/FinseFromPhoto4D/DSC03111.JPG'
output_file = './example/FinseFromPhoto4D/finse_proj_2m_FreeCalibFK1K2K3.tif'
GCP_to_img_file= './example/FinseFromPhoto4D/GCPs_to_img_FreeCalibFK1K2K3.png'

# finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    # free_param=['omega', 'phi', 'kappa',  'K1', 'K2', 'K3'],
                    # param_bounds=([-6.30, -6.30, -6.30,-1,-1,-1],
                                  # [ 6.30,  6.30,  6.30,  1, 1, 1])
                    # )

finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['omega', 'phi', 'kappa', 'X_ini','Y_ini','Z_ini', 'Foc',  'K1', 'K2', 'K3'],
                    param_bounds=([-6.30, -6.30, -6.30,-np.inf,-np.inf,-np.inf, 4600, -1,-1,-1],
                                  [ 6.30,  6.30,  6.30, np.inf, np.inf, np.inf, 5000,  1, 1, 1])
                    )


finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)
#matplotlib
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



# Input Finse Roof
#

from projimdem import resection as rs
from projimdem import projection as pi
import numpy as np

cam_file = './example/FinseDemoData/CamFinseInit.json'
GCP_file = './example/FinseDemoData/GCPs_pointagev4.csv'
dem_file = './example/FinseDemoData/time_lapse_finse_DSM_midfilled.tif'
viewshed_file = './example/FinseDemoData/viewshed.tif'
image_file = './example/FinseDemoData/2019-05-24_12-00_ori_marked.jpg'
image_file = './example/FinseDemoData/2022-07-08.jpg'
output_file = './example/FinseDemoData/FinseRoof_2019-05-24_12-00.tif'

# GCP_to_img_file= './example/FinseDemoData/GCPs_to_img_FinseRoof_2019-05-24_12-00_LOCKED.png'

# finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    # free_param=['Foc'],
                    # param_bounds=([1483],[1485])
                    # )
# finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)
# finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)
# matplotlib
# finse.project_GCPs_to_img()
# from matplotlib import pyplot as plt
# plt.savefig(GCP_to_img_file, dpi=300)


# GCP_to_img_file= './example/FinseDemoData/GCPs_to_img_FinseRoof_2019-05-24_12-00.png'
# finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    # free_param=['omega', 'phi', 'kappa', 'X_ini','Y_ini','Z_ini', 'Foc'],
                    # param_bounds=([-6.30, -6.30, -6.30,-np.inf,-np.inf,-np.inf, 1400],
                                  # [ 6.30,  6.30,  6.30, np.inf, np.inf, np.inf, 1550])
                    # )
					
GCP_to_img_file= './example/FinseDemoData/GCPs_to_img_FinseRoof_2019-05-24_12-00.png'
finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    free_param=['omega', 'phi', 'kappa', 'X_ini','Y_ini','Z_ini', 'Foc', 'DCx', 'DCy', 'K1', 'K2', 'K3', 'K4', 'K5', 'K6', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6'],
                    param_bounds=([-6.30, -6.30, -6.30,  -np.inf,-np.inf,-np.inf,  1480,   900,  500,  -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100],
                                  [ 6.30,  6.30,  6.30,   np.inf, np.inf, np.inf,  1490,  1100,  600,   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100])
                    )

# GCP_to_img_file= './example/FinseDemoData/GCPs_to_img_FinseRoof_2019-05-24_12-00.png'
# finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=',',
                    # free_param=['omega', 'phi', 'kappa', 'X_ini','Y_ini','Z_ini','DCx', 'DCy', 'K1', 'K2', 'K3'],
                    # param_bounds=([-6.30, -6.30, -6.30,-np.inf,-np.inf,-np.inf, 900,  500,-20,-20,-20],
                                  # [ 6.30,  6.30,  6.30, np.inf, np.inf, np.inf, 1100, 600, 20, 20,  20])
                    # )


finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)
#matplotlib
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
plt.imshow(finse_proj.image_undistort)
plt.scatter(finse_proj.pt_proj.X_distort,finse_proj.pt_proj.Y_distort,alpha=0.05)
plt.savefig('./example/FinseDemoData/distortK6P6.png', dpi=600)




plt.figure(figsize=(10,10))
plt.imshow(finse_proj.image_undistort)
plt.scatter(finse_proj.pt_proj.X_undistort,finse_proj.pt_proj.Y_undistort,alpha=0.2)
plt.savefig('./example/FinseDemoData/undistort.png', dpi=600)









# With R3 free calib
# First iteration with nothing free to create the object structure and limit the freedom of parameters
finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['Foc'],
                    param_bounds=([4860],[4862])
                    )
finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)
# Second iteration with some free parameters, freed with limits defined around initial values
finse.change_free_params(free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc', 'K1', 'K2', 'K3'],
                    param_bounds=(
                        [finse.cam.eop.omega-0.2, finse.cam.eop.phi-0.2, finse.cam.eop.kappa-0.2, finse.cam.eop.X_ini-100-finse.x_offset,finse.cam.eop.Y_ini-100-finse.y_offset,finse.cam.eop.Z_ini-100-finse.z_offset,finse.cam.iop.Foc*0.5,-1,-1,-1],
                        [finse.cam.eop.omega+0.2, finse.cam.eop.phi+0.2, finse.cam.eop.kappa+0.2, finse.cam.eop.X_ini+100-finse.x_offset,finse.cam.eop.Y_ini+100-finse.y_offset,finse.cam.eop.Z_ini+100-finse.z_offset,finse.cam.iop.Foc*2, 1, 1, 1]))
finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)

#matplotlib
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

# Fully locked
output_file = './example/FinseFromPhoto4D/finse_proj_2m_Locked.tif'
GCP_to_img_file= './example/FinseFromPhoto4D/GCPs_to_img_Locked.png'
finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['Foc'],
                    param_bounds=([4860],[4862])
                    )
finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)

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


# Fully free

finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['omega', 'phi', 'kappa',  'K1', 'K2', 'K3'],
                    param_bounds=([-6.30, -6.30, -6.30,-1,-1,-1],
                                  [ 6.30,  6.30,  6.30,  1, 1, 1])
                    )

finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)
#matplotlib
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




# With free Foc only (no distortions)
# First iteration with nothing free to create the object structure and limit the freedom of parameters
output_file = './example/FinseFromPhoto4D/finse_proj_2m_FreeCalibF.tif'
GCP_to_img_file= './example/FinseFromPhoto4D/GCPs_to_img_FreeCalibF.png'
finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['Foc'],
                    param_bounds=([4860],[4862])
                    )
finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)
# Second iteration with some free parameters, freed with limits defined around initial values
finse.change_free_params(free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                    param_bounds=(
                        [-6.30, -6.30, -6.30,finse.cam.eop.X_ini-30-finse.x_offset,finse.cam.eop.Y_ini-30-finse.y_offset,finse.cam.eop.Z_ini-30-finse.z_offset,finse.cam.iop.Foc*0.95],
                        [ 6.30,  6.30,  6.30,finse.cam.eop.X_ini+30-finse.x_offset,finse.cam.eop.Y_ini+30-finse.y_offset,finse.cam.eop.Z_ini+30-finse.z_offset,finse.cam.iop.Foc*1.05]))
finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)







# With fixed calib
finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini'],
                    param_bounds=(
                        [-3.15, -3.15, -3.15,-np.inf,-np.inf,-np.inf], [3.15,3.15,3.15,np.inf,np.inf,np.inf]))
finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)
#finse.proj_GCPs2img()

# With R3 free calib
finse = rs.Resection(cam_file, GCP_file, image_file, delimiter_GCP=' ',
                    free_param=['Foc'],
                    param_bounds=([4860],[4862])
                    )
finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)
finse.change_free_params(free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc', 'K1', 'K2', 'K3'],
                    param_bounds=(
                        [-3.15, -3.15, -3.15,finse.cam.eop.X_ini-30-finse.x_offset,finse.cam.eop.Y_ini-30-finse.y_offset,finse.cam.eop.Z_ini-30-finse.z_offset,finse.cam.iop.Foc*0.95,-1e-5,-1e-10,-1e-15],
                        [ 3.15,  3.15,  3.15,finse.cam.eop.X_ini+30-finse.x_offset,finse.cam.eop.Y_ini+30-finse.y_offset,finse.cam.eop.Z_ini+30-finse.z_offset,finse.cam.iop.Foc*1.05, 1e-5, 1e-10, 1e-15]))
finse.estimate_cam(method='trf', loss='soft_l1', ftol=1e-12)






# trying an iterative appraoach...not sure if weights are any good
for i in range(5):
    finse.change_free_params(free_param=['DCx', 'DCy', 'K1', 'K2', 'K3', 'P1', 'P2','P3'],
                        param_bounds=([-50,-50,-1e-5,-1e-10,-1e-15,-1e-2,-1e-2,-1e-2], 
                                      [ 50, 50, 1e-5, 1e-10, 1e-15, 1e-2, 1e-2, 1e-2]))
    
    finse.estimate_cam('trf',xtol=1e-16, loss='soft_l1', ftol=1e-12)
    
    finse.change_free_params(free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc'],
                        param_bounds=([finse.new_cam.omega-0.1, finse.new_cam.phi-0.1, finse.new_cam.kappa-0.1,finse.new_cam.X_ini-2,finse.new_cam.Y_ini-2,finse.new_cam.Z_ini-2, finse.new_cam.Foc*0.95],
                                      [finse.new_cam.omega+0.1, finse.new_cam.phi+0.1, finse.new_cam.kappa+0.1,finse.new_cam.X_ini+2,finse.new_cam.Y_ini+2,finse.new_cam.Z_ini+2, finse.new_cam.Foc*1.05]))
    
    finse.estimate_cam('trf',xtol=1e-16, loss='soft_l1', ftol=1e-12)
    finse.change_free_params(free_param=['omega', 'phi', 'kappa','X_ini','Y_ini','Z_ini', 'Foc', 'K1', 'K2', 'K3'],
                    param_bounds=([finse.new_cam.omega-0.1, finse.new_cam.phi-0.1, finse.new_cam.kappa-0.1,finse.new_cam.X_ini-2,finse.new_cam.Y_ini-2,finse.new_cam.Z_ini-2, finse.new_cam.Foc*0.95, -1e-5,-1e-10,-1e-15],
                                  [finse.new_cam.omega+0.1, finse.new_cam.phi+0.1, finse.new_cam.kappa+0.1,finse.new_cam.X_ini+2,finse.new_cam.Y_ini+2,finse.new_cam.Z_ini+2, finse.new_cam.Foc*1.05, 1e-5,1e-10,1e-15]))

#matplotlib
finse.project_GCPs_to_img()
from matplotlib import pyplot as plt
plt.savefig('./example/FinseFromPhoto4D/GCPs_to_img.png')

finse_proj = pi.Projection(dem_file=dem_file,
                          viewshed_file=viewshed_file,
                          image_file=image_file,
                          cam_param=finse.new_cam.proj_param,
                           output_file=output_file
                          )
finse_proj.project_img_to_DEM(return_raster=True, epsg=32632)

#kwargs={'alpha':.7,'cmap':pyplot.cm.terrain}
#finse_proj.plot_DEM_on_img(**kwargs)










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


