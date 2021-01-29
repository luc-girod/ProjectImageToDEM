# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 14:36:49 2021

@author: lucg
"""
# Input Finse
#
camera_file = './/FinseDemoData//CamFinseInit.inp'
point_file = './/FinseDemoData//GCPs_WebcamFinse_Centered3.inp'
Foc=1484
dem_file='.//FinseDemoData//time_lapse_finse_DSM_mid.tif'
viewshed_file='.//FinseDemoData//viewshed_mid.tif'
image_file='.//FinseDemoData//2019-05-24_12-00.jpg'
output='.//FinseDemoData//output.tif'


data = CollinearityData(camera_file, point_file)

# Perform spatial resection
x0 = np.zeros(6)
# initilaize guesses for eop as read from file
eop = data.eop
x0[0] = eop['omega']
x0[1] = eop['phi']
x0[2] = eop['kappa']
x0[3] = eop['XL']
x0[4] = eop['YL']
x0[5] = eop['ZL']

x, cov_x, info, msg, ier = leastsq(coll_func, x0, Dfun=coll_Dfunc, full_output=True)

R=RotMatrixFromAngles(x[0],x[1],x[2])
C=[x[3]+410000,x[4]+6710000,x[5]]
C=[419169.86,6718421.39,1215]
R=RotMatrixFromAngles(np.pi/2.1,-np.pi/8,-np.pi/8)
aCam=[C,R,Foc]
ProjectImage2DEM(dem_file, viewshed_file, image_file, output, aCam, dem_nan_value=1137.75)




#Input Cucza
camera_file = './/CuczaDemoData//CamCucza.inp'
point_file = './/CuczaDemoData//GCPs_Centered.inp'
Foc=2011.887
dem_file='.//CuczaDemoData//DEM.tif'
image_file='.//CuczaDemoData//Abbey-IMG_0209.jpg'
output='.//CuczaDemoData//output.ply'

data = CollinearityData(camera_file, point_file)

# Perform spatial resection
x0 = np.zeros(6)
# initilaize guesses for eop as read from file
eop = data.eop
x0[0] = eop['omega']
x0[1] = eop['phi']
x0[2] = eop['kappa']
x0[3] = eop['XL']
x0[4] = eop['YL']
x0[5] = eop['ZL']

x, cov_x, info, msg, ier = leastsq(coll_func, x0, Dfun=coll_Dfunc, full_output=True)

R=RotMatrixFromAngles(x[0],x[1],x[2])
C=[x[3],x[4],x[5]]
aCam=[C,R,Foc]

# R=[[0.993534882295323163,0.0929006109947966841,0.0652526944977123435],[0.0878277479180877285,-0.993176223631756505,0.0767285833845516713],[0.0719355569802246908,-0.0705015268583363691,-0.994914473888378059]]
# aCam=[[209.89527614679403023,91.20530793577831,107.031846453497209],R,2011.8874387887106]
ProjectImage2DEM(dem_file, image_file, output, aCam)