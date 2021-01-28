# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:27:46 2021

@author: lucg
"""
import numpy as np
#import pandas as pd
import gdal
from optparse import OptionParser
from matplotlib import pyplot
import plyfile
import time
import resection_leastsq_Dfun as resec
from scipy.optimize import leastsq

def RotMatrixFromAngles(O,P,K):
    
    RX=np.array([[1,0,0],
                 [0,np.cos(O),-np.sin(O)],
                 [0,np.sin(O),np.cos(O)]])    
    RY=np.array([[np.cos(P),0,np.sin(P)],
                 [0,1,0],    
                 [-np.sin(P),0,np.cos(P)]])
    RZ=np.array([[np.cos(K),-np.sin(K),0],
                 [np.sin(K),np.cos(K),0],
                 [0,0,1]])
    
    return RX.dot(RY.dot(RZ)).dot(np.array([[1,0,0],[0,-1,0],[0,0,-1]]))


def XYZ2Im(aPtWorld,aCam,aImSize):
    '''
    Function to project a point in world coordinate into an image

    :param aPtWorld: 3d point in world coordinates
    :param aCam: array describing a camera [position, rotation, focal]
    :param aImSize: Size of the image
    :return:  2d point in image coordinates
    '''
    # World to camera coordinate
    aPtCam=np.linalg.inv(aCam[1]).dot(aPtWorld-aCam[0])
    #print(aPtCam)
    # Test if point is behind camera (Z positive in Cam coordinates)
    if aPtCam[2]<0:
        return None
    #print("PtCam =", aPtCam)
    # Camera to 2D projected coordinate
    aPtProj=[aPtCam[0]/aPtCam[2],aPtCam[1]/aPtCam[2]]
    #print("PtProj =", aPtProj)
    # 2D projected to image coordinate
    aPtIm=[aImSize[0]/2,aImSize[1]/2]+np.array(aCam[2]).dot(aPtProj)
    if aPtIm[0]>=0 and aPtIm[1]>=0 and np.round(aPtIm[0])<aImSize[0] and np.round(aPtIm[1])<aImSize[1]:
        return aPtIm
    else:
        return None


def Raster2Array(raster_file, raster_band=1, nan_value=-9999):
    '''
    Function to convert a raster to an array

    :param raster_file: raster filepath
    :param raster_band: band to export
    :nan_value: value recognized as nan in raster. Default -9999
    :return:  array with the columns X,Y,value.
    '''
    
    # Open reference DEM
    myRaster=gdal.Open(raster_file)
    geot = myRaster.GetGeoTransform()
    Xsize=myRaster.RasterXSize
    Ysize=myRaster.RasterYSize
    data=myRaster.GetRasterBand(raster_band).ReadAsArray(0, 0, Xsize, Ysize)
    data=data.astype(float)
    data[data==nan_value] = np.nan

    # define extent and resoltuion from geotiff metadata
    extent = [geot[0], geot[0] + np.round(geot[1],3)*Xsize, geot[3], geot[3] + np.round(geot[5],3)*Ysize]

    
    # Create the X,Y coordinate meshgrid
    Xs = np.linspace(extent[0]+np.round(geot[1],3),extent[1], Xsize)
    Ys = np.linspace(extent[2]+ np.round(geot[5],3), extent[3], Ysize)
    XX, YY = np.meshgrid(Xs, Ys)
    
    XYZ = np.vstack((XX.flatten(),YY.flatten(),data.flatten())).T 
    return XYZ


def ProjectImage2DEM(dem_file, viewshed_file, image_file, output, aCam, dem_nan_value=-9999, PlyOffset=[0,0]):
    '''
    Function to project an image to a DEM

    :param dem_file: DEM raster filepath
    :param image_file: image filepath
    :param output: output point cloud filepath
    :param aCam: array describing a camera [position, rotation, focal]
    '''
    print('aCam=', aCam)
    tic = time.perf_counter()
    aDEM_as_list=Raster2Array(dem_file,nan_value=dem_nan_value)
    aViewshed_as_list=Raster2Array(viewshed_file)

    toc = time.perf_counter()
    print(f"DEM converted in {toc - tic:0.4f} seconds")

    # Load in image
    anImage=pyplot.imread(image_file)
    if len(anImage.shape)==2:
        anImage=anImage.T
        anImage = np.stack((anImage,anImage,anImage), axis=2)
    else:
            anImage = np.stack((anImage[:,:,0].T,anImage[:,:,1].T,anImage[:,:,2].T), axis=2)
            

    
    # Create output object
    aXYZinImage=np.zeros([anImage.shape[0],anImage.shape[1],3])*np.nan
    
    # Project all DEM points in the viewshed to the image
    tic = time.perf_counter()
    for p in range(aDEM_as_list.shape[0]):
        # Check if in viewshed
        if aViewshed_as_list[p][2]==255:
        # Project a the DEM pixel in image
            aWorldPt=np.array([aDEM_as_list[p][0],aDEM_as_list[p][1],aDEM_as_list[p][2]])
            aProjectedPoint=XYZ2Im(aWorldPt,aCam,anImage.shape)
            if not (aProjectedPoint is None):
                    xim=int(aProjectedPoint[0])
                    yim=int(aProjectedPoint[1])
                    aXYZinImage[xim,yim]=aWorldPt
                    
    toc = time.perf_counter()
    print(f"Position of each pixel computed in {toc - tic:0.4f} seconds")
            
    # export ply file
    vertex = np.array([],dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'),('red', 'u1'), ('green', 'u1'),('blue', 'u1')])
    # export each point X, Y, Z, R, G, B
    for i in range(anImage.shape[0]):
        for j in range(anImage.shape[1]):
            if not np.isnan(aXYZinImage[i,j][0]):
                aPoint=np.array([(aXYZinImage[i,j][0]-PlyOffset[0],aXYZinImage[i,j][1]-PlyOffset[1],aXYZinImage[i,j][2], anImage[i,j][0],anImage[i,j][1],anImage[i,j][2])],dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'),('red', 'u1'), ('green', 'u1'),('blue', 'u1')])
                vertex=np.concatenate((vertex,aPoint), axis=0)
    el = plyfile.PlyElement.describe(vertex, 'vertex')
    plyfile.PlyData([el], text=True).write(output)
    print('Total points : ', vertex.shape)
    print('PLY file extracted')
    
    return 0

def main():
    parser = OptionParser(usage="%prog [options]", version="%prog 0.1")

    # Define options
    parser.add_option(
        '-d', '--dem',
        default=None,
        help="read input DEM from FILE",
        metavar='FILE')
    
    parser.add_option(
        '-i', '--image',
        default=None,
        help="read input image from FILE",
        metavar='FILE')

    parser.add_option(
        '-o', '--output',
        default="result.ply",
        help="name of output file, default value is \"result.ply\"",
        metavar='FILE')

    parser.add_option(
        '-c', '--campos',
        type='float',
        dest='C',
        default=None,
        help="defines the camera position",
        metavar='N')
    parser.add_option(
        '-r', '--rot',
        type='float',
        dest='C',
        default=None,
        help="defines the camera rotation matrix",
        metavar='N')
    parser.add_option(
        '-f', '--foc',
        type='float',
        dest='C',
        default=None,
        help="defines the camera's focal lenght (in pixels)'",
        metavar='N')
    # Instruct optparse object
    (options, args) = parser.parse_args()
    
    aCam=[options.campos,options.rot,options.foc];
    
    ProjectImage2DEM(options.dem, options.image, options.output, aCam)

    return 0

# Input Finse
#
camera_file = './/FinseDemoData//CamFinseInit.inp'
point_file = './/FinseDemoData//GCPs_WebcamFinse_Centered3.inp'
Foc=1484
dem_file='.//FinseDemoData//time_lapse_finse_DSM_mid.tif'
viewshed_file='.//FinseDemoData//viewshed_mid.tif'
image_file='.//FinseDemoData//2019-05-24_12-00.jpg'
output='.//FinseDemoData//output.ply'


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
ProjectImage2DEM(dem_file, viewshed_file, image_file, output, aCam, dem_nan_value=1137.75, PlyOffset=[410000,6710000])




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



