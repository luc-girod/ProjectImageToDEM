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
import time
from scipy.optimize import leastsq
from osgeo import osr
    
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


def ProjectImage2DEM(dem_file, viewshed_file, image_file, output, aCam, dem_nan_value=-9999):
    '''
    Function to project an image to a DEM

    :param dem_file: DEM raster filepath
    :param viewshed_file: viewshed raster filepath, must be same geometry as DEM
    :param image_file: image filepath
    :param output: output point cloud filepath
    :param aCam: array describing a camera [position, rotation, focal]
    '''
    print('aCam=', aCam)
    
    # Load in DEM
    aDEM=gdal.Open(dem_file)
    geot = aDEM.GetGeoTransform()
    Xsize=aDEM.RasterXSize
    Ysize=aDEM.RasterYSize
    aDEMdata=aDEM.GetRasterBand(1).ReadAsArray(0, 0, Xsize, Ysize)
    aDEMdata[aDEMdata==dem_nan_value] = np.nan
    # define extent and resoluion from geotiff metadata
    extent = [geot[0], geot[0] + np.round(geot[1],3)*Xsize, geot[3], geot[3] + np.round(geot[5],3)*Ysize]
    Xs = np.linspace(extent[0]+np.round(geot[1],3),extent[1], Xsize)
    Ys = np.linspace(extent[2]+ np.round(geot[5],3), extent[3], Ysize)
    
    # Load in viewshed
    # TODO, generate it here
    aViewshed=gdal.Open(viewshed_file)
    aViewshedData=aViewshed.GetRasterBand(1).ReadAsArray(0, 0, Xsize, Ysize)

    # Load in image
    anImage=pyplot.imread(image_file)
    if len(anImage.shape)==2:
        anImage=anImage.T
        anImage = np.stack((anImage,anImage,anImage), axis=2)
    else:
            anImage = np.stack((anImage[:,:,0].T,anImage[:,:,1].T,anImage[:,:,2].T), axis=2)
            

    
    # Create output object
    anOrtho=np.zeros((Ysize, Xsize,3), dtype=np.uint8)
    
    # Project all DEM points in the viewshed to the image
    tic = time.perf_counter()
    for x in range(0,Xsize):
        for y in range(0,Ysize):
            if aViewshedData[y][x]==255:
                aWorldPt=np.array([Xs[x],Ys[y],aDEMdata[y][x]])
                aProjectedPoint=XYZ2Im(aWorldPt,aCam,anImage.shape)
                if not (aProjectedPoint is None):
                    anOrtho[y][x]=anImage[int(aProjectedPoint[0])][int(aProjectedPoint[1])]
                        
    toc = time.perf_counter()    
    
    print(f"Ortho computed in {toc - tic:0.4f} seconds")
    
    fileformat = "GTiff"
    driver = gdal.GetDriverByName(fileformat)
    metadata = driver.GetMetadata()
    if metadata.get(gdal.DCAP_CREATE) == "YES":
        print("Driver {} supports Create() method.".format(fileformat))
    
    if metadata.get(gdal.DCAP_CREATECOPY) == "YES":
        print("Driver {} supports CreateCopy() method.".format(fileformat))
    
    dst_ds = driver.Create(output, xsize=Xsize, ysize=Ysize,
                    bands=3, eType=gdal.GDT_Byte)

    dst_ds.SetGeoTransform(geot)
    srs = osr.SpatialReference()
    srs.SetUTM(32, 1)
    srs.SetWellKnownGeogCS("WGS84")
    dst_ds.SetProjection(srs.ExportToWkt())
    dst_ds.GetRasterBand(1).WriteArray(anOrtho[:,:,0])
    dst_ds.GetRasterBand(2).WriteArray(anOrtho[:,:,1])
    dst_ds.GetRasterBand(3).WriteArray(anOrtho[:,:,2])
    # Once we're done, close properly the dataset
    dst_ds = None
    print('Ortho file extracted')
    
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



