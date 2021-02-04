# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:27:46 2021

@author: lucg
"""
import numpy as np
import gdal
from matplotlib import pyplot
import time
from osgeo import osr
import os

class ProjIm2dem():
    def __init__(self, dem_file, viewshed_file, image_file, cam_param, output_file, dem_nan_value=-9999):
        
        self.output_file = output_file
        self.cam_param = cam_param
        # Load in DEM
        self.dem_raster = gdal.Open(dem_file)
        self.geot = self.dem_raster.GetGeoTransform()
        self.Xsize = self.dem_raster.RasterXSize
        self.Ysize = self.dem_raster.RasterYSize
        self.dem_data = self.dem_raster.GetRasterBand(1).ReadAsArray(0, 0, self.Xsize, self.Ysize)
        self.dem_data[self.dem_data == dem_nan_value] = np.nan
        # define extent and resoluion from geotiff metadata
        self.extent = [self.geot[0], self.geot[0] + np.round(self.geot[1],3)*self.Xsize, self.geot[3], self.geot[3] + np.round(self.geot[5],3)*self.Ysize]
        self.Xs = np.linspace(self.extent[0]+np.round(self.geot[1],3), self.extent[1], self.Xsize)
        self.Ys = np.linspace(self.extent[2]+ np.round(self.geot[5],3), self.extent[3], self.Ysize)
        # TODO, generate viewshed here 
        # Compute DEM Z at Camera XY
        col = int((self.cam_param[0][0] - self.geot[0]) / self.Xsize)
        row = int((self.geot[3] - self.cam_param[0][1] ) / self.Ysize)
        ZDEMatCamera=self.dem_data[row][col]
        ZCameraOverDEM=self.cam_param[0][2]-ZDEMatCamera
        command='gdal_viewshed -ox ' + str(self.cam_param[0][0]) +' -oy ' + str(self.cam_param[0][1])  +' -oz ' + str(ZCameraOverDEM) + ' ' + str(dem_file) + ' ' + str(viewshed_file)
        print(command)
        os.system(command)
        self.viewshed_raster = gdal.Open(viewshed_file)
        self.viewshed_data = self.viewshed_raster.GetRasterBand(1).ReadAsArray(0, 0, self.Xsize, self.Ysize)
        
        self.image = pyplot.imread(image_file)
        if len(self.image.shape)==2:
            self.image_T = self.image.T
            self.image_T = np.stack((self.image_T, self.image_T, self.image_T), axis=2)
        else:
                self.image_T = np.stack((self.image[:,:,0].T, self.image[:,:,1].T, self.image[:,:,2].T), axis=2)
                
        # Create output object
        self.ortho = np.zeros((self.Ysize, self.Xsize, 3), dtype=np.uint8)
    
    def XYZ2Im(self, aPtWorld, aImSize):
        '''
        Function to project a point in world coordinate into an image
    
        :param aPtWorld: 3d point in world coordinates
        :param aCam: array describing a camera [position, rotation, focal]
        :param aImSize: Size of the image
        :return:  2d point in image coordinates
        '''
        # World to camera coordinate
        aPtCam = np.linalg.inv(self.cam_param[1]).dot(aPtWorld - self.cam_param[0])
        #print(aPtCam)
        # Test if point is behind camera (Z positive in Cam coordinates)
        if aPtCam[2] < 0:
            return None
        #print("PtCam =", aPtCam)
        # Camera to 2D projected coordinate
        aPtProj = [aPtCam[0]/aPtCam[2], aPtCam[1]/aPtCam[2]]
        #print("PtProj =", aPtProj)
        # 2D projected to image coordinate
        aPtIm = [aImSize[0]/2, aImSize[1]/2]+np.array(self.cam_param[2]).dot(aPtProj)
        if aPtIm[0]>=0 and aPtIm[1]>=0 and np.round(aPtIm[0])<aImSize[0] and np.round(aPtIm[1])<aImSize[1]:
            return aPtIm
        else:
            return None
    
    
    def ProjectImage2DEM(self, return_raster=True, epsg=None):
        '''
        Function to project an image to a DEM
        '''
        
        # Project all DEM points in the viewshed to the image
        tic = time.perf_counter()
        for x in range(0, self.Xsize):
            for y in range(0, self.Ysize):
                if self.viewshed_data[y][x]==255:
                    aWorldPt = np.array([self.Xs[x], self.Ys[y], self.dem_data[y][x]])
                    aProjectedPoint = self.XYZ2Im(aWorldPt, self.image_T.shape)
                    if not (aProjectedPoint is None):
                        self.ortho[y][x] = self.image_T[int(aProjectedPoint[0])][int(aProjectedPoint[1])]
                            
        toc = time.perf_counter()    
        
        print(f"Ortho computed in {toc - tic:0.4f} seconds")
        
        if return_raster:
            fileformat = "GTiff"
            driver = gdal.GetDriverByName(fileformat)
            metadata = driver.GetMetadata()
            if metadata.get(gdal.DCAP_CREATE) == "YES":
                print("Driver {} supports Create() method.".format(fileformat))

            if metadata.get(gdal.DCAP_CREATECOPY) == "YES":
                print("Driver {} supports CreateCopy() method.".format(fileformat))

            dst_ds = driver.Create(self.output_file, xsize=self.Xsize, ysize=self.Ysize,
                            bands=3, eType=gdal.GDT_Byte)

            dst_ds.SetGeoTransform(self.geot)
            srs = osr.SpatialReference()
            if epsg is None:
                srs.ImportFromWkt(self.dem_raster.GetProjectionRef())
            else:
                srs.ImportFromEPSG(epsg)
            dst_ds.SetProjection(srs.ExportToWkt())
            dst_ds.GetRasterBand(1).WriteArray(self.ortho[:, :, 0])
            dst_ds.GetRasterBand(2).WriteArray(self.ortho[:, :, 1])
            dst_ds.GetRasterBand(3).WriteArray(self.ortho[:, :, 2])
            # Once we're done, close properly the dataset
            dst_ds = None
            print('Ortho computed')
            
        return 0
