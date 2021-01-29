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
    
class ProjIm2dem():
    def __init__(self, dem_file, viewshed_file, image_file, cam, output_file, dem_nan_value=-9999):
               
        # Load in DEM
        self.aDEM.file=gdal.Open(self.dem_file)
        self.aDEM.geot = self.aDEM.file.GetGeoTransform()
        self.aDEM.Xsize=self.aDEM.RasterXSize
        self.aDEM.Ysize=self.aDEM.RasterYSize
        self.aDEM.data=self.aDEM.GetRasterBand(1).ReadAsArray(0, 0, self.aDEM.Xsize, self.aDEM.Ysize)
        self.aDEM.data[self.aDEM.data==self.dem_nan_value] = np.nan
        # define extent and resoluion from geotiff metadata
        self.aDEM.extent = [self.aDEM.geot[0], self.aDEM.geot[0] + np.round(self.aDEM.geot[1],3)*self.aDEM.Xsize, self.aDEM.geot[3], self.aDEM.geot[3] + np.round(self.aDEM.geot[5],3)*self.aDEM.Ysize]
        self.aDEM.Xs = np.linspace(self.aDEM.extent[0]+np.round(self.aDEM.geot[1],3),self.aDEM.extent[1], self.aDEM.Xsize)
        self.aDEM.Ys = np.linspace(self.aDEM.extent[2]+ np.round(self.aDEM.geot[5],3), self.aDEM.extent[3], self.aDEM.Ysize)
        # TODO, generate viewshed here 
        self.aViewshed=gdal.Open(viewshed_file)
        self.aViewshed.data=self.aViewshed.GetRasterBand(1).ReadAsArray(0, 0, self.aDEMXsize, self.aDEMYsize)
        
        self.anImage=pyplot.imread(self.image_file)
        if len(self.anImage.shape)==2:
            self.anImage=self.anImage.T
            self.anImage = np.stack((self.anImage,self.anImage,self.anImage), axis=2)
        else:
                self.anImage = np.stack((self.anImage[:,:,0].T,self.anImage[:,:,1].T,self.anImage[:,:,2].T), axis=2)
                
        class ProjectImage2DEM:
            def __init__(self):
                self.success=False
        self.ProjectImage2DEM = ProjectImage2DEM
    
    def XYZ2Im(self, aPtWorld,aImSize):
        '''
        Function to project a point in world coordinate into an image
    
        :param aPtWorld: 3d point in world coordinates
        :param aCam: array describing a camera [position, rotation, focal]
        :param aImSize: Size of the image
        :return:  2d point in image coordinates
        '''
        # World to camera coordinate
        aPtCam=np.linalg.inv(self.cam[1]).dot(aPtWorld-self.cam[0])
        #print(aPtCam)
        # Test if point is behind camera (Z positive in Cam coordinates)
        if aPtCam[2]<0:
            return None
        #print("PtCam =", aPtCam)
        # Camera to 2D projected coordinate
        aPtProj=[aPtCam[0]/aPtCam[2],aPtCam[1]/aPtCam[2]]
        #print("PtProj =", aPtProj)
        # 2D projected to image coordinate
        aPtIm=[aImSize[0]/2,aImSize[1]/2]+np.array(self.cam[2]).dot(aPtProj)
        if aPtIm[0]>=0 and aPtIm[1]>=0 and np.round(aPtIm[0])<aImSize[0] and np.round(aPtIm[1])<aImSize[1]:
            return aPtIm
        else:
            return None
    
    
    def ProjectImage2DEM(self):
        '''
        Function to project an image to a DEM
        '''

        # Create output object
        anOrtho=np.zeros((self.aDEM.Ysize, self.aDEM.Xsize,3), dtype=np.uint8)
        
        # Project all DEM points in the viewshed to the image
        tic = time.perf_counter()
        for x in range(0,self.aDEM.Xsize):
            for y in range(0,self.aDEM.Ysize):
                if self.aViewshed.data[y][x]==255:
                    aWorldPt=np.array([self.aDEM.Xs[x],self.aDEM.Ys[y],self.aDEM.data[y][x]])
                    aProjectedPoint=self.XYZ2Im(aWorldPt,self.anImage.shape)
                    if not (aProjectedPoint is None):
                        anOrtho[y][x]=self.anImage[int(aProjectedPoint[0])][int(aProjectedPoint[1])]
                            
        toc = time.perf_counter()    
        
        print(f"Ortho computed in {toc - tic:0.4f} seconds")
        
        fileformat = "GTiff"
        driver = gdal.GetDriverByName(fileformat)
        metadata = driver.GetMetadata()
        if metadata.get(gdal.DCAP_CREATE) == "YES":
            print("Driver {} supports Create() method.".format(fileformat))
        
        if metadata.get(gdal.DCAP_CREATECOPY) == "YES":
            print("Driver {} supports CreateCopy() method.".format(fileformat))
        
        dst_ds = driver.Create(self.output, xsize=self.aDEM.Xsize, ysize=self.aDEM.Ysize,
                        bands=3, eType=gdal.GDT_Byte)
    
        dst_ds.SetGeoTransform(self.aDEM.geot)
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
