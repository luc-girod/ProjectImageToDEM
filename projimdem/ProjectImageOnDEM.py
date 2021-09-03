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
import pandas as pd

class ProjIm2dem():
    def __init__(self, dem_file, viewshed_file, image_file, cam_param, output_file, dem_nan_value=-9999):
        '''
        dem_file: geotif of the DEM
        viewshed_file: name for the viewshade file (will be computed by gdal based on camera parameters)
        image_file:
        cam_param:
        output_file:
        dem_nan_value:
        
        TODO: 
        - review the viewshed file. add option to not have to compute it each time
        '''
        
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
        col = int((self.cam_param[0][0] - self.geot[0]) / self.geot[1])
        row = int((self.cam_param[0][1] - self.geot[3]) / self.geot[5])
        ZDEMatCamera=self.dem_data[row][col]
        print("ZDEMatCamera", ZDEMatCamera)        
        print("ZCamera", self.cam_param[0][2])
        ZCameraOverDEM=self.cam_param[0][2]-ZDEMatCamera
        
        # Compute viewshade file using GDAL and the new camear prosition
        command='gdal_viewshed -ox ' + str(self.cam_param[0][0]) +' -oy ' + str(self.cam_param[0][1])  +' -oz ' + str(ZCameraOverDEM) + ' ' + str(dem_file) + ' ' + str(viewshed_file)
        print(command)
        os.system(command)
        self.viewshed_raster = gdal.Open(viewshed_file)
        self.viewshed_data = self.viewshed_raster.GetRasterBand(1).ReadAsArray(0, 0, self.Xsize, self.Ysize)
        self.viewshed_raster = None
        self.image = pyplot.imread(image_file)
        if len(self.image.shape)==2:
            # in case of  black and white image
            self.image_T = self.image.T
            self.image_T = np.stack((self.image_T, self.image_T, self.image_T), axis=2)
        else:
            # generak RGB case
            self.image_T = np.stack((self.image[:,:,0].T, self.image[:,:,1].T, self.image[:,:,2].T), axis=2)
                
        # Create output object
        self.ortho = np.zeros((self.Ysize, self.Xsize, 3), dtype=np.uint8)
        
        # Compute remove radial and tangential distortion based on camera parameters
        self.X_undistort, self.Y_undistort, self.image_undistort = self.img_correct_distortion()
        self.pt_proj = pd.DataFrame()
    
    def XYZ2Im_all(self, pts_world):
        imsize = self.image.shape
        
        pts_cam = pts_world.copy()
        
        # World to camera coordinate
        XYZ_world = pts_cam[['X_world','Y_world','Z_world']].values
        pts_cam[['X_cam','Y_cam','Z_cam']] = pd.DataFrame((self.cam_param[1].dot(np.subtract(XYZ_world, self.cam_param[0]).T)).T)
        
        XYZ_world = None
        
        # Removeing points behind the camera
        pts_cam.Z_cam.loc[pts_cam.Z_cam > 0] = np.nan
        pts_cam = pts_cam.dropna(axis=0)
        
        pts_cam['X_proj'] = pts_cam.X_cam.values/pts_cam.Z_cam.values
        pts_cam['Y_proj'] = pts_cam.Y_cam.values/pts_cam.Z_cam.values
        pts_cam['X_img'] = imsize[1]/2 - pts_cam['X_proj'] * self.cam_param[2]
        pts_cam['Y_img'] = imsize[0]/2 + pts_cam['Y_proj'] * self.cam_param[2]
        print('\n step 1: \n', pts_cam.head())
        
        pts_cam.X_img.loc[(pts_cam.X_img<0) | (pts_cam.X_img>imsize[1])] = np.nan
        pts_cam.Y_img.loc[(pts_cam.Y_img<0) | (pts_cam.Y_img>imsize[0])] = np.nan
        pts_cam = pts_cam.dropna(axis=0)
        print('\n step 2: \n', pts_cam.head())
        print (self.X_undistort.min())
        print (self.Y_undistort.min())
        pts_cam['X_distort'] =  self.X_undistort[pts_cam.Y_img.astype(int), pts_cam.X_img.astype(int)] - self.X_undistort.min()
        pts_cam['Y_distort'] =  self.Y_undistort[pts_cam.Y_img.astype(int), pts_cam.X_img.astype(int)] - self.Y_undistort.min()
        
        print('\n step 3: \n', pts_cam.head())
        
        
        pts_cam.X_distort.loc[(pts_cam.X_distort<0) | (pts_cam.X_distort>self.image_undistort.shape[1])] = np.nan
        pts_cam.Y_distort.loc[(pts_cam.Y_distort<0) | (pts_cam.Y_distort>self.image_undistort.shape[0])] = np.nan
        pts_cam = pts_cam.dropna(axis=0)
        
        print('\n step 4: \n', pts_cam.head())
        
        return pts_cam
        
        
        #pts_proj = [pts_cam.X/pts_cam.Z, pts_cam.Y/pts_cam.Z]
        #pts_img = np.add(np.array([imsize[0]/2, imsize[1]/2]) + pts_proj* self.cam_param[2])
        
        #pts_img.X.loc[pts_img.X<0 and pts_img.X>imsize[1]] = np.nan
        #pts_img.Y.loc[pts_img.Y<0 and pts_img.Y>imsize[0]] = np.nan
        #pts_img = pts_img.dropna()
        
        #pts_distort = np.array([self.X_undistort[pts_img.Y.astype(int),pts_img.X.astype(int)],
                                #self.Y_undistort[pts_img.Y.astype(int), pts_img.Y.astype(int)]])
    
        #pts_distort.X.loc[pts_distort.X<0 and pts_distort.X>imsize[1]] = np.nan
        #pts_distort.Y.loc[pts_distort.Y<0 and pts_distort.Y>imsize[0]] = np.nan
        #pts_distort = pts_distort.dropna()
        
        #return pts_distort
        
        
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

        # Test if point is behind camera (Z positive in Cam coordinates)
        if aPtCam[2] < 0:
            return None
        #print("PtCam =", aPtCam)
        # Camera to 2D projected coordinate
        aPtProj = [aPtCam[0]/aPtCam[2], aPtCam[1]/aPtCam[2]]
        #print("PtProj =", aPtProj)
        # 2D projected to image coordinate
        aPtIm = np.array([aImSize[0]/2, aImSize[1]/2])+np.array(self.cam_param[2]).dot(aPtProj)
        
        
        # check that points falls within the pixel domain of the image
        if aPtIm[0]>=0 and aPtIm[1]>=0 and np.round(aPtIm[0])<aImSize[0] and np.round(aPtIm[1])<aImSize[1]:
            #print("aPtIm = ", aPtIm)
            
            aPtIm_distort = np.array([self.X_undistort[aPtIm[1].astype(int),aPtIm[0].astype(int)], self.Y_undistort[aPtIm[1].astype(int), aPtIm[0].astype(int)]])
            #print("aPtIm_distort = ",aPtIm_distort)
            aPtIm_distort = (aPtIm_distort - np.array([self.X_undistort.min(), self.Y_undistort.min()])).astype(int)
            
            if aPtIm_distort[0]>=0 and aPtIm_distort[1]>=0 and np.round(aPtIm_distort[0])<aImSize[0] and np.round(aPtIm_distort[1])<aImSize[1]:
                #print("aptim_distro = ", aPtIm_distort)
                
                return aPtIm_distort
            else:
                return None
        else:
            return None
    
    def img_correct_distortion(self, return_image=True):
        """
        Function to correction radial and tangential distortion of an image
        :param img: original image
        :return 
        """
        img = self.image
        [DCx, DCy] = self.cam_param[3]
        [K1,K2,K3,K4,K5,P1,P2,P3,P4,P5,P6,P7] = self.cam_param[4]

        Xs_dis, Ys_dis = np.meshgrid(np.arange(-img.shape[1]/2,img.shape[1]/2), 
                            np.arange(-img.shape[0]/2,img.shape[0]/2)) 
        # Distortion model
        R = np.sqrt(pow(Xs_dis-DCx,2)+pow(Ys_dis-DCy,2))

        x_im_nodist = Xs_dis + (Xs_dis - DCx) * (K1*pow(R,2) + K2*pow(R,4)+ K3*pow(R,6)+ K4*pow(R,8)+ K5*pow(R,10)) + (P1 * (pow(R,2) + 2*pow((Xs_dis - DCx),2)) + 2*P2*(Xs_dis - DCx))*(Ys_dis - DCy)*(1+ P3 *pow(R,2) + P4*pow(R,4)+ P5 *pow(R,6) + P6*pow(R,8) + P7*pow(R,10))
        y_im_nodist = Ys_dis + (Ys_dis - DCy) * (K1*pow(R,2) + K2*pow(R,4)+ K3*pow(R,6)+ K4*pow(R,8)+ K5*pow(R,10)) + (P1 * (pow(R,2) + 2*pow((Ys_dis - DCy),2)) + 2*P2*(Ys_dis - DCy))*(Xs_dis - DCx)*(1+ P3 *pow(R,2) + P4*pow(R,4)+ P5 *pow(R,6) + P6*pow(R,8) + P7*pow(R,10))
        
        if return_image:

            # Project distorted image
            img_cor = np.empty([y_im_nodist.astype(int).max() - y_im_nodist.astype(int).min()+1, x_im_nodist.astype(int).max() - x_im_nodist.astype(int).min()+1,3])
            img_cor[(y_im_nodist.astype(int)-y_im_nodist.astype(int).min()), (x_im_nodist.astype(int)-x_im_nodist.astype(int).min()),:] = img
            img_cor = img_cor.astype(int)

            img = None
            return x_im_nodist, y_im_nodist, img_cor
        else:
            return x_im_nodist, y_im_nodist
    
    def ProjectImage2DEM(self, return_raster=True, epsg=None):
        '''
        Function to project an image to a DEM
        '''
        
        # Project all DEM points in the viewshed to the image
        tic = time.perf_counter()
        
        x_mesh, y_mesh = np.meshgrid(self.Xs, self.Ys)
        pt_world = pd.DataFrame()
        pt_world['X_world'] = x_mesh[self.viewshed_data==255].flatten()
        pt_world['Y_world'] = y_mesh[self.viewshed_data==255].flatten()
        pt_world['Z_world'] = self.dem_data[self.viewshed_data==255].flatten()
        pt_world['ind_xdem'] = np.argwhere(self.viewshed_data==255)[:,1]
        pt_world['ind_ydem'] = np.argwhere(self.viewshed_data==255)[:,0]
                
        self.pt_proj = self.XYZ2Im_all(pt_world)
        
        print('pt_proj:', self.pt_proj.head())
        for i, row in self.pt_proj.iterrows():
            self.ortho[row.ind_ydem.astype(int), row.ind_xdem.astype(int),:] = self.image_undistort[ row.Y_distort.astype(int), row.X_distort.astype(int), :]
        
        #for x in range(0, self.Xsize):
        #    for y in range(0, self.Ysize):
        #        if self.viewshed_data[y][x]==255:
                    
        #            aWorldPt = np.array([self.Xs[x], self.Ys[y], self.dem_data[y][x]])
                    
        #            aProjectedPoint = self.XYZ2Im(aWorldPt, self.image_T.shape)
        #            if not (aProjectedPoint is None):
        #                self.ortho[y][x] = self.image_T[int(aProjectedPoint[0])][int(aProjectedPoint[1])]
                            
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


    
#=================================================================================
#--------------------------  Practical Functions -------------------------
#=================================================================================


def img_correct_distortion(img, cam_param):
    """
    Function to correction radial and tangential distortion of an image
    :param img: original image
    :param cam_param: Object containing camera calibration coefficients. extracted from the resection.
    :return 
    """
    
    Xs_dis, Ys_dis = np.meshgrid(np.arange(-img.shape[1]/2,img.shape[1]/2), 
                        np.arange(-img.shape[0]/2,img.shape[0]/2)) 
    [DCx, DCy] = cam_param[3]
    R = np.sqrt(pow(Xs_dis-DCx,2)+pow(Ys_dis-DCy,2))
    
    [K1,K2,K3,K4,K5,P1,P2,P3,P4,P5,P6,P7] = cam_param[4]
    
    x_im_nodist = Xs_dis + (Xs_dis - DCx) * (K1*pow(R,2) + K2*pow(R,4)+ K3*pow(R,6)+ K4*pow(R,8)+ K5*pow(R,10)) + (P1 * (pow(R,2) + 2*pow((Xs_dis - DCx),2)) + 2*P2*(Xs_dis - DCx))*(Ys_dis - DCy)*(1+ P3 *pow(R,2) + P4*pow(R,4)+ P5 *pow(R,6) + P6*pow(R,8) + P7*pow(R,10))
    y_im_nodist = Ys_dis + (Ys_dis- DCy)*(K1*pow(R,2) + K2*pow(R,4)+ K3*pow(R,6)+ K4*pow(R,8)+ K5*pow(R,10)) + (P1 * (pow(R,2) + 2*pow((Ys_dis - DCy),2)) + 2*P2*(Ys_dis - DCy))*(Xs_dis - DCx)*(1+ P3 *pow(R,2) + P4*pow(R,4)+ P5 *pow(R,6) + P6*pow(R,8) + P7*pow(R,10))
    
    
    img_cor = np.empty([y_im_nodist.astype(int).max() - y_im_nodist.astype(int).min()+1, x_im_nodist.astype(int).max() - x_im_nodist.astype(int).min()+1,3])
    img_cor[(y_im_nodist.astype(int)-y_im_nodist.astype(int).min()), (x_im_nodist.astype(int)-x_im_nodist.astype(int).min()),:] = img
    img_cor = img_cor.astype(int)
    
    
    return x_im_nodist, y_im_nodist, img_cor
    