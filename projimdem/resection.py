#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Adaptation from:

> Author:  Jeffrey T. Walton, Paul Smith's College, New York
>
>   Single-photo resection - calculates the camera orientation and location
>       given camera calibration parameters, control point photo and world
>        coordinates and initial guesses for camera exterior orientation.
>
>   based on MATLAB code from:
>   Introduction to Modern Photogrammetry by Mikhail, Bethel, McGlone
>   John Wiley & Sons, Inc. 2001
>
> https://github.com/jeffwalton/photogrammetry-resection
'''


import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize 
import json
import pandas as pd
from types import SimpleNamespace
from math import sin, cos
from matplotlib import pyplot

class resection():
    def __init__(self, camera_file, GCP_file, image_file, delimiter_GCP=' ', x_offset=None, y_offset=None, z_offset=None, free_param=['omega', 'phi', 'kappa'], param_bounds=([-3.15, -3.15, -3.15], [3.15,3.15,3.15])):
        
        #Load camera parameters
        with open(camera_file, 'r') as myfile:
            self.cam = json.loads(myfile.read(), object_hook=lambda d: SimpleNamespace(**d))
        
        # Load GCP coordinates
        # headers must be: name x_img y_img x_world y_world z_world
        self.GCPs = pd.read_csv(GCP_file, delimiter=delimiter_GCP)
        
        if x_offset is None:
            self.x_offset = self.GCPs.x_world.mean()
        else:
            self.x_offset = x_offset
            
        if y_offset is None:
            self.y_offset = self.GCPs.y_world.mean()
        else:
            self.y_offset = y_offset
        
        if z_offset is None:
            self.z_offset = self.GCPs.z_world.mean()
        else:
            self.z_offset = z_offset
        
        self.free_param = free_param
        
        self.x0_dict = {}
        for param in free_param:
            if param in ['omega', 'kappa', 'phi', 'X_ini', 'Y_ini', 'Z_ini']:
                p = self.cam.eop.__getattribute__(param)
                self.x0_dict[param] = p
                if param =='X_ini':
                    self.x0_dict[param] = p - self.x_offset
                elif param =='Y_ini':
                    self.x0_dict[param] = p - self.y_offset
                elif param =='Z_ini':
                    self.x0_dict[param] = p - self.z_offset
            if param in ['Foc']:
                p = self.cam.iop.__getattribute__(param)
                self.x0_dict[param] = p
        self.x0 = list(self.x0_dict.values())
        self.param_bounds = (param_bounds)
        
        class new_cam:
            def __init__(self):
                self.center = None
                self.rotation = None
                self.proj_param = None
                self.omega = None
                self.kappa = None
                self.phi = None
        self.new_cam = new_cam
        
        
        self.GCPs['x_world_offset'] = self.GCPs.x_world - self.x_offset
        self.GCPs['y_world_offset'] = self.GCPs.y_world - self.y_offset
        self.GCPs['z_world_offset'] = self.GCPs.z_world - self.z_offset
        
        res_ini = self.collinearity_func(self.x0)
        
        idx = np.arange(0,res_ini.__len__(),2)
        self.GCPs['residual_x_lstsq'] = np.nan
        self.GCPs['residual_y_lstsq'] = np.nan
        self.GCPs['residual_x_ini'] = np.nan
        self.GCPs['residual_y_ini'] = np.nan
        self.GCPs['residual_x_ini'].loc[self.GCPs.lstsq_IO.astype(bool)] = res_ini[idx]
        self.GCPs['residual_y_ini'].loc[self.GCPs.lstsq_IO.astype(bool)] = res_ini[idx+1]
        self.image = pyplot.imread(image_file)
        
        
        
    
    def RotMatrixFromAngles(self, omega, phi, kappa):
        '''
        Rotation matrix from angles following Micmac convention
        '''
        RX = np.array([[1,0,0],
                 [0, cos(omega), -sin(omega)],
                 [0, sin(omega), cos(omega)]])    
        RY = np.array([[cos(phi), 0, sin(phi)],
                 [0,1,0],    
                 [-sin(phi), 0, cos(phi)]])
        RZ = np.array([[cos(kappa),-sin(kappa),0],
                 [sin(kappa), cos(kappa),0],
                 [0,0,1]])
        M = RX.dot(RY.dot(RZ)).dot(np.array([[1,0,0],[0,-1,0],[0,0,-1]]))

        #Mom = np.matrix([[1, 0, 0], [0, cos(omega), sin(omega)], [0, -sin(omega), cos(omega)]])
        #Mph = np.matrix([[cos(phi), 0, -sin(phi)], [0, 1, 0], [sin(phi), 0, cos(phi)]])
        #Mkp = np.matrix([[cos(kappa), sin(kappa), 0], [-sin(kappa), cos(kappa), 0], [0, 0, 1]])
        #M = Mkp * Mph * Mom    
        return M
    
    def collinearity_func(self, indep_vars):
        """
        collinearity function calculates a sum of the squared residuals of the
            collinearity equations for all of the control points
        This function is passed to scipy.optimize.minimize()

        Inputs:
            indep_vars (passed) are the exterior orientation parameters of the camera

        Returns:
            sum of squared residuals of collinearity eqns
        """
        
        omega = self.cam.eop.omega
        phi = self.cam.eop.phi
        kappa = self.cam.eop.kappa
        XL = self.cam.eop.X_ini - self.x_offset
        YL = self.cam.eop.Y_ini - self.y_offset
        ZL = self.cam.eop.Z_ini - self.z_offset
        Foc = self.cam.iop.Foc
        
        if 'omega' in self.x0_dict.keys():
            omega = indep_vars[list(self.x0_dict.keys()).index('omega')]
        if 'phi' in self.x0_dict.keys():
            phi = indep_vars[list(self.x0_dict.keys()).index('phi')]
        if 'omega' in self.x0_dict.keys():
            kappa = indep_vars[list(self.x0_dict.keys()).index('kappa')]
        if 'X_ini' in self.x0_dict.keys():
            XL = indep_vars[list(self.x0_dict.keys()).index('X_ini')]
        if 'Y_ini' in self.x0_dict.keys():
            YL = indep_vars[list(self.x0_dict.keys()).index('Y_ini')]
        if 'Z_ini' in self.x0_dict.keys():
            ZL = indep_vars[list(self.x0_dict.keys()).index('Z_ini')]
        if 'Foc' in self.x0_dict.keys():
            Foc = indep_vars[list(self.x0_dict.keys()).index('Foc')]

        Mom = np.matrix([[1, 0, 0], [0, cos(omega), sin(omega)], [0, -sin(omega), cos(omega)]])
        Mph = np.matrix([[cos(phi), 0, -sin(phi)], [0, 1, 0], [sin(phi), 0, cos(phi)]])
        Mkp = np.matrix([[cos(kappa), sin(kappa), 0], [-sin(kappa), cos(kappa), 0], [0, 0, 1]])
        M = Mkp * Mph * Mom

        tmp = self.GCPs.loc[self.GCPs.lstsq_IO.astype(bool)].reset_index(drop=True)
        
        F = np.zeros(2*tmp.shape[0])
        
        
        for i, row in tmp.iterrows():
            uvw = M * np.matrix([[row.x_world_offset - XL], [row.y_world_offset - YL], [row.z_world_offset - ZL]])
            resx = row.x_img - self.cam.iop.x0 + Foc * uvw[0,0] / uvw[2,0]
            resy = row.y_img - self.cam.iop.y0 + Foc * uvw[1,0] / uvw[2,0]
            
            F[2*i], F[2*i+1] = resx, resy

        return F
    
    def proj_GCPs2img(self):
            self.GCPs['x_img_repoj']=self.GCPs['x_img']-self.GCPs['residual_x_lstsq']+self.image.shape[1]/2
            self.GCPs['y_img_repoj']=-(self.GCPs['y_img']-self.GCPs['residual_y_lstsq']-self.image.shape[0]/2)
            fig, ax = plt.subplots()
            ax.imshow(self.image)
            pyplot.scatter(self.GCPs['x_img']+self.image.shape[1]/2,-(self.GCPs['y_img']-self.image.shape[0]/2), label='Original positions')
            for i, txt in enumerate(self.GCPs['name']):
                ax.annotate(self.GCPs['name'][i], (self.GCPs['x_img'][i]+self.image.shape[1]/2,-(self.GCPs['y_img'][i]-self.image.shape[0]/2)))
                
            
            pyplot.scatter(self.GCPs['x_img_repoj'],self.GCPs['y_img_repoj'], label='Reprojected positions')
            for i, txt in enumerate(self.GCPs['name']):
                ax.annotate(self.GCPs['name'][i], (self.GCPs['x_img_repoj'][i],self.GCPs['y_img_repoj'][i]))
            ax.legend()                


    def proj_XYZ2img(self, XYZ, omega, phi, kappa):
        
        Mom = np.matrix([[1, 0, 0], [0, cos(omega), sin(omega)], [0, -sin(omega), cos(omega)]])
        Mph = np.matrix([[cos(phi), 0, -sin(phi)], [0, 1, 0], [sin(phi), 0, cos(phi)]])
        Mkp = np.matrix([[cos(kappa), sin(kappa), 0], [-sin(kappa), cos(kappa), 0], [0, 0, 1]])
        M = Mkp * Mph * Mom
        
        
    
    def estimate_cam(self, method='dogbox', loss='cauchy', verbose=1): 
        # see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
        
        
        res = optimize.least_squares(self.collinearity_func, self.x0, 
                                     loss=loss, method=method, verbose=verbose,
                                    bounds=self.param_bounds)
        self.new_cam.RMSE = np.sqrt(np.sum(res.fun**2)/res.fun.__len__())
        self.new_cam.lstsq_results = res
        
        
        if 'omega' in self.x0_dict.keys():
            self.new_cam.omega = res.x[list(self.x0_dict.keys()).index('omega')]
        else:
            self.new_cam.omega = self.cam.eop.omega
            
        if 'phi' in self.x0_dict.keys():
            self.new_cam.phi = res.x[list(self.x0_dict.keys()).index('phi')]
        else:
            self.new_cam.phi = self.cam.eop.phi
            
        if 'kappa' in self.x0_dict.keys():
            self.new_cam.kappa = res.x[list(self.x0_dict.keys()).index('kappa')]
        else:
            self.new_cam.kappa = self.cam.eop.kappa
            
        if 'X_ini' in self.x0_dict.keys():
            self.new_cam.X_ini = res.x[list(self.x0_dict.keys()).index('X_ini')]
        else:
            self.new_cam.X_ini = self.cam.eop.X_ini - self.x_offset
            
        if 'Z_ini' in self.x0_dict.keys():
            self.new_cam.Z_ini = res.x[list(self.x0_dict.keys()).index('Z_ini')]
        else:
            self.new_cam.Z_ini = self.cam.eop.Z_ini - self.z_offset
            
        if 'Y_ini' in self.x0_dict.keys():
            self.new_cam.Y_ini = res.x[list(self.x0_dict.keys()).index('Y_ini')]
        else:
            self.new_cam.Y_ini = self.cam.eop.Y_ini - self.y_offset
            
        if 'Foc' in self.x0_dict.keys():
            self.new_cam.Foc = res.x[list(self.x0_dict.keys()).index('Foc')]
        else:
            self.new_cam.Foc = self.cam.iop.Foc
            
        self.new_cam.center = [self.new_cam.X_ini + self.x_offset, self.new_cam.Y_ini + self.y_offset, self.new_cam.Z_ini + self.z_offset]
        self.new_cam.rotation = self.RotMatrixFromAngles(self.new_cam.omega, self.new_cam.phi, self.new_cam.kappa)
        self.new_cam.proj_param = [self.new_cam.center, self.new_cam.rotation, self.new_cam.Foc]
        
        idx = np.arange(0,res.fun.__len__(),2)
        self.GCPs['residual_x_lstsq'] = np.nan
        self.GCPs['residual_y_lstsq'] = np.nan
        self.GCPs['residual_x_lstsq'].loc[self.GCPs.lstsq_IO.astype(bool)] = res.fun[idx]
        self.GCPs['residual_y_lstsq'].loc[self.GCPs.lstsq_IO.astype(bool)] = res.fun[idx+1]
        
        RMSE_ini = np.sqrt(np.sum(self.GCPs.residual_x_ini**2 + self.GCPs.residual_y_ini**2)/res.fun.__len__())
        print('RMSE initial = ', RMSE_ini)
        print('RMSE lstsq = ', self.new_cam.RMSE)
        return self.new_cam.proj_param
    
    def print_residuals(self):
        
        print(self.GCPs[['name','residual_x_ini', 'residual_y_ini', 'residual_x_lstsq',
       'residual_y_lstsq']].to_string())

    def plot_residuals(self):
        fig, ax = plt.subplots(3,2,sharex=True, sharey=True)
        sc1 = ax[0,0].scatter(self.GCPs.x_img, self.GCPs.y_img, 
                              c=self.GCPs.residual_x_ini,
                              cmap=plt.cm.RdBu,
                              vmin= -self.GCPs[['residual_x_ini', 'residual_y_ini']].abs().max().max(),
                              vmax= self.GCPs[['residual_x_ini', 'residual_y_ini']].abs().max().max())
        ax[0,0].set_title('Initial Residuals in X')
        plt.colorbar(sc1, ax=ax[0,0])
        sc2 = ax[0,1].scatter(self.GCPs.x_img, self.GCPs.y_img,
                              c=self.GCPs.residual_y_ini,
                              cmap=plt.cm.RdBu,
                              vmin= -self.GCPs[['residual_x_ini', 'residual_y_ini']].abs().max().max(),
                              vmax= self.GCPs[['residual_x_ini', 'residual_y_ini']].abs().max().max())
        ax[0,1].set_title('Initial Residuals in Y')
        plt.colorbar(sc2, ax=ax[0,1])
        sc3 = ax[1,0].scatter(self.GCPs.x_img, self.GCPs.y_img,
                              c=self.GCPs.residual_x_lstsq,
                              cmap=plt.cm.RdBu,
                              vmin= -self.GCPs[['residual_x_lstsq', 'residual_y_lstsq']].abs().max().max(),
                              vmax= self.GCPs[['residual_x_lstsq', 'residual_y_lstsq']].abs().max().max())
        ax[1,0].set_title('LstSq Residuals in X')
        plt.colorbar(sc3, ax=ax[1,0])
        sc4 = ax[1,1].scatter(self.GCPs.x_img, self.GCPs.y_img,
                              c=self.GCPs.residual_y_lstsq,
                              cmap=plt.cm.RdBu,
                              vmin= -self.GCPs[['residual_x_lstsq', 'residual_y_lstsq']].abs().max().max(),
                              vmax= self.GCPs[['residual_x_lstsq', 'residual_y_lstsq']].abs().max().max())
        ax[1,1].set_title('LstSq Residuals in Y')
        plt.colorbar(sc4, ax=ax[1,1])

        sc5 = ax[2,0].scatter(self.GCPs.x_img, 
                              self.GCPs.y_img, 
                              c=self.GCPs.residual_x_ini - self.GCPs.residual_x_lstsq)
        ax[2,0].set_title('Diff in X')
        plt.colorbar(sc5, ax=ax[2,0])
        sc6 = ax[2,1].scatter(self.GCPs.x_img, 
                              self.GCPs.y_img, 
                              c=self.GCPs.residual_y_ini - self.GCPs.residual_y_lstsq)
        ax[2,1].set_title('Diff in Y')
        plt.colorbar(sc6, ax=ax[2,1])