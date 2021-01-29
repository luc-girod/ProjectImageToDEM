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

class resection():
    def __init__(self, camera_file, GCP_file, delimiter_GCP=' ', x_offset=None, y_offset=None, z_offset=None):
        
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
        
        class estimate:
            def __init__(self):
                self.center = None
                self.rotation = None
                self.new_cam = None
        self.estimate = estimate
        
        self.x0 = [self.cam.eop.omega,
                  self.cam.eop.phi,
                  self.cam.eop.kappa,
                  self.cam.eop.X_ini - self.x_offset,
                  self.cam.eop.Y_ini - self.y_offset,
                  self.cam.eop.Z_ini - self.z_offset]
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

        omega = indep_vars[0]
        phi = indep_vars[1]
        kappa = indep_vars[2]
        XL = indep_vars[3]
        YL = indep_vars[4]
        ZL = indep_vars[5]

        Mom = np.matrix([[1, 0, 0], [0, cos(omega), sin(omega)], [0, -sin(omega), cos(omega)]])
        Mph = np.matrix([[cos(phi), 0, -sin(phi)], [0, 1, 0], [sin(phi), 0, cos(phi)]])
        Mkp = np.matrix([[cos(kappa), sin(kappa), 0], [-sin(kappa), cos(kappa), 0], [0, 0, 1]])
        M = Mkp * Mph * Mom

        tmp = self.GCPs.loc[self.GCPs.lstsq_IO.astype(bool)].reset_index(drop=True)
        
        F = np.zeros(2*tmp.shape[0])
        
        
        for i, row in tmp.iterrows():
            uvw = M * np.matrix([[row.x_world_offset - XL], [row.y_world_offset - YL], [row.z_world_offset - ZL]])
            resx = row.x_img - self.cam.iop.x0 + self.cam.iop.Foc * uvw[0,0] / uvw[2,0]
            resy = row.y_img - self.cam.iop.y0 + self.cam.iop.Foc * uvw[1,0] / uvw[2,0]
            
            F[2*i], F[2*i+1] = resx, resy

        return F
    
    
    def estimate_cam(self, method='dogbox', loss='cauchy', verbose=1): 
        # see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
        res = optimize.least_squares(self.collinearity_func, self.x0, loss=loss, method=method, verbose=verbose)
        self.estimate.RMSE = np.sqrt(np.sum(res.fun**2)/res.fun.__len__())
        self.estimate.lstsq_results = res
        self.estimate.center = [res.x[3] + self.x_offset, res.x[4] + self.y_offset, res.x[5] + self.z_offset]
        self.estimate.rotation = self.RotMatrixFromAngles(res.x[0], res.x[1], res.x[2])
        self.estimate.new_cam = [self.estimate.center, self.estimate.rotation, self.cam.iop.Foc]
        
        idx = np.arange(0,res.fun.__len__(),2)
        self.GCPs['residual_x_lstsq'] = np.nan
        self.GCPs['residual_y_lstsq'] = np.nan
        self.GCPs['residual_x_lstsq'].loc[self.GCPs.lstsq_IO.astype(bool)] = res.fun[idx]
        self.GCPs['residual_y_lstsq'].loc[self.GCPs.lstsq_IO.astype(bool)] = res.fun[idx+1]
        
        RMSE_ini = np.sqrt(np.sum(self.GCPs.residual_x_ini**2 + self.GCPs.residual_y_ini**2)/res.fun.__len__())
        print('RMSE initial = ', RMSE_ini)
        print('RMSE lstsq = ', self.estimate.RMSE)
        return self.estimate.new_cam
    
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