# -*- coding: utf-8 -*-
#!/usr/bin/python


import sys
from sys import platform
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from pylab import *
from pylab import *
import scipy
from scipy import ndimage,optimize, signal
import scipy.ndimage.interpolation as spni
import scipy.fftpack as spf
from scipy.signal import *
from skimage.feature import match_template
from PIL import Image
import os
from os.path import isfile, join
import string
import pyqtgraph as pg
from PySide import QtGui, QtCore
from pyqtgraph import QtGui, QtCore
import h5py
import tomopy
import time
from numpy import linalg as LA


'''
needs sequence of tiff files and theta file.
saves shift file in ./shift
and other files in ./sav (projection and reprojection
'''

class Example(QtGui.QMainWindow):
    
        def __init__(self):
                super(Example, self).__init__()
                self.initUI()

        def initUI(self):
                runAction = QtGui.QAction("Run",self)
                runAction.triggered.connect(self.run)

                menubar = self.menuBar()
                fileMenu = menubar.addMenu('&File')
                fileMenu.addAction(runAction)

                self.setWindowTitle("iterative")

                self.show()

        def run(self):
                self.openFolder()
                self.openThetaFile()
                self.readTheta()
                self.convert2array()
                self.refineShift(self.data,saveok=True)
                
        def openFolder(self):
                try:
                        folderName = QtGui.QFileDialog.getExistingDirectory(self,"Open Folder",
                                                           QtCore.QDir.currentPath())
                        global RH,fileName
                        RH=folderName
                        folderName=str(folderName)
                        file_name_array = [ f for f in os.listdir(folderName) if isfile(join(folderName,f))]
                        file_name_array = [ f for f in os.listdir(folderName) if string.find(f,"tif")!=-1]
                        file_name_array = [folderName+"/"+f for f in file_name_array]

                        self.directory=1
                        self.fileNames=file_name_array
                        fileName=file_name_array
                        self.projections=len(self.fileNames)
                        self.convert2array()
                        
                except IndexError:
                        print "no folder has been selected"
                except OSError:
                        print "no folder has been selected"
        def convert2array(self):
                self.y, self.x=np.asarray(Image.open(str(self.fileNames[0])),dtype=float32).shape
                self.data=zeros([self.projections,self.y,self.x],dtype=float32)
                for i in arange(self.projections):
                        self.data[i,:,:]=np.asarray(Image.open(str(self.fileNames[i])),dtype=float32)[...]
                
        def openThetaFile(self):
                try:
                        self.thetaFileName = QtGui.QFileDialog.getOpenFileName(self, "Open Theta File",
                                                                          QtCore.QDir.currentPath(),filter="txt (*.txt)")
                        
                except IndexError:
                        print "choose file!"
                except OSError:
                        print "choose file!"
        def readTheta(self):
                f=open(str(self.thetaFileName),'r')
                read=f.readlines()
                self.thetaArray = zeros(len(read),dtype=float32)
                for i,theta in enumerate(read):
                        self.thetaArray[i]=theta[theta.rfind(",")+1:]
        
        def xcor(self,img1, img2):
                fa=fft2(img1)
                fb=fft2(img2)
                c=abs(ifft2(fa.conj()*fb))#/(abs(fa)*abs(fb)))
                t0,t1=np.unravel_index(np.argmax(c),img1.shape)
                if t0 > img1.shape[0]//2:
                        t0 -=img1.shape[0]
                if t1 > img1.shape[1]//2:
                        t1-=img1.shape[1]
                return t0,t1

        def diff_write_txt2(self,mat, filename):
                f = open(filename,"w")
                for i in arange(mat.shape[1]):
                        f.writelines("shift for " + str(i).zfill(3)+ " is y: " + str(mat[0,i]) + "x: " + str(mat[1,i]) + "\n")

        def x_rotation(self,matrix,ang):
                temp = zeros([len(ang),matrix.shape[0],matrix.shape[2]])
                print temp.shape
                for i in arange(len(ang)):
                        b=ndimage.rotate(matrix,-ang[i],axes=(1,2),reshape=False)
                        b1=sum(b,axis=1)
                        print b1.shape, ang[i]
                        temp[i,:,:]=b1
                return temp
        def xcor_data_proj_onlyx(self,data,proj,printok=False):
                copy=zeros(data.shape,dtype=float32)
                diff = zeros([2,data.shape[0]],dtype=float32)
                for i in arange(data.shape[0]):
                        t0,t1=self.xcor(data[i,:,:],proj[i,:,:])
                        copy[i,:,:]=np.roll(data[i,:,:],-t1,axis=1)
                        diff[0,i]=t0
                        diff[1,i]=t1
                        if printok:
                                print t0,t1
                return copy, diff
        def xcor_data_proj(self,data,proj):
                copy=zeros(data.shape,dtype=float32)
                diff = zeros([2,data.shape[0]],dtype=float32)
                for i in arange(data.shape[0]):
                        t0,t1=self.xcor(data[i,:,:],proj[i,:,:])
                        copy[i,:,:]=np.roll(np.roll(data[i,:,:],-t0,axis=0),-t1,axis=1)
                        diff[0,i]=t0
                        diff[1,i]=t1
                        if printok:
                                print t0,t1
                return copy, diff
        def refineShift(self,orig,num_iter=1,fname="./2ms_full/",saveok=False):
                num_loop=0
                mat=zeros(orig.shape,dtype=np.float32)
                mat[...]=orig[...]
                while num_loop < num_iter:
 
                        recon=tomopy.recon(mat,self.thetaArray*np.pi/180,algorithm="mlem",
                                           reg_par=np.array([1,0.01],dtype=np.float32),num_iter=30)
                        print "reconstruction for "+str(num_loop) +" is done"
                        temp=self.x_rotation(recon,self.thetaArray)
                        temp2=zeros(temp.shape,dtype=float32)
                        temp2[...]=wiener(temp)[...]
                        #temp2[np.where(temp2<temp2.max()*0.01)]=0.1**5
                        #print temp2
                        mat,diff=self.xcor_data_proj_onlyx(mat,temp2)
                        if saveok:
                                tomopy.write_tiff_stack(mat, fname=fname+str(num_loop)+"/"+"stack_",dtype=float32,digit=3)
                                #tomopy.write_tiff_stack(temp, fname=fname+str(num_loop)+"_bw/stack_",dtype=float32,digit=3)
                                tomopy.write_tiff_stack(temp2, fname=fname+str(num_loop)+"_reprojection/stack_",dtype=float32,digit=3)
                                tomopy.write_tiff_stack(recon, fname=fname+str(num_loop)+"_recon/stack_",dtype=float32,digit=3)
                                self.diff_write_txt2(diff,"./shift/shift_"+str(num_loop+1)+".txt")
                                print diff
                        num_loop+=1
        
def main():
    
            app = QtGui.QApplication(sys.argv)
            ex = Example()
            sys.exit(app.exec_())

if __name__ == '__main__':
            main()
            




