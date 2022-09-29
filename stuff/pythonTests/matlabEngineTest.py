#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 07:44:26 2021

@author: marlon
"""
import numpy as np

import matlab.engine
eng = matlab.engine.start_matlab()

tf = eng.isprime(37)
print(tf)

# temp = eng.solve([1 3;4 5],[1;1])


eng.eval('tspan = [0 5];',nargout = 0)
eng.eval('y0 = 0;',nargout = 0)
t,y=eng.eval('ode45(@(t,y) 2*t, tspan, y0)',nargout = 2)
eng.plot(t,y,'-o')

# patch with matlab data
x = eng.eval('[0 1 1 0];')
y = eng.eval('[0 0 1 1];')
eng.patch(x,y,'red')

# patch with python data
xP = np.array([[0,1,1,0],[0,1,1,0]])
yP = np.array([[0,0,1,1],[1,1,2,2]])
xP = matlab.double(xP.tolist())
yP = matlab.double(yP.tolist())
eng.patch(matlab.double(xP),matlab.double(yP),'red')
# eng.patch('XData',matlab.double(xP),'YData',matlab.double(yP),'red')

##############################################################
nodes = np.array([[0. , 0. , 0. ],
       [0.5, 0. , 0. ],
       [1. , 0. , 0. ],
       [0. , 0.5, 0. ],
       [0.5, 0.5, 0. ],
       [1. , 0.5, 0. ],
       [0. , 1. , 0. ],
       [0.5, 1. , 0. ],
       [1. , 1. , 0. ],
       [0. , 0. , 0.5],
       [0.5, 0. , 0.5],
       [1. , 0. , 0.5],
       [0. , 0.5, 0.5],
       [0.5, 0.5, 0.5],
       [1. , 0.5, 0.5],
       [0. , 1. , 0.5],
       [0.5, 1. , 0.5],
       [1. , 1. , 0.5],
       [0. , 0. , 1. ],
       [0.5, 0. , 1. ],
       [1. , 0. , 1. ],
       [0. , 0.5, 1. ],
       [0.5, 0.5, 1. ],
       [1. , 0.5, 1. ],
       [0. , 1. , 1. ],
       [0.5, 1. , 1. ],
       [1. , 1. , 1. ]])
edof = np.array([[ 1.,  2.,  5.,  4., 10., 11., 14., 13.],
       [ 2.,  3.,  6.,  5., 11., 12., 15., 14.],
       [ 4.,  5.,  8.,  7., 13., 14., 17., 16.],
       [ 5.,  6.,  9.,  8., 14., 15., 18., 17.],
       [10., 11., 14., 13., 19., 20., 23., 22.],
       [11., 12., 15., 14., 20., 21., 24., 23.],
       [13., 14., 17., 16., 22., 23., 26., 25.],
       [14., 15., 18., 17., 23., 24., 27., 26.]],dtype=int)
# edof = edof - np.ones(edof.shape, dtype="i")
SX = np.array([[1, 2, 3, 4],
      [5, 8, 7, 6],
      [1, 5, 6, 2],
      [2, 6, 7, 3],
      [3, 7, 8, 4],
      [4, 8, 5, 1]],dtype=int)
SX = SX - np.ones(SX.shape, dtype="i")
Faces = np.concatenate([edof[:,SX[0,:]], edof[:,SX[1,:]], edof[:,SX[2,:]], edof[:,SX[3,:]], edof[:,SX[4,:]], edof[:,SX[5,:]]])
FaceColor = np.array([0,1,0])
FaceAlpha = 0.5
eng.patch('Vertices',matlab.double(nodes.tolist()),'Faces',matlab.uint32(Faces.tolist()),'FaceColor',matlab.double(FaceColor.tolist()),'FaceAlpha',FaceAlpha)
eng.eval('view(3)')
