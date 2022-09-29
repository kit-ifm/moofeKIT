#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 11:03:34 2021

@author: marlon
"""

import matlab.engine
eng = matlab.engine.start_matlab()
eng.startScriptSimple(nargout=0)
# setupObject = eng.startScriptSimple
# eng.getfield(setupObject, 'fileName')
# eng.subsref(setupObject, {'type':'.','subs':'fileName'})
# eng.edit('startScriptSimple',nargout=0)
