# -*- coding: utf-8 -*-
import os, sys


class codes:
	def __init__(self,curdir):
		self.curdir=os.path.abspath(curdir)
		self.tooldir=None
	def update(self,input):
		print "Tool not found"
	def set_directory(self,directory):
		self.tooldir = os.path.abspath(directory)





result = sketch.run(idata)
sketch