# -*- coding: utf-8 -*-
from os.path import join
from ctypes import *

decor = lambda dirdll, dllname: join(dirdll, dllname + ".dll")

class libhandle:
	def __init__(self,dirdll=None, dlls=None):
		"""Load and upload libraries
		dirdll : directory with libraries
		dlls   : list of libraries for load 
		"""		
		self.dlls = {}
		self.libHandle = {}
		if type(dlls) == list:
			for dll in dlls:
				self.libHandle[dll] = windll.kernel32.LoadLibraryA(decor(dirdll,dll))
				self.dlls[dll] = WinDLL(None, handle=self.libHandle[dll])
		if type(dlls) == str:
			self.libHandle[dlls] = windll.kernel32.LoadLibraryA(decor(dirdll,dlls))
			self.dlls[dlls] = WinDLL(None, handle=self.libHandle[dlls])

	def upload(self):
		""" Unload library
		"""
		for key in self.libHandle:
			del self.dlls[key] # clean up by removing reference to the ctypes library object
			windll.kernel32.FreeLibrary(self.libHandle[key]) # unload the DLL
