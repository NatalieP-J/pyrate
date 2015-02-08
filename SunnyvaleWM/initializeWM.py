import multiprocessing
import subprocess
from numpy import *
import os
import sys
from loadWM import LoadData

def makeWMfiles(rangeval):
	for i in range(rangeval):
		f = open('WMrateget{0}.py'.format(i+1),'wb')
		f.write('i = {0}\n'.format(i))
		f.close()
		os.system('cat templateWMrateget.py >> WMrateget{0}.py'.format(i+1))

