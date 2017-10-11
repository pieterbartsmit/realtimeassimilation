#!/Library/Frameworks/Python.framework/Versions/3.6/bin/IPython

import sys
import os

from importlib import reload

# Add the model directory to the searchpath
cwd = os.getcwd()
sys.path.append(cwd)
