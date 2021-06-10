import pandas as pd
import numpy as np

D=np.load('pcangsd.selection.npy')
np.savetxt("out.selection.2.txt", D)
