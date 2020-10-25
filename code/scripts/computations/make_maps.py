import sys
sys.path.append('../../')
import numpy as np
import pandas as pd
from processing import generate_ref_maps_v2, generate_dust_map_1998
from array import array

full_path_ref = "../../../data/processed_maps/ref_maps/"
full_path_dust = "../../../data/processed_maps/dust_maps/"

for nside in [128, 256, 512]:
    maps = generate_ref_maps_v2(nside)
    for m,z1,z2 in zip(maps, [0.1,0.3,0.5,1.2], [0.2,0.4,0.6,1.3]):
        filename = str(z1) + 'z' + str(z2) + '_' + str(nside) + '_' + 'overdensity_field.bin'
        with open(full_path_ref + filename, "wb") as f:
            array('d', m).tofile(f)
            
    
    dust_map = generate_dust_map_1998(nside)
    filename = 'ebv_' + str(nside) + '_1998.bin'
    with open(full_path_dust + filename, "wb") as f:
            array('d', m).tofile(f) #potential mistake here, should be dust_map instead of m!
        
        
       
