import os
import datetime
import re
import numpy as np
from pre_processing_functions import DirStructure
from pre_processing_functions import LoadInputData
from pre_processing_functions import write_matrix_case

# inputs for running
folder_path = os.path.join(os.environ["HOMEPATH"], "Desktop")
code_folder_path = "MarkovESValuation"
start = datetime.datetime.strptime("01-01-2019", "%m-%d-%Y")  # day case starts on, (needs to be greater than 03-03-2011, due to NYISO has two days missing in MAR 2021)
end = datetime.datetime.strptime("12-31-2019", "%m-%d-%Y")  # day case ends on
RTP_file = "RTP_WEST_2010_2019.mat"
DAP_file = "DAP_WEST_2010_2019.mat"
DABias = True # Train Markov process model for real-time model (Flase) or DAP-RTP bias model (True)

# create a file structure object
f = DirStructure(
    folder_path, 
    code_folder=code_folder_path, 
    RTP_folder="RTP_data", 
    DAP_folder="DAP_data", 
    results_folder="transition_matrix", 
    RTP_file = RTP_file,
    t1 = start, 
    t2 = end, 
) # the first arg should be the local directory you put code folder in
# the second is whatever you named the folder where you put the code
f.make_directories()

# load input data
data_class = LoadInputData(f)
kw_dict = data_class.load_input_data(RTP_file, DAP_file)

# write transition matrice for cases
optional_kwargs = {
    "summer_start": 125,
    "summer_end": 285,
    "time_step": 12, # slices in one hour
    "price_bar": 200, # price greater than price bar is price spike
    "state_gap": 10, # state gap
    "bias_bar": 50, # bias greater than 50 is gournded to positive spike node (smaller than -50 is grounded to negative spike node)
    "bias_gap": 10, # bias state gap
    "matrix_num": 24,
}

write_matrix_case(kw_dict, start, end, f, RTP_file, DABias, **optional_kwargs)
