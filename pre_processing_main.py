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
start = datetime.datetime.strptime("01-01-2018", "%m-%d-%Y")  # day case starts on, (needs to be greater than 02-29-2016, due to the leap year)
end = datetime.datetime.strptime("12-31-2018", "%m-%d-%Y")  # day case ends on
RTP_file = "RTP_WEST_2010_2019.mat"

# create a file structure object
f = DirStructure(
    folder_path, 
    code_folder=code_folder_path, 
    input_folder="RTP_data", 
    results_folder="transition_matrix", 
    RTP_file = RTP_file,
    t1 = start, 
    t2 = end, 
) # the first arg should be the local directory you put code folder in
# the second is whatever you named the folder where you put the code
f.make_directories()

# load input data
data_class = LoadInputData(f)
kw_dict = data_class.load_input_data(RTP_file)

# write transition matrice for cases
optional_kwargs = {
    "summer_start": 125,
    "summer_end": 285,
    "time_step": 12, # slices in one hour
    "price_bar": 200, # price greater than price bar is price spike
    "state_gap": 10, # state gap
    "matrix_num": 24,
}

write_matrix_case(kw_dict, start, end, f, RTP_file, **optional_kwargs)
