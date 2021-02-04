import os
import datetime
import re
import numpy as np
import pandas as pd
from scipy.io import loadmat

# directory structure for inputs/outputs
class DirStructure(object):
    """Creates directory structure for inputs (from RTP) and 
    outputs (transition matrices)"""

    def __init__(
        self,
        code_directory,
        code_folder="MarkovESValuation",
        input_folder="RTP_data",
        results_folder="transition_matrix",
        RTP_file = "RTP_NYC_2010_2019.mat",
        t1 = datetime.datetime.strptime("01-01-2010", "%m-%d-%Y"),
        t2 = datetime.datetime.strptime("01-01-2011", "%m-%d-%Y"),
    ):
        self.BASE_DIRECTORY = code_directory
        self.MESV_DIRECTORY = os.path.join(self.BASE_DIRECTORY, code_folder)
        self.RTP_DIRECTORY = os.path.join(self.MESV_DIRECTORY, input_folder)
        self.RESULTS_DIRECTORY = os.path.join(
            self.MESV_DIRECTORY, results_folder + "\\" + re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(t1.year) + "_" + str(t2.year)
        )
        self.RESULTS_DIRECTORY_S = os.path.join(self.RESULTS_DIRECTORY + "\\" + "Summer")
        self.RESULTS_DIRECTORY_W = os.path.join(self.RESULTS_DIRECTORY + "\\" + "Winter")
        self.RESULTS_DIRECTORY_SF = os.path.join(self.RESULTS_DIRECTORY + "\\" + "Spring_Fall")
        self.RESULTS_DIRECTORY_Y = os.path.join(self.RESULTS_DIRECTORY + "\\" + "Year")

    def make_directories(self):
        """create any directories that don't already exist (applies only to output directories)
        """
        if not os.path.exists(self.RESULTS_DIRECTORY):
            os.mkdir(self.RESULTS_DIRECTORY)
        if not os.path.exists(self.RESULTS_DIRECTORY_S):
            os.mkdir(self.RESULTS_DIRECTORY_S)
        if not os.path.exists(self.RESULTS_DIRECTORY_W):
            os.mkdir(self.RESULTS_DIRECTORY_W)
        if not os.path.exists(self.RESULTS_DIRECTORY_SF):
            os.mkdir(self.RESULTS_DIRECTORY_SF)
        if not os.path.exists(self.RESULTS_DIRECTORY_Y):
            os.mkdir(self.RESULTS_DIRECTORY_Y)

# class for loading RTP case data

class LoadInputData(object):
    def __init__(self, f):
        self.f = f
        self.input_dict = {}

    def load_input_data(self, RTP_file):
        """
        Arguments:
            f {class(DirStructure)} -- a folder directory for the case. Needs pointer to input data
        Returns:
            [dict] -- dictionary containing dataframes with loaded data
        """
        ## NREL SourceData folder imports ##
        self.input_dict["RTP"] = pd.DataFrame(loadmat(os.path.join(self.f.RTP_DIRECTORY, RTP_file))['RTP']).T
        print("...completed load data")
        return self.input_dict

def transition_matrix(df, state_num):
    n = state_num
    M = [[0]*n for _ in range(n)]
    for x in range(df.shape[0]):
        transitions = df.iloc[x]
        for (i,j) in zip(transitions,transitions[1:]):
            M[i][j] += 1

    #now convert to probabilities:
    for row in M:
        s = sum(row)
        if s > 0:
            row[:] = [f/s for f in row]
    return M

def write_matrix_case(kw_dict, start, end, dir_structure, **kwargs):
    zero_day = datetime.datetime.strptime("01-01-2010", "%m-%d-%Y")
    max_day = datetime.datetime.strptime("01-01-2020", "%m-%d-%Y")
    # some checks
    assert (
        start - zero_day
    ).days >= 0.0  # check you start after the first day for which data is available
    assert (
        max_day - end
    ).days >= 0.0  # check you end before the last day for which data is available
    assert (
        end - start
    ).days >= 0.0  # check you end after the day you start
    
    try:
        wm = kwargs["winter_month"]
    except KeyError:
        print("NOTE: no winter_month, add winter_month based on default behavior")
        wm = [1,2]
    try:
        sm = kwargs["summer_month"]
    except KeyError:
        print("NOTE: no summer_month, add summer_month based on default behavior")
        sm = [7,8]
    try:
        ts = kwargs["time_step"]
    except KeyError:
        print("NOTE: no time_step, add time_step based on default behavior")
        ts = 12
    try:
        pb = kwargs["price_bar"]
    except KeyError:
        print("NOTE: no price_bar, add price_bar based on default behavior")
        pb = 200
    try:
        sg = kwargs["state_gap"]
    except KeyError:
        print("NOTE: no state_gap, add state_gap based on default behavior")
        sg = 10
    try:
        mn = kwargs["matrix_num"]
    except KeyError:
        print("NOTE: no matrix_num, add matrix_num based on default behavior")
        mn = 24

    sn = pb//sg + 2 # state number
    rtp = kw_dict['RTP'].loc[(start - zero_day).days : (end - zero_day).days - 1,:].copy() #slicing RTP data
    rtp = (rtp//sg+1).astype(int)
    rtp[rtp < 0] = 0
    rtp[rtp > sn - 2] = sn - 1
    rtp['date'] = pd.date_range(start, end-datetime.timedelta(days = 1))
    rtp_w = pd.DataFrame(columns=list(rtp))
    rtp_s = pd.DataFrame(columns=list(rtp))
    rtp_sf = pd.DataFrame(columns=list(rtp))
    for m in list(range(1,13)):
        if m in wm:
            rtp_w = rtp_w.append(rtp[rtp['date'].map(lambda x: x.month) == m])
        elif m in sm:
            rtp_s = rtp_s.append(rtp[rtp['date'].map(lambda x: x.month) == m])
        else:
            rtp_sf = rtp_sf.append(rtp[rtp['date'].map(lambda x: x.month) == m])

    for i in range(mn):
        t1 = i*ts
        t2 = (i+1)*ts
        slice = rtp.iloc[:,t1:t2]
        slice_w = rtp_w.iloc[:,t1:t2]
        slice_s = rtp_s.iloc[:,t1:t2]
        slice_sf = rtp_sf.iloc[:,t1:t2]

        m = transition_matrix(slice, sn)
        m_w = transition_matrix(slice_w, sn)
        m_s = transition_matrix(slice_s, sn)
        m_sf = transition_matrix(slice_sf, sn)
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_Y,'matrix'+str(i)+".csv"), m, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_W,'matrix'+str(i)+".csv"), m_w, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_S,'matrix'+str(i)+".csv"), m_s, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_SF,'matrix'+str(i)+".csv"), m_sf, delimiter = ',')

    print("...completed creating case transition matrices!")
    return None


