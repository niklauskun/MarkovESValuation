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
        self.RESULTS_DIRECTORY_SEASON = os.path.join(self.RESULTS_DIRECTORY + "\\" + "Season")
        self.RESULTS_DIRECTORY_SEASON_S = os.path.join(self.RESULTS_DIRECTORY_SEASON + "\\" + "Summer")
        self.RESULTS_DIRECTORY_SEASON_NS = os.path.join(self.RESULTS_DIRECTORY_SEASON + "\\" + "Non_Summer")
        self.RESULTS_DIRECTORY_WEEK = os.path.join(self.RESULTS_DIRECTORY + "\\" + "Week")
        self.RESULTS_DIRECTORY_WEEK_WD = os.path.join(self.RESULTS_DIRECTORY_WEEK + "\\" + "Weekday")
        self.RESULTS_DIRECTORY_WEEK_WE = os.path.join(self.RESULTS_DIRECTORY_WEEK + "\\" + "Weekend")
        self.RESULTS_DIRECTORY_SW = os.path.join(self.RESULTS_DIRECTORY + "\\" + "Season_Week")
        self.RESULTS_DIRECTORY_SW_SWD = os.path.join(self.RESULTS_DIRECTORY_SW + "\\" + "Summer_Weekday")
        self.RESULTS_DIRECTORY_SW_SWE = os.path.join(self.RESULTS_DIRECTORY_SW + "\\" + "Summer_Weekend")
        self.RESULTS_DIRECTORY_SW_NSWD = os.path.join(self.RESULTS_DIRECTORY_SW + "\\" + "Non_Summer_Weekday")
        self.RESULTS_DIRECTORY_SW_NSWE = os.path.join(self.RESULTS_DIRECTORY_SW + "\\" + "Non_Summer_Weekend")
        self.RESULTS_DIRECTORY_Y = os.path.join(self.RESULTS_DIRECTORY + "\\" + "Year")
        self.RESULTS_DIRECTORY_N = os.path.join(self.RESULTS_DIRECTORY + "\\" + "Null")

    def make_directories(self):
        """create any directories that don't already exist (applies only to output directories)
        """
        if not os.path.exists(self.RESULTS_DIRECTORY):
            os.mkdir(self.RESULTS_DIRECTORY)
        if not os.path.exists(self.RESULTS_DIRECTORY_SEASON):
            os.mkdir(self.RESULTS_DIRECTORY_SEASON)
        if not os.path.exists(self.RESULTS_DIRECTORY_SEASON_S):
            os.mkdir(self.RESULTS_DIRECTORY_SEASON_S)
        if not os.path.exists(self.RESULTS_DIRECTORY_SEASON_NS):
            os.mkdir(self.RESULTS_DIRECTORY_SEASON_NS)
        if not os.path.exists(self.RESULTS_DIRECTORY_WEEK):
            os.mkdir(self.RESULTS_DIRECTORY_WEEK)
        if not os.path.exists(self.RESULTS_DIRECTORY_WEEK_WD):
            os.mkdir(self.RESULTS_DIRECTORY_WEEK_WD)
        if not os.path.exists(self.RESULTS_DIRECTORY_WEEK_WE):
            os.mkdir(self.RESULTS_DIRECTORY_WEEK_WE)
        if not os.path.exists(self.RESULTS_DIRECTORY_SW):
            os.mkdir(self.RESULTS_DIRECTORY_SW)
        if not os.path.exists(self.RESULTS_DIRECTORY_SW_SWD):
            os.mkdir(self.RESULTS_DIRECTORY_SW_SWD)
        if not os.path.exists(self.RESULTS_DIRECTORY_SW_SWE):
            os.mkdir(self.RESULTS_DIRECTORY_SW_SWE)
        if not os.path.exists(self.RESULTS_DIRECTORY_SW_NSWD):
            os.mkdir(self.RESULTS_DIRECTORY_SW_NSWD)
        if not os.path.exists(self.RESULTS_DIRECTORY_SW_NSWE):
            os.mkdir(self.RESULTS_DIRECTORY_SW_NSWE)
        if not os.path.exists(self.RESULTS_DIRECTORY_Y):
            os.mkdir(self.RESULTS_DIRECTORY_Y)
        if not os.path.exists(self.RESULTS_DIRECTORY_N):
            os.mkdir(self.RESULTS_DIRECTORY_N)

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

def null_transition_matrix(df, state_num):
    num = df._get_numeric_data()
    n = state_num
    s = num.shape[0]*num.shape[1]
    M = [[0]*n for _ in range(n)]
    for i, j in ((a,b) for a in range(n) for b in range(n)):
        f = num[num==i].count().sum()
        M[j][i] = f/s
    return M

def write_matrix_case(kw_dict, start, end, dir_structure, RTP_file, **kwargs):
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
        ss = kwargs["summer_start"]
    except KeyError:
        print("NOTE: no summer_start, add summer_start based on default behavior")
        ss = 125
    try:
        se = kwargs["summer_end"]
    except KeyError:
        print("NOTE: no summer_end, add summer_end based on default behavior")
        se = 285
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
    rtp = kw_dict['RTP'].loc[(start - zero_day).days - 2 : (end - zero_day).days - 2,:].copy() #slicing RTP data, -2 because 2 leap days in 2012 and 2016 were ignored
    rtp = (rtp//sg+1).astype(int)
    rtp[rtp < 0] = 0
    rtp[rtp > sn - 2] = sn - 1
    rtp['date'] = pd.date_range(start, end)
    rtp_summer = pd.DataFrame(columns=list(rtp)) # initialize summer RTP data
    rtp_nsummer = pd.DataFrame(columns=list(rtp)) # initialize non-summer RTP data
    rtp_weekday = pd.DataFrame(columns=list(rtp)) # initialize weekday RTP data
    rtp_weekend = pd.DataFrame(columns=list(rtp)) # initialize weekend RTP data
    rtp_sweekday = pd.DataFrame(columns=list(rtp)) # initialize summer weekday RTP data
    rtp_sweekend = pd.DataFrame(columns=list(rtp)) # initialize summer weekend RTP data
    rtp_nsweekday = pd.DataFrame(columns=list(rtp)) # initialize non-summer weekday RTP data
    rtp_nsweekend = pd.DataFrame(columns=list(rtp)) # initialize non-summer weekend RTP data

    for n in range(rtp.shape[0]):
        date = rtp.iloc[n]['date']
        weekday = date.weekday() in [0,1,2,3,4]
        year_start = datetime.datetime(date.year,1,1)
        summer = ss - 1 <= (date-year_start).days <= se - 1
        if weekday == True:
            rtp_weekday = rtp_weekday.append(rtp.iloc[n])
            if summer == True:
                rtp_summer = rtp_summer.append(rtp.iloc[n])
                rtp_sweekday = rtp_sweekday.append(rtp.iloc[n])
            else:
                rtp_nsummer = rtp_nsummer.append(rtp.iloc[n])
                rtp_nsweekday = rtp_nsweekday.append(rtp.iloc[n])
        else:
            rtp_weekend = rtp_weekend.append(rtp.iloc[n])
            if summer == True:
                rtp_summer = rtp_summer.append(rtp.iloc[n])
                rtp_sweekend = rtp_sweekend.append(rtp.iloc[n])
            else:
                rtp_nsummer = rtp_nsummer.append(rtp.iloc[n])
                rtp_nsweekend = rtp_nsweekend.append(rtp.iloc[n])
        # or use lambda? rtp_s = rtp_s.append(rtp[rtp['date'].map(lambda x: x.month) == m])
    m_n = null_transition_matrix(rtp, sn) #create null transition matrices, which means price independent from previous timespoint
    
    print("...RTP subset created")

    for i in range(mn):
        t1 = i*ts
        t2 = (i+1)*ts
        slice = rtp.iloc[:,t1:t2]
        slice_summer = rtp_summer.iloc[:,t1:t2]
        slice_nsummer = rtp_nsummer.iloc[:,t1:t2]
        slice_weekday = rtp_weekday.iloc[:,t1:t2]
        slice_weekend = rtp_weekend.iloc[:,t1:t2]
        slice_sweekday = rtp_sweekday.iloc[:,t1:t2]
        slice_sweekend = rtp_sweekend.iloc[:,t1:t2]
        slice_nsweekday = rtp_nsweekday.iloc[:,t1:t2]
        slice_nsweekend = rtp_nsweekend.iloc[:,t1:t2]

        m = transition_matrix(slice, sn)
        m_s = transition_matrix(slice_summer, sn)
        m_ns = transition_matrix(slice_nsummer, sn)
        m_wd = transition_matrix(slice_weekday, sn)
        m_we = transition_matrix(slice_weekend, sn)
        m_swd = transition_matrix(slice_sweekday, sn)
        m_swe = transition_matrix(slice_sweekend, sn)
        m_nswd = transition_matrix(slice_nsweekday, sn)
        m_nswe = transition_matrix(slice_nsweekend, sn)

        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_Y, re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(start.year) + "_" + str(end.year) + "_" + 'matrix_year_'+str(i)+".csv"), m, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_SEASON_S, re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(start.year) + "_" + str(end.year) + "_" + 'matrix_summer_'+str(i)+".csv"), m_s, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_SEASON_NS, re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(start.year) + "_" + str(end.year) + "_" + 'matrix_nonsummer_'+str(i)+".csv"), m_ns, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_WEEK_WD, re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(start.year) + "_" + str(end.year) + "_" + 'matrix_weekday_'+str(i)+".csv"), m_wd, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_WEEK_WE, re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(start.year) + "_" + str(end.year) + "_" + 'matrix_weekend_'+str(i)+".csv"), m_we, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_SW_SWD, re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(start.year) + "_" + str(end.year) + "_" + 'matrix_summer_weekday_'+str(i)+".csv"), m_swd, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_SW_SWE, re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(start.year) + "_" + str(end.year) + "_" + 'matrix_summer_weekend_'+str(i)+".csv"), m_swe, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_SW_NSWD, re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(start.year) + "_" + str(end.year) + "_" + 'matrix_nonsummer_weekday_'+str(i)+".csv"), m_nswd, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_SW_NSWE, re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(start.year) + "_" + str(end.year) + "_" + 'matrix_nonsummer_weekend_'+str(i)+".csv"), m_nswe, delimiter = ',')
        np.savetxt(os.path.join(dir_structure.RESULTS_DIRECTORY_N, re.findall(r'_(.+?)_', RTP_file )[0] + "_" + str(start.year) + "_" + str(end.year) + "_" + 'matrix_null_'+str(i)+".csv"), m_n, delimiter = ',')

    print("...completed creating case transition matrices!")
    return None


