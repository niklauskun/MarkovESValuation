import os
import scipy.io

class DirStructure(object):
    """Creates directory structure for inputs and outputs"""

    def __init__(
        self,
        code_directory,
        data_folder="data",
        code_folder="",
        results_folder="",
        results_case="",
    ):
        self.BASE_DIRECTORY = code_directory

        self.DATA_DIRECTORY = os.path.join(
            self.BASE_DIRECTORY, data_folder
        )
        self.CODE_DIRECTORY = os.path.join(self.BASE_DIRECTORY, code_folder)