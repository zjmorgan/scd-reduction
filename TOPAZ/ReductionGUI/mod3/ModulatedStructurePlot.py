import os
import sys
import ReduceDictionary
from SatellitePlot import SatellitePlot

#
# Get the config file name and the run number to process from the command line
#
if (len(sys.argv) < 2):
  print("You MUST give the config file name on the command line")
  exit(0)

config_file_name = sys.argv[1]

#
# Load the parameter names and values from the specified configuration file 
# into a dictionary and set all the required parameters from the dictionary.
#
config_file_name = os.path.dirname(os.path.realpath(__file__)) +'/' + config_file_name 


params_dictionary = ReduceDictionary.LoadDictionary(config_file_name )

output_directory      = params_dictionary[ "output_directory" ]
exp_name              = params_dictionary[ "exp_name" ]
cell_type             = params_dictionary[ "cell_type" ] 
centering             = params_dictionary[ "centering" ]
run_nums              = params_dictionary[ "run_nums" ]

test = SatellitePlot()
for run in run_nums:
    output_name = output_directory + "/" + str(run) + "_" + exp_name + "_" + cell_type + "_" + centering
    output_name = output_directory + "/" + str(run) + "_Niggli"
    integrate_file = output_name + ".integrate"
    if run == run_nums[0]:
        matrix_file = output_name + ".mat"
    test.load_peaks(integrate_file,matrix_file)
    test.plot_Qpeaks()

if len(run_nums) < 2:
    sys.exit(0)
#output_name = output_directory + "/" + exp_name + "_" + cell_type + "_" + centering
output_name = output_directory + "/" + exp_name + "_Niggli"
integrate_file = output_name + ".integrate"
#matrix_file = output_name + ".mat"
test.load_peaks(integrate_file,matrix_file)
test.plot_Qpeaks()

#
# Try to get this to terminate when run by ReduceSCD_Parallel.py, from NX session
#
sys.exit(0)

