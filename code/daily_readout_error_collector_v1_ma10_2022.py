#import qiskit
from qiskit import IBMQ
import pandas as pd
from datetime import datetime, timedelta


IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q-ornl', group='ornl', project='csc406')

DEVICE = 'ibm_washington' # Enter the device name here
backend = provider.get_backend(DEVICE)

(START_YEAR, START_MONTH, START_DAY) = (2021, 12, 1)
(STOP_YEAR, STOP_MONTH, STOP_DAY)=(2022, 3, 9) 
start_dt = datetime(year = START_YEAR, month = START_MONTH, day = START_DAY)
stop_dt = datetime(year=STOP_YEAR, month = STOP_MONTH, day = STOP_DAY)

columns = ['last_update_date', 'input_date']
for i in range(127):
    columns = columns + ['q'+str(i)+'_readout_err']+['q'+str(i)+'_readout_err_cal_time']

device_data = pd.DataFrame(columns = columns)
df = pd.DataFrame(columns = columns)

i=0
dt = stop_dt
while (dt >= start_dt):
    print(dt.strftime('%Y - %m - %d - %H'))
    prop = backend.properties(datetime=dt)

    prop_dict=prop.to_dict()

    device_data.loc[i,'last_update_date'] = prop_dict['last_update_date']
    device_data.loc[i,'input_date'] = dt.strftime('%Y - %m - %d')
        
    ### Readout Error
    for qubit_num in range(127):
        my_dict_array = prop_dict['qubits'][qubit_num]
        for dict_item in my_dict_array:
            if dict_item['name']=='readout_error':
                device_data.loc[i,'q'+str(qubit_num)+'_readout_err'] = dict_item['value']
                device_data.loc[i,'q'+str(qubit_num)+'_readout_err_cal_time'] = dict_item['date']
    i=i+1
    dt = dt -  timedelta(days=1)

reversed_data = device_data.iloc[::-1].reset_index(drop=True)

#reversed_data.to_csv('test.csv')