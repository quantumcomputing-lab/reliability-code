
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 16:26:57 2022

@author: samudra
"""

import pandas as pd
from qiskit import *
import numpy as np
from qiskit import IBMQ
from qiskit.providers.aer.noise import NoiseModel
from sklearn import metrics
from datetime import datetime
import time

IBMQ.load_account()
#QHUB = IBMQ.get_provider(hub='ibm-q-ornl', group='ornl', project='sys-reserve')
QHUB = IBMQ.get_provider(hub='ibm-q-ornl', group='ornl', project='csc406')

device_list = (
#'ibm_washington',
#'ibmq_brooklyn',
# 'ibmq_montreal',
# 'ibmq_mumbai',
# 'ibm_cairo',
# 'ibm_auckland',
# 'ibm_hanoi',
# 'ibmq_toronto',
# 'ibmq_guadalupe',
# 'ibm_perth',
# 'ibm_lagos',
# 'ibmq_jakarta',
# 'ibmq_manila',
# 'ibmq_bogota',
# 'ibmq_santiago',
# 'ibmq_quito',
# 'ibmq_belem', # belem offline on Mar 22 2022
 'ibmq_lima',
'ibmq_armonk',)

# get parameters
for my_real_device in device_list:
    device = QHUB.get_backend(my_real_device)
    print(device.name())
    properties = device.properties()
    #noise_model=NoiseModel.from_backend(device)
    #simulator = Aer.get_backend('qasm_simulator')
    #coupling_map = device.configuration().coupling_map
    
    nqubits = len(properties._qubits.keys())
    print(nqubits) 
    qr= QuantumRegister(nqubits,'qr')
    cr = ClassicalRegister(nqubits,'cr')
    qcirc = QuantumCircuit(qr, cr)
    
    for i in range(nqubits):   
        qcirc.measure(qr[i],cr[i])
    
    '''
    #rc = execute(qcirc, backend = simulator, memory = True, shots = 8192).result()
    rc= execute(qcirc, backend  = simulator,\
                  #initial_layout = layout,\
                  noise_model    = noise_model,\
                  coupling_map   = coupling_map,\
                  #basis_gates    = test_config['basis_gates'],\
                  #seed_simulator = seed,
                  optimization_level = 0,\
                  memory = True,\
                  shots          = 8192).result()
    '''
    
    # Real Device
    #rc= 
    execute(qcirc, backend  = device,\
            #initial_layout = layout,\
            #noise_model    = test_config['noise_model'],\
            #coupling_map   = coupling_map,\
            #basis_gates    = test_config['basis_gates'],\
            #seed_simulator = seed,
            optimization_level = 0,\
            memory = True,\
            shots          = 8192)
            #.result()

'''
memory = rc.get_memory()
q_ts={}
for qubit_counter in range(nqubits):
    s=''.join([y[qubit_counter] for y in memory])
    q_ts[qubit_counter] = [int(s[i]) for i in range(len(s))]
    
df = pd.DataFrame(q_ts)
#metrics.normalized_mutual_info_score(df.loc[:,0], df.loc[:,1])

timestamp_string = str(datetime.now().strftime("%Y-%m-%d %H%M%S"))
filename = 'register_' + 'time_' + timestamp_string + '.csv'
df.to_csv(filename , index=False)