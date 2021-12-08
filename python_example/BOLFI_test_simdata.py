import numpy as np
import elfi
from scipy.io import loadmat
import elfi.clients.multiprocessing

from simulator_test_simdata import simulator, percentile_values

# This file is drawing samples from the approximated posterior distribution by utilising BOLFI.

# Written by: Anttoni Ojanperä, Nitin Williams
# Based on code written by: Antti Karvanen

num_processes = 2
seed = 23112021
np.random.seed(seed)

# Loading observed summary statistics
mat_model = loadmat('WCnet_noisy_300s_MEG_simulatednetwork_sub1.mat')
observed_model = mat_model['wPLImags_vals']


# Set up multiprocessing
elfi.set_client(elfi.clients.multiprocessing.Client(num_processes=num_processes))


# Specifying prior distributions
wee = elfi.Prior('norm',20,5)
wei = elfi.Prior('norm',18,6)
wie = elfi.Prior('norm',18,6)
wii = elfi.Prior('norm',1,0.2)
be = elfi.Prior('norm',3,1)
bi = elfi.Prior('norm',5,1)
taue = elfi.Prior('norm',18.6,3.8)
taui = elfi.Prior('norm',15.1,4.7)
k = elfi.Prior('norm',1.5,0.5)
IH = elfi.Prior('norm',2.5,0.5)
psi_sigma = elfi.Prior('norm',0.15,0.05)

sim = elfi.Simulator(simulator, wee, wei, wie, wii, be, bi, taue, taui, k, IH, psi_sigma, observed=observed_model)

# Summary statistics

S1 = elfi.Summary(percentile_values,sim,0)
S2 = elfi.Summary(percentile_values,sim,10)
S3 = elfi.Summary(percentile_values,sim,20)
S4 = elfi.Summary(percentile_values,sim,30)
S5 = elfi.Summary(percentile_values,sim,40)
S6 = elfi.Summary(percentile_values,sim,50)
S7 = elfi.Summary(percentile_values,sim,60)
S8 = elfi.Summary(percentile_values,sim,70)
S9 = elfi.Summary(percentile_values,sim,80)
S10 = elfi.Summary(percentile_values,sim,90)
S11 = elfi.Summary(percentile_values,sim,100)

# Distance

d = elfi.Distance('euclidean', S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11)

################################
# BOLFI implementation:

# log_d = elfi.Operation(np.log, sim.model['d'])

output_pool = elfi.OutputPool(outputs=['wee','wei','wie','wii','be','bi','taue','taui','k','IH','psi_sigma','sim','d',
               'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11'])

# pool=output_pool
bolfi = elfi.BOLFI(d, batch_size=1, initial_evidence=4, update_interval=4,
                   bounds={'wee': (10, 30), 'wei': (6, 30), 'wie': (6, 30), 'wii': (0.6, 1.4), 'be':(1,5),
                           'bi':(3,7), 'taue':(11,26.2), 'taui':(5.7,24.5), 'k':(0.5,2.5), 'IH':(1.5,3.5),
                           'psi_sigma':(0.05,0.25)},
                   acq_noise_var=0.1, seed=seed)

# np.save('bolfi_test_23112021_sim_data', bolfi)

post = bolfi.fit(n_evidence=6)

bolfi.plot_discrepancy()

result_BOLFI = bolfi.sample(6)

result_BOLFI.plot_marginals();
result_BOLFI.plot_pairs();

np.save('bolfi_test_19112021_sim_data_bolfisample', result_BOLFI.samples_array)

