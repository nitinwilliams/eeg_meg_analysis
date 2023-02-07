import matlab.engine
import numpy as np
import elfi
import GPy
from scipy.io import loadmat

eng = matlab.engine.start_matlab()

seed = 9092022
np.random.seed(seed)
test_indicator = 0

# Loading observed summary statistics 
mat_model = loadmat('Z:\\Documents\\MEGMOD\\MATLAB\\data\\WCnetnoisy_distancedependentdelays_GPfitting_expMEGdata_75subs_12params_connstrengths_refdata.mat')
observed_model = mat_model['wPLImags_connstrengths_ref']

# Loading IO combinations
model_simdata = loadmat('Z:\\Documents\\MEGMOD\\MATLAB\\data\\WCnetnoisy_distancedependentdelays_SSI_logdiscrepancies_connstrengths_expMEGdata_75subs_12params_radius1pt5_9092sims.mat')
Y = model_simdata['Y']
X = model_simdata['X']

# Specifying prior distributions
wee = elfi.Prior('norm',0,1)
wei = elfi.Prior('norm',0,1)
wie = elfi.Prior('norm',0,1)
wii = elfi.Prior('norm',0,1)
be = elfi.Prior('norm',0,1)
bi = elfi.Prior('norm',0,1)
taue = elfi.Prior('norm',0,1)
taui = elfi.Prior('norm',0,1)
k = elfi.Prior('norm',0,1)
IH = elfi.Prior('norm',0,1)
psi_sigma = elfi.Prior('norm',0,1)
velocity = elfi.Prior('norm',0,1)

# Specifiying simulator function
def simulator(wee, wei, wie, wii, be, bi, taue, taui, k, IH, psi_sigma, velocity, batch_size=1, random_state=None):
    global test_indicator
    test_indicator = test_indicator + 1
    print(f"\n\nSTARTING NEW SIMULATION: {test_indicator}")
    print(f"With parameters: {wee.item(),wei.item(),wie.item(),wii.item(),be.item(),bi.item(),taue.item(),taui.item(),k.item(),IH.item(),psi_sigma.item(),velocity.item()}")
    ret = eng.simulate_WCnetnoisy_distancedependentdelays(wee.item(),wei.item(),wie.item(),wii.item(),be.item(),bi.item(),taue.item(),taui.item(),k.item(),IH.item(),psi_sigma.item(),velocity.item())
    return np.array(ret)

sim = elfi.Simulator(simulator,wee,wei,wie,wii,be,bi,taue,taui,k,IH,psi_sigma,velocity,observed=observed_model)

# Specifying summary function
def same_values(x):
    C=np.transpose(x);
    return C

# Summary
S = elfi.Summary(same_values,sim)

# Specifying distance function
def ssi_est(x,y):
    d = eng.ssim_calc(matlab.double(x.tolist()),matlab.double(y.tolist()))
    return np.atleast_1d(np.array(d))

# Distance
d = elfi.Distance(ssi_est,S)

# Log of distance
log_d = elfi.Operation(np.log, sim.model['d'])

# Setting BOLFI input parameters
output_pool = elfi.OutputPool(outputs=['wee','wei','wie','wii','be','bi','taue','taui','k','IH','psi_sigma','velocity','d','log_d'])

final_kernel=GPy.kern.RBF(12, ARD=True, name='rbf') + GPy.kern.Bias(12, name='bias')
rbf_kernel_details = elfi.GPyRegression(kernel=final_kernel, parameter_names=sim.model.parameter_names,
  bounds={'wee':(-2,2),'wei':(-2,2),'wie':(-2,2),'wii':(-2,2),'be':(-2,2),
          'bi':(-2,2),'taue':(-2,2),'taui':(-2,2),'k':(-2,2),'IH':(-2,2),
          'psi_sigma':(-2,2),'velocity':(-2,2)},noise_var=2e-02)

# Creating dictionary
init_evidence_dict={'wee':X[:,0],'wei':X[:,1],'wie':X[:,2],'wii':X[:,3],'be':X[:,4],'bi':X[:,5],'taue':X[:,6],'taui':X[:,7],'k':X[:,8],'IH':X[:,9],'psi_sigma':X[:,10],'velocity':X[:,11],'log_d':Y[:,0]}

bolfi = elfi.BOLFI(log_d, initial_evidence=init_evidence_dict, update_interval=1, target_model=rbf_kernel_details,
                   exploration_rate=10, batch_size=1, async_acq=False, seed=seed)

# Fixing noise variance
bolfi.target_model._gp.Gaussian_noise.variance.fix

# Fitting BOLFI surrogate model
post = bolfi.fit(n_evidence=9093, bar=True)

bolfi.plot_discrepancy()
bolfi.plot_gp()

eng.quit()

np.save('Z:\\Documents\\MEGMOD\\Python\\data\\bolfi_12params_distancedependentdelaysmodel_expMEGdata_09092022_9093sims_paramarray',bolfi.target_model._gp.param_array)
np.save('Z:\\Documents\\MEGMOD\\Python\\data\\bolfi_12params_distancedependentdelaysmodel_expMEGdata_09092022_9093sims_Ypred',bolfi.target_model.predict(bolfi.target_model.X, noiseless=True)[0].ravel())
np.save('Z:\\Documents\\MEGMOD\\Python\\data\\bolfi_12params_distancedependentdelaysmodel_expMEGdata_09092022_9093sims_Yobs',bolfi.target_model.Y.ravel())

# Posterior sampling
result_BOLFI = bolfi.sample(1000, warmup=500, n_chains=4, algorithm='nuts', target_prob=0.8)
print(result_BOLFI)

result_BOLFI.plot_marginals(bins=30)
result_BOLFI.plot_pairs(bins=30)

np.save('Z:\\Documents\\MEGMOD\\Python\\data\\bolfi_12params_distancedependentdelaysmodel_expMEGdata_09092022_9093sims_defthresh_SSIlogdisc_targetprob0pt8_bolfisample', result_BOLFI.samples_array)
# elfi.OutputPool(output_pool)