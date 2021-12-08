import matlab.engine
import numpy as np
import multiprocessing as mp

eng = None


def get_matlab():
    global eng
    if eng is None:
        print(mp.current_process(), "Creating new matlab")
        eng = matlab.engine.start_matlab()
    return eng


# Function to estimate percentile values
def percentile_values(x, pctval=0):
    C = np.percentile(x, pctval, axis=0)
    return C


# Defining characteristics of simulator
def simulator(wee, wei, wie, wii, be, bi, taue, taui, k, IH, psi_sigma, batch_size=1, random_state=None):
    eng = get_matlab()

    print(f"\n\nSTARTING NEW SIMULATION")
    print(f"With parameters: {wee.item(), wei.item(), wie.item(), wii.item(), be.item(), bi.item(),taue.item(),taui.item(), k.item(), IH.item(), psi_sigma.item()}")

    ret = eng.simulate_WCnetnoisy(wee.item(), wei.item(), wie.item(), wii.item(), be.item(), bi.item(), taue.item(), taui.item(), k.item(), IH.item(), psi_sigma.item())
    return np.array(ret)
