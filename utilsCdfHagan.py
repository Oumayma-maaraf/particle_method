import numpy as np
import scipy.stats as stats

norm  = stats.norm
def bs_call_price(s_0, matur, strikes, sigma, r):
    forwrd = s_0 * np.exp(r * matur)
    discount = np.exp(- r * matur)
    d_1 = (np.log(forwrd / strikes) + 0.5 * sigma**2 * matur) / (sigma * np.sqrt(matur))
    d_2 = d_1 - sigma * np.sqrt(matur)
    return discount * (forwrd * norm.cdf(d_1) - strikes *  norm.cdf(d_2))

def bs_put_price(s_0, matur, strikes, sigma, r):
    forwrd = s_0 * np.exp(r * matur)
    discount = np.exp(- r * matur)
    d_1 = (np.log(forwrd / strikes) + 0.5 * sigma**2 * matur) / (sigma * np.sqrt(matur))
    d_2 = d_1 - sigma * np.sqrt(matur)
    return  - discount * forwrd * norm.cdf( - d_1) + strikes * discount *  norm.cdf(-d_2)

def get_forwards(maturities, s_0 = 1., r = 3e-2):
    return s_0 * np.exp(r * maturities)

# market implied vol Hagan
def implied_vol_hagan(maturities, strikes, forwards, dtype_ = np.float64, gamma = .2, rho = -0.9, alpha = 0.2, beta = 0.5):
    out = np.empty_like(strikes, dtype = dtype_)
    
    for n, mat in enumerate(maturities):
        lamb = (gamma / alpha) * ((forwards[n])**(1 - beta))
        log_m = np.log(strikes[:, n] / forwards[n])
        out[:, n] = alpha * (1. - .5 * (1. - beta - rho * lamb) * log_m + (((1. - beta)**2 + (2. - 3 * rho**2) * lamb**2) * log_m**2) / 12.) / (forwards[n]**(1. - beta))
    return out 

def get_price_surface(s_0, maturities, strikes, vol_surface, r = 3e-2):
    out = np.empty_like(vol_surface)
    for n, matur in enumerate(maturities):
        out[:, n] = bs_call_price(s_0, matur, strikes[:, n], vol_surface[:, n], r)
    return out

def get_price_surface_put(s_0, maturities, strikes, vol_surface, r = 3e-2):
    out = np.empty_like(vol_surface)
    for n, matur in enumerate(maturities):
        out[:, n] = bs_put_price(s_0, matur, strikes[:, n], vol_surface[:, n], r)
    return out

def get_market_cdf(s_0, strikes, maturities, forwards, eps = 1e-6, deriv = 'centered', r = 3e-2):
    #dtype_ = 'complex_' if deriv == 'complex' else np.float64
    out = np.empty_like(strikes)
    
    if (deriv == 'centered'):
            vols_surf_peps = implied_vol_hagan(maturities, strikes +  eps, forwards = forwards)
            vols_surf_meps= implied_vol_hagan(maturities, strikes - eps, forwards = forwards)

    if (deriv =='complex'): 
            vols_surf_comp  = implied_vol_hagan(maturities, strikes +  1j * eps, forwards = forwards, dtype_= 'complex_')
    

    for n, mat in enumerate(maturities):
        if (deriv =="centered"):
                price_peps = bs_call_price(s_0, mat, strikes[:, n] + eps, vols_surf_peps[:, n], r)
                price_meps = bs_call_price(s_0, mat, strikes[:, n] - eps, vols_surf_meps[:, n], r)
                derivative = 0.5 * (price_peps - price_meps) / eps
        if (deriv =='complex'):
                price_comp = bs_call_price(s_0, mat, strikes[:, n] + 1j * eps, vols_surf_comp[:, n], r)
                derivative = np.imag(price_comp) / eps

        out[:, n] = 1. + np.exp(r * mat) * derivative
    return out

def get_market_density(s_0, strikes, maturities, forwards, vols, eps = 1e-6, deriv = 'centered', r = 3e-2):
    #dtype_ = 'complex_' if deriv == 'complex' else np.float64
    out = np.empty_like(strikes)
    
    if (deriv == 'centered'):
            vols_surf_peps = implied_vol_hagan(maturities, strikes +  eps, forwards = forwards)
            vols_surf_meps= implied_vol_hagan(maturities, strikes - eps, forwards = forwards)

    

    for n, mat in enumerate(maturities):
        if (deriv =="centered"):
                price_peps = bs_call_price(s_0, mat, strikes[:, n] + eps, vols_surf_peps[:, n], r)
                price_meps = bs_call_price(s_0, mat, strikes[:, n] - eps, vols_surf_meps[:, n], r)
                price_seps = bs_call_price(s_0, mat, strikes[:, n] , vols[:, n], r)
                derivative =  (price_peps - 2 * price_seps + price_meps) / (eps * eps)
        
        out[:, n] = np.exp(r * mat) * derivative
    return out