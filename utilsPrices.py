import numpy as np

def call_MC(matur, s_t, r, strikes):
    price = np.empty_like(strikes)
    for k, strike in enumerate(strikes):
        price[k] = np.maximum(s_t - strike, 0.).mean()
    return np.exp(- r * matur) * price

def put_MC(matur, s_t, r, strikes):
    price = np.empty_like(strikes)
    for k, strike in enumerate(strikes):
        price[k] = np.maximum(strike - s_t, 0.).mean()
    return np.exp(- r * matur) * price

def payoff_digitale(s_t, strike):
    n = s_t.shape[0] + 0.0
    sum = 0.0
    for i in range(s_t.shape[0]):
        if s_t[i] >= strike:
            sum += 1.0
    return sum / n

def option_digitale(matur, s_t, strikes, r):
    price = np.empty_like(strikes)
    for k, strike in enumerate(strikes):
        price[k] = payoff_digitale(s_t, strikes[k])
    return np.exp(- r * matur) * price

def price_spread(s_t1, s_t2, strikes, matur, r):
    price = np.empty_like(strikes)
    for k, strike in enumerate(strikes):
        price[k] = np.maximum(s_t1 - s_t2 - strike, 0.).mean()
    return np.exp(- r * matur) * price

def forward_start(t_2, s_t_1, s_t_2, strikes, r):
    price = np.empty_like(strikes)
    for k, strike in enumerate(strikes):
        price[k] = np.maximum((s_t_2 / s_t_1) - strike, 0.).mean()
    return np.exp(- r * t_2) * price


def price_asian(s_t_w, strikes, T, r):
    price = np.maximum(np.subtract.outer(s_t_w.mean(axis = 1), strikes), 0.).mean(0)
    return np.exp(- r * T) * price

def LSMC_BermudanOption(exercice_dates, spot_at_dates, payoff, r, deg_ = 2):
    
    H = payoff(spot_at_dates)
    V = np.zeros_like(H)               
    V[:, -1] = H[:, -1]
    N = exercice_dates.shape[0]
    for t in range(N - 2, 0, -1):
        validPaths = H[:, t] > 0
        
        df_ = np.exp(- r * (exercice_dates[t + 1] - exercice_dates[t]))
        rg = np.polyfit(spot_at_dates[validPaths, t], V[validPaths, t + 1] * df_, deg_)
        C = np.polyval(rg, spot_at_dates[validPaths, t])

        exercise = np.zeros(len(validPaths), dtype=bool)     # initialize
        exercise[validPaths] = H[validPaths, t] > C          # paths where it is optimal to exercise

        V[exercise, t] = H[exercise, t]                        # set V equal to H where it is optimal to exercise 
        V[exercise, t+1:] = 0.                                 # set future cash flows, for that path, equal to zero  
        discount_path = (V[:, t] == 0)                        # paths where we didn't exercise 
        V[discount_path, t] = V[discount_path, t + 1] * df_
    
    return V.mean(0)

def bermudian_puts(spot_samples, exercice_dates, strikes, df):
    result =  np.empty((strikes.shape[0], spot_samples.shape[1]))
    for k in range(strikes.shape[0]):
        payoff_ = lambda s: np.maximum(strikes[k] - s, 0.)
        result[k, :] =  df * LSMC_BermudanOption(exercice_dates, spot_samples, payoff_, 8)
    return result

def if_1(X, B):
    result = []
    for i in range(X.shape[0]):
        if X[i] < B:
            result.append(1.)
        else:
            result.append(0.)
    return np.array(result)

def payoff_barriere(matur, s_t, B, strikes, r):
    max_s_t = s_t.max(axis = 1)
    max_cond = if_1(max_s_t, B)
    s_T = s_t[:, -1]
    price = np.empty_like(strikes)
    for k in range(strikes.shape[0]):
        price[k] = (np.maximum(s_T - strikes[k], 0.) * max_cond).mean()
    return np.exp(- r * matur) * price