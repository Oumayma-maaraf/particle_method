import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize
from scipy.optimize import fsolve, newton, bisect
from scipy.interpolate import interp1d, splrep, sproot, PPoly, lagrange
from numpy.polynomial.polynomial import Polynomial


norm = stats.norm
dt = 1. / 252.
identity = lambda x: x


class utilsClv:
    def __init__(self, mdl_params, rng)-> None:
        self.mdl_params = mdl_params
        self.rng = rng

    def update_mdl_params(self, mdl_params_):
        self.mdl_params = mdl_params_
    
    def generate_traj(self, x_0, matur, n_samples = 10_000, freq = 1.):
        kappa, theta, sigma = self.mdl_params.values()
        time_step = freq * dt
        disc_steps = int(matur / time_step)

        x = np.empty((n_samples, disc_steps + 1))
        x[:, 0] = x_0

        for n in range(1, disc_steps + 1):
            z = self.rng(n_samples)
            x[:, n] = x[:, n - 1] + kappa * (theta - x[:, n - 1]) * time_step + sigma * np.sqrt(time_step) * z
        
        return x
    
    def get_gauss_collocations(self, deg):
        return np.sqrt(2) * np.polynomial.hermite.hermgauss(deg)[0]
    
    def get_colloc_points(self, x_0, time_grid, order):
        kappa, theta, sigma = self.mdl_params.values()
        mean_t = x_0 * np.exp(-kappa * time_grid) + theta * (1. - np.exp(-kappa * time_grid))
        sd_t = np.sqrt(.5 * sigma**2 * (1. - np.exp(-2 * kappa * time_grid)) / kappa)

        gauss_pnts = self.get_gauss_collocations(order)
        return mean_t[:, None] + np.outer(sd_t, gauss_pnts)
    
    def get_colloc_points_x_t(self, t, x_0, order):
        kappa, theta, sigma = self.mdl_params.values()
        mean_t = x_0 * np.exp(-kappa * t) + theta * (1. - np.exp(-kappa * t))
        sd_t = np.sqrt(.5 * sigma**2 * (1. - np.exp(-2 * kappa * t)) / kappa)

        gauss_pnts = self.get_gauss_collocations(order)
        return mean_t.mean() + sd_t * gauss_pnts

    def get_spots(self, 
                  x_0,
                  maturities, 
                  strikes, 
                  colloc_points_x,
                  market_cdf, 
                  order = 4, 
                  solver_kwargs = {"x0": .205},
                  solver_name = 'brent'):
        
        out = np.empty_like(colloc_points_x)
        kappa, theta, sigma = self.mdl_params.values()
        if (solver_name == 'bisect'):
            wrap = identity
            solver_ = lambda func: optimize.bisect(func, a = 1e-3, b = 2., maxiter = 150, xtol = 1e-10)
        if (solver_name =='brent'):
            wrap = identity
            solver_ = lambda func: fsolve(func = func, x0 = 0.205)
        if (solver_name =='minimize'):
            wrap = np.abs
            bnds = optimize.Bounds(1e-3, 2.5)
            solver_ = lambda func : optimize.minimize(func, x0 = 0.1, method = 'L-BFGS-B', 
                                                      bounds=bnds, options={'maxiter': 200, 
                                                                            'ftol': 1e-10, 
                                                                            'gtol': 1e-15}).x
            
        for n, maturitie in enumerate(maturities):
            cdf_func = interp1d(strikes[:, n], market_cdf[:, n], kind = 'cubic', fill_value='extrapolate')

            mean_t = x_0 * np.exp(-kappa * maturitie) + theta * (1. - np.exp(-kappa * maturitie))
            sd_t = np.sqrt(.5 * sigma**2 * (1. - np.exp(-2 * kappa * maturitie)) / kappa)

            cdf_X = lambda x: norm.cdf(x, loc = mean_t.mean(), scale = sd_t)
        
            for i in range(order):
                xi_n = cdf_X(colloc_points_x[n][i])
                func_ = lambda x: wrap(cdf_func(x) - xi_n)
                out[n, i] = solver_(func_)
        return out            

    def search_index(self, y, X):
        X = np.array(X)
        for i in range(X.shape[0] - 1):
            if (y == X[i]):
                return i
            if (y < X[i + 1] and y > X[i]):
                return i
        if (y == X[-1]):
            return X.shape[0] - 1
        if (y < X[-1] and y > X[-2]):
            return X.shape[0] - 2

    def get_s_t(self, t, maturities, s):    
        index_i = self.search_index(t, maturities)
        if (t in maturities):
            return s[index_i]
        else:
            return (s[index_i] + ((s[index_i + 1] - s[index_i]) * (t - maturities[index_i]) / (maturities[index_i + 1] - maturities[index_i])))

    def get_g_func(self, t, maturities, s, x_0, order):
        s_t = self.get_s_t(t, maturities, s)
        x_t = self.get_colloc_points_x_t(t, x_0, order)
        poly = lagrange(x_t, s_t)
        return Polynomial(poly.coef[::-1])
    
    def auto_clv(self, t, t_grid, maturities, s, x_0, x_t, order):
        idx_t = np.where(t_grid == t)[0][0]
        g_t = self.get_g_func(t, maturities, s, x_0, order)
        s_t = g_t(x_t[:, idx_t])
        return s_t