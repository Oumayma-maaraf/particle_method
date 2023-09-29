import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize
from scipy.optimize import fsolve, newton, bisect
from scipy.interpolate import interp1d, splrep, sproot, PPoly

class utilsDef:
              
        def __init__(self, nbStrike, sizeSample) -> None:
            self.nbStrike = nbStrike
            self.sizeSample = sizeSample
        
        

        def findStrikes(self, n, cdfInverse):
            strikes = []
            strikes.append(0.0)
            for i in range(1, n - 1):
                value = i / (n - 1.)
                if ((value < 1.) & (value > 0.0)):
                    strikes.append(cdfInverse[i - 1])
                if (value >= 1.):
                    strikes.append(np.inf)
                    break
                if (value == 0.):
                    strikes.append(0.0)
            strikes.append(np.inf)
            return np.array(strikes)
        
        
        def findStrikes_norm(self, n, average, stdDev):
            strikes = []
            strikes.append(0.0)
            for i in range(1, n - 1):
                value = i / (n - 1.)
                if ((value < 1.) & (value > 0.0)):
                    strikes.append(np.exp(stats.norm.ppf(value, loc = average, scale = stdDev)))
                if (value >= 1.):
                    strikes.append(np.inf)
                    break
                if (value == 0.):
                    strikes.append(0.0)
            strikes.append(np.inf)
            return np.array(strikes)

        def firstDeformationBucket(self, s_t, strikes, j, CardBucket):
            nbStrike = self.nbStrike
            result = np.empty(CardBucket)
            if (j != nbStrike - 2):
                result[0] = strikes[j]
                for i in range(1, CardBucket - 1):
                    result[i] = result[0] + ((strikes[j + 1] - strikes[j]) / (s_t[(j + 1) * CardBucket - 1] - s_t[j * CardBucket])) * (s_t[i + j * CardBucket] - s_t[j * CardBucket])
                result[CardBucket - 1] = strikes[j + 1]
            else:
                for i in range(j * CardBucket, (j + 1) * CardBucket):
                    result[i - j * CardBucket] = np.maximum(strikes[j], s_t[i])
            return result
        
        def firstDeformation(self, s_t, strikes, CardBucket):
            sizeSample = self.sizeSample
            nbStrike = self.nbStrike
            result = np.empty(sizeSample)
            for j in range(0, nbStrike - 1):
                deformation = self.firstDeformationBucket(s_t, strikes, j, CardBucket)
                result[j * CardBucket: (j + 1) * CardBucket] = deformation
            return result
        
        def secondDeformationBucket(self, s_tilde, alpha_i, pricetheo, pricetilde, strikes, i):
            cardBucket = s_tilde.shape[0]
            result = np.empty(cardBucket)
            if (pricetilde < pricetheo):
                for j in range(0, cardBucket):
                    result[j] = alpha_i * s_tilde[j] + (1 - alpha_i) * strikes[i + 1]
            else:
                for j in range(0, cardBucket):
                    result[j] = alpha_i * s_tilde[j] + (1 - alpha_i) * strikes[i]
            return result
        
        def secondDeformation(self, s_t, a_i, pricetheo, pricetilde, strikes, cardBucket):
            result = np.empty(s_t.shape[0])
            nbstrike = self.nbStrike
            m = a_i.shape[0]

            s_tilde_n = self.firstDeformationBucket(s_t, strikes, nbstrike - 2, cardBucket)
            strikes[nbstrike - 1] = np.max(s_tilde_n)
            
            for i in range(0, m):
                s_tilde = self.firstDeformationBucket(s_t, strikes, i, cardBucket)
                deformation = self.secondDeformationBucket(s_tilde, a_i[i], pricetheo[i], pricetilde[i], strikes, i)
                result[i * cardBucket: (i + 1) * cardBucket] = deformation
            
            return result
        
        def alpha_i(self, s_t, strikes, pricetheo, r, Maturity):
            sizeSample = self.sizeSample
            nbStrike = self.nbStrike
            cardBucket = int(sizeSample / (nbStrike - 1))
            alpha_i = np.empty(nbStrike - 1)

            s_tilde_n = self.firstDeformationBucket(s_t, strikes, nbStrike - 2, cardBucket)
            strikes[nbStrike - 1] = np.maximum(pricetheo[nbStrike - 2] * np.exp(r * Maturity) * sizeSample + strikes[nbStrike - 2], np.max(s_tilde_n))

            pricetilde = np.empty(nbStrike - 1)
            sumterm = 0.0

            for i in range(nbStrike - 2, -1, -1):
                s_tilde = self.firstDeformationBucket(s_t, strikes, i, cardBucket)
                pricetilde[i] = np.exp(-r * Maturity) * (s_tilde.sum() + sumterm - cardBucket * (nbStrike - 1 - i) * strikes[i]) / sizeSample

                if (pricetilde[i] < pricetheo[i]):
                    alpha_i[i] = (sizeSample * pricetheo[i] * np.exp(r * Maturity) + cardBucket * ((nbStrike - 1 - i) * strikes[i] - strikes[i + 1]) - sumterm) / (s_tilde.sum() - cardBucket * strikes[i + 1])

                else:
                    alpha_i[i] = (sizeSample * pricetheo[i] * np.exp(r * Maturity) + cardBucket * ((nbStrike - 1 - i) * strikes[i] - strikes[i]) - sumterm) / (s_tilde.sum() - cardBucket * strikes[i])

                sumterm += self.secondDeformationBucket(s_tilde, alpha_i[i], pricetheo[i], pricetilde[i], strikes, i).sum()

            return alpha_i, pricetilde  
        
        
        def inversecdf(self, n, solver_name, cdf_f, solver_kwargs):
            identity = lambda x: x
            if solver_name =='bisect':
                wrap = identity
                solver_ = lambda func: optimize.bisect(func, **solver_kwargs)
            if solver_name =='brent':
                wrap = identity
                solver_ = lambda func: fsolve(func = func, **solver_kwargs)
            if solver_name == 'minimize':
                wrap = np.abs
                bnds = optimize.Bounds(1e-3, 2.5)
                solver_ = lambda func : optimize.minimize(func, bounds = bnds,**solver_kwargs).x

            values = np.empty(n - 2)
            for i in range(1, n - 1):
                values[i- 1] = i / (n - 1.)
            out = np.empty(n - 2)
            for i in range(0, n - 2):
                func_ = lambda x: wrap(cdf_f(x) - values[i])
                out[i] = solver_(func_)
            return out 
      