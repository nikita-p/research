import numpy as np
import pandas as pd
import iminuit
from iminuit import Minuit
import matplotlib.pyplot as plt
from scipy.special import erfc, expit
import scipy.integrate as integrate
import uproot
import re

class PhotonEff():
    def __init__(self):
        self.df = None
        self.histos = None
        self.fit_name = None
        self.fit_results = dict()
        self.good_fits = dict()
        self.par_names = ['$\mu$', '$\sigma$', 'c', 'N']
    def __variance(self, k, n):
        return (k+1)*(k+2)/(n+2)/(n+3) - (k+1)**2/(n+2)**2
    def __chi2(self, mu, s, c, N):
        xedges, points, errs = self.get_histo_by_name(self.fit_name)
        return np.sum( np.square((points - self.sigFunc(xedges, mu, s, c, N))/errs) )
    def sigFunc(self, x, mu, s, c, N):
        return N*( 1 - expit((x-mu)/s) + c)
    def get_signal_by_name(self, name, x):
        return self.sigFunc( x, *(self.fit_results[name][0]) )
    def get_histo_by_name(self, name):
        xedges = self.histos[name][0]
        points = self.histos[name][1]
        errs = self.histos[name][2]
        return (xedges, points, errs)
    def fit_histo_by_name(self, name):
        self.fit_name = name
        parameters = {
            'mu': 0.16, #float(name)*2e-3*0.157-0.1547, #
            's': 1/250, #float(name)*2e-3*0.0487-0.0493 , #
            'c':0.01,
            'N':0.35,
        }
        limits = {
            'limit_mu' : [0.005, 1], #[parameters['mu'], parameters['mu']], #
            'limit_s' : [0,1], #[parameters['s'], parameters['s']], #
            'limit_c' : [0,1],
            'limit_N' : [0,1],
        }
        errors = {
            'error_mu' : 0.001,
            'error_s' : 1/250,
            'error_c' : 0.01,
            'error_N' : 0.05,
        }
        m = Minuit( self.__chi2, **parameters, **limits, **errors, errordef=1)
        m.migrad()
        ndf = len(self.get_histo_by_name(name)[0]) - len(m.np_values())
        print(m.migrad_ok())
        self.good_fits[name] = m.migrad_ok()
        return ( m.np_values(), m.np_errors(), m.fval, ndf )
    def show_histo_by_name(self, name):
        xedges, points, errs = self.get_histo_by_name(name)
        plt.errorbar(xedges, points, yerr=errs, fmt='o', ms=10, lw=2)
        x = np.linspace(0, xedges.max(), 100)
        text_string = "$\\sqrt{s}$: " + "{0:.2f} GeV\n".format( float(re.findall(r'_(\d+.?\d*)_',name)[0])*2e-3)
        if name in self.fit_results:
            plt.plot(x,  self.sigFunc(x, *(self.fit_results[name][0]) ), lw=2 )
            chi2, ndf = self.fit_results[name][2], self.fit_results[name][3]
            text_string += f'$\chi^2$/ndf = {chi2:.2f} / {ndf}\n'
            for i in range(len(self.par_names)):
                par, val, err = self.par_names[i], self.fit_results[name][0][i], self.fit_results[name][1][i]
                text_string += f'{par} = {val:.2e} $\pm$ {err:.2e}\n'
        plt.text(0.5*xedges.max(), 0.55*points.max(), text_string)
        plt.xlabel('Photon energy ($E_{\gamma}$), GeV')
        plt.ylabel('$\\varepsilon_{reg} (x)$')
        plt.xlim(0, xedges.max())
        return     
    def open(self, *files, min_entries=None): #min_entries drops files with small events
        f = [_ for _ in files ]
        it = uproot.pandas.iterate(f, treepath='mc_passed', reportpath=True)
        self.df = pd.concat([j.assign(name=re.findall(r'/([^/]*.root)',i)[0]) for i,j in it])
        if min_entries is not None:
            self.df = self.df.set_index('name').\
            loc[self.df.groupby('name').agg({'passed_cuts':'count'}).\
                query('passed_cuts>@min_entries').index].reset_index()
    def get_histo(self, df, bins=81):
        xedges = np.linspace(0, df['ph_energy'].max()/2., bins)
        hpassed, _ = np.histogram(df.query('passed_cuts==True')['ph_energy'], bins=xedges)
        hall, _ = np.histogram(df['ph_energy'], bins=xedges)
        xedges = (xedges[1:] + xedges[:-1])/2
        points = (hpassed+1)/(hall+2)
        errs = np.sqrt( self.__variance(hpassed, hall) )
        return np.array( [xedges, points, errs] )
    def get_histos(self):
        df = self.df
        self.histos = df.groupby('name').apply(lambda x: self.get_histo(x) )
        return self.histos
    def fit_histos(self):
        for name in self.histos.index:
            self.fit_results[name] = self.fit_histo_by_name(name)
        return
    def get_names(self):
        return self.histos.index
    def get_efficiency(self):
        x, y, yerr = [], [], []
        for name in self.histos.index:
            x.append( float(re.findall(r'_(\d+.?\d*)_',name)[0])*2e-3 )
            fr = self.fit_results[name]
            y.append( self.sigFunc(0, *fr[0]) )
            m = ( self.sigFunc(0, *(fr[0]+fr[1])) - self.sigFunc(0, *(fr[0]-fr[1])) )/2
            yerr.append( m )
        return x, y, yerr

class RadCor:
    alpha = 7.297e-3
    me = 0.511 #MeV
    
    def __init__(self):
        return
    def l(self, E, dE):
        return np.log(E/dE)
    def L(self, s):
        me = self.me
        return np.log(s/(me**2))
    def beta(self, s):
        a = self.alpha
        p = np.pi
        return (2*a/p)*(self.L(s) - 1)
    
    def soft_terms_rad(self, E, dE):
        a = self.alpha
        p = np.pi
        s = 4*(E**2)
        l = self.l(E, dE)
        L = self.L(s)
        
        s1 = -2*l*(L-1) + 3*L/2 + (p**2)/3 - 2
        s2 = (1/2)*((-2*l*(L-1))**2)
        s3 = (3*L/2 + (p**2)/3 - 2)*(-2*l*(L-1))
        s4 = (L**2)*( -l/3 + 11/8 - (p**2)/3 )
        s5 = L*(2*(l**2)/3 + 10*l/9 - 203/48 + 11*(p**2)/12 + 3*1.202)
        s6 = -(4/9)*(l**3) - (10/9)*(l**2) - (2/9)*(28/3 - p**2)*l
        
        result = 1 + (a/p)*s1 + ((a/p)**2)*(s2 + s3 + s4 + s5 + s6)
        return result
    def soft_terms_cs(self, E, cs, dE):
        return self.soft_terms_rad(E, dE)/cs
        
    def F(self, x, s):
        a = self.alpha
        b = self.beta(s)
        p = np.pi
        L = self.L(s)
        m = self.me
        E = np.sqrt(s)/4
        
        s1 = (a/p)*( (p**2)/3 - 1/2 ) + 3*b/4
        s2 = (-(b**2)/24)*(L/3 + 2*(p**2) - 37/4 )
        s3 = -b*(1-x/2)
        s4 = 4*(2-x)*np.log(1/x)
        s5 = (1/x)*(1+3*((1-x)**2))*np.log(1/(1-x))
        s6 = - 6 + x
        
        s7 = 0 if x<(2*m/E) else (1/(6*x))*((x - 2*m/E)**b)*\
            ((np.log(s*(x**2)/(m**2)) - 5/3)**2)*\
            (2 - 2*x + x**2 + (b/3)*(np.log(s*(x**2)/(m**2)) - 5/3))
        s8 = 0 if x<(2*m/E) else ((L**2)/2)*((2/3)*((1-(1-x)**3)/(1-x)) -\
            (2-x)*np.log(1/(1-x)) + x/2 )
        
        
        result = b*(x**(b-1))*( 1 + s1 + s2 ) + s3 + \
        (1/8)*(b**2)*(s4 + s5 + s6) + \
        ((a/p)**2)*(s7 + s8)
        return result
    
    def F_Integral(self, e, cs, e_beam, params, Xmax=1):
        p = PhotonEff()
        s = 4*(e_beam**2)
        s_cs = 4*(e**2)
        if not( np.all(np.diff(s_cs) > 0) ):
            print('Problem')
        return integrate.quad( lambda x: self.F(x, s)*np.interp(s*(1-x), s_cs, cs)*\
                              p.sigFunc(np.sqrt(s)*(1 - np.sqrt(1-x))*1e-3, *params), 
                              0., Xmax, points=[0, 1], 
                              limit=5000, epsrel=0.0001)
    def F_Radcor(self, e, cs, e_beam, params, Xmax=1):
        integral = self.F_Integral( e, cs, e_beam, params, Xmax)
        return ( integral[0]/np.interp(e_beam, e, cs), integral[1]/np.interp(e_beam, e, cs) )

class MDVM():
    def __init__(self):
        self.ALPHA = 7.297352e-3
        self.C = 0.3893793656e12 #(MeV)^2 * nb

        self.mPhi = 1019.464
        self.mRho = 775.26
        self.mOmg = 782.65

        self.mK0 = 497.611
        self.mP0 = 135.
        self.mKC = 493.677
        self.mPC = 139.57
        self.mKstar = 891.76

        self.w0Phi = 4.247
        self.w0Rho = 149.1
        self.w0Omg = 8.49
    
    def BETA(self, s, M_K): #бета
        E = np.sqrt(s)/2.
        P = np.where( E<M_K, 0, np.sqrt( E**2 - M_K**2 ) )
        return P/E
    def PV2_diff(self, s, M, mi, mj): #фазовый объём распада на 2 разные частицы https://arxiv.org/pdf/hep-ph/9609216.pdf (2.6)
        q0 = np.sqrt( (M**2 - (mi - mj)**2)*(M**2 - (mi + mj)**2) )/(2*M)
        E = np.sqrt(s)
        q_temp = np.where( E <= (mi + mj), 0, (E**2 - (mi - mj)**2)*(E**2 - (mi + mj)**2) )
        q1 = np.sqrt( q_temp )/(2*E)
        return (q1/q0)**3
    def PV2(self, s, M, Mn): #фазовый объём распада на 2 одинаковые частицы
        w = np.where(s <= 4*Mn*Mn, 0, np.power( (s - 4*Mn**2)/(M**2 - 4*Mn**2), 3./2 )*(M**2)/s )
        return w
#     def FAS3(self, s):
#         twoe = (np.sqrt(np.array(s))*1e-3).reshape((-1,1)) #GeV        
#         A = np.array([-6.12394622, 25.0341405, -34.1311022, 15.5413717])
#         B = np.array([5.29354148, -7.90990714, -2.26007613, 5.21453902])
#         LM = 1.1
#         E = np.array([twoe**i for i in range(A.shape[0])])
#         POL = (np.where(twoe >= LM, B, A)*E).sum(axis=(0,2))
#         Fval = (POL/0.393728*(0.00749/0.0361478))
#         return Fval
    def FAS3(self, s):
        e = np.sqrt(s)*1e-3
        f1 = 5.196 + 59.17*(e-1) + 227.7 * (e-1)**2 + 147 * (e-1)**3 - 998 * (e-1)**4 - 1712 * (e-1)**5
        f2 = 5.196 + 80*(e-1) + 200 * (e-1)**2 + 590 * (e-1)**3 - 510 * (e-1)**4 + 220 * (e-1)**5
        return np.where(e<1, f1, f2)
    def FAS3_RhoPiPi(self, s):
        e = np.sqrt(s)*1e-3
        f1 = 1.9e-3 + 2.68e-2*(e-1.2) + 2.446e-1 * (e-1.2)**2 + 3.1487 * (e-1.2)**3 + 23.3131 * (e-1.2)**4 + 59.7669 *  (e-1.2)**5
        f2 = 1.9e-3 + 2.33e-2*(e-1.2) + 6.65e-2 * (e-1.2)**2 - 3.84e-2 * (e-1.2)**3 + 2.36e-2 * (e-1.2)**4 - 6.5e-3 * (e-1.2)**5
        return np.where(e<1.2, f1, f2)
    def FAS3_OmgPiPi(self, s):
        e = np.sqrt(s)*1e-3
        f1 = 1.7e-3 + 2.61e-2*(e-1.2) + 2.67e-1 * (e-1.2)**2 + 3.61199 * (e-1.2)**3 + 27.6 * (e-1.2)**4 + 73.6433 * (e-1.2)**5
        f2 = 1.7e-3 + 2.18e-2*(e-1.2) + 6.89e-2 * (e-1.2)**2 - 4.52e-2 * (e-1.2)**3 + 3.25e-2 * (e-1.2)**4 - 1.09e-2 * (e-1.2)**5
        return np.where(e<1.2, f1, f2)
    def PV3(self, s, MX, m1, m2, m3): #фазовый объём (почти: он типа нормирован) распада на 3 частицы
        #WARNING. TRY TO FIND THE RIGHT VARIANT
        pv = self.FAS3(s)
        pv0 = self.FAS3(MX**2)
        return pv/pv0
    def PVG(self, s, MX, Mn): #распад частицы M на фотон и частицу Mn
        pv  = ((s - Mn**2)/(2*np.sqrt(s)))**3
        pv0 = ((MX**2 - Mn**2)/(2*MX))**3
        return pv/pv0
    def P_RhoPiPi(self, s, MX):
        return np.where( np.sqrt(s)>self.mRho+2*self.mP0, self.FAS3_RhoPiPi(s)/self.FAS3_RhoPiPi(MX**2), 0 )
    def P_OmgPiPi(self, s, MX):
        return np.where( np.sqrt(s)>self.mOmg+2*self.mP0, self.FAS3_OmgPiPi(s)/self.FAS3_OmgPiPi(MX**2), 0 )
    
    def WOmg(self, s, W0, MX):    
        Br_3Pi = 0.892    
        Br_Pi0G = 0.084
        Br_2Pi = 0.0153
        ost = 1 - Br_Pi0G - Br_2Pi - Br_3Pi;
        
        mPC = self.mPC
        mP0 = self.mP0

        W = W0 * ((Br_3Pi + ost) * self.PV3(s, MX, mPC, mPC, mP0) + \
                  Br_Pi0G * self.PVG(s, MX, mP0) + Br_2Pi * self.PV2(s, MX, mPC));
        return W
    def WPhi(self, s, W0, MX):
        Br_KC = 0.492
        Br_KN = 0.34
        Br_3Pi = 0.1524
        Br_EG = 0.01303
        ost = 1 - Br_KC - Br_KN - Br_3Pi - Br_EG

        mEta = 547.862
        mKC = self.mKC
        mK0 = self.mK0
        mPC = self.mPC
        mP0 = self.mP0

        W = W0*((Br_KC + ost)* self.PV2(s, MX, mKC) + Br_KN* self.PV2(s, MX, mK0) + \
                Br_3Pi * self.PV3(s, MX, mPC, mPC, mP0) + Br_EG * self.PVG(s, MX, mEta));
        return W
    def WPhi1680(self, s, W0, MX): #fully KK* decay
        W = W0*self.PV2_diff(s, MX, self.mKC, self.mKstar)
        return W
    
    def WRho(self, s, W0, MX):
        return W0 * self.PV2(s, MX, self.mPC)
    def WRhoX(self, s, W0, MX):
        return W0 * self.PV2_diff(s, MX, self.mP0, self.mOmg)
    def WOmgX(self, s, W0, MX):
        return W0 * self.PV2_diff(s, MX, self.mRho, self.mP0) #self.PV3(s, MX, self.mPC, self.mPC, self.mP0)
    def WRhoN(self, s, W0, MX):
        return W0 * self.P_RhoPiPi(s, MX)
    def WOmgN(self, s, W0, MX):
        return W0 * self.P_OmgPiPi(s, MX)
    def WPhiX(self, s, W0, MX):
        return W0 * self.PV2(s, MX, self.mKC)
    
    def BW(self, s, MX, WX0, WX): #функция Брейта-Вигнера
        bw = (MX**2)/( MX**2 - s - 1j*np.sqrt(s)*WX(s, WX0, MX) )#MX
        return bw
    def BW_RhoX(self, s, mX, wX):
        return self.BW(s, mX, wX, self.WRhoX)
    def BW_OmgX(self, s, mX, wX):
        return self.BW(s, mX, wX, self.WOmgX)
    def BW_PhiX(self, s, mX, wX):
        return self.BW(s, mX, wX, self.WPhiX)
    def BW_RhoN(self, s, mX, wX):
        return self.BW(s, mX, wX, self.WRhoN)
    def BW_OmgN(self, s, mX, wX):
        return self.BW(s, mX, wX, self.WOmgN)
    
    def BW_Rho(self, s):
        return self.BW(s, self.mRho, self.w0Rho, self.WRho)
    def BW_Omg(self, s):
        return self.BW(s, self.mOmg, self.w0Omg, self.WOmg)
    def BW_Phi(self, s):
        return self.BW(s, self.mPhi, self.w0Phi, self.WPhi)
    def BW_Rho1(self, s, m=1465, g=400):
        return self.BW_RhoX(s, m, g)
    def BW_Omg1(self, s, m=1420, g=220):
        return self.BW_OmgX(s, m, g)
    def BW_Rho2(self, s): #no found in PDG
        return self.BW_RhoX(s, 1574, 234)
    def BW_Omg2(self, s, m=1670, g=315):
        return self.BW_OmgN(s, m, g)
    def BW_Phi1(self, s, m=1673, g=182):
        return self.BW(s, m, g, self.WPhi1680)
    def BW_Rho3(self, s, m=1720, g=250):
        return self.BW_RhoN(s, m, g)
    def BW_Rho4(self, s, m=1880, g=160): #no found in PDG
        return self.BW_RhoX(s, m, g)
    def BW_Rho5(self, s, m=2134, g=343): #no found in PDG
        return self.BW_RhoX(s, m, g)
    def BW_Phi2(self, s, mPhi2=2198, gPhi2=71):
        return self.BW_PhiX(s, mPhi2, gPhi2)#BESIII:2239.2, 139.8 #PDG: 2198, 71 
    
    #PUBLIC
    def F0(self, x, KR, KO, KP, n, mode=False): #формфактор, нулевое приближение; mode: 0 - short/long; 1 - charged;
        s = (x*1e3)**2 
        F = KR * self.BW_Rho(s) + KO * self.BW_Omg(s) + n * KP * self.BW_Phi(s)
        return F
    def F1(self, x, par, mode=False): #формфактор с учётом omega(1400)
        if(len(par)<2):
            par = [0.9802, 1.1787, 1.1659, 1.0397, -0.0699, -0.0841, -0.1617, -0.0851, -0.1474, 0.0151, 0.0334, 0.2793, 1515.0, 248.7102, 1449.0932, 280.1714, 1671.4159, 181.8241, 1700.0, 449.9986, 1610.0001, 366.2361, 1958.3253, 219.2568, 1996.7143, 336.7568, 2210.1137, 107.5767]
        n = 1 if mode else par[0] #1.027
        s = (x*1e3)**2
        CR = np.array([par[1], par[4], par[7], par[10]])
        CO = np.array([par[2], par[5], par[8], par[11]])
        CP = np.array([par[3], par[6], par[9]])
        
        KR = CR/2. if mode else -CR/2.
        KO = CO/6.
        KP = CP/3. #здесь уже нет дополнительного параметра n, как в F0 (в соотв. с моделью)
        
        F1 = self.F0(x, KR[0], KO[0], KP[0], n, mode)
        F1 += KR[1] * self.BW_Rho1(s, par[12], par[13]) + KR[2] * self.BW_Rho3(s, par[18], par[19]) + KR[3] * self.BW_Rho3(s, par[22], par[23])
        F1 += KO[1] * self.BW_Omg1(s, par[14], par[15]) + KO[2] * self.BW_Omg2(s, par[20], par[21]) + KO[3] * self.BW_Omg2(s, par[24], par[25]) 
        F1 += KP[1] * self.BW_Phi1(s, par[16], par[17]) + KP[2] * self.BW_Phi2(s, par[26], par[27])
        return F1
        
    def Cross_Section(self, x, par, mode=False):
        s = (x*1e3)**2
        fabs2 = np.abs( self.F1(x, par, mode) )**2
        M_K = self.mKC if mode else self.mK0
        constant = (np.pi/3.) * (self.ALPHA**2) * self.C
        cs = np.where( x<0.4976*2, 0, constant * (self.BETA(s, M_K)**3) * fabs2 / s  )
        return cs
    def Cross_Section_Neutral(self, x, par):
        return self.Cross_Section(x, par, False)
    def Cross_Section_Charged(self, x, par):
        return self.Cross_Section(x, par, True)    
    def Cross_Section2Formfacror(self, x, cs, mode=False):
        s = (x*1e3)**2
        M_K = self.mKC if mode else self.mK0
        constant = (np.pi/3.) * (self.ALPHA**2) * self.C
        return cs*s/(constant * ( self.BETA(s, M_K)**3 ) )
        
    
def get_KSKL_phi(file='../../formfactor/data/k0k0_koz.dat'):
    #KSKL на phi(1020)
    cs_k0k0_phi = pd.read_csv(file, delimiter='\t', \
                              names=['index', 'energy', 'energy_err', 'cs', 'cs_err'])
    cs_k0k0_phi['energy'] *= 1e-3
    cs_k0k0_phi['energy_err'] *= 1e-3
    cs_k0k0_phi['cs_err'] = np.sqrt( cs_k0k0_phi.cs_err**2 + (0.018*cs_k0k0_phi.cs)**2 ) #учёт систематики: 1.8% из статьи
    cs_k0k0_phi.drop('index', axis=1, inplace=True)
    return cs_k0k0_phi

def get_KPKM_phi(file='../../formfactor/data/k+k-_koz.dat'):
    #K+K- на phi(1020)
    cs_kpkm_phi = pd.read_csv(file, delimiter='\t', \
                              names=['energy', 'energy_err', 'cs', 'cs_err'])
    cs_kpkm_phi['energy'] *= 1e-3
    cs_kpkm_phi['cs_err'] = np.sqrt( cs_kpkm_phi.cs_err**2 + (0.02*cs_kpkm_phi.cs)**2 )#учёт систематики: мин.2% из статьи
    cs_kpkm_phi['energy_err'] *= 1e-3
    return cs_kpkm_phi

def get_KPKM_up(file='../../formfactor/data/k+k-.dat'):
    #K+K- на больших энергиях
    cs_kpkm = pd.read_csv(file, \
                              names=['energy', 'cs', 'stat', 'syst'])
    cs_kpkm['cs_err'] = np.sqrt(cs_kpkm.stat**2 + cs_kpkm.syst**2)
    cs_kpkm.drop(['stat','syst'], axis=1, inplace=True)
    cs_kpkm['energy_err'] = 0.
    return cs_kpkm

def get_KPKM_up2(file='../../formfactor/data/k+k-_bes.dat'):
    cs_kpkm = pd.read_csv(file)
    cs_kpkm['cs_err'] = np.sqrt( cs_kpkm.cs_stat**2 + cs_kpkm.cs_sys**2 )
    cs_kpkm['cs'] *= 1e-3
    cs_kpkm['cs_err'] *= 1e-3
    cs_kpkm['energy_err'] = 0.
    return cs_kpkm[['energy','energy_err','cs','cs_err']]

def get_KSKL_babar(file='../../formfactor/data/kskl_babar.dat'):
    cs_kpkm = pd.read_csv(file)
    cs_kpkm['energy_err'] = 0.
    return cs_kpkm[['energy','energy_err','cs','cs_err']]

def get_KPKM_babar(file='../../formfactor/data/k+k-_babar.dat'):
    cs_kpkm = pd.read_csv(file)
    cs_kpkm['energy_err'] = 0.
    cs_kpkm['cs_err'] = cs_kpkm['cs_err_rel']*cs_kpkm['cs']
    return cs_kpkm[['energy','energy_err','cs','cs_err','formfac']].query('energy>1.05&energy<2.1')

def get_KSKL_up(file='./data/fit_frame_19.csv', radcorsfile=None):
    cs_kskl = pd.read_csv(file).drop_duplicates(subset=['label', 'lum'])
    if radcorsfile is not None:
        rads = pd.read_csv(radcorsfile, index_col='name').sort_index()
        rads['rad_err_rel'] = rads['err']/rads['eff'] if 'eff' in rads else 0
        cs_kskl = cs_kskl.merge(rads, on='mcname', how='inner')
#         cs_kskl['rad'] = np.interp(cs_kskl.label, rads.index, rads.rad)
#         cs_kskl['rad_err_rel'] = np.interp(cs_kskl.label, rads.index, rads.rel_err)
    else:
        cs_kskl['rad'] = 1
        cs_kskl['rad_err_rel'] = 0
    cs_kskl['cs'] = cs_kskl.Ns/(cs_kskl.lum*cs_kskl.TrigEff*cs_kskl.rad)
    cs_kskl['energy'] = cs_kskl.Emean*2e-3
    cs_kskl['cs_err'] = cs_kskl.cs*np.sqrt( (cs_kskl.Ns_err/cs_kskl.Ns)**2 + 
                                            (cs_kskl.TrigErr/cs_kskl.TrigEff)**2 +
                                            cs_kskl.rad_err_rel**2 )
    cs_kskl['dEmin'] *= 2e-3
    cs_kskl['dEmax'] *= 2e-3
    cs_kskl = cs_kskl.rename({'dEmin':'energy_err_min', 'dEmax':'energy_err_max'}, axis=1)
    return cs_kskl[['energy','energy_err_min', 'energy_err_max','cs','cs_err','rad']]

def plot_cs(fit_func, params, css, labels):
    xmin, xmax = 10, 0
    for cs,lab in zip(css, labels):
        x_err = cs.energy_err if 'energy_err' in cs.columns else cs[['energy_err_min', 'energy_err_max']].values.T
        plt.errorbar(data=cs, x='energy', y='cs', xerr=x_err, yerr='cs_err', fmt='o', label=lab, ms=5, lw=1)
        xmin, xmax = min(xmin, cs.energy.min()), max(xmax, cs.energy.max())
    if fit_func is not None:
        x1 = np.linspace(xmin, xmax, 1000)
        plt.plot(x1, fit_func(x1, params), lw=2, c='tomato' )
    plt.legend()
    plt.yscale('log')