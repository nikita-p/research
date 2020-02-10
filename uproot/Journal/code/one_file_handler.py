
import uproot
import numpy as np
import pandas as pd

class HandleTree:
    mks = 497.611
    
    nt_cut = 2
    pmin_cut = 40
    z_cut = 13
    theta_cut = 0.6
    theta2_cut = np.pi - theta_cut
    hit_cut = 6
    chi2r_cut = 30
    chi2z_cut = 25
    rho_cut = 0.1
    align_cut = 0.8
    dedx_cut = 2000
    
    def pidedx_cut(self, P):
        return 5.58030e+9 / np.power(P + 40., 3) + 2.21228e+3 - 3.77103e-1 * P;

    def p_cut(self, E):
        return 2 * (0.0869 * E - 36.53)

    def p_ideal(self, E):
        return np.sqrt(E**2 - self.mks**2)

    def __init__(self, path, MC_soft=False):
        self.path = path
        self.ks = None
        self.tracks = None
        self.result = None
        self.MC_soft = MC_soft

    def open_file(self):
        df = uproot.open(self.path, xrootdsource={'parallel': True})['tr_ph']
        print('Opening file...', self.path)
        self.tracks = df.pandas.df(branches=['tptot', 'nt', 'tdedx', 'tz', 'tth', 'tphi', 'tnhit', \
                                              'tchi2r', 'tchi2z', 'trho', 'emeas'])
        print('Tracks block is ready.')
        self.ks = df.pandas.df(branches=['nks', 'ksptot', 'ksminv', 'emeas', 'ksalign', 'ksvind', 'trigbits', 'ebeam'])
        print('KS block is ready.')
        self.sim = df.pandas.df(branches=['emeas', 'simtype', 'nsim', 'simorig', 'simmom']) if self.MC_soft else None
        print('Generator block is ready.')

        print('File has been opened.')
        return
        
    def two_good_tracks_events(self):
        print('First cuts...')
        
        self.tracks = self.tracks.query('\
                    nt>=@self.nt_cut&\
                    tptot>@self.pmin_cut&\
                    tptot<1.1*emeas&\
                    abs(tz)<@self.z_cut&\
                    tchi2r<@self.chi2r_cut&\
                    tchi2z<@self.chi2z_cut&\
                    tth>@self.theta_cut&tth<@self.theta2_cut&\
                    tnhit>@self.hit_cut&\
                    abs(trho)>@self.rho_cut').drop(['nt', 'tz', 'tchi2r', 'tchi2z', 'tnhit', 'trho'], axis=1)
        self.tracks = self.tracks[ np.abs( self.pidedx_cut(self.tracks.tdedx) ) < self.dedx_cut ]
        print('Cuts have been applied')
        
        df_two_good_tracks = self.tracks.groupby('entry').agg({'tptot':'count'}).query('tptot==2')
        self.tracks = self.tracks.loc[ df_two_good_tracks.index ]
        print('Two good track events have been selected')

        return
    
    def good_ks(self):
        if 'ksvind[2]' in self.ks:
            ksvind_drops = [f'ksvind[{i}]' for i in range(2,20)]
            self.ks.drop(ksvind_drops, axis=1, inplace=True)
        
        print('Second cuts...')

        self.ks = self.ks.assign(difmass = np.abs(self.ks['ksminv'] - self.mks) )
        min_difmasses = self.ks.groupby('entry').agg({'difmass':np.min})
        self.ks = pd.merge(min_difmasses, self.ks, on=['difmass', 'entry']).drop(['difmass'], axis=1)
        
        print('Only the best kaon in every event has been selected')

        self.ks = self.ks.query('ksalign>@self.align_cut')
        print('Align cut has been done')
        
        self.ks = self.ks[ np.abs( self.ks['ksptot'] - self.p_ideal(self.ks['emeas']) ) < self.p_cut(self.ks['emeas']) ]
        print('Momentum cut has been done')
        
        self.ks = self.ks.rename({'ksvind[0]':'ksvind_0', 'ksvind[1]':'ksvind_1'}, axis=1)

        self.ks.drop(['nks'], axis=1, inplace=True)
        return 
        
    def merge(self):
        print('Two tables merging...')
        
        print('Tracks block len', self.tracks.shape[0])       
        print('KS block len', self.ks.shape[0])
        
        good_tracks_indexes = pd.DataFrame( self.tracks.to_records() ).set_index('entry').\
                                        groupby('entry').agg(ksvind_0 = ('subentry', 'min'),\
                                                             ksvind_1 = ('subentry','max') )
        self.result = pd.merge(good_tracks_indexes, self.ks, on=['entry', 'ksvind_0', 'ksvind_1'])
        self.result.drop(['ksvind_0', 'ksvind_1', 'ksptot', 'ksalign'], axis=1, inplace=True)
        
        return
    
    def soft_photon_events(self):
        self.sim = self.sim.query('simtype==310').copy()
        
        if self.sim.simmom.mean() < 1:
            self.sim.simmom = self.sim.simmom*1e3
        self.sim = self.sim[ np.abs( self.sim['simmom'] - self.p_ideal(self.sim['emeas']) ) < \
                            1.1*self.p_cut(self.sim['emeas']) ]
        
        only_one_ks_per_event = self.sim.groupby('entry').agg({'simmom':'count'}).query('simmom==1')
        self.sim = self.sim.loc[only_one_ks_per_event.index]
        
        self.sim = self.sim.reset_index().set_index('entry').drop(['subentry'], axis=1)
        return
        
    def analyze(self):
        self.open_file()
        self.two_good_tracks_events()
        self.good_ks()
        self.merge()
        
        if self.MC_soft:
            self.soft_photon_events()
            self.result = self.result.loc[self.sim.index & self.result.index]
            self.result = self.result.assign(total_soft_num = self.sim.shape[0])
        return
        
    def save_csv(self, folder):
        filename = folder + f"{self.result.emeas.mean():.2f}.csv"
        self.result.to_csv(filename, index=False)
        print('File has been saved in', filename)
        return
        
def pool_func(filename):
    a = HandleTree(filename, False)
    a.analyze()
    a.save_csv('outputs/2017/')
    return
