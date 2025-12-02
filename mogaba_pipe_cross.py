
from __main__ import *
from sklearn.cluster import DBSCAN as dbs
import pyclass as p

c_cs = p.comm
g_cs = p.gdict

c_cs("set variable general")
c_cs("set variable calibration")
c_cs("set variable position")

class CrossScan:
    def __init__(self,
                 path_p=None , path_c=None , path_dir=None  , saveplot=False,
                 file=None   , azel=None   , cs_source=None , station=None  ,
                 telname=None, snr=3       , pipe_log=None
                 ):
        self.path_p    = path_p         # path used in python
        self.path_c    = path_c         # path used in class
        self.path_dir  = path_dir       # path to load Tb excel file
        self.file      = file           # sdd file name
        self.station   = station        # station name (e.g., KYS)
        self.saveplot  = saveplot       # path to save plots (Tsys, Tau, AzEl, ..)
        self.azel      = azel           # Az or El info.
        self.telname   = telname        # telescope name in format of GILDAS/CLASS (e.g., KTN21M22L)
        self.cs_source = None           # source for cross-scan fitting
        self.snr       = snr            # signal-to-noise ratio used in Gaussian model-fitting of cross-scan profile
        self.pipe_log  = pipe_log
        self.nwalker   = None
        self.nstep     = None

    def set_init(self):
        if self.station is None : self.station=self.file.split('.sdd')[0].split('_')[-1]

        if   self.station == 'KPC':
            self.ant, self.lat, self.lon, self.height = name_PC, lat_PC, lon_PC,  height_PC
        elif self.station == 'KYS':
            self.ant, self.lat, self.lon, self.height = name_YS, lat_YS, lon_YS,  height_YS
        elif self.station == 'KUS':
            self.ant, self.lat, self.lon, self.height = name_US, lat_US, lon_US,  height_US
        elif self.station == 'KTN':
            self.ant, self.lat, self.lon, self.height = name_TN, lat_TN, lon_TN,  height_TN

        path_fig = self.path_dir + 'Figures/'
        path_dat = self.path_dir + 'data_cs/'
        mkdir('%s/'   %(path_fig))
        mkdir('%s/cs/'%(path_fig))

        mkdir('%s/'    %(path_dat))
        mkdir('%s/22/' %(path_dat))
        mkdir('%s/43/' %(path_dat))
        mkdir('%s/86/' %(path_dat))
        mkdir('%s/129/'%(path_dat))

        try:
            c_cs("sic directory %s"%(self.path_c))
        except:
            print('    !!! Given path "%s" does not exist.'%(self.path_c))
            abort()
        try:
            c_cs("file in %s"%(self.file))
        except:
            print('    !!! Given sdd file "%s" does not exist.'%(self.file))
            abort()
        c_cs("set default")
        c_cs("set ty c")
        c_cs("find")

        c_cs("get first")
        c_cs("define character*12 sourdate")
        c_cs("let sourdate 'r%head%gen%cdobs'")
        date = format_date(str(g_cs.sourdate.__sicdata__))
        delattr(g_cs, "sourdate")

        freqs = np.array([]).astype(int)
        bands = self.file.split('_')[-3].upper()
        if 'K' in bands : freqs=np.append(freqs, 22 )
        if 'Q' in bands : freqs=np.append(freqs, 43 )
        if 'W' in bands : freqs=np.append(freqs, 86 )
        if 'D' in bands : freqs=np.append(freqs, 129)
        self.bands = bands
        self.freqs = freqs
        self.date  = date

    def set_cs_default(self):
        if self.cs_source is None:
            print('    !!! Source for cross-scan fitting is not given.')
            abort()
        c_cs("set default")
        c_cs("set ty c")

    def get_planet(self):
        if   self.planet.upper() == 'VENUS'  : r=12104 *u.km
        elif self.planet.upper() == 'MARS'   : r=6792  *u.km
        elif self.planet.upper() == 'JUPITER': r=142984*u.km
        elif self.planet.upper() == 'SATURN' : r=120536*u.km
        time   = Ati(self.cs_log_sour[self.cs_log_sour.Source==self.planet]['Date'].values[0])
        ant    = EarthLocation(lat=self.lat, lon=self.lon, height=self.height)
        planet = get_body_barycentric(self.planet, time, ephemeris='builtin')
        earth  = get_body_barycentric('earth', time, ephemeris='builtin')
        x_p, y_p, z_p = planet.x.to(u.km), planet.y.to(u.km), planet.z.to(u.km)
        x_e, y_e, z_e = earth .x.to(u.km), earth .y.to(u.km), earth .z.to(u.km)
        d      = sqrt((x_p-x_e)**2 + (y_p-y_e)**2 + (z_p-z_e)**2)
        s      = ((r/d) * u.rad).to(u.arcsec)
        self.planet_dist = d    # units of km
        self.planet_size = s    # units of arc-second

    def load_source_info(self, time_thresh=60):
        Nscans  = int(g_cs.found)
        sources = np.array([], dtype="U32")
        scans   = np.array([], dtype="i4")
        dtime   = np.array([], dtype="U32")
        tsys    = np.array([], dtype="f8")
        tau     = np.array([], dtype="f8")
        az      = np.array([], dtype="f8")
        el      = np.array([], dtype="f8")
        freq    = np.array([], dtype="i4")
        line    = np.array([], dtype="U32")
        lr_pol  = np.array([], dtype="U32")
        for n in range(Nscans):
            if n==0: c_cs("get first")
            else   : c_cs("get next")
            c_cs("define character*12 sour")
            c_cs("define character*12 sourdate")
            c_cs("let sour     'r%head%pos%sourc'")
            c_cs("let sourdate 'r%head%gen%cdobs'")
            date     = format_date(str(g_cs.sourdate.__sicdata__))
            datetime = format_time(date, g_cs.ut.astype(float))
            # source = str(g_cs.sour.__sicdata__ ).replace('b','').replace("'",'').replace(' ','')
            # teles = str(g_cs.teles.astype(str)).replace('b','').replace("'",'').replace(' ','')
            source = np.char.decode(g_cs.sour.__sicdata__, "utf-8")
            teles = np.char.decode(g_cs.teles, "utf-8")
            source = str(np.char.rstrip(source))
            teles = str(np.char.rstrip(teles))
            rl = teles[-1]
            f  = str(teles.split('%s21M'%(self.station))[1].split('%s'%(rl))[0])

            sources = np.append(sources, source)
            scans   = np.append(scans  , g_cs.scan.astype(int))
            dtime   = np.append(dtime  , str(datetime))
            freq    = np.append(freq   , f)
            tsys    = np.append(tsys   , g_cs.tsys.astype(float))
            tau     = np.append(tau    , g_cs.tau .astype(float))
            az      = np.append(az     , g_cs.az  .astype(float)*u.rad.to(u.deg))
            el      = np.append(el     , g_cs.el  .astype(float)*u.rad.to(u.deg))
            line    = np.append(line   , str(g_cs.line.astype(str)).replace(' ',''))
            lr_pol  = np.append(lr_pol , rl)
            delattr(g_cs, "sour")
            delattr(g_cs, "sourdate")
        self.npol = len(np.unique(freq))
        cs_log_all = pd.DataFrame([sources, scans, dtime, tsys, tau, az, el, freq, line, lr_pol],
                                  index=['Source', 'ScanNum', 'Date', 'Tsys', 'Tau', 'Az', 'El', 'Freq', 'Line', 'LR_pol']).transpose()
        cs_log_all['Year'] = Ati(np.array(cs_log_all['Date']).astype(str), format='iso').byear
        cs_log_all['MJD']  = Ati(np.array(cs_log_all['Date']).astype(str), format='iso').mjd
        cs_log_all = cs_log_all[['Date', 'Year', 'MJD', 'Source', 'ScanNum', 'Tsys', 'Tau', 'Az', 'El', 'Freq', 'Line', 'LR_pol']]
        cs_log_all = cs_log_all.sort_values(by='ScanNum').reset_index(drop=True)
        self.cs_log_all  = cs_log_all

        log_mjd  =np.array([]).astype(float)
        log_sour =np.array([]).astype(str  )
        log_scan1=np.array([]).astype(int  )
        log_scan2=np.array([]).astype(int  )
        log_tsys =np.array([]).astype(float) ; log_dtsys=np.array([]).astype(float)
        log_tau  =np.array([]).astype(float) ; log_dtau =np.array([]).astype(float)
        log_az   =np.array([]).astype(float)
        log_el   =np.array([]).astype(float)

        map_time = lambda x : x.split(' ')[1]
        date  = cs_log_all['Date'][0].split(' ')[0]
        nchan = int(self.npol*2)
        cs_log_all['Time'] = np.round((cs_log_all['MJD']-Ati(date, format='iso').mjd) * u.day.to(u.s),0).astype(int)

        path_dat_cs = self.path_dir + 'data_cs/'
        mkdir(path_dat_cs)
        if cs_log_all.shape[0]<nchan:
            cols = ['Source', 'Date', 'Year', 'MJD', 'ScanNum', 'Nseq', 'Nscan', 'Tsys_1', 'dTsys_1', 'Tsys_2', 'dTsys_2', 'Tau_1', 'dTau_1', 'Tau_2', 'dTau_2', 'Az', 'El', 'Scan1', 'Scan2']
            self.cs_log_sour = pd.DataFrame([np.nan for i in range(len(cols))], index=cols).transpose()
        else:
            cs_log_all['Time'][nchan:]  = np.array(cs_log_all['Time'][nchan:]) - np.array(cs_log_all['Time'][0:-nchan])
            cs_log_all['Time'][0:nchan] = np.array(cs_log_all['Time'][0:nchan]) - cs_log_all['Time'][0]

            sep = np.append(np.array([0]), np.array(np.where(cs_log_all['Time'] > time_thresh)[0][0::nchan]))
            cs_log_sour = cs_log_all.iloc[sep].reset_index(drop=True)
            cs_log_sour['Nscan'] = np.zeros(cs_log_sour.shape[0]).astype(int)
            cs_log_sour['Nscan'][0:-1] = sep[1:] - sep[0:-1]
            cs_log_sour['Nscan'][cs_log_sour.shape[0]-1]  = int(cs_log_all.shape[0]-sep[-1])
            cs_log_sour['Nseq'] = np.ones(cs_log_sour.shape[0]).astype(int)
            sources=np.array([])
            for i in range(cs_log_sour.shape[0]):
                sour_name = cs_log_sour['Source'][i]
                if not sour_name  in sources:
                    sources = np.append(sources, sour_name)
                else:
                    cs_log_sour['Nseq'][i] = len(sources[sources==sour_name])+1
                    sources = np.append(sources, sour_name)
            cs_log_sour['Scan1' ]  = cs_log_sour['ScanNum']
            cs_log_sour['Scan2' ]  = cs_log_sour['ScanNum']+cs_log_sour['Nscan']-1
            cs_log_sour['Tsys_1']  = np.zeros(cs_log_sour.shape[0])
            cs_log_sour['Tau_1' ]  = np.zeros(cs_log_sour.shape[0])
            cs_log_sour['Tsys_2']  = np.zeros(cs_log_sour.shape[0])
            cs_log_sour['Tau_2' ]  = np.zeros(cs_log_sour.shape[0])
            cs_log_sour['dTsys_1'] = np.zeros(cs_log_sour.shape[0])
            cs_log_sour['dTau_1' ] = np.zeros(cs_log_sour.shape[0])
            cs_log_sour['dTsys_2'] = np.zeros(cs_log_sour.shape[0])
            cs_log_sour['dTau_2' ] = np.zeros(cs_log_sour.shape[0])
            endrows = [0]
            endrow = 0
            for Nsour in range(cs_log_sour.shape[0]):
                endrow += cs_log_sour['Nscan'][Nsour]
                endrows.append(endrow)

            for Nsour in range(cs_log_sour.shape[0]):
                cs_log_all_  = cs_log_all.iloc[endrows[Nsour]:endrows[Nsour+1]]
                cs_log_all_["Freq"] = cs_log_all_["Freq"].astype(int)
                cs_log_all_1 = cs_log_all_[cs_log_all_.Freq==int(self.freqs[0])]
                cs_log_sour['Tsys_1' ][Nsour] = np.nanmean(cs_log_all_1['Tsys'])
                cs_log_sour['dTsys_1'][Nsour] = np.nanstd (cs_log_all_1['Tsys'])
                cs_log_sour['Tau_1'  ][Nsour] = np.nanmean(cs_log_all_1['Tau' ])
                cs_log_sour['dTau_1' ][Nsour] = np.nanstd (cs_log_all_1['Tau' ])
                try:
                    cs_log_all_2 = cs_log_all_[cs_log_all_.Freq==int(self.freqs[1])]
                    cs_log_sour['Tsys_2' ][Nsour] = np.nanmean(cs_log_all_2['Tsys'])
                    cs_log_sour['dTsys_2'][Nsour] = np.nanstd (cs_log_all_2['Tsys'])
                    cs_log_sour['Tau_2'  ][Nsour] = np.nanmean(cs_log_all_2['Tau' ])
                    cs_log_sour['dTau_2' ][Nsour] = np.nanstd (cs_log_all_2['Tau' ])
                except:
                    pass

            cs_log_sour  = cs_log_sour[['Source', 'Date', 'Year', 'MJD', 'ScanNum', 'Nseq', 'Nscan', 'Tsys_1', 'dTsys_1', 'Tsys_2', 'dTsys_2', 'Tau_1', 'dTau_1', 'Tau_2', 'dTau_2', 'Az', 'El', 'Scan1', 'Scan2']]
            drop_lowscan = np.where(cs_log_sour['Nscan']<int(16*self.npol/4))[0]
            cs_log_sour  = cs_log_sour.drop(drop_lowscan, axis=0).reset_index(drop=True)

            self.cs_log_sour = cs_log_sour

        self.cs_log_all .to_excel(path_dat_cs + '%s_Log_CS_%s_all.xlsx'%(self.station, self.project))
        self.cs_log_sour.to_excel(path_dat_cs + '%s_Log_CS_%s_avg.xlsx'%(self.station, self.project))


    def fit_base(self):
        lr_pol = self.telname[-1]

        log_scans = self.cs_log_sour.copy()
        log_all   = self.cs_log_all.copy()
        log_scans = log_scans[log_scans.Source==self.cs_source]
        log_scans = log_scans[log_scans.Nseq  ==self.nseq     ].reset_index(drop=True)
        scan1, scan2  = log_scans['Scan1'][0], log_scans['Scan2'][0]

        log_   = self.cs_log_all.copy()
        log_   = log_[np.array(log_['Source']).astype(str)==str(self.cs_source)]
        log_   = log_[np.array(log_['Line']  ).astype(str)==str(self.azel)     ]
        log_   = log_[np.array(log_['Freq']  ).astype(str)==str(self.freq)     ]
        log_   = log_[np.array(log_['LR_pol']).astype(str)==str(lr_pol)        ].reset_index(drop=True)
        log_   = log_[np.logical_and(int(scan1)<=np.array(log_['ScanNum']).astype(int),
                                     np.array(log_['ScanNum']).astype(int)<=int(scan2))].reset_index(drop=True)
        if log_.shape[0]>=2:
            self.set_cs_default()
            c_cs("set tel %s" %(self.telname))
            c_cs("set sou %s" %(self.cs_source))
            c_cs("set line %s"%(self.azel))
            c_cs("set scan %s %s"%(log_['ScanNum'][0], log_['ScanNum'][log_.shape[0]-1]))
            c_cs("find")
            c_cs("set angle s")
            c_cs("aver /nocheck")

            self.x = g_cs.rx * u.rad.to(u.arcsec)
            self.y = g_cs.ry
            self.fit_cs()

            if self.savecsfit:
                self.SaveCSFit()
            self.peak=self.fit_p[0] ; self.dpeak=self.fit_c[0]
            self.offs=self.fit_p[1] ; self.doffs=self.fit_c[1]
            self.fwhm=self.fit_p[2] ; self.dfwhm=self.fit_c[2]

        else:
            self.dy=1
            self.peak =0.1    ; self.dpeak =np.nan
            self.offs =np.nan ; self.doffs =np.nan
            self.fwhm =np.nan ; self.dfwhm =np.nan

    def fit_cs(self):
        x  = self.x
        y  = self.y
        fwhm  = self.fwhm_b
        sigma = fwhm/(2*sqrt(2*np.log(2)))

        data_fit = pd.DataFrame([x, y*1e3], index=['x', 'y']).transpose()
        data_fit = data_fit.dropna(axis=0).reset_index(drop=True)
        y_min, y_max = abs(np.min(y)), abs(np.max(y))
        try:
            p, c  = curve_fit(Gaussian, data_fit['x'], data_fit['y'],
                            bounds=[[-0.8*y_max*1e3, -1.0*fwhm, 0.6*fwhm, -1*y_min*1e3],
                                    [+1.2*y_max*1e3, +1.0*fwhm, 1.4*fwhm, +1*y_min*1e3]]
                            ,maxfev=10000)
        except:
           p = [np.max(y)*1e3, 0, fwhm, 0]
        if np.logical_and(self.freq==129, self.cs_source in ['Venus', 'Jupiter']):
            fwhm = fwhm*2
        cmc = CrossMCMC(x=x, y=y*1e3, yerr=np.zeros(len(y)), p=p, fwhm=fwhm)
        cmc.nwalker = self.nwalker
        cmc.nstep   = self.nstep
        cmc.freq    = self.freq
        cmc.run_cmcmc()
        p_mcmc = cmc.p_mcmc
        u_mcmc = cmc.u_mcmc
        f_smpl = cmc.f_smpl

        y_model = Gaussian(x, *p_mcmc)
        y_residual = y-y_model
        dy=np.std(y_residual)
        self.fit_p, self.fit_c, self.f_smpl  = p_mcmc, u_mcmc, f_smpl
        self.x, self.y, self.dy = x, y, dy

    def run_cs(self):
        cs_log_all  = self.cs_log_all
        print(self.cs_log_sour.columns)
        cs_log_sour = self.cs_log_sour[np.logical_and(self.cs_log_sour.Tsys_1>0, self.cs_log_sour.Tsys_2>0)].reset_index(drop=True)
        cs_log_sour = check_neg_tsys(cs_log_all, cs_log_sour)
        cs_log_all['Freq'] = np.array(cs_log_all['Freq']).astype(int).astype(str)
        self.npol = len(np.unique(cs_log_all['Freq']))

        sources = cs_log_sour['Source']
        lr_pols = np.unique(cs_log_all['LR_pol'])
        azels   = np.unique(cs_log_all['Line'  ])
        freqs   = np.sort(np.unique(cs_log_all['Freq']).astype(int))
        self.freqs = freqs

        self.path_fig_cs = self.path_dir + 'Figures/cs/%s_%s/'%(self.station, self.date)
        mkdir(self.path_fig_cs)
        mkdir(self.path_fig_cs+'../CS_Logs/')
        if self.saveplot:
            self.SaveCSLog()

        for Nfreq, freq in enumerate(freqs):
            fwhm  = ((C.c/(freq*u.GHz)).to(u.m)/(21*u.m)) * u.rad.to(u.arcsecond)
            self.fwhm_b = fwhm
            if Nfreq==0:
                self.date = self.cs_log_all['Date'][0].split(' ')[0]
            self.freq = freq
            self.path_dat_cs = self.path_dir + 'data_cs/%s/'%(self.freq)
            sours = []
            mjds  = []
            nseqs = []
            x_l, y_l, dy_l   = [], [], []
            x_r, y_r, dy_r   = [], [], []
            peak_l , peak_r  = [], []
            dpeak_l, dpeak_r = [], []
            fwhm_l, fwhm_r = [], []
            offs_l, offs_r = [], []
            els = []
            for Nsour, source in enumerate(sources):
                scan1 = cs_log_sour.iloc[Nsour]['Scan1']
                scan2 = cs_log_sour.iloc[Nsour]['Scan2']
                nseq  = cs_log_sour.iloc[Nsour]['Nseq' ]
                mjd   = cs_log_sour[np.logical_and(cs_log_sour.Source==source, cs_log_sour.Nseq==nseq)]['MJD'].values[0]
                el    = cs_log_sour[np.logical_and(cs_log_sour.Source==source, cs_log_sour.Nseq==nseq)].reset_index(drop=True)['El'][0]
                els.append(el)
                for lr_pol in lr_pols:  # [L R]
                    self.lr_pol = lr_pol
                    for azel in azels:  # [Az, El]
                        self.saveplot  = True
                        self.cs_source = source
                        self.azel      = azel
                        self.telname   = '%s21M%s%s'%(self.station, freq, lr_pol)
                        self.scan1     = scan1
                        self.scan2     = scan2
                        self.nseq      = nseq
                        self.fit_base()
                        peak, dpeak = self.peak, self.dpeak
                        offs, doffs = self.offs, self.doffs
                        fwhm, dfwhm = self.fwhm, self.dfwhm
                        fwhm_b      = self.fwhm_b

                        dy = self.dy
                        #if np.logical_or(peak/dy<self.snr, np.isnan(peak)   ) or\
                        if np.logical_or(peak/dpeak<self.snr, np.isnan(peak)) or\
                           np.logical_or(offs>=0.7*fwhm_b, doffs>=0.2*fwhm_b) or\
                           np.logical_or(0.5*fwhm_b>fwhm , fwhm>2.0*fwhm_b  ) or\
                           fwhm<=3*dfwhm                                      :
                            if azel=='AZ':
                                peak_az=ufloat(np.nan, np.nan)
                                fwhm_az=ufloat(np.nan, np.nan)
                                offs_az=ufloat(np.nan, np.nan)
                            elif azel=='EL':
                                peak_el=ufloat(np.nan, np.nan)
                                fwhm_el=ufloat(np.nan, np.nan)
                                offs_el=ufloat(np.nan, np.nan)
                        else:
                            if azel=='AZ':
                                peak_az=ufloat(peak, dpeak)
                                fwhm_az=ufloat(fwhm, dfwhm)
                                offs_az=ufloat(offs, doffs)
                            elif azel=='EL':
                                peak_el=ufloat(peak, dpeak)
                                fwhm_el=ufloat(fwhm, dfwhm)
                                offs_el=ufloat(offs, doffs)
                    peak_az = peak_az * e**( 4*ln(2) * (offs_el/fwhm_az)**2 )
                    peak_el = peak_el * e**( 4*ln(2) * (offs_az/fwhm_el)**2 )
                    peak_T  = (peak_az+peak_el)/2

                    fwhm_ = (fwhm_az + fwhm_el)/2
                    offs_ = (offs_az**2 + offs_el**2)**0.5

                    gain = cal_gain_curve(self.station, int(self.date[:4]), self.freq, el)
                    peak_T /= gain
                    if unp.nominal_values(peak_T)/unp.std_devs(peak_T)<=self.snr:
                        if lr_pol=='L':
                            print(source)
                            sours  .append(source)
                            mjds   .append(mjd   )
                            nseqs  .append(nseq  )
                            peak_l .append(np.nan)
                            dpeak_l.append(np.nan)
                            fwhm_l .append(np.nan)
                            offs_l .append(np.nan)
                            dy_l   .append(np.nan)
                        elif lr_pol=='R':
                            peak_r .append(np.nan)
                            dpeak_r.append(np.nan)
                            fwhm_r .append(np.nan)
                            offs_r .append(np.nan)
                            dy_r   .append(np.nan)
                    else:
                        if lr_pol=='L':
                            sours  .append(source)
                            mjds   .append(mjd   )
                            nseqs  .append(nseq  )
                            peak_l .append(np.round(unp.nominal_values(peak_T), 5))
                            dpeak_l.append(np.round(unp.std_devs(peak_T)      , 5))
                            fwhm_l .append(np.round(unp.nominal_values(fwhm_) , 5))
                            offs_l .append(np.round(unp.nominal_values(offs_) , 5))
                            dy_l   .append(np.round(self.dy                   , 5))
                        elif lr_pol=='R':
                            peak_r .append(np.round(unp.nominal_values(peak_T), 5))
                            dpeak_r.append(np.round(unp.std_devs(peak_T)      , 5))
                            fwhm_r .append(np.round(unp.nominal_values(fwhm_) , 5))
                            offs_r .append(np.round(unp.nominal_values(offs_) , 5))
                            dy_r   .append(np.round(self.dy                   , 5))

            cs_fit_info = pd.DataFrame([sours, mjds, nseqs, els, peak_l, dpeak_l, dy_l, peak_r, dpeak_r, dy_r, fwhm_l, fwhm_r, offs_l, offs_r],
                                       index=['Source', 'MJD', 'Nseq', 'El', 'Peak_L', 'dPeak_L', 'stdT_L', 'Peak_R', 'dPeak_R', 'stdT_R', "fwhm_L", "fwhm_R", "offset_L", "offset_R"]).transpose()
            cs_fit_info.dropna(axis=0, inplace=True)
            cs_fit_info.reset_index(drop=True, inplace=True)

            planets = np.array([])
            aeff_planets = []
            if str(self.freq) in ['22', '43', '86', '129'] : aeff_planets.append('VENUS')
            if str(self.freq) in ['43', '86', '129']       : aeff_planets.append('MARS')
            if str(self.freq) in ['22', '43', '86']        : aeff_planets.append('JUPITER')
            for Nsour in range(cs_fit_info.shape[0]):
                source = cs_fit_info['Source'][Nsour]
                if source.upper() in aeff_planets:
                    planets = np.append(planets, source)

            etas_l   = np.array([])
            etas_r   = np.array([])
            planets  = np.unique(planets)
            mjd_date = int(Ati(self.date, format='iso').mjd)
            template = pd.DataFrame([])
            if len(planets)!=0:
                for planet in planets:
                    self.planet=planet
                    self.get_planet()
                    cs_fit_sour = cs_fit_info[cs_fit_info.Source==self.planet].reset_index(drop=True)
                    cs_fit_sour["eta_l" ] = np.full(cs_fit_sour.shape[0], np.nan)
                    cs_fit_sour["eta_r" ] = np.full(cs_fit_sour.shape[0], np.nan)
                    cs_fit_sour["deta_l"] = np.full(cs_fit_sour.shape[0], np.nan)
                    cs_fit_sour["deta_r"] = np.full(cs_fit_sour.shape[0], np.nan)
                    if planet=='Mars':
                        Tb_mars  = pd.read_excel(self.path_dir+'Tb_Mars.xlsx', index_col=[0])
                        Tb=Tb_mars[Tb_mars.MJD==mjd_date]['Tb%s'%(self.freq)].values
                    if planet=='Jupiter' : Tb=get_Tb(self.path_dir, 'Tb_Jupiter.xlsx', self.planet, self.freq)
                    if planet=='Venus'   : Tb=get_Tb(self.path_dir, 'Tb_Venus.xlsx'  , self.planet, self.freq)
                    for i in range(cs_fit_sour.shape[0]):
                        Ta_l   = ufloat(cs_fit_sour['Peak_L' ][i], cs_fit_sour['dPeak_L'][i])
                        Ta_r   = ufloat(cs_fit_sour['Peak_R' ][i], cs_fit_sour['dPeak_R'][i])
                        eta_l  = cal_aeff(nu=self.freq*u.GHz, s_p=self.planet_size, Tb=Tb, Ta=Ta_l).value
                        eta_r  = cal_aeff(nu=self.freq*u.GHz, s_p=self.planet_size, Tb=Tb, Ta=Ta_r).value
                        etas_l = np.append(etas_l, eta_l)
                        etas_r = np.append(etas_r, eta_r)
                        cs_fit_sour["eta_l" ][i] = unp.nominal_values(eta_l)
                        cs_fit_sour["deta_l"][i] = unp.std_devs(eta_l)
                        cs_fit_sour["eta_r" ][i] = unp.nominal_values(eta_r)
                        cs_fit_sour["deta_r"][i] = unp.std_devs(eta_r)
                    template = pd.concat([template, cs_fit_sour], axis=0).reset_index(drop=True)
                template["Tsys" ] = np.zeros(template.shape[0])
                template["Tau"  ] = np.zeros(template.shape[0])
                template["dTsys"] = np.zeros(template.shape[0])
                template["dTau" ] = np.zeros(template.shape[0])
                for i in range(template.shape[0]):
                    if freq==freqs[0]:
                        template["Tsys" ][i] = cs_log_sour[cs_log_sour.MJD==template["MJD"][i]]["Tsys_1" ]
                        template["Tau"  ][i] = cs_log_sour[cs_log_sour.MJD==template["MJD"][i]]["Tau_1"  ]
                        template["dTsys"][i] = cs_log_sour[cs_log_sour.MJD==template["MJD"][i]]["dTsys_1"]
                        template["dTau" ][i] = cs_log_sour[cs_log_sour.MJD==template["MJD"][i]]["dTau_1" ]
                    elif freq==freqs[1]:
                        template["Tsys" ][i] = cs_log_sour[cs_log_sour.MJD==template["MJD"][i]]["Tsys_2" ]
                        template["Tau"  ][i] = cs_log_sour[cs_log_sour.MJD==template["MJD"][i]]["Tau_2"  ]
                        template["dTsys"][i] = cs_log_sour[cs_log_sour.MJD==template["MJD"][i]]["dTsys_2"]
                        template["dTau" ][i] = cs_log_sour[cs_log_sour.MJD==template["MJD"][i]]["dTau_2" ]
                rndidx = ['MJD', 'El', 'Peak_L', 'dPeak_L', 'stdT_L', 'Peak_R', 'dPeak_R', 'stdT_R', "fwhm_L", "fwhm_R", "offset_L", "offset_R", "Tsys", "dTsys", "Tau", "dTau", "eta_l", "deta_l", "eta_r", "deta_r"]
                for i in range(len(rndidx)):
                    template[rndidx[i]] = np.round(np.array(template[rndidx[i]], dtype="f8"), 3)
                template = template[["Source", "MJD", "Nseq", "Peak_L", "dPeak_L", "stdT_L", "Peak_R", "dPeak_R", "stdT_R", "fwhm_L", "fwhm_R", "offset_L", "offset_R", "Tsys", "dTsys", "Tau", "dTau", "El", "eta_l", "deta_l", "eta_r", "deta_r"]]
                template = template.sort_values(by="MJD").reset_index(drop=True)
                print("\n # {0} GHz #\n".format(freq), template, "\n")
                dropidx = input("Enter index numbers to drop (if not, press 'enter'): ")
                if dropidx:
                    dropidx = dropidx.replace(",", " ")
                    dropidx = dropidx.split(" ")
                    dropidx = list(filter(None, dropidx))
                    dropidx = list(map(int, dropidx))
                    if dropidx:
                        dropidx.sort()
                        template = template.drop(dropidx, axis=0).reset_index(drop=True)

                print(template)

                if template.shape[0] != 0:
                    eta_l, deta_l = np.average(template["eta_l"], weights=1/template["deta_l"]**2, returned=True)
                    eta_r, deta_r = np.average(template["eta_r"], weights=1/template["deta_r"]**2, returned=True)
                    cs_fit_info['eta_L'] = np.array([np.round(eta_l, 5) for i in range(cs_fit_info.shape[0])])
                    cs_fit_info['eta_R'] = np.array([np.round(eta_r, 5) for i in range(cs_fit_info.shape[0])])
                else:
                    eta_l, deta_l = np.nan, np.nan
                    eta_r, deta_r = np.nan, np.nan
                    cs_fit_info['eta_L'] = np.array([np.nan for i in range(cs_fit_info.shape[0])])
                    cs_fit_info['eta_R'] = np.array([np.nan for i in range(cs_fit_info.shape[0])])

            else:
                cs_fit_info["eta_L"] = np.full(cs_fit_info.shape[0], np.nan)
                cs_fit_info["eta_R"] = np.full(cs_fit_info.shape[0], np.nan)
            peak_l = unp.uarray(cs_fit_info['Peak_L'], cs_fit_info['dPeak_L'])
            peak_r = unp.uarray(cs_fit_info['Peak_R'], cs_fit_info['dPeak_R'])
            cs_fit_info[ 'S_L'] = np.round(unp.nominal_values(8*C.k_B.value/(pi*21**2*cs_fit_info['eta_L'])*peak_l*(u.J/(u.m)**2).to(u.Jy)) , 5)
            cs_fit_info['dS_L'] = np.round(unp.std_devs      (8*C.k_B.value/(pi*21**2*cs_fit_info['eta_L'])*peak_l*(u.J/(u.m)**2).to(u.Jy)) , 5)
            cs_fit_info[ 'S_R'] = np.round(unp.nominal_values(8*C.k_B.value/(pi*21**2*cs_fit_info['eta_R'])*peak_r*(u.J/(u.m)**2).to(u.Jy)) , 5)
            cs_fit_info['dS_R'] = np.round(unp.std_devs      (8*C.k_B.value/(pi*21**2*cs_fit_info['eta_R'])*peak_r*(u.J/(u.m)**2).to(u.Jy)) , 5)

            if Nfreq==0:
                self.cs_log_all_1  = self.cs_log_all
                self.cs_log_sour_1 = self.cs_log_sour
                self.cs_fit_info_1 = cs_fit_info
                if np.logical_or(etas_l.shape[0]!=0, etas_r.shape[0]!=0):
                    self.eta_l_1 = eta_l
                    self.eta_r_1 = eta_r
                else:
                    self.eta_l_1 = np.nan
                    self.eta_r_1 = np.nan
            if Nfreq==1:
                self.cs_log_all_2  = self.cs_log_all
                self.cs_log_sour_2 = self.cs_log_sour
                self.cs_fit_info_2 = cs_fit_info
                if np.logical_or(etas_l.shape[0]!=0, etas_r.shape[0]!=0):
                    self.eta_l_2       = eta_l
                    self.eta_r_2       = eta_r
                else:
                    self.eta_l_2 = np.nan
                    self.eta_r_2 = np.nan

            cs_fit_info.to_excel(self.path_dat_cs + './%s_Log_CS_fit_%s_%s.xlsx' %(self.station, self.date, self.freq))


    def SaveCSFit(self):
        x, y, p, c, f_smpl = self.x, np.array(self.y,dtype="f8"), self.fit_p, self.fit_c, self.f_smpl
        model = Gaussian(x, *p)
        path_dir = self.path_fig_cs
        if not os.path.isdir(path_dir): os.system('mkdir %s'%(path_dir))

        labels = [r'$T_{\rm a,peak}$', r'$\Delta \delta$', r'${\rm FWHM}$', r'$T_{\rm base}$']
        fig_corner = corner.corner(f_smpl, labels=labels)
        ax_cs = fig_corner.add_axes([0.65, 0.57,0.33,0.33])
        ax_cs.plot(x, y    , c='black', lw=2.0)
        ax_cs.plot(x, model, c='red'  , lw=2.0)
        fig_corner.suptitle(
            'Cross-Scan\n%s | %s%s | %s\nT:%.2f (%.2f) | Off:%.2f (%.2f) | FWHM:%.2f (%.2f)'%(
                self.cs_source, self.freq, self.lr_pol, self.azel, p[0], c[0], p[1], c[1], p[2], c[2]),
            fontsize=14, x=0.63)
        ax_cs.set_xlabel(r'$\delta \, {\rm(arcsec)}$', fontsize=15)
        ax_cs.set_ylabel(r'$T_{\rm a} \, {\rm (K)}$' , fontsize=15)
        fig_corner.savefig(self.path_fig_cs + '%s_%s_CSFit_%s_%s%s_%s_%s.png'%(self.date, self.station, self.cs_source, self.freq, self.lr_pol, self.azel, self.nseq))
        close_figure(fig_corner)

    def SaveCSLog(self):
        log_all = self.cs_log_all
        log_avg = self.cs_log_sour
        date  = self.date
        freqs = np.sort(np.unique(log_all['Freq']).astype(int))
        LRs   = ['L', 'R']
        self.freqs = freqs

        fsize = 16
        fig_cslog, ax_cslog = plt.subplots(3, 4, figsize=(fsize, fsize*9/16), sharex=True)
        tsys_1_l, tsys_1_r, tsys_2_l, tsys_2_r = ax_cslog[0,0], ax_cslog[0,1], ax_cslog[0,2], ax_cslog[0,3]
        tau_1_l , tau_1_r , tau_2_l , tau_2_r  = ax_cslog[1,0], ax_cslog[1,1], ax_cslog[1,2], ax_cslog[1,3]
        el_1_l  , el_1_r  , el_2_l  , el_2_r   = ax_cslog[2,0], ax_cslog[2,1], ax_cslog[2,2], ax_cslog[2,3]

        tsys_1_l.set_title('%sL'%(freqs[0])) ; tsys_1_r.set_title('%sR'%(freqs[0]))
        tsys_1_l.set_ylabel(r'$T_{\rm sys}$'+r'$\rm \,(K)$'    , fontsize=15)
        tau_1_l .set_ylabel(r'$\tau$'                          , fontsize=15)
        el_1_l  .set_ylabel(r'$\rm Elevation$'+r'$\rm \,(deg)$', fontsize=12)
        if self.npol==2:
            tsys_2_l.set_title('%sL'%(freqs[1])) ; tsys_2_r.set_title('%sR'%(freqs[1]))

        axes_xlabel = [tsys_1_l, tsys_1_r, tsys_2_l, tsys_2_r, tau_1_l, tau_1_r, tau_2_l, tau_2_r]
        axes_ylabel = [tsys_1_r, tsys_2_r, tau_1_r , tau_2_l, tau_2_r, el_1_r , el_2_l , el_2_r]
        axes_ylim   = [el_1_l  , el_1_r  , el_2_l  , el_2_r]
        axes_all    = [tsys_1_l, tsys_1_r, tsys_2_l, tsys_2_r,
                       tau_1_l , tau_1_r , tau_2_l , tau_2_r ,
                       el_1_l  , el_1_r  , el_2_l  , el_2_r  ]

        for n in range(len(axes_xlabel)) : axes_xlabel[n].tick_params(labelbottom=False)
        for n in range(len(axes_ylabel)) : axes_ylabel[n].tick_params(labelleft=False)
        for n in range(len(axes_ylim  )) : axes_ylim[n].set_ylim(0,90)
        for n in range(len(axes_all   )) :
            axes_all[n].xaxis.set_major_locator(MultipleLocator(3))
            axes_all[n].xaxis.set_minor_locator(MultipleLocator(1))

        project = self.file.split('_%s'%(self.station))[0]
        fig_cslog.text(0.5, 0.02, '(%s, %s) | UTC (hour)'%(project, date), ha='center', fontsize=15)

        colors = ['pink'     , 'coral', 'red' , 'lime'   , 'green' ,
                  'lightblue', 'aqua' , 'blue', 'magenta', 'purple']

        log_all["Tsys"][np.isinf(log_all.Tsys.values.astype("f8"))] = np.nan
        log_all["Tau" ][np.isinf(log_all.Tau .values.astype("f8"))] = np.nan
        log_all["Time"] = (log_all["MJD"] - Ati(date, format='iso').mjd) * u.day.to(u.hr)

        date_all = log_all["Date"].values.astype(str)
        mjd_all = Ati(date_all, format="iso").mjd
        time_sec = ((mjd_all - int(np.min(mjd_all))) * 24 * 3600)

        db = dbs(eps=60, min_samples=1).fit((time_sec).reshape(-1, 1))
        scannums = db.labels_
        uscan = np.unique(scannums)

        pols = ["L", "R"]
        for nscan, scan in enumerate(uscan):
            log_lr = log_all[scannums == scan].reset_index(drop=True)

            if log_lr.shape[0] == 0:
                continue

            for npol, pol in enumerate(pols):
                mask_pol = log_lr["LR_pol"] == pol
                log_pol = log_lr[mask_pol].reset_index(drop=True)

                if log_pol.shape[0] == 0:
                    continue

                ax_tsys1 = tsys_1_l if npol == 0 else tsys_1_r
                ax_tsys2 = tsys_2_l if npol == 0 else tsys_2_r
                ax_tau1  = tau_1_l  if npol == 0 else tau_1_r
                ax_tau2  = tau_2_l  if npol == 0 else tau_2_r
                ax_el1   = el_1_l   if npol == 0 else el_1_r
                ax_el2   = el_2_l   if npol == 0 else el_2_r

                for nfreq, freq in enumerate(freqs):
                    val_freq = log_pol["Freq"].values.astype(int)
                    mask_freq = val_freq == freq
                    log_freq = log_pol.loc[mask_freq].reset_index(drop=True)

                    if log_freq.shape[0] == 0:
                        continue

                    if nfreq == 0:
                        ax_tsys, ax_tau, ax_el = ax_tsys1, ax_tau1, ax_el1
                    elif nfreq == 1:
                        ax_tsys, ax_tau, ax_el = ax_tsys2, ax_tau2, ax_el2

                    if npol == 0 and nfreq == 0:
                        ax_tsys.scatter(log_freq["Time"], log_freq["Tsys"], marker="x", s=30, color=colors[nscan % len(colors)], label=log_freq["Source"][0])
                    else:
                        ax_tsys.scatter(log_freq["Time"], log_freq["Tsys"], marker="x", s=30, color=colors[nscan % len(colors)])
                    ax_tau .scatter(log_freq["Time"], log_freq["Tau" ], marker="x", s=30, color=colors[nscan % len(colors)])
                    ax_el  .scatter(log_freq["Time"], log_freq["El"  ], marker="x", s=30, color=colors[nscan % len(colors)])

        tsys_min = np.nanmin(log_all["Tsys"])
        tsys_max = np.nanmax(log_all["Tsys"])
        tau_min = np.nanmin(log_all["Tau"])
        tau_max = np.nanmax(log_all["Tau"])
        tsys_1_l.set_yscale("log")
        tsys_1_r.set_yscale("log")
        tsys_2_l.set_yscale("log")
        tsys_2_r.set_yscale("log")
        tsys_1_l.set_ylim(0.9*tsys_min, 1.1*tsys_max)
        tsys_1_r.set_ylim(0.9*tsys_min, 1.1*tsys_max)
        tsys_2_l.set_ylim(0.9*tsys_min, 1.1*tsys_max)
        tsys_2_r.set_ylim(0.9*tsys_min, 1.1*tsys_max)
        tau_1_l.set_ylim(0.0, 1.1*tau_max)
        tau_1_r.set_ylim(0.0, 1.1*tau_max)
        tau_2_l.set_ylim(0.0, 1.1*tau_max)
        tau_2_r.set_ylim(0.0, 1.1*tau_max)
        tsys_1_l.grid(True)
        tsys_1_r.grid(True)
        tsys_2_l.grid(True)
        tsys_2_r.grid(True)
        tau_1_l.grid(True)
        tau_1_r.grid(True)
        tau_2_l.grid(True)
        tau_2_r.grid(True)
        el_1_l.grid(True)
        el_1_r.grid(True)
        el_2_l.grid(True)
        el_2_r.grid(True)
        tsys_1_l.legend(ncol=10, fontsize=10, bbox_to_anchor=(4.50, 1.50), fancybox=True)
        fig_cslog.savefig(self.path_fig_cs + '../CS_Logs/%s_%s_CSLog.png'%(self.station, self.date))
        close_figure(fig_cslog)

class CrossMCMC:
    def __init__(self, x=None, y=None, yerr=None, p=None, fwhm=None):
        self.x       = x
        self.y       = y
        self.yerr    = yerr
        self.p       = p        # Initial Gaussian-fit parameters, including 'log_f' value
        self.fwhm    = fwhm     # Full-width half maximum of the (Gaussian) beam at a given frequency
        self.nwalker = 4 * 2
        self.nstep   = 2000
        self.nburn   = int(0.5*self.nstep)

    def run_cmcmc(self):
        x, y, yerr = self.x, self.y, self.yerr
        p, fwhm = self.p, self.fwhm
        freq = self.freq
        model= Gaussian(x, *p)
        resi = y - model
        rms  = np.std(resi)
        yerr = np.array([rms for i in range(len(y))])

        nll = lambda *args: -log_likelihood(*args)
        p0 = p
        nstep   = self.nstep
        nburn   = int(0.35*nstep)
        ndim    = len(p0)
        nwalker = self.nwalker
        ymin, ymax = abs(np.min(y)), abs(np.max(y))
        p_ml = minimize(nll, p0, args=(x, y, yerr),
                    bounds=[[+0.8*ymax, +1.2*ymax],
                            [-0.1*fwhm, +0.1*fwhm],
                            [+0.6*fwhm, +1.4*fwhm],
                            [-3.0*ymin, +3.0*ymin]])

        pos = p_ml.x + 1e-4 * np.random.randn(int(nwalker), ndim)
        sampler = emcee.EnsembleSampler(nwalker, len(p0), log_probability, args=(x, y, yerr, ymin, ymax, fwhm))
        sampler.run_mcmc(pos, nstep, progress=True)
        flat_samples = sampler.get_chain(discard=self.nburn, thin=1, flat=True)

        p_mcmc, u_mcmc = np.array([]), np.array([])
        for i in range(ndim):
            mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
            q    = np.diff(mcmc)
            v, err = mcmc[1], (q[0]+q[1])/2
            p_mcmc = np.append(p_mcmc, v)
            u_mcmc = np.append(u_mcmc, err)

        p_mcmc[0]/=1e3 ; u_mcmc[0]/=1e3
        p_mcmc[3]/=1e3 ; u_mcmc[3]/=1e3

        self.p_mcmc = p_mcmc
        self.u_mcmc = u_mcmc
        self.f_smpl = flat_samples


def log_likelihood(theta, x, y, yerr):
    peak, offs, fwhm, base = theta
    model  = Gaussian(x, peak, offs, fwhm, base)
    return -0.5 * np.sum((y - model)**2 / yerr**2 + np.log(2 * np.pi * yerr**2))

def log_prior(theta, ymin, ymax, fwhm):
    peak, offs, fwgm_f, base = theta
    if +0.6*ymax < peak   < +1.5*ymax and\
       -20       < offs   < +20       and\
       +0.4*fwhm < fwgm_f < +1.6*fwhm and\
       -3.0*ymin < base   < +3.0*ymin :
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr, ymin, ymax, fwhm):
    lp = log_prior(theta, ymin, ymax, fwhm)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

def Gaussian(x, peak, offs, fwhm, baseline):
    model = peak * exp(-4*ln(2)/fwhm**2 * (x-offs)**2) + baseline
    return model

def get_cslog(cs_log, source, freq, line, lr_pol):
    cs_log_ = cs_log.copy()
    cs_log_ = cs_log_[cs_log_.Source == source]
    cs_log_ = cs_log_[cs_log_.Freq   == freq  ]
    cs_log_ = cs_log_[cs_log_.LR_pol == lr_pol]
    cs_log_ = cs_log_.reset_index().drop('index', axis=1)
    return cs_log_

def get_Tb(path, file, planet, freq):
    openTb = pd.read_excel(path+file)
    freq   = float(freq)
    if freq==22 and planet=="Jupiter":
        tb = ufloat(openTb['Tb'][0], openTb['dTb'][0])
    else:
        N=0
        while openTb['Freq'][N]<=freq:
            N+=1
        N -= 1
        f1 = openTb['Freq'][N+0]
        f2 = openTb['Freq'][N+1]
        t1 = ufloat(openTb['Tb'][N+0], openTb['dTb'][N+0])
        t2 = ufloat(openTb['Tb'][N+1], openTb['dTb'][N+1])
        tb = t1 + (freq-f1) * (t2-t1)/(f2-f1)
    return tb

def cal_aeff(Ta, Tb, s_p, nu):
    """
    w       : observing wavelength              in [m]
    nu      : observing frequency               in [GHz]
    D       : physical antenna diameter         in [m]
    Ta      : observed antenna temperature      in [K]
    Tb      : observed antenna temperature      in [K]
    s_p     : angular size of a source (planet) in [arc-second]
    s_m     : FWHM of the main beam             in [arc-second]
    # See (Lee et al. 2011, PASP, 123, 1398L) for more details
    """
    D   = 21*u.m
    w   = (C.c/nu).to(u.m)
    s_m = 0.89 * ((w/D) * u.rad).to(u.arcsec)
    s_p = 1.00 * s_p

    o_m = (1.133 * s_m**2).to(u.sr)
    o_s = (o_m * (1 - exp(-ln(2)*(s_p/s_m)**2))).to(u.sr)

    aeff = (w**2 * Ta / (np.pi*(D/2)**2) / Tb / o_s * u.sr)
    return aeff

def check_neg_tsys(log_all, log_sour):
    drop_index = []
    for Nsour, source in enumerate(log_sour['Source']):
        scan1, scan2 = log_sour['Scan1'][Nsour], log_sour['Scan2'][Nsour]
        tsys_sour = log_all[np.logical_and(scan1<=log_all.ScanNum, log_all.ScanNum<=scan2)].reset_index(drop=True)
        if np.min(tsys_sour['Tsys'])<0:
            drop_index.append(Nsour)
    log_sour = log_sour.drop(drop_index, axis=0).reset_index(drop=True)
    return log_sour

name_PC, lat_PC, lon_PC,  height_PC = 'PC', (37*u.deg+32*u.arcmin+00.1*u.arcsec).to(u.deg), (128*u.deg+26*u.arcmin+55.1*u.arcsec).to(u.deg), 541*u.m
name_YS, lat_YS, lon_YS,  height_YS = 'YS', (37*u.deg+33*u.arcmin+54.9*u.arcsec).to(u.deg), (126*u.deg+56*u.arcmin+27.4*u.arcsec).to(u.deg), 139*u.m
name_US, lat_US, lon_US,  height_US = 'US', (35*u.deg+32*u.arcmin+44.2*u.arcsec).to(u.deg), (129*u.deg+14*u.arcmin+59.3*u.arcsec).to(u.deg), 170*u.m
name_TN, lat_TN, lon_TN,  height_TN = 'TN', (33*u.deg+17*u.arcmin+20.9*u.arcsec).to(u.deg), (126*u.deg+27*u.arcmin+34.4*u.arcsec).to(u.deg), 452*u.m
