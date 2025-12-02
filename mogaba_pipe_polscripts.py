
from __main__ import *
from scipy import optimize
import pyclass as p

class PosPolScan:
    def __init__(self,
                 mode=None, scannum=None, npol=None, delay_fit=False,
                 delay=0, c_p=None, g_p=None, station=None,
                 bchan=500, echan=3500
                 ):
        self.c_p        = c_p
        self.g_p        = g_p
        self.station    = station
        self.binnum     = 128
        self.nswitch    = 16
        self.mode       = mode
        self.npol       = npol
        self.delay_fit  = delay_fit
        self.delay      = delay
        self.scannum    = scannum
        self.scannums   = []
        self.version    = 0
        self.scan_tbl   = None
        self.az         = []
        self.el         = []
        self.ra         = 0
        self.dec        = 0
        self.d          = []        # correlation output
        self.c          = []        # Q+iU | Linearly polarized intensity
        self.si         = None      # Stokes I
        self.sv         = None      # Stokes V
        self.cc         = None      # Stokes Q+iU
        self.d_vfc      = []        # VFC data
        self.Tvfc1      = None      # VFC temperature 1 (LL)
        self.Tvfc2      = None      # VFC temperature 2 (RR)
        self.tsys1      = []
        self.tsys2      = []
        self.v          = []        # correlation output of Vane
        self.vc         = []
        self.bad_index  = []
        self.leak       = None
        self.data_unpol = None
        self.data_aref  = None
        self.sideband   = 1
        self.bchan      = bchan
        self.echan      = echan

    def read_scan(self, flag_subchan=None):
        c_p = self.c_p
        g_p = self.g_p
        binnum    = self.binnum
        nswitch   = self.nswitch
        scannum   = self.scannum
        npol      = self.npol
        delay     = self.delay
        delay_fit = self.delay_fit
        # scannums : total 8 sets of scans for 8 repeats
        # scannum  : one scan set in scannums

        """
        get vane infomation on initial scan (On-source position)
        v0 : vane info. on the first scan (LL, RR, Q, U)
        """
        self.v0, self.v0_vfc, c_p, g_p = read_vane(c_p, g_p, scannum)

        nch = len(self.v0[0])
        chans = np.arange(nch)
        Lpol = self.v0[2] + 1j * self.v0[3]  # Linearly polarized intensity

        c_p("get %s"%(scannum + 4 * npol))
        tsys1 = float(g_p.tsys) # LL

        c_p("get %s"%(scannum + 4 * npol + 1))
        tsys2 = float(g_p.tsys) # RR

        if g_p.image > g_p.frequency:
            self.sideband = -1   # lower side band
        else:
            self.sideband = +1   # upper side band

        if delay_fit:
            delay = fit_delay(Lpol)

        if delay != 0:
            phs = np.arange(nch) * delay
            Lpol = Lpol * exp(-1j*phs)            # rotate phase as amount of delay
            self.v0[2] = real(Lpol)
            self.v0[3] = imag(Lpol)

        if binnum > 0:
            self.v = conv(self.v0, binnum)[:, binnum: -binnum]
        else:
            self.v = self.v0
        self.vc = self.v[2] + 1j * self.v[3]

        self.m0, self.r0, self.m0_vfc, self.r0_vfc, self.az, self.el, c_p, g_p = read_wps(c_p, g_p, scannum+4*npol, npol, nswitch, self.mode)

        if flag_subchan == "all" or isinstance(flag_subchan, list):
            if flag_subchan == "all":
                flag_subchan = np.arange(self.m0.shape[0])

            mask_flag = np.zeros(self.m0.shape[0], dtype=bool)
            mask_flag[flag_subchan] = True
            self.m0[mask_flag,:,:] *= np.nan
            self.r0[mask_flag,:,:] *= np.nan
            self.m0_vfc[mask_flag,:,:] *= np.nan
            self.r0_vfc[mask_flag,:,:] *= np.nan

        self.m0 = np.nanmean(self.m0, axis=0)
        self.r0 = np.nanmean(self.r0, axis=0)
        self.m0_vfc = self.m0_vfc
        self.r0_vfc = self.r0_vfc
        self.az = np.nanmean(self.az)
        self.el = np.nanmean(self.el)

        self.d0 = self.m0 - self.r0
        # m0     : main_mean      (RR, LL, Q, U)
        # r0     : refs_mean      (RR, LL, Q, U)
        # d0     : source_mean    (RR, LL, Q, U)
        # m0_vfc : main_vfc power (RR, LL, Q, U)
        # v0_vfc : refs_vfc power (RR, LL, Q, U)
        # az     : az list
        # el     : el list

        Lpol = self.d0[2] + 1j*self.d0[3]
        if delay != 0:
            phs = np.arange(nch) * delay
            Lpol = Lpol*exp(-1j*phs)
            self.d0[2] = real(Lpol)
            self.d0[3] = imag(Lpol)

        self.d_vfc = self.m0_vfc - self.r0_vfc

        """
        Scaling using Auto Power Spectrum
        """
        if binnum > 0:
            self.d = conv(self.d0, binnum)[:, binnum: -binnum]
        else:
            self.d = self.d0
        self.c = self.d[2] + 1.j*self.d[3]  # definition of linear polarization // L = Q + j*U

        self.ra, self.dec = float(g_p.lam), float(g_p.bet)  # in radian
        self.bin_num = binnum
        self.nswitch = nswitch
        self.scannum = scannum
        self.tau0 = float(g_p.tau_signal)
        self.tchop = float(g_p.tchop)
        self.tcold = float(g_p.tcold)
        self.delay = delay
        self.tsys1 = tsys1
        self.tsys2 = tsys2


class PosPolData:
    def __init__(self,
                 mode=None, scannum=None, npol=None, delay_fit=False,
                 delay=0, c_p=None, g_p=None, station =None,
                 bchan=500, echan=3500
                 ):
        self.c_p        = c_p
        self.g_p        = g_p
        self.station    = station
        self.binnum     = 128
        self.nswitch    = 16
        self.mode       = mode
        self.npol       = npol
        self.delay_fit  = delay_fit
        self.delay      = delay
        self.scannum    = scannum
        self.scannums   = []
        self.version    = 0
        self.scan_tbl   = None
        self.az         = []
        self.el         = []
        self.ra         = 0
        self.dec        = 0
        self.d          = []        # correlation output
        self.c          = []        # Q+iU | Linearly polarized intensity
        self.si         = None      # Stokes I?
        self.sv         = None      # Stokes V
        self.cc         = None      # Stokes Q+iU
        self.d_vfc      = []        # VFC data
        self.Tvfc1      = None      # VFC temperature 1 (LL)
        self.Tvfc2      = None      # VFC temperature 2 (RR)
        self.tsys1      = []        # system temperature at LL
        self.tsys2      = []        # system temperature at RR
        self.v          = []        # correlation output of Vane
        self.vc         = []
        self.bad_index  = []
        self.leak       = None
        self.data_unpol = None
        self.data_aref  = None
        self.sideband   = 1
        self.bchan      = bchan
        self.echan      = echan

        if   self.station=='KPC' : self.ant_n, self.ant_lat, self.ant_lon,  self.ant_height = name_PC, lat_PC, lon_PC,  height_PC
        elif self.station=='KYS' : self.ant_n, self.ant_lat, self.ant_lon,  self.ant_height = name_YS, lat_YS, lon_YS,  height_YS
        elif self.station=='KUS' : self.ant_n, self.ant_lat, self.ant_lon,  self.ant_height = name_US, lat_US, lon_US,  height_US
        elif self.station=='KTN' : self.ant_n, self.ant_lat, self.ant_lon,  self.ant_height = name_TN, lat_TN, lon_TN,  height_TN


    def cal_leak(self):
        if self.mode != 'unpol':
            raise Exception("Mode mismatch in calculating leakage term")
        self.leak = self.c / (self.d[:,0] * self.d[:,1])**0.5

    def cal_stokes(self):
        if self.mode == 'unpol':
            raise Exception("Mode mismatch!")

        if self.data_unpol == None:
            raise Exception("Data for a unpolarized source do not exit!")

        unpol   = self.data_unpol
        nscans  = len(self.scannums)
        unpol_v = np.nanmean(unpol.v, axis=0)
        uv0 = np.nansum(unpol_v[0])
        uv1 = np.nansum(unpol_v[1])

        unpol_d = np.nanmean(unpol.d, axis=0)

        v = np.nansum(self.v, axis=2)
        r0 = (v[:,0]/uv0)**0.5
        r1 = (v[:,1]/uv1)**0.5

        c_phase = np.exp(1.j*(np.angle(np.nanmean(unpol.vc, axis=0)) - np.angle(self.vc)))
        c_phase_mean = np.angle(np.nanmean(c_phase[:,self.bchan:self.echan]))

        self._r0 = r0
        self._r1 = r1
        self._c_phase = c_phase

        si = []
        sv = []
        tr = []
        tl = []
        norm1 = []
        norm2 = []
        Tvfc1_unpol = np.abs(np.nanmean(unpol.d_vfc[:,:,0,0]))
        Tvfc2_unpol = np.abs(np.nanmean(unpol.d_vfc[:,:,1,0]))
        Tvfc_mean_unpol = (Tvfc1_unpol * Tvfc2_unpol)**0.5
        for i in range(nscans):
            # norm1_ = 1
            # norm2_ = 1
            # _si = 0.5*(1./r0[i]**2*self.d[i,0]/unpol_d[0]+1./r1[i]**2*self.d[i,1]/unpol_d[1])
            # _sv = 0.5*(1./r0[i]**2*self.d[i,0]/unpol_d[0]-1./r1[i]**2*self.d[i,1]/unpol_d[1])

            norm1_ = np.abs(np.nanmean(1./r0[i]**2*self.d[i,0]/unpol_d[0]) * Tvfc_mean_unpol / np.nanmean(self.d_vfc[i,:,0,0]))
            norm2_ = np.abs(np.nanmean(1./r1[i]**2*self.d[i,1]/unpol_d[1]) * Tvfc_mean_unpol / np.nanmean(self.d_vfc[i,:,1,0]))
            _si = 0.5*(1./r0[i]**2*self.d[i,0]/unpol_d[0]/norm1_ + 1./r1[i]**2*self.d[i,1]/unpol_d[1]/norm2_)
            _sv = 0.5*(1./r0[i]**2*self.d[i,0]/unpol_d[0]/norm1_ - 1./r1[i]**2*self.d[i,1]/unpol_d[1]/norm2_)

            _tr = 1./r0[i]**2*self.d[i,0]/unpol_d[0]
            _tl = 1./r1[i]**2*self.d[i,1]/unpol_d[1]
            si.append(_si)
            sv.append(_sv)
            tr.append(_tr)
            tl.append(_tl)
            norm1.append(norm1_)
            norm2.append(norm2_)



        self.si = np.array(si)
        self.sv = np.array(sv)
        self.tr = np.array(tr)
        self.tl = np.array(tl)

        leak = np.nanmean(unpol.leak, axis=0)
        unpol_c = np.nanmean(unpol.c, axis=0)

        cc = []
        xamp = 1. / (np.nanmean(unpol.d[:,0], axis=0) * np.nanmean(unpol.d[:,1], axis=0))**0.5
        for i in range(nscans):
            _cc = xamp *\
                    (
                        self.c[i] / r0[i] / r1[i] / np.sqrt(np.abs(norm1[i] * norm2[i]))
                        * np.exp(+1j*c_phase_mean) - self.si[i] * unpol_c
                    ) * np.exp(-1j*np.angle(leak))
            cc.append(_cc)
        self.cc = np.array(cc)

    def cal_pangle(self):
        lat = self.ant_lat.value*pi/180
        ra  = self.ra
        dec = self.dec
        az  = self.az
        el  = self.el
        sin_pa = cos(lat)/cos(dec)*sin(az)
        cos_pa = (sin(lat)-sin(dec)*sin(el)) / (cos(dec)*cos(el))
        self.p_angle = angle(cos_pa+1.j*sin_pa)
        self.out_pang = pd.DataFrame([self.p_angle * u.rad.to(u.deg)], columns=["pang{0}".format(i+1) for i in range(len(self.p_angle))])

    def correct_rot(self):
        self.rotation = exp(-1.j*self.p_angle*2)
        for i in range(self.cc.shape[0]):
            ri = self.rotation[i]
            self.cc[i] *= ri


class dsmCal:
    def __init__(self):
        self.NFFT     = 2**14   # the total number of FFT points is set at 4096 (2**14) channels per stream | 512 MHz -> 128 MHz * 4 streams | 4 * 4096 = 2**14
        self.N_sample = 1e+8
        self.sigma    = 1.0
        self.w        = 4

    def cal_Pwr(self, bds00):
        s = self.sigma
        w = self.w
        psi = (0.5-bds00)*2
        p   = float(s)**2/2 / pow(ss.erfinv(psi), 2)
        return p

    def cal_CF(self, v):
        s   = self.sigma
        w   = self.w
        psi = ss.erf(1/np.sqrt(2.)/v)
        E   = np.exp(-0.5/v**2)
        cf  = 2*((w-1)*E+1)**2/np.pi/(psi+w**2*(1-psi))
        return psi, E, cf

    def cal_bds(self, iarea):
        s = self.sigma
        w = self.w
        N_sample=self.N_sample,
        NFFT    =self.NFFT
        b = ((iarea/2/N_sample*NFFT/8) - 0.5/w**2)*w**2/(w**2-1)
        return b

    def Acal_Pwr(self, iarea):
        s = self.sigma
        w = self.w
        b = self.cal_bds(iarea)[0]
        p = self.cal_Pwr(b)
        return p

    def Acal_CF(self, iarea):
        s = self.sigma
        w = self.w
        p = self.Acal_Pwr(iarea)
        psi, E, cf = self.cal_CF(p**0.5)
        return cf

    def Auto_correct(self, ps):
        s     = self.sigma
        w     = self.w
        n     = len(ps)
        iarea = ps.sum()
        p     = self.Acal_Pwr(iarea)
        cf    = self.Acal_CF(iarea)
        scale_factor = 1/cf/iarea*2*n*p
        offset = (1-1/cf)*p
        return ps*scale_factor+offset


dc = dsmCal()

def fit_delay(lpol):
    # pfft = np.fft.fft(lpol)
    # pfft_a = np.absolute(pfft)
    # n = len(lpol)
    # argmax = pfft_a.argmax()

    # if (1 <= argmax) & (argmax <= 4094):
    #     y = pfft_a[argmax-1:argmax+2]
    #     r = np.polyfit(np.arange(-1,2), y, 2)   # 2nd order polynomial fitting
    #     x0 = -r[1]/r[0]/2
    #     delay = argmax+x0
    # else:
    #     if argmax < 1:
    #         delay_offset = 3
    #     else:
    #         delay_offset = -3
    #     phs = 2*np.pi*np.arange(n)/n*delay_offset
    #     Lpol3 = Lpol*exp(-1.j*phs)
    #     delay = fit_delay(Lpol3)+delay_offset

    # if delay > n/2:
    #     delay = delay - n

    lpol = lpol[500:3500]
    n = len(lpol)
    nzp = 100 * n
    zp_lpol = np.pad(lpol, (0, nzp - n), mode="constant")
    fft = np.fft.fft(zp_lpol)
    xout = np.fft.fftfreq(nzp, d=1)
    delay = xout[np.argmax(np.abs(fft))] * 2 * np.pi

    # print(soln.x)
    # fitphs = delay * np.arange(n)
    # fitpol = np.exp(1j * fitphs)
    # fig, axes = plt.subplots(4, 1, figsize=(16, 10))
    # axes[0].plot(np.arange(3000), np.unwrap(np.angle(Lpol[500:3500]),period=np.pi), c="black")
    # axes[0].plot(np.arange(3000), fitphs[500:3500], c="red")
    # axes[1].plot(np.arange(3000), np.unwrap(np.angle(Lpol[500:3500]), period=np.pi) - fitphs[500:3500], c="red")
    # axes[2].plot(np.arange(n), Lpol.real/np.abs(Lpol), c="black")
    # axes[2].plot(np.arange(n), Lpol.imag/np.abs(Lpol), c="red")
    # axes[3].loglog(np.arange(n), np.abs(np.fft.fft(Lpol.real/np.abs(Lpol))), c="black")
    # axes[3].loglog(np.arange(n), np.abs(np.fft.fft(Lpol.imag/np.abs(Lpol))), c="red")
    # fig.tight_layout()
    # plt.show()
    # sys.exit()
    return delay

def conv(pol, bin_num):
    r = []
    for i in range(4):
        r.append(np.convolve(pol[i], [1./bin_num]*bin_num, mode='same'))
    return np.array(r)

def walsh(i):
    bool = 0
    mask = 0x8000
    while mask:
        if mask & i:
            bool += 1
        mask = mask >> 1
    return bool%2

def switching_pattern(n):
    pattern = []
    for i in range(2**n):
        pattern.append(walsh(i))
    return pattern

def read_vane(c_p, g_p, scannum):
    pol     = []
    pol_vfc = []
    bds     = []
    p_auto  = []
    iareas  = []
    for i in range(4):
        c_p("get %s"%(scannum+i)) # las('get %d'%(scannum)) | read LL, RR, Q, U
        bw = fabs(g_p.freq_step)*g_p.channels
        y  = np.array(g_p.ry)*512./bw
        if i < 2:
            n     = len(y)
            iarea = y.sum()
            bds00 = dc.cal_bds(iarea)

            iareas.append(iarea)
            bds   .append(bds00)
            p_auto.append(dc.cal_Pwr(bds00))

            y = dc.Auto_correct(y)
        else:
            p = (p_auto[0]*p_auto[1])**0.5
            psi, E, cf = dc.cal_CF(p**0.5)
            y = y*p/cf/(iareas[0]*iareas[1])**0.5*2*n
        pol    .append(y)
        pol_vfc.append(np.array(g_p.count))
    return np.array(pol), np.array(pol_vfc), c_p, g_p

def read_pol(c_p, g_p, scannum, mode):
    bw      = 512
    pol     = []
    pol_vfc = []
    bds     = []
    p_auto  = []
    iareas  = []
    for i in range(4):
        c_p("get %s"%(scannum+i))
        bw = fabs(g_p.freq_step) * g_p.channels
        y  = np.array(g_p.ry)*512./bw
        if i == 0:
            az  = float(g_p.azimuth)
            el  = float(g_p.elevation)
            tau = float(g_p.tau_signal)/sin(el)
        gain  = (g_p.tchop - g_p.tcold) / (g_p.count[1] - g_p.count[2]) # Kelvin per Voltage (KPV) | count[1], count[2] = hot load Volatage, cold load Voltage
        sigma = 1.0
        if i < 2:
            n     = len(y)
            iarea = y.sum()
            bds00 = dc.cal_bds(iarea)

            iareas.append(iarea)
            bds   .append(bds00)
            p_auto.append(dc.cal_Pwr(bds00))

            y = dc.Auto_correct(y)
        else:
            p = (p_auto[0] * p_auto[1])**0.5
            psi, E, cf = dc.cal_CF(p**0.5)
            y = y / cf / (iareas[0]*iareas[1])**0.5 * 2 * n * p
        pol.append(y*exp(tau))
        pol_vfc.append(np.array(g_p.count)*exp(tau)*gain)
    return np.array(pol), np.array(pol_vfc), az, el, c_p, g_p

def read_wps(c_p, g_p, init_scannum, npol, nswitch, mode):
    norder              = int(np.log2(nswitch))+1
    pos_pattern         = switching_pattern(norder)
    refs_dsm = []
    mains_dsm = []
    refs_vfc = []
    mains_vfc = []
    az_list = []
    el_list = []
    for i,w in enumerate(pos_pattern):
        scannum = init_scannum +i*4*npol
        dsm, p_vfc, az, el, c_p, g_p = read_pol(c_p, g_p, scannum, mode)
        if w == 0:  # on sky
            refs_dsm.append(dsm)
            refs_vfc.append(p_vfc)
        else:       # on source
            mains_dsm.append(dsm)
            mains_vfc.append(p_vfc)
            az_list  .append(az)
            el_list  .append(el)
    mains_dsm = np.array(mains_dsm)
    refs_dsm = np.array(refs_dsm)
    mains_vfc = np.array(mains_vfc)
    refs_vfc = np.array(refs_vfc)
    az_list = np.array(az_list)
    el_list = np.array(el_list)
    return mains_dsm, refs_dsm, mains_vfc, refs_vfc, az_list, el_list, c_p, g_p


name_PC, lat_PC, lon_PC, height_PC = 'PC', (37*u.deg+32*u.arcmin+00.1*u.arcsec).to(u.deg), (128*u.deg+26*u.arcmin+55.1*u.arcsec).to(u.deg), 541*u.m
name_YS, lat_YS, lon_YS, height_YS = 'YS', (37*u.deg+33*u.arcmin+54.9*u.arcsec).to(u.deg), (126*u.deg+56*u.arcmin+27.4*u.arcsec).to(u.deg), 139*u.m
name_US, lat_US, lon_US, height_US = 'US', (35*u.deg+32*u.arcmin+44.2*u.arcsec).to(u.deg), (129*u.deg+14*u.arcmin+59.3*u.arcsec).to(u.deg), 170*u.m
name_TN, lat_TN, lon_TN, height_TN = 'TN', (33*u.deg+17*u.arcmin+20.9*u.arcsec).to(u.deg), (126*u.deg+27*u.arcmin+34.4*u.arcsec).to(u.deg), 452*u.m
