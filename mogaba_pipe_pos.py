
from __main__ import *
import copy
from sklearn.cluster import DBSCAN as dbs
import pyclass as p

from mogaba_pipe_imports import cal_gain_curve

c_ps = p.comm
g_ps = p.gdict

c_ps("set variable general")
c_ps("set variable calibration")
c_ps("set variable position")

# pps = PosPolScan()    : temporary class to read scan infomation
# ppd = PosPolData()    : fianl class to assign polarization data
# Main stream           : scan-by-scan data                       (dsmCal & read_wps)
#                           -> subset of scans                    (pps)
#                           -> overall data for individual source (ppd)
#                           -> load to PoSwitch

class PoSwitch:
    def __init__(self,
                 path_p=None, path_c=None, path_dir=None, saveplot=False,
                 file=None, station=None, telname=None, unpol=None,
                 unpol_lst=None, unpol_all=None, polnum=None, aref_n=None,
                 mode=None, lr_swap=None, autoflag=False, pipe_log=None,
                 bchan=500, echan=3500
                 ):
        self.path_p   = path_p          # path used in python
        self.path_c   = path_c          # path used in class
        self.path_dir = path_dir        # path to load Tb excel file
        self.file     = file            # sdd file name
        self.station  = station         # station name (e.g., KYS)
        self.saveplot = saveplot        # path to save plots (Tsys, Tau, AzEl, ..)
        self.unpol    = unpol           # unpol source name
        self.unpol_lst= unpol_lst       # list of unpol sources in given data
        self.unpol_all= unpol_all       # list of available sources as unpol
        self.aref_n   = aref_n          # name of angle reference source
        self.polnum   = polnum          # polnum (e.g., 1:22 GHz , 2:43 GHz for KQ data)
        self.mode     = mode            # data reduction mode ("unpol", "aref", "target")
        self.binnum   = 128
        self.nsetup   = 4 * 2           # N(Stokes | RR, LL, Q, U) * Nfreq (e.g., 22, 43)
        self.nswitch  = 16              # number of switching (off-on-on- ... | refer to Kang et al., 2015, JKAS, 48, 257)
        self.nonoff   = 2               # on- or off-source position
        self.pos_sour = None            # source name for data processing
        self.lr_swap  = lr_swap         # toggle whether to swap L-R pol into R-L order
        self.autoflag = autoflag
        self.pipe_log = pipe_log
        self.errmsg   = None
        self.out_pang = pd.DataFrame([])
        self.bchan    = bchan
        self.echan    = echan

    def set_init(self):
        if self.station is None : self.station=self.file.split(".sdd")[0].split("_")[-1]
        self.proj = self.file.split("_%s"%(self.station))[0]

        path_fig = self.path_dir + "Figures/"
        path_dat = self.path_dir + "data_pos/"
        mkdir("%s"         %(path_fig))
        mkdir("%s/pos/"    %(path_fig))
        mkdir("%s/pos/22/" %(path_fig))
        mkdir("%s/pos/43/" %(path_fig))
        mkdir("%s/pos/86/" %(path_fig))
        mkdir("%s/pos/94/" %(path_fig))
        mkdir("%s/pos/129/"%(path_fig))
        mkdir("%s/pos/141/"%(path_fig))

        mkdir("%s"     %(path_dat))
        mkdir("%s/22/" %(path_dat))
        mkdir("%s/43/" %(path_dat))
        mkdir("%s/86/" %(path_dat))
        mkdir("%s/94/"%(path_dat))
        mkdir("%s/129/"%(path_dat))
        mkdir("%s/141/"%(path_dat))

        try:
            c_ps("sic directory %s"%(self.path_c))
        except:
            self.errmsg = "    !!! Given path '%s' does not exist."%(self.path_c)
            raise Exception(self.errmsg)
        try:
            c_ps("file in %s"%(self.file))
        except:
            self.errmsg = "    !!! Given sdd file '%s' does not exist."%(self.file)
            raise Exception(self.errmsg)
        c_ps("set default")
        c_ps("set variable general")
        c_ps("set ty l")
        c_ps("find")

    def load_source_info(self):
        Nscans  = int(g_ps.found)
        sources = np.array([]).astype(str)
        scans   = np.array([]).astype(int)
        dtime   = np.array([]).astype(str)
        tsys    = np.array([]).astype(float)
        tau     = np.array([]).astype(float)
        az      = np.array([]).astype(float)
        el      = np.array([]).astype(float)
        freq    = np.array([]).astype(int)
        stokes  = np.array([]).astype(str)
        for n in range(Nscans):
            if n==0: c_ps("get first")
            else   : c_ps("get next")
            c_ps("define character*12 sour")
            c_ps("define character*12 sourdate")
            c_ps("let sour     'r%head%pos%sourc'")
            c_ps("let sourdate 'r%head%gen%cdobs'")
            date     = format_date(str(g_ps.sourdate.__sicdata__))
            datetime = format_time(date, g_ps.ut.astype(float))
            source   = str(g_ps.sour.__sicdata__ ).replace("b","").replace("'","").replace(" ","")
            teles    = str(g_ps.teles.astype(str)).replace("b","").replace("'","").replace(" ","")
            S = teles[-1]   # Stokes parameter(LL, RR, Q, U)
            f = int(g_ps.line.astype(float)/1000)
            if source in ["IK_TAU", "V1111_OPH", "R_LEO", "U_ORI", "ORION_KL", "TX_CAM", "OMI_CET", "VY_CMA", "R_HYA", "W_HYA", "U_HER", "VX_SGR", "R_AQL", "R_AQR", "R_CAS", "T_CEP"]:
                delattr(g_ps, "sour")
                delattr(g_ps, "sourdate")
                continue

            sources = np.append(sources, source)
            scans   = np.append(scans  , g_ps.scan.astype(int))
            dtime   = np.append(dtime  , str(datetime))
            freq    = np.append(freq   , f)
            tsys    = np.append(tsys   , g_ps.tsys.astype(float))
            tau     = np.append(tau    , g_ps.tau .astype(float))
            az      = np.append(az     , g_ps.az  .astype(float)*u.rad.to(u.deg))
            el      = np.append(el     , g_ps.el  .astype(float)*u.rad.to(u.deg))
            stokes  = np.append(stokes , S)
            delattr(g_ps, "sour")
            delattr(g_ps, "sourdate")
        npol = len(np.unique(freq))
        pos_log_all = pd.DataFrame([sources, scans, dtime, tsys, tau, az, el, freq, stokes],
                                  index=["Source", "ScanNum", "Date", "Tsys", "Tau", "Az", "El", "Freq", "Stokes"]).transpose()
        if not pos_log_all.shape[0]<2112:
            pos_log_all["Year"] = Ati(np.array(pos_log_all["Date"]).astype(str), format="iso").byear
            pos_log_all["MJD"]  = Ati(np.array(pos_log_all["Date"]).astype(str), format="iso").mjd
            pos_log_all = pos_log_all[["Date", "Year", "MJD", "Source", "ScanNum", "Tsys", "Tau", "Az", "El", "Freq", "Stokes"]]

            self.pos_log_all  = pos_log_all
            mjd0  = int(np.min(self.pos_log_all["MJD"]))
            nchan = int(npol*4)
            self.pos_log_all["Time"]         =np.round((self.pos_log_all["MJD"]-mjd0) * u.day.to(u.s),0).astype(int)
            self.pos_log_all["Time"][nchan:] =np.array(self.pos_log_all["Time"][nchan:]) -np.array(self.pos_log_all["Time"][0:-nchan])
            self.pos_log_all["Time"][0:nchan]=np.array(self.pos_log_all["Time"][0:nchan])-self.pos_log_all["Time"][0]
            self.pos_log_all = correct_time(self.pos_log_all, nchan)

            time_error = self.pos_log_all[np.abs(self.pos_log_all.Time)<120]["Time"]
            thresh_avg, thresh_err = np.mean(time_error), np.std(time_error)

            gap_ = np.where(self.pos_log_all["Time"]>thresh_avg+10*thresh_err)[0][0::nchan]
            gap_ = np.append(np.array([0]), gap_)
            pos_log_sour = self.pos_log_all.iloc[gap_].reset_index()

            check_sour = lambda x, series : x in series
            list_avg   = np.array(pos_log_sour["Source"]).astype(str)
            list_all   = np.array(self.pos_log_all["Source"]).astype(str)
            CRABs    = ["CRAB", "CRAB1", "CRAB2", "3C279"]
            row_crab = []
            for crab in CRABs:
                if not np.logical_and(check_sour(crab, list_avg), check_sour(crab, list_all)):
                    row_ = np.where(self.pos_log_all.Source==crab)[0]
                    if row_.shape[0]!=0:
                        row_crab.append(row_[0])
            if len(row_crab)!=0:
                crab_dat = pd.DataFrame(self.pos_log_all.iloc[row_crab]).reset_index()
                pos_log_sour = pd.concat([pos_log_sour, crab_dat], axis=0).sort_values(by="index").reset_index(drop=True)
            index1 = np.append(np.array(pos_log_sour["index"][1:]), np.array([self.pos_log_all.shape[0]]))
            index2 = np.array(pos_log_sour["index"])
            pos_log_sour["Nscan"] = np.array(index1-index2).astype(int)
            pos_log_sour = pos_log_sour[["Source", "Date", "Year", "MJD", "ScanNum", "Nscan", "Az", "El"]]
            pos_log_sour["Tsys_1"] = np.full(pos_log_sour.shape[0], np.nan)
            pos_log_sour["Tsys_2"] = np.full(pos_log_sour.shape[0], np.nan)
            pos_log_sour["Tau_1" ] = np.full(pos_log_sour.shape[0], np.nan)
            pos_log_sour["Tau_2" ] = np.full(pos_log_sour.shape[0], np.nan)
            ufreq = np.unique(freq)
            for nsour, source in enumerate(pos_log_sour["Source"]):
                scannum1, scannum2 = pos_log_sour["ScanNum"][nsour], pos_log_sour["ScanNum"][nsour]+pos_log_sour["Nscan"][nsour]
                pla  = pos_log_all[pos_log_all.Source==source].reset_index(drop=True)
                pla  = pla[np.logical_and(scannum1<=pla.ScanNum, pla.ScanNum<scannum2)].reset_index(drop=True)
                pla1 = pla[pla.Freq==ufreq[0]].reset_index(drop=True)
                pos_log_sour["Tsys_1"][nsour] = np.median(pla1["Tsys"])
                pos_log_sour["Tau_1" ][nsour] = np.median(pla1["Tau" ])
                if len(ufreq)==2:
                    pla2 = pla[pla.Freq==ufreq[1]].reset_index(drop=True)
                    pos_log_sour["Tsys_2"][nsour] = np.median(pla2["Tsys"])
                    pos_log_sour["Tau_2" ][nsour] = np.median(pla2["Tau" ])

            pos_log_sour["Nrep"] = np.array(pos_log_sour["Nscan"]//(self.nsetup * (self.nswitch*self.nonoff+1))).astype(int)
            lowrep_      = np.where(pos_log_sour["Nrep"] < 4)[0]
            pos_log_sour = pos_log_sour.drop(lowrep_, axis=0).reset_index(drop=True)
            pos_log_sour["Nswitch"] = [self.nswitch for i in range(pos_log_sour.shape[0])]

            self.pos_log_sour = pos_log_sour

            log_all = self.file.split(".sdd")[0] + "_All.xlsx"
            self.pos_log_all .to_excel(self.path_p + log_all )
            self.pos_log_sour.to_excel(self.path_p + self.log)
        else:
            self.pos_log_all = pos_log_all

    def remake_log_target(self):
        log_target = self.pos_log_sour.copy()
        log_target = log_target[np.logical_and(log_target.ScanNum!=self.scan_aref,
                                               log_target.ScanNum!=self.scan_unpol)].reset_index(drop=True)
        self.log_target = log_target

    def get_info_calib(self):
        self.date = self.pos_log_all["Date"][0].split(" ")[0]
        Nsour     = self.pos_log_sour.shape[0]
        Sources   = np.array(self.pos_log_sour["Source"]).astype(str)
        Sources   = np.char.upper(Sources)
        Unpol_all = np.char.upper(np.array(self.unpol_all))
        unpols_n  = []  # unpol sources name
        for N in range(Nsour):
            if Sources[N] in Unpol_all:
                unpols_n.append(self.pos_log_sour["Source"][N])
        self.unpols_n = list(np.unique(unpols_n))

        if self.aref_n in Sources:
            self.aref_s = self.pos_log_sour[self.pos_log_sour==self.aref_n]["ScanNum"]
            self.aref_r = self.pos_log_sour[self.pos_log_sour==self.aref_n]["Nrep"]
        else:
            if "CRAB2" in Sources:
                self.aref_n = "CRAB2"
                self.aref_s = self.pos_log_sour[self.pos_log_sour==self.aref_n]["ScanNum"]
                self.aref_r = self.pos_log_sour[self.pos_log_sour==self.aref_n]["Nrep"]
            elif "CRAB1" in Sources:
                self.aref_n = "CRAB1"
                self.aref_s = self.pos_log_sour[self.pos_log_sour==self.aref_n]["ScanNum"]
                self.aref_r = self.pos_log_sour[self.pos_log_sour==self.aref_n]["Nrep"]
            elif "3C279" in Sources:
                self.aref_n = "3C279"
                self.aref_s = self.pos_log_sour[self.pos_log_sour==self.aref_n]["ScanNum"]
                self.aref_r = self.pos_log_sour[self.pos_log_sour==self.aref_n]["Nrep"]
            else:
                self.errmsg = "There is no Pol. Angle reference source (CRAB, CRAB1, CRAB2)."
                print(self.errmsg)
                print("File name : %s"%(self.file))
                print("End Process")
                writelog(self.path_dir, self.pipe_log, "No Pol. Angle reference source in %s"%(self.file), "a")

    def get_info_scan(self):
        source  = self.pos_sour
        nswitch = self.nswitch
        nsetup  = self.npol*4 ; self.nsetup = nsetup
        nonoff  = self.nonoff
        scan1r  = nsetup*(nswitch*nonoff+1)
        polnum  = self.polnum
        if self.npol > 2:
            out_txt =\
                "File Name: %s"%(self.file) \
                + "\nNumber of frequency exceeds two.\n"
            raise ValueError(out_txt)
        if self.mode != "target":
            scan_ = self.pos_log_sour.copy()

            if self.mode == "unpol":
                mask = scan_.ScanNum == self.scan_unpol

            if self.mode == "aref" :
                mask = scan_.ScanNum == self.scan_aref

            scani = scan_.loc[mask].reset_index(drop=True)
            nrep = scani["Nrep"].values[0]
            scannum = scani["ScanNum"].values[0]
            self.scannums  = [scannum + scan1r * n + 4 * (polnum - 1) for n in range(nrep)]
        elif self.mode =="target":
            scan_ = self.log_target.copy()
            scani = scan_[scan_.Source == source].reset_index(drop=True)
            scani = scani.iloc[self.pos_nseq]
            nrep  = scani["Nrep"]
            self.scannums  = [scani["ScanNum"] + scan1r * n + 4 * (polnum - 1) for n in range(nrep)]
        self.scan_info = scani
        self.nrepeat   = nrep


    def get_bad_scans(self):
        bad_scans = self.bad_chans
        nflag = len(bad_scans)

        scannums = []
        chans = []
        subchans = []
        for Nscan, flag_info in enumerate(bad_scans):
            scannum = flag_info[0]
            chan = list(flag_info[1].keys())
            subchan = list(flag_info[1].values())
            scannums.append(scannum)
            chans.append(chan)
            subchans.append(subchan)
        self.flag_info = pd.DataFrame(
            [scannums, chans, subchans],
            index=["scan", "chan", "subchan"]
        ).transpose()

    def get_poldata(self):
        ppd = self.ppd
        mode        = self.mode
        scannums    = self.scannums
        nscans      = len(scannums)
        sqrtn       = sqrt(nscans)
        ch1, ch2    = self.bchan, self.echan
        stat_method = np.nanmean
        elevation   = self.elevation
        year        = int(self.date[:4])
        self.gain   = cal_gain_curve(self.station, year, self.freq, elevation)

        if mode != "unpol":
            cc = ppd.cc
            si = ppd.si
            sv = ppd.sv

            Tvfc1, Tvfc2 = self.data_unpol.Tvfc1, self.data_unpol.Tvfc2
            v1, v2 = stat_method(Tvfc1), stat_method(Tvfc2)
            v_mean = (v1*v2)**0.5
            sa     = stat_method(si[:, ch1:ch2], axis=1)*v_mean # si  : I/I_0
            va     = stat_method(sv[:, ch1:ch2], axis=1)        # sv  : V/I

            sa_mean = stat_method(sa)
            si_mean = sa_mean/v_mean

            cc_angle   = angle(cc)
            cc_phs     = stat_method(cc[:, ch1:ch2], axis=1)

            angle_mean = angle(stat_method(cc_phs))
            angle_std  = np.nanstd(angle(cc_phs*exp(-1j*angle_mean)))

            if ppd.data_aref:
                data_aref    = self.data_aref
                data_aref_cc = data_aref.cc
                aref_phs     = exp(1.j*angle(stat_method(data_aref_cc, axis=0)))*exp(-1.j*angle(stat_method(data_aref_cc[:,ch1:ch2])))

                cc_phs_mean = exp(1.j*angle_mean)

                cca   = cc/cc_phs_mean/aref_phs
                cca_r = real(cca)[:, ch1:ch2]
                cca_i = imag(cca)[:, ch1:ch2]
                pm    = np.abs(stat_method(cca[:, ch1:ch2]))/si_mean
                dpm   = sqrt(np.nanstd(stat_method(cca_r, axis=1))**2+np.nanstd(stat_method(cca_i, axis=1))**2)/sqrtn/np.abs(si_mean)/sqrt(2)
                tp    = np.abs(stat_method(stat_method(cca[:, ch1:ch2]*v_mean, axis=1)))
                dtp   = np.nanstd(stat_method(cca[:, ch1:ch2]*v_mean, axis=1))/sqrtn

                cc_phs = stat_method(cca[:, ch1:ch2] * cc_phs_mean, axis=1)
                self.cal_cc = cca * cc_phs_mean
            else:
                ca  = stat_method(np.abs(cc)[:,ch1:ch2], axis=1)
                pm  = stat_method(ca/si_mean)
                dpm = np.nanstd(ca/si_mean)/sqrtn
                tp  = np.abs(stat_method(ca*v_mean))
                dtp = np.nanstd(ca*v_mean)/sqrtn
                self.cal_cc = cc

            out_angle_mean = angle(stat_method(cc_phs))
            out_angle_std  = np.nanstd(angle(cc_phs*exp(-1j*out_angle_mean)))


            vm = va/sa_mean*100
            out_angle_mean *= 180/pi/2
            out_angle_std  *= 180/pi/2

            self.getp_ti, self.getp_dti = sa_mean/self.gain   , np.nanstd(sa)/sqrtn/self.gain
            self.getp_tp, self.getp_dtp = tp/self.gain        , dtp/self.gain
            self.getp_tv, self.getp_dtv = stat_method(vm), np.nanstd(vm)/sqrtn
            self.getp_pm, self.getp_dpm = pm             , dpm
            self.getp_pa, self.getp_dpa = out_angle_mean , out_angle_std/sqrtn
        self.getp_t1, self.getp_dt1 = stat_method(ppd.Tvfc1)/self.gain, np.nanstd(ppd.Tvfc1)/self.gain
        self.getp_t2, self.getp_dt2 = stat_method(ppd.Tvfc2)/self.gain, np.nanstd(ppd.Tvfc2)/self.gain

    def rd_polscans(self):
        self.get_info_scan()
        polnum       = self.polnum
        binnum       = self.binnum
        mode         = self.mode
        nswitch      = self.nswitch
        npol         = self.npol
        scani        = self.scan_info
        scannums     = self.scannums
        sour_scannum = scani["ScanNum"]

        ppd = PosPolData(station=self.station, bchan=self.bchan, echan=self.echan)
        ppd.scannums = scannums
        ppd.c_p      = c_ps
        ppd.g_p      = g_ps
        ppd.station  = self.station
        ppd.binnum   = binnum
        ppd.polnum   = polnum
        ppd.mode     = mode
        ppd.nswitch  = nswitch
        ppd.npol     = npol

        tsys1, tsys2 = [],[]

        self.get_bad_scans()
        flag_info = self.flag_info.copy()
        flag_info = flag_info[flag_info.scan == scannums[0]].reset_index(drop=True)

        for iscan, scannum in enumerate(ppd.scannums):
            if self.autoflag:
                flag_subchan = "none"
            else:
                if flag_info.empty:
                    flag_subchan = "none"
                else:
                    flag_chan_ = flag_info["chan"][0]
                    if iscan in flag_chan_:
                        loc = np.where(np.array(flag_chan_) == iscan)[0][0]
                        flag_subchan = flag_info["subchan"][0][loc]
                    else:
                        flag_subchan = "none"

            if np.logical_and(mode=="unpol", iscan==0):
                pps = PosPolScan(mode=mode, scannum=scannum, npol=npol, delay_fit=True, delay=0, station=self.station, bchan=self.bchan, echan=self.echan)
                pps.c_p        = c_ps
                pps.g_p        = g_ps
                pps.read_scan(flag_subchan=flag_subchan)
                ppd.delay      = pps.delay
                ppd.sideband   = pps.sideband
            elif np.logical_and(mode=="unpol", iscan!=0):
                pps.scannum    = scannum
                pps.delay      = ppd.delay
                pps.delay_fit  = False
                pps.read_scan(flag_subchan=flag_subchan)

            if np.logical_and(mode!="unpol", iscan==0):
                pps = PosPolScan(mode=mode, scannum=scannum, npol=npol, delay_fit=False, delay=self.ppd.delay, station=self.station, bchan=self.bchan, echan=self.echan)
                pps.c_p        = c_ps
                pps.g_p        = g_ps
                pps.read_scan(flag_subchan=flag_subchan)
                ppd.sideband   = pps.sideband
            elif np.logical_and(mode!="unpol", iscan!=0):
                pps.scannum    = scannum
                pps.read_scan(flag_subchan=flag_subchan)

            ppd.az   .append(pps.az)
            ppd.el   .append(pps.el)
            ppd.d    .append(pps.d)
            ppd.c    .append(pps.c)
            ppd.v    .append(pps.v)
            ppd.vc   .append(pps.vc)    # convolved vane information
            ppd.d_vfc.append(pps.d_vfc)
            ppd.tsys1.append(pps.tsys1)
            ppd.tsys2.append(pps.tsys2)

        ppd.ra, ppd.dec = pps.ra, pps.dec
        ppd.az    = np.array(ppd.az)
        ppd.el    = np.array(ppd.el)
        ppd.d     = np.array(ppd.d)
        ppd.c     = np.array(ppd.c)
        ppd.v     = np.array(ppd.v)
        ppd.vc    = np.array(ppd.vc)
        ppd.d_vfc = np.array(ppd.d_vfc)
        ppd.Tvfc1 = ppd.d_vfc[:,:,0,0]
        ppd.Tvfc2 = ppd.d_vfc[:,:,1,0]
        ppd.tsys1 = np.array(ppd.tsys1)
        ppd.tsys2 = np.array(ppd.tsys2)

        if self.autoflag:
            Tvfc1_90l = np.percentile(ppd.Tvfc1, self.perc_tvfc)
            Tvfc1_90h = np.percentile(ppd.Tvfc1, 100 - self.perc_tvfc)
            Tvfc2_90l = np.percentile(ppd.Tvfc2, self.perc_tvfc)
            Tvfc2_90h = np.percentile(ppd.Tvfc2, 100 - self.perc_tvfc)
            mask1_l = ppd.Tvfc1 < Tvfc1_90l
            mask1_h = ppd.Tvfc1 > Tvfc1_90h
            mask2_l = ppd.Tvfc2 < Tvfc2_90l
            mask2_h = ppd.Tvfc2 > Tvfc2_90h
            mask1 = (mask1_l) | (mask1_h)
            mask2 = (mask2_l) | (mask2_h)

            ppd.az    = []
            ppd.el    = []
            ppd.d     = []
            ppd.c     = []
            ppd.v     = []
            ppd.vc    = []
            ppd.d_vfc = []
            ppd.Tvfc1 = None
            ppd.Tvfc2 = None
            ppd.tsys1 = []
            ppd.tsys2 = []
            chans     = []
            subchans  = []

            for iscan, scannum in enumerate(ppd.scannums):
                mask = (mask1[iscan]) | (mask2[iscan])
                flag_subchan = np.where(mask)[0].tolist()

                if not flag_info.empty:
                    flag_chan_ = flag_info["chan"][0]
                    if iscan in flag_chan_:
                        loc = np.where(np.array(flag_chan_) == iscan)[0][0]
                        flag_subchan_ = flag_info["subchan"][0][loc]
                        if flag_subchan_ == "all":
                            flag_subchan_ = [isubchan for isubchan in range(nswitch)]
                        flag_subchan = flag_subchan + flag_subchan_
                        flag_subchan = list(set(flag_subchan))
                        flag_subchan.sort()

                if not flag_subchan:
                    flag_subchan = "none"

                if np.logical_and(mode=="unpol", iscan==0):
                    pps = PosPolScan(
                        mode=mode, scannum=scannum, npol=npol, delay_fit=True, delay=0,
                        station=self.station, bchan=self.bchan, echan=self.echan
                    )
                    pps.c_p        = c_ps
                    pps.g_p        = g_ps
                    pps.read_scan(flag_subchan=flag_subchan)
                    ppd.delay      = pps.delay
                    ppd.sideband   = pps.sideband
                elif np.logical_and(mode=="unpol", iscan!=0):
                    pps.scannum    = scannum
                    pps.delay      = ppd.delay
                    pps.delay_fit  = False
                    pps.read_scan(flag_subchan=flag_subchan)

                if np.logical_and(mode!="unpol", iscan==0):
                    pps = PosPolScan(
                        mode=mode, scannum=scannum, npol=npol, delay_fit=False, delay=self.ppd.delay,
                        station=self.station, bchan=self.bchan, echan=self.echan
                    )
                    pps.c_p        = c_ps
                    pps.g_p        = g_ps
                    pps.read_scan(flag_subchan=flag_subchan)
                    ppd.sideband   = pps.sideband
                elif np.logical_and(mode!="unpol", iscan!=0):
                    pps.scannum    = scannum
                    pps.read_scan(flag_subchan=flag_subchan)

                ppd.az   .append(pps.az)
                ppd.el   .append(pps.el)
                ppd.d    .append(pps.d)
                ppd.c    .append(pps.c)
                ppd.v    .append(pps.v)
                ppd.vc   .append(pps.vc)    # convolved vane information
                ppd.d_vfc.append(pps.d_vfc)
                ppd.tsys1.append(pps.tsys1)
                ppd.tsys2.append(pps.tsys2)

            ppd.ra, ppd.dec = pps.ra, pps.dec
            ppd.az    = np.array(ppd.az)
            ppd.el    = np.array(ppd.el)
            ppd.d     = np.array(ppd.d)
            ppd.c     = np.array(ppd.c)
            ppd.v     = np.array(ppd.v)
            ppd.vc    = np.array(ppd.vc)
            ppd.d_vfc = np.array(ppd.d_vfc)
            ppd.Tvfc1 = ppd.d_vfc[:,:,0,0]
            ppd.Tvfc2 = ppd.d_vfc[:,:,1,0]
            ppd.tsys1 = np.array(ppd.tsys1)
            ppd.tsys2 = np.array(ppd.tsys2)

        if mode == "unpol":
            self.data_unpol = ppd
        else:
            ppd.unpol_data = self.data_unpol
            self.pol_data = ppd
            if mode == "aref": self.data_aref=ppd
            else             : ppd.data_aref =self.data_aref

        if pps.delay != ppd.delay:
            ppd.delay = pps.delay

        self.pps = pps
        self.ppd = ppd

    def run_pos(self):
        self.mode = "unpol"
        self.pos_sour = self.unpol
        log_unpol = self.pos_log_sour[self.pos_log_sour.Source == self.pos_sour].reset_index(drop=True)
        mask_nrep = log_unpol["Nrep"] >= 4
        log_unpol = log_unpol.loc[mask_nrep].reset_index(drop=True)
        if self.polnum==1:
            mask_tsys = log_unpol.Tsys_1 == np.min(log_unpol.Tsys_1)
            mask_tau = log_unpol.Tau_1 == np.min(log_unpol.Tau_1)
        if self.polnum==2:
            mask_tsys = log_unpol.Tsys_2 == np.min(log_unpol.Tsys_2)
            mask_tau = log_unpol.Tau_2 == np.min(log_unpol.Tau_2)
        mask_unpol = mask_tau
        self.mask_unpol = mask_unpol
        self.mjd = log_unpol.loc[mask_unpol]["MJD"].values[0]
        self.nswitch = log_unpol.loc[mask_unpol]["Nswitch"].values[0]
        self.scannum = log_unpol.loc[mask_unpol]["ScanNum"].values[0]
        self.scan_unpol = log_unpol.loc[mask_unpol]["ScanNum"].values[0]
        self.elevation = log_unpol.loc[mask_unpol]["El"].values[0]
        self.rd_polscans()
        self.ppd.mode     = self.mode
        self.ppd.pos_sour = self.pos_sour
        self.ppd.cal_pangle()
        self.out_pang = pd.concat([self.out_pang, self.ppd.out_pang], axis=0).reset_index(drop=True)
        self.ppd.cal_leak()
        self.get_poldata()
        self.PlotPolDat()
        self.SavePolDat()

        self.mode = "aref"
        self.pos_sour = self.aref_n
        log_aref = self.pos_log_sour[self.pos_log_sour.Source==self.pos_sour].reset_index(drop=True)
        mask_nrep = log_aref["Nrep"] > 4
        if self.polnum==1:
            mask_tsys = log_aref.Tsys_1 == np.min(log_aref.Tsys_1)
            mask_tau = log_aref.Tau_1 == np.min(log_aref.Tau_1)
        if self.polnum==2:
            mask_tsys = log_aref.Tsys_2 == np.min(log_aref.Tsys_2)
            mask_tau = log_aref.Tau_2 == np.min(log_aref.Tau_2)
        mask_aref = mask_tau
        self.mask_aref = mask_aref
        self.mjd = log_aref.loc[mask_aref]["MJD"].values[0]
        self.nswitch = log_aref.loc[mask_aref]["Nswitch"].values[0]
        self.scannum = log_aref.loc[mask_aref]["ScanNum"].values[0]
        self.scan_aref = log_aref.loc[mask_aref]["ScanNum"].values[0]
        self.elevation = log_aref.loc[mask_aref]["El"].values[0]
        self.rd_polscans()
        self.ppd.data_unpol = self.data_unpol
        self.ppd.mode       = self.mode
        self.ppd.pos_sour   = self.pos_sour
        self.ppd.cal_pangle()
        self.out_pang = pd.concat([self.out_pang, self.ppd.out_pang], axis=0).reset_index(drop=True)
        self.ppd.cal_stokes()
        sign = 1
        if self.lr_swap:
            sign = -1
        self.ppd.sv *= sign
        sideband = self.ppd.sideband
        if sign*sideband==-1 : self.ppd.cc = np.conj(self.ppd.cc)
        self.ppd.correct_rot()
        self.get_poldata()
        self.PlotPolDat()
        self.SavePolDat()

        self.remake_log_target()
        count = np.array([])
        for ntarget in range(self.log_target.shape[0]):
            self.mode      = "target"
            self.pos_sour  = self.log_target["Source" ][ntarget]
            self.mjd       = self.log_target["MJD"    ][ntarget]
            self.nswitch   = self.log_target["Nswitch"][ntarget]
            self.scannum   = self.log_target["ScanNum"][ntarget]
            self.elevation = self.log_target["El"     ][ntarget]
            self.pos_nseq  = count[count==self.pos_sour].shape[0]
            count = np.append(count, self.pos_sour)
            self.rd_polscans()
            self.ppd.data_unpol = self.data_unpol
            self.ppd.mode       = self.mode
            self.ppd.pos_sour   = self.pos_sour
            self.ppd.cal_pangle()
            self.out_pang = pd.concat([self.out_pang, self.ppd.out_pang], axis=0).reset_index(drop=True)
            self.ppd.cal_stokes()

            sign = 1
            if self.lr_swap : sign = -1
            self.ppd.sv *= sign
            sideband = self.ppd.sideband
            if sign*sideband==-1 : self.ppd.cc = np.conj(self.ppd.cc)
            self.ppd.correct_rot()

            self.get_poldata()
            self.PlotPolDat()
            self.SavePolDat()

        freq = self.freq
        path_dat = self.path_dir + "data_pos/%s/"%(self.freq)
        mkdir(path_dat)
        self.poldat.to_excel(path_dat + "%s_%s_%s_%s_%s.xlsx"%(self.station, self.date, self.freq, self.unpol, self.proj))


    def PlotPolDat(self):
        if   self.mode !="target": log_sour = self.pos_log_sour
        elif self.mode =="target": log_sour = self.log_target
        mode     = self.mode
        ppd      = self.ppd
        scannums = ppd.scannums
        log_scan = log_sour[log_sour.ScanNum == self.scannum].reset_index(drop=True)
        date_time= log_scan["Date"][0]
        nrep     = log_scan["Nrep"][0]
        az       = log_scan["Az"][0]
        el       = log_scan["El"][0]
        nswitch  = self.nswitch
        Tintg    = 3.0
        ch1, ch2 = self.bchan, self.echan

        if mode!="unpol":
            si_mean = np.nanmean(ppd.si[:, ch1:ch2])

        if self.SaveACPlot:
            mkdir(self.path_dir + "Figures/pos_ac/")
            mkdir(self.path_dir + "Figures/pos_ac/%s/"         %(self.freq))
            mkdir(self.path_dir + "Figures/pos_ac/%s/%s_%s/"   %(self.freq, self.station, self.date))
            mkdir(self.path_dir + "Figures/pos_ac/%s/%s_%s/%s/"%(self.freq, self.station, self.date, self.unpol))
            fig_ph, ax_ph = plt.subplots(2,2, figsize=(10,10))
            for i in range(nrep):
                ax_ph[0,0].plot(ppd.d[i,0])
                ax_ph[0,1].plot(ppd.d[i,1])
                ax_ph[1,0].plot(ppd.d[i,2])
                ax_ph[1,1].plot(ppd.d[i,3])
            fig_ph.suptitle("Power Spectrum\n%s (%s GHz)"%(self.pos_sour, self.freq), fontsize=20)
            ax_ph[0,0].set_title("%sRR"%(self.freq), fontsize=16)
            ax_ph[0,1].set_title("%sLL"%(self.freq), fontsize=16)
            ax_ph[1,0].set_title("%sQ"%(self.freq), fontsize=16)
            ax_ph[1,1].set_title("%sU"%(self.freq), fontsize=16)
            fig_ph.savefig(self.path_dir + "Figures/pos_ac/%s/%s_%s/%s/%s_%s.png"%(self.freq, self.station, self.date, self.unpol, self.freq, self.pos_sour))
            close_figure(fig_ph)

        fsize  = 14
        fs = 10
        lines  = ["-", "--", ":"]
        marks  = ["D", "o" , "s"]
        fig_pol, ax_pol = plt.subplots(5, 2, figsize=(fsize, fsize*11/16))
        cpwr_amp, cpwr_phs = ax_pol[0,0], ax_pol[0,1]   # cross power spectrum
        lpol_amp, lpol_phs = ax_pol[1,0], ax_pol[1,1]   # (leakage for "unpol" | cross-correlation  for "pol")
        ccpl_1  , ccpl_2   = ax_pol[2,0], ax_pol[2,1]   # (LL, RR  for "unpol" | Stokes I, Stokes V for "pol")
        tvfc_1  , tvfc_2   = ax_pol[3,0], ax_pol[3,1]   # VFC temperature      at LL and RR
        tvfc_1_ , tvfc_2_  = ax_pol[4,0], ax_pol[4,1]   # mean VFC temperature at LL and RR

        label_cpwr       = r"$(Q+iU)~({\rm uncor.})$"
        label_lpol_unpol = r"$d_{\rm L} - d_{\rm R}$"
        label_lpol_pol   = r"$(Q+iU)/I$"
        label_ccpl_unpol_1, label_ccpl_unpol_2 = r"$v_{\rm LL}$", r"$v_{\rm RR}$"
        label_ccpl_pol_1  , label_ccpl_pol_2   = r"$I/I_{0}$"   , r"$V/I$"
        label_tvfc_1      , label_tvfc_2       = r"$T_{\rm VFC, LL}$"  , r"$T_{\rm VFC, RR}$"
        label_tvfc_1_     , label_tvfc_2_      = r"<$T_{\rm VFC, LL}>$", r"$<T_{\rm VFC, RR}>$"

        if self.lr_swap : LRpol = ["R", "L"]
        else            : LRpol = ["L", "R"]
        title1 = "%s"%(date_time) + " "*40 \
        +"%s %s%s %s%s %s %s %s %s %.2f %.2f %.2f"%(self.pos_sour, self.freq, LRpol[0], self.freq, LRpol[1], int(scannums[0]), nrep, nswitch, Tintg, az, el, self.gain)
        title2 = "%.3f (%.3f) | %.3f (%.3f)"%(self.getp_t1, self.getp_dt1, self.getp_t2, self.getp_dt2)
        if mode=="unpol":
            label_lpol = label_lpol_unpol
            label_ccpl_1, label_ccpl_2 = label_ccpl_unpol_1, label_ccpl_unpol_2
        else:
            s = " | %.3f (%.3f) | %.3f (%.3f) | %.3f (%.3f) | %.3f (%.3f) | %.3f (%.3f)"%(
                 self.getp_ti, self.getp_dti, self.getp_tp, self.getp_dtp, self.getp_tv, self.getp_dtv, self.getp_pm*100, self.getp_dpm*100, self.getp_pa, self.getp_dpa)
            title2 += s
            label_lpol = label_lpol_pol
            label_ccpl_1, label_ccpl_2 = label_ccpl_pol_1, label_ccpl_pol_2

        cpwr_amp.set_ylabel(label_cpwr, fontsize=fs)
        cpwr_phs.set_ylabel(label_cpwr, fontsize=fs)
        lpol_amp.set_ylabel(label_lpol, fontsize=fs)
        lpol_phs.set_ylabel(label_lpol, fontsize=fs)
        ccpl_1.set_ylabel(label_ccpl_1, fontsize=fs)
        ccpl_2.set_ylabel(label_ccpl_2, fontsize=fs)
        tvfc_1.set_ylabel(label_tvfc_1, fontsize=fs)
        tvfc_2.set_ylabel(label_tvfc_2, fontsize=fs)
        tvfc_1_.set_ylabel(label_tvfc_1_, fontsize=fs)
        tvfc_2_.set_ylabel(label_tvfc_2_, fontsize=fs)

        index_list = np.arange(ppd.c.shape[0])
        n=0

        for nidx, idx in enumerate(index_list):
            o = n // 10
            ls, mk = lines[o], marks[o]
            n+=1

            cpow_dat_amp = np.abs(ppd.c[idx])
            cpow_dat_phs = np.angle(ppd.c[idx], deg=True)
            if mode=="unpol":
                lpol_dat_amp = np.abs(ppd.leak[idx])
                lpol_dat_phs = np.angle(ppd.leak[idx], deg=True)
                crspol_1 = np.abs(ppd.d[idx,0])
                crspol_2 = np.abs(ppd.d[idx,1])
            else:
                # lpol_dat_amp = np.abs(ppd.cc[idx])/si_mean
                # lpol_dat_phs = np.angle(ppd.cc[idx], deg=True)
                lpol_dat_amp = np.abs(self.cal_cc[idx])/si_mean
                lpol_dat_phs = np.angle(self.cal_cc[idx], deg=True)
                crspol_1 = ppd.si[idx]
                crspol_2 = ppd.sv[idx]

            Tvfc1 = ppd.Tvfc1[idx]
            Tvfc2 = ppd.Tvfc2[idx]
            Tx1 = np.arange(len(Tvfc1))
            Tx2 = np.arange(len(Tvfc2))
            mask_nan1 = ~np.isnan(Tvfc1)
            mask_nan2 = ~np.isnan(Tvfc2)
            Tvfc1 = Tvfc1[mask_nan1]
            Tvfc2 = Tvfc2[mask_nan2]
            Tx1 = Tx1[mask_nan1]
            Tx2 = Tx2[mask_nan2]
            cpwr_amp.plot(cpow_dat_amp, ls=ls)
            cpwr_phs.plot(cpow_dat_phs, ls=ls)
            lpol_amp.plot(lpol_dat_amp, ls=ls)
            lpol_phs.plot(lpol_dat_phs, ls=ls)
            ccpl_1.plot(crspol_1, ls=ls)
            ccpl_2.plot(crspol_2, ls=ls)
            tvfc_1.plot(Tx1, Tvfc1, ls=ls, marker="o", markersize=6)
            tvfc_2.plot(Tx2, Tvfc2, ls=ls, marker="o", markersize=6)
            tvfc_1_.plot([idx], [np.nanmean(Tvfc1)], marker=mk)
            tvfc_2_.plot([idx], [np.nanmean(Tvfc2)], marker=mk)

        fig_pol.suptitle("%s \n %s"%(title1, title2))

        freq = self.freq
        if freq == 21:
            freq = 22

        path_fig = self.path_dir + "Figures/pos/%s/%s_%s_%s/"%(freq, self.station, self.date, self.proj)
        mkdir(path_fig)
        path_fig += "%s/"%(self.unpol)
        mkdir(path_fig)

        fig_pol.savefig(path_fig + "%s_%s_%s_%s.png"%(self.mode, int(self.scannums[0]), self.pos_sour, freq))
        close_figure(fig_pol)

    def SavePolDat(self):
        t1, dt1 = self.getp_t1, self.getp_dt1
        t2, dt2 = self.getp_t2, self.getp_dt2
        el = self.pos_log_sour[np.logical_and(self.pos_log_sour.Source==self.pos_sour, self.pos_log_sour.ScanNum==self.scannum)].reset_index(drop=True)["El"][0]
        if self.polnum==1: eta = self.eta1
        if self.polnum==2: eta = self.eta2
        if self.mode == "unpol":
            t1, t2= ufloat(t1, dt1), ufloat(t2, dt2)
            vt1t2 = unp.nominal_values(t1*t2)
            if   np.logical_or(vt1t2<0, np.isnan(vt1t2)) : uti = ufloat(np.nan, np.nan)
            elif vt1t2>=0                             : uti = unp.sqrt(t1*t2)
            usi = 8*C.k_B.value/(pi*21**2*eta)*uti * (u.J/u.m**2).to(u.Jy)
            eta = np.round(unp.nominal_values(eta), 5)
            ti, dti = np.round(unp.nominal_values(uti), 5), np.round(unp.std_devs(uti), 5)
            si, dsi = np.round(unp.nominal_values(usi), 5), np.round(unp.std_devs(usi), 5)
            poldat = pd.DataFrame([self.pos_sour, self.mjd, el, ti ,dti, 0 ,0, 0, 0, 0, 0, 0, 0, si ,dsi, 0 ,0, eta],
                                  index=["Source", "MJD", "El", "Ti", "dTi", "Tp", "dTp", "PM", "dPM", "PA", "dPA", "PA_c", "dPA_c", "Si", "dSi", "Sp", "dSp", "eta"]).transpose()
            self.poldat = poldat
        else:
            ti, dti = self.getp_ti, self.getp_dti
            tp, dtp = self.getp_tp, self.getp_dtp
            tv, dtv = self.getp_tv, self.getp_dtv
            pm, dpm = self.getp_pm, self.getp_dpm
            pa, dpa = self.getp_pa, self.getp_dpa
            uti, utp = ufloat(ti, dti), ufloat(tp, dtp)
            usi = 8*C.k_B.value/(pi*21**2*eta)*uti * (u.J/u.m**2).to(u.Jy)
            usp = 8*C.k_B.value/(pi*21**2*eta)*utp * (u.J/u.m**2).to(u.Jy)
            eta = np.round(unp.nominal_values(eta), 5)
            ti, dti = np.round(unp.nominal_values(uti), 5), np.round(unp.std_devs(uti), 5)
            tp, dtp = np.round(unp.nominal_values(utp), 5), np.round(unp.std_devs(utp), 5)
            si, dsi = np.round(unp.nominal_values(usi), 5), np.round(unp.std_devs(usi), 5)
            sp, dsp = np.round(unp.nominal_values(usp), 5), np.round(unp.std_devs(usp), 5)
            pm, dpm = np.round(pm, 5), np.round(dpm, 5)

            if self.mode == "aref":
                self.pa_aref = ufloat(pa, dpa)
                pa = np.round(pa, 5)
                dpa = np.round(dpa, 5)

                mjd_data = Ati(self.date, format="iso").mjd
                mjd_crit = Ati("2019-03-31", format="iso").mjd

                mask_pagan = "TRIPEE" in self.file.upper()
                mask_date  = mjd_data < mjd_crit

                if self.pos_sour == "3C286":
                    if 20 < int(self.freq) < 30:
                        pa_c, dpa_c  = 32.323, 0
                    elif 40 < int(self.freq) < 50:
                        pa_c, dpa_c  = 33.473, 0
                    elif 80 < int(self.freq) < 100:
                        pa_c, dpa_c  = 36.876, 0
                    elif 120 < int(self.freq) < 150:
                        pa_c, dpa_c  = 38.618, 0
                elif self.pos_sour in ["CRAB", "CRAB1", "CRAB2", "CRABP"]:
                    if mask_pagan:
                        if self.pos_sour == "CRABP":
                            pa_c, dpa_c  = 152, 0
                        else:
                            if 20 < int(self.freq) < 30:
                                pa_c, dpa_c  = 151.19, 0
                            elif 40 < int(self.freq) < 50:
                                pa_c, dpa_c  = 148.90, 0
                            elif 80 < int(self.freq) < 100:
                                pa_c, dpa_c  = 146.57, 0
                            elif 120 < int(self.freq) < 150:
                                pa_c, dpa_c  = 145.69, 0
                    else:
                        if mask_date:
                            if 20 < int(self.freq) < 30:
                                pa_c, dpa_c  = 151.19, 0
                            elif 40 < int(self.freq) < 50:
                                pa_c, dpa_c  = 148.90, 0
                            elif 80 < int(self.freq) < 100:
                                pa_c, dpa_c  = 146.57, 0
                            elif 120 < int(self.freq) < 150:
                                pa_c, dpa_c  = 145.69, 0
                        else:
                            if self.pos_sour == "CRAB":
                                pa_c, dpa_c  = 152, 0
                            else:
                                if 20 < int(self.freq) < 30:
                                    pa_c, dpa_c  = 151.19, 0
                                elif 40 < int(self.freq) < 50:
                                    pa_c, dpa_c  = 148.90, 0
                                elif 80 < int(self.freq) < 100:
                                    pa_c, dpa_c  = 146.57, 0
                                elif 120 < int(self.freq) < 150:
                                    pa_c, dpa_c  = 145.69, 0
                else:
                    pa_c = np.nan
                    dpa_c = np.nan
                self.pa_corr = pa_c
            else:
                pa_source = ufloat(pa, dpa)
                pa_source = pa_source - self.pa_aref + self.pa_corr
                pa  , dpa   = np.round(pa, 5), np.round(dpa, 5)
                pa_c, dpa_c = np.round(unp.nominal_values(pa_source), 5), np.round(unp.std_devs(pa_source), 5)
                if pa_c > 180:
                    pa_c -= 180
                if pa_c < 0:
                    pa_c += 180

            sourdat = pd.DataFrame([self.pos_sour, self.mjd, el, ti, dti, tp, dtp, pm*100, dpm*100, pa, dpa, pa_c, dpa_c, si, dsi, sp, dsp, eta],
                                   index=["Source", "MJD", "El", "Ti", "dTi", "Tp", "dTp", "PM", "dPM", "PA", "dPA", "PA_c", "dPA_c", "Si", "dSi", "Sp", "dSp", "eta"]).transpose()
            self.poldat = pd.concat([self.poldat, sourdat], axis=0)

    def SavePSLog(self):
        pos_log_all = self.pos_log_all
        pos_log_avg = self.pos_log_sour
        date  = self.date
        freqs = np.sort(np.unique(pos_log_all["Freq"]).astype(int))
        LRs   = ["L" , "R" ]
        self.freqs = freqs

        fsize = 16
        fig_pslog, ax_pslog = plt.subplots(3, 4, figsize=(fsize, fsize*9/16), sharex=True)
        tsys_1_l, tsys_1_r, tsys_2_l, tsys_2_r = ax_pslog[0,0], ax_pslog[0,1], ax_pslog[0,2], ax_pslog[0,3]
        tau_1_l , tau_1_r , tau_2_l , tau_2_r  = ax_pslog[1,0], ax_pslog[1,1], ax_pslog[1,2], ax_pslog[1,3]
        el_1_l  , el_1_r  , el_2_l  , el_2_r   = ax_pslog[2,0], ax_pslog[2,1], ax_pslog[2,2], ax_pslog[2,3]

        tsys_1_l.set_title("%sL"%(freqs[0])) ; tsys_1_r.set_title("%sR"%(freqs[0]))
        tsys_1_l.set_ylabel(r"$T_{\rm sys}$"+r"$\rm \,(K)$"    , fontsize=15)
        tau_1_l .set_ylabel(r"$\tau$"                          , fontsize=15)
        el_1_l  .set_ylabel(r"$\rm Elevation$"+r"$\rm \,(deg)$", fontsize=12)
        if self.npol==2:
            tsys_2_l.set_title("%sL"%(freqs[1])) ; tsys_2_r.set_title("%sR"%(freqs[1]))

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

        project = self.proj
        fig_pslog.text(0.5, 0.02, "(%s, %s) | UTC (hour)"%(project, date), ha="center", fontsize=15)

        colors = ["pink"     , "coral", "red" , "lime"   , "green" ,
                  "lightblue", "aqua" , "blue", "magenta", "purple"]
        log_all = copy.deepcopy(pos_log_all)
        log_avg = copy.deepcopy(pos_log_avg)
        log_all["Tsys"][np.isinf(log_all.Tsys.values.astype("f8"))] = np.nan
        log_all["Tau" ][np.isinf(log_all.Tau .values.astype("f8"))] = np.nan
        idx_chop = []

        log_all["Time"] = (log_all["MJD"] - Ati(date, format="iso").mjd) * u.day.to(u.hr)

        for nscan in range(log_avg.shape[0]):
            scannum = int(log_avg["ScanNum"][nscan])
            nswitch = log_avg["Nswitch"][nscan]
            nonoff = 2
            nrep = int(log_avg["Nrep"][nscan])
            ncount = self.npol * 4 * (nswitch * nonoff + 1)
            for nrep_ in range(nrep):
                idx_chop_ = [scannum + ncount * nrep_ + i for i in range(self.npol * 4)]
                idx_chop = idx_chop + idx_chop_

        mask_chop = np.array(list(map(lambda x: int(x) in idx_chop, log_all["ScanNum"])))
        log_all = log_all[~mask_chop].reset_index(drop=True)

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
                mask_pol = log_lr["Stokes"] == pol
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
                    mask_freq = log_pol["Freq"] == freq
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
        if np.isfinite(tsys_max):
            tsys_1_l.set_yscale("log")
            tsys_1_r.set_yscale("log")
            tsys_2_l.set_yscale("log")
            tsys_2_r.set_yscale("log")
            tsys_1_l.set_ylim(0.9*tsys_min, 1.1*tsys_max)
            tsys_1_r.set_ylim(0.9*tsys_min, 1.1*tsys_max)
            tsys_2_l.set_ylim(0.9*tsys_min, 1.1*tsys_max)
            tsys_2_r.set_ylim(0.9*tsys_min, 1.1*tsys_max)
            tsys_1_l.grid(True)
            tsys_1_r.grid(True)
            tsys_2_l.grid(True)
            tsys_2_r.grid(True)
        if np.isfinite(tau_max):
            tau_1_l.set_ylim(0.0, 1.1*tau_max)
            tau_1_r.set_ylim(0.0, 1.1*tau_max)
            tau_2_l.set_ylim(0.0, 1.1*tau_max)
            tau_2_r.set_ylim(0.0, 1.1*tau_max)
            tau_1_l.grid(True)
            tau_1_r.grid(True)
            tau_2_l.grid(True)
            tau_2_r.grid(True)
        el_1_l.grid(True)
        el_1_r.grid(True)
        el_2_l.grid(True)
        el_2_r.grid(True)
        tsys_1_l.legend(ncol=10, fontsize=10, bbox_to_anchor=(4.50, 1.50), fancybox=True)

        path_fig = self.path_dir + "Figures/pos/PS_Logs/"
        mkdir(path_fig)
        fig_pslog.savefig(path_fig + "/%s_%s_%s_PSLog.png"%(self.station, self.date, self.proj))
        close_figure(fig_pslog)


def correct_time(data, nchan):
    time_err = data[data.Time > 3600*20].reset_index()
    daysec   = 3600*24
    for nrow in range(time_err.shape[0]):
        time = time_err["MJD"][nrow]-int(time_err["MJD"][nrow]) * u.day.to(u.s)
        if np.logical_or(time<60, time>daysec-60):
           data["MJD"][time_err["index"][nrow]] -= 1
           data["Date"] = Ati(data["MJD"], format="mjd").iso
           data["Year"] = Ati(data["MJD"], format="mjd").byear
           mjd0 = int(np.min(data["MJD"]))
           data["Time"]     =np.round((data["MJD"]-mjd0) * u.day.to(u.s),0).astype(int)
           data["Time"][nchan:] =np.array(data["Time"][nchan:]) -np.array(data["Time"][0:-nchan])
           data["Time"][0:nchan]=np.array(data["Time"][0:nchan])-data["Time"][0]
    return data

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
