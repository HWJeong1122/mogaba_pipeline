from imports import dropbox_path
from mogaba_pipe_imports     import *
from mogaba_pipe_cross       import *
from mogaba_pipe_polscripts  import *
from mogaba_pipe_pos         import *

import warnings
warnings.filterwarnings("ignore")

""" MOGABA pipeline script """
SaveCSFit  = True           # if 'True', save cross-scan Gaussian fitting plots
SaveCSLog  = True           # if 'True', save cross-scan log info
SavePSLog  = True           # if 'True', save position-switching log info
SaveACPlot = False          # if 'True', auto-correlation plots will be saved
Auto_Flag  = False          # if 'True', auto-flagging mode is applied in position-switching data
Run_CSFit  = True           # if 'True', cross-scan fit will be performed using the MCMC ; elsewhere, skip cs-fit
LR_Swap    = False          # if 'True', LR rx-pol switches to RL rx-pol
FlagChannel= False
perc_tvfc  = 16
                            # !!! Please note that LR-swapping is forecd at 129 GHz (@ line 353) !!!
antenna = f"tn".upper()
station = f"K{antenna}"     # should be in format of 'KYS' / 'KUS' / 'KPC' / 'KTN'
nw, nr  = 5*2, 2000         # the number of walkers & total step of MCMC in cs-profile fitting
Polnum  = 0                 # 0:all(1&2) / 1:1-only / 2:2-only
                            #   e.g., 1 and 2 mean K- and Q-band for KQ data

path_p     = "absolute/path/to/your/sdd/files/"             # sdd directory (python)
path_c     = "relative/path/to/your/sdd/files/"             # sdd directory (GILDAS/CLASS)
path_dir   = "path/to/your/working/directory/"              # working directory
path_cslog = "path/to/your/working/directory/data_cs/"
path_pslog = path_p

# path_p     = f"{dropbox_path}/ForShare/MOGABA_pipe/FITS/" # sdd directory (python)
# path_c     = f"~/Dropbox/ForShare/MOGABA_pipe/FITS/"      # sdd directory (GILDAS/CLASS)
# path_dir   = f"{dropbox_path}/ForShare/MOGABA_pipe/"      # working directory
# path_cslog = f"{dropbox_path}/ForShare/MOGABA_pipe/data_cs/"
# path_pslog = path_p

pipe_log = f"mogaba_pipelog_{antenna}_SE.log"

files = [
# files to run
# you may need to download the data from the server! (e.g., 'scp' command)
"MOGABA_K_11_KTN.sdd"
]

# some data contains different epochs of recording into one '.sdd' file.
# each number indicates the corresponding epoch by order.
version=[1, 2, 3]

flag_file =[
"files having some issues (optional)"
]

# unpolarized sources (mainly planets)
Unpols = ["JUPITER", "MARS", "VENUS", "SATURN"]

# angle reference source
Aref   = "CRAB"

flag_scan1 = [
# Pol_1
# [scannum, {channum:[subchans]}],

# epoch 1 (2011-11-21)
[775918, {0:[4,10], 1:[2,6,7], 2:[2], 3:[0,1,2,3,8]}],
[769400, {0:[3,4,5,6,10,11,13,14,15], 1:[1,2,3,9,10,11,12], 2:[5,6,7,8,9,10], 3:[6,12,15]}],
[770472, {0:[0,9,12,13,14,15], 1:[12]}],
[771562, {0:[0,1,2,3], 2:[0,1,2,3,4,9], 3:[5]}],
[772634, {0:[4,5,8,9], 1:[13,14,15], 2:[7,8], 3:[5,10]}],
[774828, {0:[2,3,4,12,13,14,15], 2:[0,1,2,3]}],
[777038, {0:[0,1,14,15], 1:[6], 3:[10,11,12,13,14,15]}],
[778128, {0:[5,6,7,8,9,10,11,12,13,14,15], 1:[3,4,5,6,7], 2:[0,1]}],
[779200, {0:[0,1,11,12,13,14,15], 1:[6,7,8,9], 3:[0,1]}],
[780322, {0:[3,4,8,11], 1:[0,1], 2:[0,5,6,7], 3:[6,7,8,9,10,11,12,13]}],
[781442, {0:[3], 1:[0,3], 3:[0,1,2,3,4]}],
[782532, {1:[1,11,13,14], 2:[0,1,2,3,4,5,6]}],
[783652, {1:[4,5,6,7,8,13], 2:[0,2,4,8,12]}],
[784742, {2:[3,4,5,6,7,8,9,10,11,12]}],
[785862, {0:[10,11,12,13,14,15], 1:[4,8], 3:[0,1]}],
[786936, {0:[0,1,2,3,4,5], 1:[0,1,2,3,4], 3:[13]}],
[789146, {0:[3,10], 1:[13], 2:[0,1,9,12,13,14,15]}],
[790234, {1:[0,1,2,3,4], 2:[0,1,2,3,4]}],
[791308, {0:[5,6,14,15], 2:[0,7,8,9,10,11,12], 1:[10,14,15]}],
[792380, {2:[1,3,4,5], 3:[9,10,11,12,13]}],
[793502, {0:[0,1,2,3,6], 2:[0,1,2,3,4,7]}],
[794622, {0:[7,8,10,12,13], 1:[4], 2:[5]}],
[775918, {1:[0,1,2,3,11,13,14,15]}],

# epoch 2 (2011-11-26)
[827322, {0:[0,1,6,7,8,9], 1:[0,1,2,3,4], 2:[0,1,2,3,4,7,8,10,12], 3:[14]}],
[803348, {0:[0,3], 1:[5], 2:[2,3,4,5,6,12], 3:[0,1,2]}],
[804422, {0:[0,6], 1:[0,1,2,3,4], 3:[0,1,2,3,11,13,14]}],
[805542, {0:[9,12], 1:[0,1,2], 2:[4,5,12], 3:[0,1,2,10,13]}],
[806664, {1:[0,4,5,6,7,8,11,12,13,14,15], 2:[6,10], 3:[4,11,12]}],
[807784, {0:[12,13,14,15], 1:[0,1,3,4,6], 2:[3,10,11], 3:[1,3,8]}],
[808906, {0:[2,5,8,10], 3:[8]}],
[810026, {0:[0,11,12,13,14], 1:[0,15], 2:[0], 3:[0]}],
[811148, {0:[3,4,5,13,14,15], 1:[0,1,2,10,12], 3:[0,1,15]}],
[812220, {0:[2,12,13], 3:[14,15]}],
[813294, {1:[0,1,2,3,4,5,6,13,14,15], 2:[0,1], 3:[0,1,2,3]}],
[814366, {0:[15], 3:[0,1,9,10,11,12,13]}],
[815456, {0:[0,1,4,5,6,11,12,13], 1:[0,1,5], 2:[0]}],
[816528, {0:[0,3], 3:[15]}],
[817618, {0:[0,13], 1:[0,3,12], 2:[0], 3:[0,1,2]}],
[818690, {0:[0], 1:[0,11,12], 3:[1,2]}],
[819780, {0:[0,1,6,9,10,15], 1:[0], 2:[0,1,2,3,4,5], 3:[0,1,2,3,4,5,6]}],
[820852, {1:[7,8], 3:[0]}],
[821926, {0:[12,13,14,15], 1:[4,11], 3:[11]}],
[823014, {1:[0], 3:[0]}],
[824088, {0:[0,1], 1:[0,1], 2:[0,1], 3:[0,1,3]}],
[825176, {0:[3,4], 1:[8,9,10], 3:[14,15]}],
[826250, {0:[0,1,4,8,12], 1:[0,5,6], 2:[0,1], 3:[0,1,2,9]}],
[828444, {0:[0,1,2,5,8], 1:[8,9,10], 2:[10]}],
[830654, {0:[7,8,9,10], 1:[8], 2:[8], 3:[13,14,15]}],
[831742, {0:[4], 1:[0,11,12], 2:[0], 3:[0,1,4,5,6]}],
[832816, {1:[0,3,4], 2:[0,3,4], 3:[0]}],
[833936, {0:[7], 2:[0,3,5], 3:[10,11,14,15]}],
[835058, {0:[1,9,11], 1:[10,11,12,13,14,15], 2:[3,4], 3:[7,8]}],
[836146, {0:[0], 1:[0], 2:[0], 3:[0]}],

# epoch 3 (2011-12-05)
[861292, {2:[14,15], 3:[0,1,2]}],
[856984, {0:[0], 1:[6,9,12], 2:[5,11,12], 3:[]}],
[858056, {0:[7,8], 1:[0,1,7,13,14,15], 2:[0,1,2], 3:[5,11,12,13,14,15]}],
[859146, {2:[2,3,4,5], 3:[8,9,10,11]}],
[860218, {0:[0], 1:[8,9], 3:[0,1,2]}],
[862412, {1:[8,9,10,11,12,13,14,15]}],
[865712, {0:[3], 1:[8], 3:[8,9,12]}],
[866784, {0:[10,11,12], 1:[15], 3:[12,13,14,15]}],
[869026, {0:[0,1,2,6]}],
[870116, {1:"all"}],
[871236, {0:[0,7,10,14], 1:[0], 2:[0,1,2], 3:[0,1]}],
[872326, {1:[0,1]}],
[873446, {1:[13,14,15], 2:[4,5,6,7]}],
[874520, {0:[8], 1:[7], 2:[8,9,13,14,15]}],
[875640, {1:[9,10,11], 1:[15], 2:[0,10,11,12,13,14,15], 3:[12,13,14,15]}],
[876730, {0:[0,10], 1:[0,1], 2:[0,1,2,3,4,6,14,15], 3:[0]}],
[878892, {1:[6,7,8], 2:[1,8], 3:[11,12,13,14,15]}],
[879964, {0:[12]}],
[881086, {0:[0,14,15], 1:[0], 2:[0], 3:[0]}],
[882206, {3:[9,10,11,12,13,14,15]}],
[863502, {1:[10,11], 3:[7,15]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
[111111, {0:[], 1:[], 2:[], 3:[]}],
]

flag_scan2 = [
# Pol_2
# [scannum, {channum:[subchans]}],
]

flag_scans = [flag_scan1, flag_scan2]
""""""""""""""""""""
""" Run Pipeline """
""""""""""""""""""""

mkpipelog(path_dir, pipe_log)
for Nfile, file in enumerate(files):
    if file in flag_file : continue
    """
    Cross-Scan
    """
    project = file.split(f"_{station}.sdd")[0]

    filter_cslogs = lambda x:np.logical_and(f"{station}_Log_CS_{project}_" in x, ".xlsx" in x)
    filter_pslogs = lambda x:np.logical_and(f"{project}_"                  in x, ".xlsx" in x)
    filter_cs_all = lambda x:"all" in x
    filter_ps_all = lambda x:"All" in x
    filter_cs_sour = lambda x:not "all" in x
    filter_ps_sour = lambda x:not "All" in x

    cslog_all_ = list(filter(filter_cslogs, os.listdir(path_cslog)))
    cslog_all  = list(filter(filter_cs_all , cslog_all_))
    cslog_sour = list(filter(filter_cs_sour, cslog_all_))
    cslog_all .sort()
    cslog_sour.sort()

    pslog_all_ = list(filter(filter_pslogs, os.listdir(path_pslog)))
    pslog_all  = list(filter(filter_ps_all , pslog_all_))
    pslog_sour = list(filter(filter_ps_sour, pslog_all_))
    pslog_all .sort()
    pslog_sour.sort()

    writelog(path_dir, pipe_log, f"***** Open SDD File : {file} (File Number:{int(Nfile+1)}) *****", "a")
    pipe_cs = CrossScan(
        path_p=path_p, path_c=path_c, path_dir=path_dir, file=file,
        station=station, pipe_log=pipe_log
    )
    pipe_cs.project = project
    pipe_cs.savecsfit = SaveCSFit
    pipe_cs.nwalker = nw
    pipe_cs.nstep = nr

    writelog(path_dir, pipe_log, "Initial CS-Setting is Done.", "a")
    pipe_cs.set_init()
    pipe_cs.saveplot = SaveCSLog
    Ncslog = len(cslog_sour)
    if Ncslog == 0:
        Ncslog = 1

    for n in range(Ncslog):
        if n+1 not in version : continue
        if Run_CSFit:
            try:
                writelog(path_dir, pipe_log, " ", "a")
                pipe_cs.cs_log_all  = pd.read_excel(path_cslog+cslog_all [n], index_col=[0]).dropna(axis=0).reset_index(drop=True)
                pipe_cs.cs_log_sour = pd.read_excel(path_cslog+cslog_sour[n], index_col=[0]).dropna(axis=0).reset_index(drop=True)
                writelog(path_dir, pipe_log, f"Load CS-Log Information from .xlsx Files (N={n+1}).", "a")
            except:
                writelog(path_dir, pipe_log, f"Make CS-Log Information from SDD Fil.", "a")
                pipe_cs.load_source_info(time_thresh=180)

            if np.logical_or(
                pipe_cs.cs_log_all.shape[0]<64,
                pipe_cs.cs_log_sour.shape[0]<1
            ):
                writelog(path_dir, pipe_log, f"Not Enough Data to proceed CS.", "a")
                writelog(path_dir, pipe_log, f"End Processes for the SDD File : {file} !!!", "a")
                writelog(path_dir, pipe_log, " ", "a")
                writelog(path_dir, pipe_log, " ", "a")
                writelog(path_dir, pipe_log, " ", "a")
                writelog(path_dir, pipe_log, "--------------------------------------------------", "a")
                continue

            pipe_cs.cs_log_all ["ScanNum"]=pipe_cs.cs_log_all ["ScanNum"].astype(int)
            pipe_cs.cs_log_all ["Freq"   ]=pipe_cs.cs_log_all ["Freq"   ].astype(int)
            pipe_cs.cs_log_sour["ScanNum"]=pipe_cs.cs_log_sour["ScanNum"].astype(int)
            pipe_cs.cs_log_sour["Nseq"   ]=pipe_cs.cs_log_sour["Nseq"   ].astype(int)
            pipe_cs.cs_log_sour["Nscan"  ]=pipe_cs.cs_log_sour["Nscan"  ].astype(int)
            pipe_cs.cs_log_sour["Scan1"  ]=pipe_cs.cs_log_sour["Scan1"  ].astype(int)
            pipe_cs.cs_log_sour["Scan2"  ]=pipe_cs.cs_log_sour["Scan2"  ].astype(int)

            date1 = np.array(pipe_cs.cs_log_sour["Date"])[0]
            date2 = np.array(pipe_cs.cs_log_sour["Date"])[-1]
            mjd_diff=Ati(date2, format="iso").mjd-Ati(date1, format="iso").mjd
            writelog(path_dir, pipe_log, f"CS-Observing Date : {date1} to {date2} ({mjd_diff:.2f} day)", "a")
            writelog(path_dir, pipe_log, f"Run Cross-Scan Fitting.", "a")
            pipe_cs.run_cs()

            writelog(path_dir, pipe_log, f"End Cross-Scan Fitting.", "a")
            writelog(path_dir, pipe_log, " ", "a")

            npol = pipe_cs.npol


        """
        Position-Switching
        """
        writelog(path_dir, pipe_log, "Start Position-Switching", "a")
        pipe_pos = PoSwitch(
            path_p=path_p, path_c=path_c, path_dir=path_dir,
            file=file, station=station, unpol_all=Unpols,
            aref_n=Aref, lr_swap=LR_Swap, pipe_log=pipe_log
        )
        try:
            pipe_pos.eta1 = unp.sqrt(pipe_cs.eta_l_1*pipe_cs.eta_r_1)
            pipe_pos.eta2 = unp.sqrt(pipe_cs.eta_l_2*pipe_cs.eta_r_2)
        except:
            try:
                pipe_pos.eta1 = unp.sqrt(pipe_cs.eta_l_1*pipe_cs.eta_r_1)
                pipe_pos.eta2 = np.nan
            except:
                pipe_pos.eta1 = np.nan
                pipe_pos.eta2 = np.nan
        writelog(path_dir, pipe_log, "Initial PS-Setting is Done.", "a")
        pipe_pos.set_init()
        pipe_pos.SaveACPlot = SaveACPlot
        pipe_pos.FlagChannel = FlagChannel
        pipe_pos.perc_tvfc = perc_tvfc

        pslog = file.replace(".sdd", ".xlsx")     # log file name
        pipe_pos.savepslog = SavePSLog
        pipe_pos.log = pslog
        try:
            writelog(path_dir, pipe_log, f"Load PS-Log Information from .xlsx Files (N={n+1}).", "a")
            pipe_pos.pos_log_all  = pd.read_excel(path_pslog+pslog_all [n]).dropna(axis=0) .reset_index(drop=True)
            pipe_pos.pos_log_sour = pd.read_excel(path_pslog+pslog_sour[n] , index_col=[0]).reset_index(drop=True)
        except:
            writelog(path_dir, pipe_log, f"Make PS-Log Information from SDD File.", "a")
            pipe_pos.load_source_info()

        if pipe_pos.pos_log_all.shape[0] < 2112:
            writelog(path_dir, pipe_log, f"Not Enough Data to proceed PS.", "a")
            writelog(path_dir, pipe_log, f"End Processes for the SDD File : {file} !!!", "a")
            writelog(path_dir, pipe_log, " ", "a")
            writelog(path_dir, pipe_log, " ", "a")
            writelog(path_dir, pipe_log, " ", "a")
            writelog(path_dir, pipe_log, "--------------------------------------------------", "a")
            continue


        pipe_pos.pos_log_all ["Freq"   ] = pipe_pos.pos_log_all ["Freq"   ].astype(int)
        pipe_pos.pos_log_all ["ScanNum"] = pipe_pos.pos_log_all ["ScanNum"].astype(int)
        pipe_pos.pos_log_sour["Nscan"  ] = pipe_pos.pos_log_sour["Nscan"  ].astype(int)
        pipe_pos.pos_log_sour["Nrep"   ] = pipe_pos.pos_log_sour["Nrep"   ].astype(int)
        pipe_pos.pos_log_sour["ScanNum"] = pipe_pos.pos_log_sour["ScanNum"].astype(int)
        pipe_pos.npol = len(np.unique(np.array(pipe_pos.pos_log_all ["Freq"])))

        date1, date2 = np.array(pipe_pos.pos_log_sour["Date"])[0], np.array(pipe_pos.pos_log_sour["Date"])[-1]
        mjd_diff = Ati(date2, format="iso").mjd-Ati(date1, format="iso").mjd
        writelog(path_dir, pipe_log, f"PS-Observing Date : {date1} to {date2} ({mjd_diff:.2f} day)", "a")
        writelog(path_dir, pipe_log, f"Get Calibrator Info.", "a")
        pipe_pos.get_info_calib()

        if not pipe_pos.errmsg is None:
            writelog(path_dir, pipe_log, pipe_pos.errmsg, "a")
            writelog(path_dir, pipe_log, f"End Processes for the SDD File : {file} !!!", "a")
            writelog(path_dir, pipe_log, " ", "a")
            writelog(path_dir, pipe_log, " ", "a")
            writelog(path_dir, pipe_log, " ", "a")
            writelog(path_dir, pipe_log, "--------------------------------------------------", "a")
            continue

        if SavePSLog:
            writelog(path_dir, pipe_log, "Save PS-Log (Tsys, tau, Elevation)", "a")
            pipe_pos.SavePSLog()

        if   Polnum == 0:
            pol_range=range(1,3)
        elif Polnum == 1:
            pol_range=range(1,2)
        elif Polnum == 2:
            pol_range=range(2,3)

        if pipe_pos.npol == 1:
            pol_range = range(1,2)

        sour_lst = pipe_pos.pos_log_sour["Source"].tolist()

        for Npn in pol_range:
            pipe_pos.polnum = Npn
            pipe_pos.freq = pipe_pos.freqs[Npn-1]

            if np.logical_and(
                not pipe_pos.unpols_n,
                "3C84" in sour_lst
            ):
                pipe_pos.unpols_n = ["3C84"]

            for Nunp, unpol in enumerate(pipe_pos.unpols_n):
                pipe_pos.unpol = unpol
                pipe_pos.bad_chans = flag_scans[Npn-1]

                if Auto_Flag:
                    pipe_pos.autoflag = Auto_Flag

                if str(pipe_pos.freq) == "129":
                    pipe_pos.lr_swap=True

                writelog(path_dir, pipe_log, f"Run Position-Switching (Unpol:{unpol} | Freq:{pipe_pos.freq} GHz).", "a")
                pipe_pos.out_pang = pd.DataFrame([])
                pipe_pos.run_pos()

    writelog(path_dir, pipe_log, f"End Processes for the SDD File : {file} !!!", "a")
    writelog(path_dir, pipe_log, " ", "a")
    writelog(path_dir, pipe_log, " ", "a")
    writelog(path_dir, pipe_log, " ", "a")
    writelog(path_dir, pipe_log, "--------------------------------------------------", "a")
