import re, csv, os, sys, time, gc, emcee, corner
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt
import scipy.special     as ss
from datetime            import datetime
from scipy.optimize      import minimize
from scipy.optimize      import curve_fit
from uncertainties       import ufloat
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body
from uncertainties       import unumpy    as unp
from astropy             import constants as C
from astropy             import units     as u
from astropy.time        import Time      as Ati
from matplotlib.ticker   import (MultipleLocator, AutoMinorLocator)
from imports import dropbox_path

abort = sys.exit

fabs=np.fabs ; nconc=np.concatenate
sin =np.sin  ; cos  =np.cos    ; tan  =np.tan
exp =np.exp  ; e    =np.e      ; pi   =np.pi
ln  =np.log  ; log  =np.log10  ; angle=np.angle
sqrt=np.sqrt ; sqr  =np.square
real=np.real ; imag =np.imag

gc_kys = pd.read_excel(f"{dropbox_path}/AGN/KVN_SD/MOGABA_pipe/KVN_GC.xlsx", sheet_name="KYS")
gc_kus = pd.read_excel(f"{dropbox_path}/AGN/KVN_SD/MOGABA_pipe/KVN_GC.xlsx", sheet_name="KUS")
gc_ktn = pd.read_excel(f"{dropbox_path}/AGN/KVN_SD/MOGABA_pipe/KVN_GC.xlsx", sheet_name="KTN")

string_month = {'JAN':'01', 'FEB':'02', 'MAR':'03', 'APR':'04', 'MAY':'05', 'JUN':'06',
                'JUL':'07', 'AUG':'08', 'SEP':'09', 'OCT':'10', 'NOV':'11', 'DEC':'12'}

def cte(covariance):
    return np.sqrt(np.diag(covariance))

def format_date(date_str):
    date_str = date_str
    patterns = [r'\d{4}/\d{1,2}/\d{1,2}',       # yyyy/mm/dd
                r'\d{1,2}/\d{1,2}/\d{4}',       # mm/dd/yyyy
                r'\d{4}-\d{1,2}-\d{1,2}',       # yyyy-mm-dd
                r'\d{1,2}-\d{1,2}-\d{4}',       # mm-dd-yyyy
                r'\d{1,2}[A-Za-z]{3}\d{2}',     # ddMonyy (e.g., 1feb20)
                r'\d{1,2}-[A-Za-z]{3}-\d{2,4}'  # dd-Mon-yyyy (e.g., 06-MAY-2020)
               ]
    for pattern in patterns:
        match = re.search(pattern, date_str)
        if match:
            date = match.group()
            if '/' in date:
                date = re.sub('/', '-', date)
                date = str(date)
            else:
                if '-' in date:
                    #date = datetime.strptime(date, '%d-%b-%Y').strftime('%Y-%m-%d')
                    date = '%s-%s-%s'%(date.split('-')[2], string_month[date.split('-')[1]], date.split('-')[0])
                else:
                    date = datetime.strptime(date, '%d%b%y').strftime('%Y-%m-%d')
            return date
    return None

def format_time(date_str, time_float):
    h = int(time_float//3600)
    m = int((time_float-h*3600)//60)
    s = int(round(float((time_float-h*3600)%60),0))
    if h<10: h='0%s'%(h)
    if m<10: m='0%s'%(m)
    if s<10: s='0%s'%(s)
    datetime = '%s %s:%s:%s'%(date_str, h, m, s)
    return datetime

def close_figure(fig):
    plt.close(fig)
    plt.close('all')
    plt.clf()
    fig.clear()
    gc.collect()

def mkdir(path):
    if not os.path.isdir(path) : os.system('mkdir %s'%(path))

def pow(value, power):
    return value**power

def check_dir(path):
    return os.path.isdir(path)

def check_file(file):
    return os.path.isfile(file)

def mkpipelog(path, file):
    if os.path.isfile(path+file):
        os.system('rm %s%s'%(path,file))
    openfile = open(path+file, 'w')
    openfile.close()

def writelog(path, file, text, mode):
    get_time = time.strftime("%Y-%m-%d %X", time.localtime(time.time()))
    openlog  = open(path+file, mode=mode)
    openlog.write('(%s) %s\n'%(get_time, text))
    openlog.close()

def iwavg(value, sigma):
    """
    Inverse-variance Weighting Average
    """
    iwavg_sigma = (1/np.sum(1/sigma**2))**0.5
    iwavg_value = np.sum(value/sigma**2) * iwavg_sigma**2
    return iwavg_value, iwavg_sigma

def cal_gain_curve(station, year, freq, el):
    if   station=="KYS":
        gcurve = gc_kys
    elif station=="KUS":
        gcurve = gc_kus
    elif station=="KTN":
        gcurve = gc_ktn
    if int(freq)==25:
        freq=22
    if int(freq)==94:
        freq=86
    if int(freq)==141:
        freq=129
    gcurve_ = gcurve[gcurve.Season==int(year)]

    if freq == 21:
        a0 = gcurve_["A0_22"]
        a1 = gcurve_["A1_22"]
        a2 = gcurve_["A2_22"]
    else:
        a0 = gcurve_["A0_%s"%(freq)]
        a1 = gcurve_["A1_%s"%(freq)]
        a2 = gcurve_["A2_%s"%(freq)]
    gain = a0*el**2 + a1*el + a2
    return gain.values[0]
