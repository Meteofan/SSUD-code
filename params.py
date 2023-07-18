import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy import ndimage
import sys
import pdb

def export_legend(legend, filename="legend.png"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

# Function filtering data
def apply_filter(raw,filter_type,sigma):
    filter=filter_type.split('-')[1]
    filtered=ndimage.gaussian_filter(raw,sigma=sigma,mode=filter)
    return filtered

# Function identifying spatial spinup distance from filtered data
def get_spatial_su_distance(array,max_ssd,n_std):
    ssu_distance=[]
    mm=np.mean(array[max_ssd:])
    std=np.std(array[max_ssd:])
    min=np.min(array[max_ssd:])
    print('mm,std,min: ',mm,std,min)
    le=array<=(mm-n_std*std)
    #le=ratio<=(min)
    count=0
    for ii in range(len(le)):
        if np.sum(le[:ii])==ii:
            count=ii
    ssu_distance.append(count) # find the minimum (concave down)
    ssu_distance.append(mm) # find the minimum (concave down)
    ssu_distance.append(std) # find the minimum (concave down)
    return ssu_distance

seasons=['DJF','MAM','JJA','SON']
#seasons=['JJA']
sim_names={'ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (ERA5)','ECan_2.5km_NAM11mCS_noQC_noVS_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (GEM12_SU)','ECan_2.5km_NAM11mCS_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (GEM12_SU-W)','ECan_2.5km_NAM11mP3_noVS_noICE_noRAIN_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (GEM12_P3-C)','ECan_2.5km_NAM11mP3_noVS_noICE_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (GEM12_P3-CR)', 'ECan_2.5km_NAM11mP3_noVS_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (GEM12_P3-CRI)','ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (GEM12_P3-WCRI)'}
#sim_names={'ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (ERA5)','ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (GEM12_P3-WCRI)'}

sim_cpu={'ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon':1.993,'ECan_2.5km_NAM11mCS_noQC_noVS_newP3_CLASS_DEEPoff_SHALon':2.078,'ECan_2.5km_NAM11mCS_newP3_CLASS_DEEPoff_SHALon':2.087,'ECan_2.5km_NAM11mP3_noVS_noICE_noRAIN_newP3_CLASS_DEEPoff_SHALon':2.090,'ECan_2.5km_NAM11mP3_noVS_noICE_newP3_CLASS_DEEPoff_SHALon':2.096, 'ECan_2.5km_NAM11mP3_noVS_newP3_CLASS_DEEPoff_SHALon':2.141,'ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon':2.145}

pilot_cpu={'ERA5':0,'NAM-11m_ERA5_GEM5_CLASS_NV_NA_CONSUN_SN8_20yrs_interp':0.153,'NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp':0.231}

sim_diskp={'ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon':3.6,'ECan_2.5km_NAM11mCS_noQC_noVS_newP3_CLASS_DEEPoff_SHALon':2.9,'ECan_2.5km_NAM11mCS_newP3_CLASS_DEEPoff_SHALon':4.6,'ECan_2.5km_NAM11mP3_noVS_noICE_noRAIN_newP3_CLASS_DEEPoff_SHALon':3.1,'ECan_2.5km_NAM11mP3_noVS_noICE_newP3_CLASS_DEEPoff_SHALon':3.4, 'ECan_2.5km_NAM11mP3_noVS_newP3_CLASS_DEEPoff_SHALon':4.0,'ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon':5.6}

pilot_names={'ERA5':'ERA5','NAM-11m_ERA5_GEM5_CLASS_NV_NA_CONSUN_SN8_20yrs_interp':'GEM12_SUN','ECan_2.5km_NAM11mCS_noQC_noVS_newP3_CLASS_DEEPoff_SHALon':'GEM12_SUN','NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp':'GEM12_P3'}

pilot_markers={'ERA5':'s','NAM-11m_ERA5_GEM5_CLASS_NV_NA_CONSUN_SN8_20yrs_interp':'v','NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp':'o'}
sim_pilots={'ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon':'ERA5','ECan_2.5km_NAM11mCS_noQC_noVS_newP3_CLASS_DEEPoff_SHALon':'NAM-11m_ERA5_GEM5_CLASS_NV_NA_CONSUN_SN8_20yrs_interp','ECan_2.5km_NAM11mCS_newP3_CLASS_DEEPoff_SHALon':'NAM-11m_ERA5_GEM5_CLASS_NV_NA_CONSUN_SN8_20yrs_interp', 'ECan_2.5km_NAM11mP3_noVS_newP3_CLASS_DEEPoff_SHALon':'NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp','ECan_2.5km_NAM11mP3_noVS_noICE_newP3_CLASS_DEEPoff_SHALon':'NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp','ECan_2.5km_NAM11mP3_noVS_noICE_noRAIN_newP3_CLASS_DEEPoff_SHALon':'NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp','ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon':'NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp'}

sim_colors={'ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon':'magenta','ECan_2.5km_NAM11mCS_noQC_noVS_newP3_CLASS_DEEPoff_SHALon':'r','ECan_2.5km_NAM11mCS_newP3_CLASS_DEEPoff_SHALon':'orange','ECan_2.5km_NAM11mP3_noVS_noICE_noRAIN_newP3_CLASS_DEEPoff_SHALon':'yellow','ECan_2.5km_NAM11mP3_noVS_noICE_newP3_CLASS_DEEPoff_SHALon':'green', 'ECan_2.5km_NAM11mP3_noVS_newP3_CLASS_DEEPoff_SHALon':'lightblue','ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon':'blue'}
pilot_colors={'ERA5':'k','NAM-11m_ERA5_GEM5_CLASS_NV_NA_CONSUN_SN8_20yrs_interp':'grey','NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp':'lightgrey'}

sim_lines={'ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon':'dashed','ECan_2.5km_NAM11mCS_noQC_noVS_newP3_CLASS_DEEPoff_SHALon':'dashed','ECan_2.5km_NAM11mCS_newP3_CLASS_DEEPoff_SHALon':'dashed','ECan_2.5km_NAM11mP3_noVS_noICE_noRAIN_newP3_CLASS_DEEPoff_SHALon':'dashed','ECan_2.5km_NAM11mP3_noVS_noICE_newP3_CLASS_DEEPoff_SHALon':'dashed', 'ECan_2.5km_NAM11mP3_noVS_newP3_CLASS_DEEPoff_SHALon':'dashed','ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon':'dashed'}

model_res={'dx':2.5,'Nx':1330,'Ny':1060}
