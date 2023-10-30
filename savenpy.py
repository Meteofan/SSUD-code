import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy import ndimage
import matplotlib as matplotlib
import sys

import netCDF4
import copy

import seaborn as sns
import pandas as pd

liste_sim=['ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon','ECan_2.5km_NAM11mCS_noQC_noVS_newP3_CLASS_DEEPoff_SHALon','ECan_2.5km_NAM11mCS_newP3_CLASS_DEEPoff_SHALon','ECan_2.5km_NAM11mP3_noVS_noICE_noRAIN_newP3_CLASS_DEEPoff_SHALon','ECan_2.5km_NAM11mP3_noVS_noICE_newP3_CLASS_DEEPoff_SHALon','ECan_2.5km_NAM11mP3_noVS_newP3_CLASS_DEEPoff_SHALon','ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon']
legende_sim=['GEM2.5 (ERA5)','GEM2.5 (GEM12_SUN)','GEM2.5 (GEM12_SUN-W)','GEM2.5 (GEM12_P3-C)','GEM2.5 (GEM12_P3-CR)','GEM2.5 (GEM12_P3-CRI)','GEM2.5 (GEM12_P3-WCRI)']

liste_pil=['ERA5','NAM-11m_ERA5_GEM5_CLASS_NV_NA_CONSUN_SN8_20yrs_interp','NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp']
#legende_pil=['ERA5','GEM12_SUN (ERA5)','GEM12_P3 (ERA5)']

start_month=12
start_year=2015

end_month=11
end_year=2017

liste_saisons     = ['SON','DJF','MAM','JJA']
#liste_saisons     = ['SON']

for saison in liste_saisons:
    raw_signal_x_sim = []
    raw_signal_y_sim = []

    raw_signal_x_pil = []
    raw_signal_y_pil = []

    scd_sims         = []
    ssu_sims         = []

    for pil in liste_pil:
        if (pil == 'ERA5'):
            filename = '/pampa/roberge/REANALYSE/ERA5/Saisons/ERA5_0p0225_1330x1060_cas_4ans_GLHBQC_' + str(start_year) + str(start_month).zfill(2) + '-' + str(end_year) + str(end_month).zfill(2) + '/' + saison + '/' + 'var_pr_' + saison + '_' +  str(start_year) + str(start_month).zfill(2) + '-' + str(end_year) + str(end_month).zfill(2) + '_1h_timmean.nc4'
            varname  = 'pr'
        else:
            filename = '/pampa2/roberge/Output/GEM5/Cascades_CORDEX_NetCDF/Saisons/' + pil  + '_saisons_' + str(start_year) + str(start_month).zfill(2) + '-' + str(end_year) + str(end_month).zfill(2) + '_NetCDF/' + saison + '/' + 'var_PR_' + saison + '_' +  str(start_year) + str(start_month).zfill(2) + '-' + str(end_year) + str(end_month).zfill(2) + '_1h_timmean.nc4'
            varname  = 'pr'


        chemin_p = filename
        fichier_netcdf_p = netCDF4.Dataset(chemin_p)
        champ_p = fichier_netcdf_p.variables[varname][0,:,:]
        lat = fichier_netcdf_p.variables['lat'][:]
        lon = fichier_netcdf_p.variables['lon'][:]

        r = copy.deepcopy(champ_p)
        fichier_netcdf_p.close()

        # Afficher les dimensions de la matrice
        nbj, nbi = champ_p.shape

        print(champ_p.shape)
        print(nbj)
        print(nbi)

        # Séparer la matrice en 0.25 *np à 0.75 *np pour exclure la contamination de l'autre direction
        najd=int(0.25*nbj)
        najf=nbj-najd
        
        A = copy.deepcopy(najd)
        print("A = " + str(A))
    
        naid=int(0.25*nbi)
        
        B = copy.deepcopy(naid)
        print("B = " + str(B))

        naif=nbi-naid

        print("najd = " + str(najd))
        print("najf = " + str(najf))
        print("naid = " + str(naid))
        print("naif = " + str(naif))


        champ_slice_xfixe = copy.deepcopy(r[najd:najf,:])
        champ_slice_yfixe = copy.deepcopy(r[:,naid:naif])
    
        # Calcul de la moyenne de la variance sur la slice
        #raw_x = np.mean(champ_slice_xfixe,axis=0)
        #raw_y = np.mean(champ_slice_yfixe,axis=1)
        raw_signal_x_pil.append(np.mean(champ_slice_xfixe,axis=0)/np.mean(r))
        raw_signal_y_pil.append(np.mean(champ_slice_yfixe,axis=1)/np.mean(r))

        # Sauvegarde pour alejandro
        np.save("npy/raw_signal_x_" + pil + "_" + saison + "_normalized.npy",(np.mean(champ_slice_xfixe,axis=0)/np.mean(r)).filled())
        np.save("npy/raw_signal_y_" + pil + "_" + saison + "_normalized.npy",(np.mean(champ_slice_yfixe,axis=1)/np.mean(r)).filled())
        np.save("npy/raw_signal_x_" + pil + "_" + saison + "_notnormalized.npy",(np.mean(champ_slice_xfixe,axis=0)).filled())
        np.save("npy/raw_signal_y_" + pil + "_" + saison + "_notnormalized.npy",(np.mean(champ_slice_yfixe,axis=1)).filled())

    for sim in liste_sim:
        filename = '/pampa2/roberge/Output/GEM5/Cascades_CORDEX_NetCDF/Saisons/' + sim  + '_saisons_' + str(start_year) + str(start_month).zfill(2) + '-' + str(end_year) + str(end_month).zfill(2) + '_NetCDF/' + saison + '/' + 'var_PR_' + saison + '_' +  str(start_year) + str(start_month).zfill(2) + '-' + str(end_year) + str(end_month).zfill(2) + '_1h_timmean.nc4'
        varname  = 'pr'

        chemin_p = filename
        fichier_netcdf_p = netCDF4.Dataset(chemin_p)
        champ_p = fichier_netcdf_p.variables[varname][0,:,:]
        lat = fichier_netcdf_p.variables['lat'][:]
        lon = fichier_netcdf_p.variables['lon'][:]

        r = copy.deepcopy(champ_p)
        fichier_netcdf_p.close()

        # Afficher les dimensions de la matrice
        nbj, nbi = champ_p.shape

        print(champ_p.shape)
        print(nbj)
        print(nbi)

        # Séparer la matrice en 0.25 *np à 0.75 *np pour exclure la contamination de l'autre direction
        najd=int(0.25*nbj)
        najf=nbj-najd
        
        A = copy.deepcopy(najd)
        print("A = " + str(A))
    
        naid=int(0.25*nbi)
        
        B = copy.deepcopy(naid)
        print("B = " + str(B))

        naif=nbi-naid

        print("najd = " + str(najd))
        print("najf = " + str(najf))
        print("naid = " + str(naid))
        print("naif = " + str(naif))


        champ_slice_xfixe = copy.deepcopy(r[najd:najf,:])
        champ_slice_yfixe = copy.deepcopy(r[:,naid:naif])
    
        # Sauvegarde
        np.save("npy/raw_signal_x_" + sim + "_" + saison + "_normalized.npy",(np.mean(champ_slice_xfixe,axis=0)/np.mean(r)).filled())
        np.save("npy/raw_signal_y_" + sim + "_" + saison + "_normalized.npy",(np.mean(champ_slice_yfixe,axis=1)/np.mean(r)).filled())
        np.save("npy/raw_signal_x_" + sim + "_" + saison + "_notnormalized.npy",(np.mean(champ_slice_xfixe,axis=0)).filled())
        np.save("npy/raw_signal_y_" + sim + "_" + saison + "_notnormalized.npy",(np.mean(champ_slice_yfixe,axis=1)).filled())


            
                


