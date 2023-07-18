import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy import ndimage
import matplotlib as matplotlib
import sys
import pdb

font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)
dpi_num=300
seasons=['DJF','MAM','JJA','SON']
#seasons=['JJA']
sim_names={'ERA5':'ERA5','ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (ERA5)','NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp':'GEM12_P3','ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon':'GEM2.5 (GEM12_P3-WCRI)'}
sim_colors={'ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon':'orange','ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon':'blue','NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp':'orange','ERA5':'blue'}
sim_lines={'ECan_2.5km_NAM11mP3_newP3_CLASS_DEEPoff_SHALon':'solid','ECan_2.5km_ERA5_newP3_CLASS_DEEPoff_SHALon':'solid','NAM-11m_ERA5_GEM5_CLASS_NV_NA_newP3-SCPF_SN8_20yrs_interp':'dashed','ERA5':'dashed'}
simulations=list(sim_names.keys())
filter_types=['gaussian-reflect','gaussian-nearest']
filter_type=filter_types[0]
sigmas=[10]
sigma=10
axes=['x','y']
tot_gps=[560] # total number of grid points to consider
max_ssds=[120] # maximum spin up
n_stds=[3]
dx=2.5

# Function filtering data
def apply_filter(raw,filter_type,sigma):
    if filter_type=='gaussian-nearest':
        filtered=ndimage.gaussian_filter(raw,sigma=sigma,mode='nearest')
    if filter_type=='gaussian-reflect':
        filtered=ndimage.gaussian_filter(raw,sigma=sigma,mode='reflect')
    return filtered

# Function identifying spatial spinup distance from filtered data
def get_spatial_su_distance(ratio,max_ssd,n_std):
    ssu_distance=[]
    mm=np.mean(ratio[max_ssd:])
    std=np.std(ratio[max_ssd:])
    min=np.min(ratio[max_ssd:])
    print('mm,std,min: ',mm,std,min)
    le=ratio<=(mm-n_std*std)
    le=ratio<=min
    #le=ratio<=(min)
    count=0
    for ii in range(len(le)):
        if np.sum(le[:ii])==ii:
            count=ii
    ssu_distance.append(count) # find the minimum (concave down)
    ssu_distance.append(mm) # find the minimum (concave down)
    ssu_distance.append(std) # find the minimum (concave down)
    return ssu_distance

# Loop of both sides of the longitudinal or meridional mean values of precipitation
dist_dic=np.zeros((int(len(simulations)/2),len(seasons),2*len(axes),len(n_stds),len(tot_gps),len(max_ssds)))

for nnstd,n_std in enumerate(n_stds):
    for ntotgp,tot_gp in enumerate(tot_gps):
        for nmaxssd,max_ssd in enumerate(max_ssds):
            border_name=[]
            for iaxe,axe in enumerate(axes):
                if axe=='x':
                    borders={0:'west',1:'east'}
                    label_t=r'$| \overline{p_{i}}^{j} |$'
                else:
                    borders={0:'south',1:'north'}
                    label_t=r'$| \overline{p_{j}}^{i} |$'
                for iborder,border in enumerate(borders):
                    print('Processing border: ',borders[iborder])
                    border_name.append(borders[border])
                    ibb=2*iaxe+iborder
                    for isea,sea in enumerate(seasons):
                        print('Processing season: ',sea)
                        sim_curve={}
                        rat_curve={}
                        count2p5=0
                        for isim,sim in enumerate(simulations):
                            print('Processing simulation: ',sim_names[sim])
                            file_name='npy/raw_signal_'+axe+'_'+sim+'_'+sea+'.npy'
                            darray=np.load(file_name) # read array
                            add_file='_n-std'+str(n_std)+'_tot-gp'+str(tot_gp)+'_max-ssd'+str(max_ssd)+'_'+filter_type+str(sigma)+'_'+axe+'_'+borders[iborder]+'_'+sea
                            if iborder==0: # for west and south keep as it is
                                raw=darray[:tot_gp]
                            else:
                                raw=darray[-tot_gp:][::-1] # for east and north mirror it
                            raw=raw/np.mean(raw) # normalise precipitation
                            x = np.arange(raw.shape[0]) # number of grid points from the border
                            filtered=apply_filter(raw,filter_type,sigma) #filtered data
                            sim_curve[sim]=filtered
                            #plt.plot(x,raw,label='raw '+sim_names[sim],color=sim_colors[sim],linewidth=1,linestyle=sim_lines[sim])
                            plt.plot(x,filtered,label=sim_names[sim],color=sim_colors[sim],linestyle=sim_lines[sim],linewidth=1.5)
                            if sim_names[sim][:6]=='GEM2.5':
                                #ratio=filtered/sim_curve[simulations[isim-1]]
                                ratio=100.*(filtered-sim_curve[simulations[isim-1]])/filtered
                                rat_curve[sim]=ratio
                                ssu_distance=get_spatial_su_distance(ratio,max_ssd,n_std) # spatial spinup distance
                                dist_dic[count2p5,isea,ibb,nnstd,ntotgp,nmaxssd]=ssu_distance[0]
                                count2p5=count2p5+1
                                        
                        #Adding text inside a rectangular box by using the keyword 'bbox'
                        plt.text(x[-1]-50,1.4,sea,fontsize = 18,bbox = dict(facecolor = 'grey', alpha = 0.5))
                        plt.xlabel('$\#$ of grid points from '+borders[iborder]+' border')
                        plt.ylabel(label_t)
                        plt.xlim([0,x[-1]+1])
                        plt.ylim([0.3,1.5])
                        if borders[iborder]=='west':
                            plt.legend(loc='lower right')
                        plt.savefig('png/raw-filtered_data'+add_file+'.png', bbox_inches="tight",dpi=dpi_num)
                        plt.close()

                        # Ratio figures
                        for sim in simulations:
                            if sim_names[sim][:6]=='GEM2.5':
                                plt.plot(x,rat_curve[sim],label=sim_names[sim],color=sim_colors[sim],linestyle='dotted',linewidth=1.5)
                                ssu_distance=get_spatial_su_distance(rat_curve[sim],max_ssd,n_std) # spatial spinup distance
                                plt.plot(ssu_distance[0],rat_curve[sim][ssu_distance[0]],label=sim_names[sim]+' SSUD',color=sim_colors[sim],marker='o',ms=10)
    
                                plt.plot(x[max_ssd:],(ssu_distance[1]+n_std*ssu_distance[2])*np.ones(x[max_ssd:].shape[0]),color=sim_colors[sim],linestyle='dashed',linewidth=.75)
                                plt.plot(x[max_ssd:],(ssu_distance[1]-n_std*ssu_distance[2])*np.ones(x[max_ssd:].shape[0]),color=sim_colors[sim],linestyle='dashed',linewidth=.75)

                        plt.text(x[-1]-50,1.32,sea, fontsize = 18,bbox = dict(facecolor = 'grey', alpha = 0.5))
                        plt.xlabel('$\#$ of grid points from '+borders[iborder]+' border')
                        plt.ylabel(label_t+' ratio')
                        plt.xlim([0,x[-1]+1])
                        #plt.ylim([0.5,1.4])
                        plt.ylim([-100,30])
                        if borders[iborder]=='west':
                            plt.legend(loc='lower right')
                        plt.savefig('png/ratio-filtered_data'+add_file+'.png', bbox_inches="tight",dpi=dpi_num)
                        plt.close()

x1=np.arange(4)+1.05
x2=np.arange(4)+.95
for nnstd,n_std in enumerate(n_stds):
    for ntotgp,tot_gp in enumerate(tot_gps):
        for nmaxssd,max_ssd in enumerate(max_ssds):
            pil=np.sum(dist_dic[0,:,:,nnstd,ntotgp,nmaxssd]>=dist_dic[1,:,:,nnstd,ntotgp,nmaxssd])
            tot=np.sum(dist_dic[0,:,:,nnstd,ntotgp,nmaxssd]>=0)
            print('Percentage of times that ssd is larger or equal without the intermediate pilot: ',100*pil/tot,' %')
            print('mean ssd with pilot ',simulations[0],': ',np.mean(dist_dic[0,:,nnstd,ntotgp,nmaxssd]),' # of points')
            print('mean ssd with pilot ',simulations[2],': ',np.mean(dist_dic[1,:,nnstd,ntotgp,nmaxssd]),' # of points')
            print('max ssd with pilot ',simulations[0],': ',np.max(dist_dic[0,:,nnstd,ntotgp,nmaxssd]),' # of points')
            print('max ssd with pilot ',simulations[2],': ',np.max(dist_dic[1,:,nnstd,ntotgp,nmaxssd]),' # of points \n')
            maxy=300
            for isea,sea in enumerate(seasons):
                add_file='_n-std'+str(n_std)+'_'+str(tot_gp)+'-'+str(max_ssd)+'_'+str(sigma)+'_'+sea
                plt.plot(x1,dx*dist_dic[0,isea,:,nnstd,ntotgp,nmaxssd].flatten(),color=sim_colors[simulations[1]],label=sim_names[simulations[1]], linestyle = 'None',marker='o',ms=10)
                plt.plot(x2,dx*dist_dic[1,isea,:,nnstd,ntotgp,nmaxssd].flatten(),color=sim_colors[simulations[3]],label=sim_names[simulations[3]], linestyle = 'None',marker='o',ms=10)
                if isea==0:
                    plt.legend()
                yy=plt.gca().get_ylim()
                xx=plt.gca().get_xlim()
                plt.text(4.2,maxy-30,sea,fontsize = 18,bbox = dict(facecolor = 'grey', alpha = 0.5))
                plt.xlabel('border')
                plt.ylabel('spatial spin-up distance (km)')
                plt.xlim([0,5])
                plt.ylim([0,maxy])
                plt.xticks(np.arange(4)+1,border_name)
                plt.savefig('png/ssd-mean_'+add_file+'.png', bbox_inches="tight",dpi=dpi_num)
                plt.close()

            maxy=200
            # Border mean figure
            add_file='borders_n-std'+str(n_std)+'_'+str(tot_gp)+'-'+str(max_ssd)+'_'+str(sigma)
            plt.plot(x1,dx*np.mean(dist_dic[0,:,:,nnstd,ntotgp,nmaxssd],axis=0).flatten(),color=sim_colors[simulations[1]],label=sim_names[simulations[1]], linestyle = 'None',marker='o',ms=10)
            plt.plot(x2,dx*np.mean(dist_dic[1,:,:,nnstd,ntotgp,nmaxssd],axis=0).flatten(),color=sim_colors[simulations[3]],label=sim_names[simulations[3]], linestyle = 'None',marker='o',ms=10)
            
            mm0=np.mean(dist_dic[0,:,:,nnstd,ntotgp,nmaxssd])
            mm1=np.mean(dist_dic[1,:,:,nnstd,ntotgp,nmaxssd])
            plt.ylabel('spatial spin-up distance (km)')
            plt.xlim([0,5])
            plt.ylim([0,maxy])
            plt.xticks(np.arange(4)+1,border_name)
            #plt.title(add_file)
            plt.savefig('png/ssd-mean_'+add_file+'.png', bbox_inches="tight",dpi=dpi_num)
            plt.close()
            
            # Season mean figure
            add_file='seasons_n-std'+str(n_std)+'_'+str(tot_gp)+'-'+str(max_ssd)+'_'+str(sigma)
            plt.plot(x1,dx*np.mean(dist_dic[0,:,:,nnstd,ntotgp,nmaxssd],axis=1).flatten(),color=sim_colors[simulations[1]],label=sim_names[simulations[1]], linestyle = 'None',marker='o',ms=10)
            plt.plot(x2,dx*np.mean(dist_dic[1,:,:,nnstd,ntotgp,nmaxssd],axis=1).flatten(),color=sim_colors[simulations[3]],label=sim_names[simulations[3]], linestyle = 'None',marker='o',ms=10)
            plt.ylabel('spatial spin-up distance (km)')
            plt.xlim([0,5])
            plt.ylim([0,maxy])
            plt.xticks(np.arange(4)+1,seasons)
            plt.savefig('png/ssd-mean_'+add_file+'.png', bbox_inches="tight",dpi=dpi_num)
            plt.close()
