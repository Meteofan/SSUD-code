import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy import ndimage
import matplotlib as matplotlib
import sys
import pdb
import params as pr

font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)
dpi_num=300
# Read params file
seasons=pr.seasons
sim_names=pr.sim_names
sim_cpu=pr.sim_cpu
sim_diskp=pr.sim_diskp
pilot_names=pr.pilot_names
sim_pilots=pr.sim_pilots
sim_colors=pr.sim_colors
pilot_colors=pr.pilot_colors
sim_lines=pr.sim_lines
pilot_markers=pr.pilot_markers
#dx=pr.model_res['dx']
Nx=pr.model_res['Nx']
Ny=pr.model_res['Ny']

simulations=list(sim_names.keys())
filter_types=['gaussian-reflect','gaussian-nearest']
filter_type=filter_types[0]
sigmas=[5]
sigma=5
axes=['x','y']
tot_gps=[300,350,400,450,500] # total number of grid points to consider
max_ssds=[120] # maximum spin up
n_stds=[1.5,2.5,3.5]
#dx=2.5

# Loop of both sides of the longitudinal or meridional mean values of precipitation
dist_dic=np.zeros((int(len(simulations)),len(seasons),2*len(axes),len(n_stds),len(tot_gps),2))
percent_dic=np.zeros((len(n_stds),len(tot_gps),2))

for nnstd,n_std in enumerate(n_stds):
    for ntotgp,tot_gp in enumerate(tot_gps):
        max_ssds=np.asarray([int(0.3*tot_gp),int(0.25*tot_gp)])
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
                            file_name=pr.path_in+'raw_signal_'+axe+'_'+sim+'_'+sea+'_normalized.npy'
                            print('Reading: '+file_name)
                            darray=np.load(file_name) # read array
                            file_name=pr.path_in+'raw_signal_'+axe+'_'+sim_pilots[sim]+'_'+sea+'_normalized.npy'
                            print('Reading: '+file_name)
                            darray_p=np.load(file_name) # read array
                            add_file='_n-std'+str(n_std)+'_tot-gp'+str(tot_gp)+'_max-ssd'+str(max_ssd)+'_'+filter_type+str(sigma)+'_'+axe+'_'+borders[iborder]+'_'+sea
                            if iborder==0: # for west and south keep as it is
                                raw=darray[:tot_gp]
                                rawp=darray_p[:tot_gp]
                            else:
                                raw=darray[-tot_gp:][::-1] # for east and north mirror it
                                rawp=darray_p[-tot_gp:][::-1]
                            raw=raw/np.mean(raw) # normalise precipitation
                            rawp=rawp/np.mean(rawp) # normalise precipitation
                            x = np.arange(raw.shape[0]) # number of grid points from the border
                            filtered=pr.apply_filter(raw,filter_type,sigma) #filtered data
                            filteredp=pr.apply_filter(rawp,filter_type,sigma) #filtered data
                            sim_curve[sim]=filtered
                            sim_curve[sim_pilots[sim]]=filteredp
                            plt.plot(x,filtered,label=sim_names[sim],color=sim_colors[sim],linestyle=sim_lines[sim],linewidth=1.5)
                            plt.plot(x,filteredp,label=pilot_names[sim_pilots[sim]],color=pilot_colors[sim_pilots[sim]],linewidth=1.5)
                            RD=(filtered-filteredp)/(filtered+filteredp)
                            rat_curve[sim]=RD
                            ssu_distance=pr.get_spatial_su_distance(RD,max_ssd,n_std) # spatial spinup distance
                            dist_dic[isim,isea,ibb,nnstd,ntotgp,nmaxssd]=ssu_distance[0]
                                        
            pil=np.sum(dist_dic[1,:,:,nnstd,ntotgp,nmaxssd]>=dist_dic[4,:,:,nnstd,ntotgp,nmaxssd])
            tot=np.sum(dist_dic[1,:,:,nnstd,ntotgp,nmaxssd]>=0)
            percent_dic[nnstd,ntotgp,nmaxssd]=100*pil/tot
            
x1=np.arange(len(tot_gps))+1.025
x2=np.arange(len(tot_gps))+.975
dist_dic_m=np.mean(dist_dic[0,:],axis=(0,1))
dist_dic_s=np.std(dist_dic[0,:],axis=(0,1))
colors=['b','r','g']
markers=['o','x','s']
fig = plt.figure()
ax = plt.subplot(111)
for nmaxssd,max_ssd in enumerate(max_ssds):
    for nnstd,n_std in enumerate(n_stds):
        if nmaxssd==0:
            plt.errorbar(x1+nnstd*0.1-0.1,dist_dic_m[nnstd,:,nmaxssd],yerr=dist_dic_s[nnstd,:,nmaxssd],color=colors[nnstd],label=str(n_std)+r'$\sigma$; 33'+r'$\%$', linestyle = 'None',marker=markers[nmaxssd],ms=10)
        else:
            plt.errorbar(x2+nnstd*0.1-0.1,dist_dic_m[nnstd,:,nmaxssd],yerr=dist_dic_s[nnstd,:,nmaxssd],color=colors[nnstd],label=str(n_std)+r'$\sigma$; 25'+r'$\%$', linestyle = 'None',marker=markers[nmaxssd],ms=10)
#plt.legend()
#plt.text(4.2,maxy-30,sea,fontsize = 18,bbox = dict(facecolor = 'grey', alpha = 0.5))
plt.xlabel('total number of grid points')
plt.ylabel('spatial spin-up distance ($\#$ of gps)')
plt.xlim([0,len(tot_gps)+1])
plt.ylim([0,120])
plt.xticks(np.arange(len(tot_gps))+1,tot_gps)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(pr.path_out+'sensitivity_SSUD_ERA5_normalized.png', bbox_inches="tight",dpi=dpi_num)

legende = plt.legend(loc="lower right", ncol=2, framealpha=1,prop={'size': 10})
frame = legende.get_frame()
frame.set_color('white')
pr.export_legend(legende,pr.path_out+'sensitivity_SSUD_ERA5_legend_normalized.png')
plt.legend().set_visible(False)

plt.close("all")

x1=np.arange(len(tot_gps))+1.025
x2=np.arange(len(tot_gps))+.975
dist_dic1=np.mean(dist_dic[0,:],axis=(0,1))
dist_dic2=np.mean(dist_dic[-1,:],axis=(0,1))
dist_dic_m=dist_dic2/dist_dic1
colors=['b','r','g']
markers=['o','x','s']
fig = plt.figure()
ax = plt.subplot(111)
for nmaxssd,max_ssd in enumerate(max_ssds):
    for nnstd,n_std in enumerate(n_stds):
        if nmaxssd==0:
            plt.plot(x1+nnstd*0.1-0.1,dist_dic_m[nnstd,:,nmaxssd],color=colors[nnstd],label=str(n_std)+r'$\sigma$; 33'+r'$\%$', linestyle = 'None',marker=markers[nmaxssd],ms=10)
        else:
            plt.plot(x2+nnstd*0.1-0.1,dist_dic_m[nnstd,:,nmaxssd],color=colors[nnstd],label=str(n_std)+r'$\sigma$; 25'+r'$\%$', linestyle = 'None',marker=markers[nmaxssd],ms=10)
#plt.legend()
#plt.text(4.2,maxy-30,sea,fontsize = 18,bbox = dict(facecolor = 'grey', alpha = 0.5))
plt.xlabel('total number of grid points')
plt.ylabel('SSUD ratio')
plt.xlim([0,len(tot_gps)+1])
plt.ylim([0,.5])
plt.xticks(np.arange(len(tot_gps))+1,tot_gps)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(pr.path_out+'sensitivity_ratio_normalized.png', bbox_inches="tight",dpi=dpi_num)

legende = plt.legend(loc="lower right", ncol=2, framealpha=1,prop={'size': 10})
frame = legende.get_frame()
frame.set_color('white')
pr.export_legend(legende,pr.path_out+'sensitivity_ratio_legend_normalized.png')
plt.legend().set_visible(False)

plt.close("all")
