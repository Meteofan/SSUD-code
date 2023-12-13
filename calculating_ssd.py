import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
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

#print(str(enumerate(seasons)))
#input("un moment")

sim_names=pr.sim_names
sim_cpu=pr.sim_cpu
pilot_cpu=pr.pilot_cpu
sim_diskp=pr.sim_diskp
pilot_names=pr.pilot_names
sim_pilots=pr.sim_pilots
sim_colors=pr.sim_colors
pilot_colors=pr.pilot_colors
sim_lines=pr.sim_lines
pilot_markers=pr.pilot_markers
dx=pr.model_res['dx']
Nx=pr.model_res['Nx']
Ny=pr.model_res['Ny']
normalise='norm-no'
#units= '(' + r'mm $day^{-1}$' + ')'
units=""
simulations=list(sim_names.keys())
filter_types=['gaussian-reflect','gaussian-nearest','gaussian-wrap']
filter_type=filter_types[0]
sigma=5
axes=['x','y']
tot_gps=[400] # total number of grid points to consider
max_ssds=[133] # maximum spin up
n_stds=[2.5]
lwidth=1.8
                    
#print(simulations)

# Loop of both sides of the longitudinal or meridional mean values of precipitation
dist_dic=np.zeros((int(len(simulations)),len(seasons),2*len(axes),len(n_stds),len(tot_gps),len(max_ssds)))
cpu_dic=np.zeros((int(len(simulations))))

for nnstd,n_std in enumerate(n_stds):
    for ntotgp,tot_gp in enumerate(tot_gps):
        for nmaxssd,max_ssd in enumerate(max_ssds):
            border_name=[]
            for iaxe,axe in enumerate(axes):
                if axe=='x':
                    borders={0:'west',1:'east'}
                    if normalise=='norm-yes':
                        label_t=r'$| \overline{p_{i}}^{j} |$'
                    else:
                        label_t=r'$\overline{p^{j}}_{i}$ ' + units

                    label_t=r'$\langle \overline{p_{i}}^{j} \rangle$'
                else:
                    borders={0:'south',1:'north'}
                    if normalise=='norm-yes':
                        label_t=r'$| \overline{p_{j}}^{i} |$'
                    else:
                        label_t=r'$\overline{p^{i}}_{j}$ ' + units

                    label_t=r'$\langle \overline{p_{i}}^{j} \rangle$'

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
                            add_file='_n-std'+str(n_std)+'_tot-gp'+str(tot_gp)+'_max-ssd'+str(max_ssd)+'_'+filter_type+str(sigma)+'_'+axe+'_'+borders[iborder]+'_'+sea+'_'+normalise
                            if iborder==0: # for west and south keep as it is
                                raw=darray[:tot_gp]
                                rawp=darray_p[:tot_gp]
                                if normalise=='norm-yes':
                                    raw=darray[:tot_gp]/np.mean(darray[:])
                                    rawp=darray_p[:tot_gp]/np.mean(darray_p[:])
                            else:
                                raw=darray[-tot_gp:][::-1] # for east and north mirror it
                                rawp=darray_p[-tot_gp:][::-1]
                                if normalise=='norm-yes':
                                    raw=darray[-tot_gp:][::-1]/np.mean(darray[:])
                                    rawp=darray_p[-tot_gp:][::-1]/np.mean(darray_p[:])

                            x = np.arange(raw.shape[0]) # number of grid points from the border
                            filtered=pr.apply_filter(raw,filter_type,sigma) #filtered data
                            filteredp=pr.apply_filter(rawp,filter_type,sigma) #filtered data
                            sim_curve[sim]=filtered
                            sim_curve[sim_pilots[sim]]=filteredp
                            plt.plot(x,filtered,label=sim_names[sim],color=sim_colors[sim],linestyle=sim_lines[sim],linewidth=lwidth)
                            plt.plot(x,filteredp,label=pilot_names[sim_pilots[sim]],color=pilot_colors[sim_pilots[sim]],linewidth=lwidth)
                            RD=(filtered-filteredp)/(filtered+filteredp)
                            rat_curve[sim]=RD
                            ssu_distance=pr.get_spatial_su_distance(RD,max_ssd,n_std) # spatial spinup distance
                            dist_dic[isim,isea,ibb,nnstd,ntotgp,nmaxssd]=ssu_distance[0]
                            cpu_dic[isim]=sim_cpu[sim]+pilot_cpu[sim_pilots[sim]]

                        #Adding text inside a rectangular box by using the keyword 'bbox'
                        plt.xlabel('$\#$ of grid points from '+borders[iborder]+' border')
                        plt.ylabel(label_t)
                        plt.xlim([0,x[-1]+1])
                        if normalise=='norm-yes':
                            plt.ylim([0.3,1.7])
                        else:
                            ylim=[0,3.3]
                                
                        plt.ylim(ylim)
                        plt.text(x[-1]-50,ylim[1]-0.4,sea,fontsize = 18,bbox = dict(facecolor = 'grey', alpha = 0.5))

                        """
                        "legende = plt.legend(loc="lower right", framealpha=1,prop={'size': 10})
                        frame = legende.get_frame()
                        frame.set_color('white')
                        pr.export_legend(legende,'png/raw-filtered_data_legend.png')
                        plt.legend().set_visible(False)
                        #plt.legend(loc='lower right')
                        """
                        if borders[iborder]=='west' and len(simulations)==2:
                            plt.legend(loc='upper left')
                        plt.savefig(pr.path_out+'raw-filtered_data'+add_file+'_normalized.png', bbox_inches="tight",dpi=dpi_num)
                        plt.close("all")

                        # Ratio figures
                        for sim in simulations:
                            if sim_names[sim][:6]=='GEM2.5':
                                error = n_std*ssu_distance[2]
                                plt.plot(x,rat_curve[sim],label=sim_names[sim],color=sim_colors[sim],linestyle='dashed',linewidth=lwidth)
                                plt.fill_between(x, rat_curve[sim]-error, rat_curve[sim]+error,color=sim_colors[sim],alpha=.5)

                        for sim in simulations:
                            if sim_names[sim][:6]=='GEM2.5':
                                ssu_distance=pr.get_spatial_su_distance(rat_curve[sim],max_ssd,n_std) # spatial spinup distance
                                plt.plot(ssu_distance[0],rat_curve[sim][ssu_distance[0]],color=sim_colors[sim],marker='o',ms=10,alpha=1,markeredgecolor='black',markeredgewidth=2.)

                        plt.text(x[-1]-50,-.9,sea, fontsize = 18,bbox = dict(facecolor = 'grey', alpha = 0.5))
                        plt.xlabel('$\#$ of grid points from '+borders[iborder]+' border')
                        plt.ylabel('Relative Difference (RD)')
                        plt.xlim([0,x[-1]+1])
                        #plt.ylim([0.5,1.4])
                        plt.ylim([-1,.3])
                        #legende = plt.legend(loc="lower right", framealpha=1,prop={'size': 10})
                        #frame = legende.get_frame()
                        #frame.set_color('white')
                        #pr.export_legend(legende,'png/ratio-filtered_data_legend.png')
                        if borders[iborder]=='west' and len(simulations)==2:
                            plt.legend(loc='lower left')
                        #plt.legend().set_visible(False)
                        plt.savefig(pr.path_out+'ratio-filtered_data'+add_file+'_normalized.png', bbox_inches="tight",dpi=dpi_num)
                        plt.close("all")

x1=np.arange(4)+1
for nnstd,n_std in enumerate(n_stds):
    for ntotgp,tot_gp in enumerate(tot_gps):
        for nmaxssd,max_ssd in enumerate(max_ssds):
            pil=np.sum(dist_dic[0,:,:,nnstd,ntotgp,nmaxssd]>=dist_dic[-1,:,:,nnstd,ntotgp,nmaxssd])
            tot=np.sum(dist_dic[0,:,:,nnstd,ntotgp,nmaxssd]>=0)
            print('Percentage of times that ssd is larger or equal without the intermediate pilot: ',100*pil/tot,' %')
            print('mean ssd with pilot ',simulations[0],': ',np.mean(dist_dic[0,:,nnstd,ntotgp,nmaxssd]),' # of points')
            print('mean ssd with pilot ',simulations[-1],': ',np.mean(dist_dic[-1,:,nnstd,ntotgp,nmaxssd]),' # of points')
            print('max ssd with pilot ',simulations[0],': ',np.max(dist_dic[0,:,nnstd,ntotgp,nmaxssd]),' # of points')
            print('max ssd with pilot ',simulations[-1],': ',np.max(dist_dic[-1,:,nnstd,ntotgp,nmaxssd]),' # of points \n')
            maxy=140
            for isea,sea in enumerate(seasons):
                add_file='_n-std'+str(n_std)+'_tot-gp'+str(tot_gp)+'_max-ssd'+str   (max_ssd)+'_'+filter_type+str(sigma)+'_'+sea+'_'+normalise
                for isim,sim in enumerate(simulations):
                    plt.plot(x1-0.1+0.05*isim,dist_dic[isim,isea,:,nnstd,ntotgp,nmaxssd].flatten(),marker=pilot_markers[sim_pilots[sim]],color=sim_colors[sim],label=sim_names[sim], linestyle = 'None',ms=10)
                #if isea==0:
                #    plt.legend()
                yy=plt.gca().get_ylim()
                xx=plt.gca().get_xlim()
                #plt.text(4.2,maxy-30,sea,fontsize = 18,bbox = dict(facecolor = 'grey', alpha = 0.5))
                plt.text(4.2,maxy-12,sea,fontsize = 18,bbox = dict(facecolor = 'grey', alpha = 0.5))
                plt.xlabel('border')
                plt.ylabel('spatial spin-up distance ($\#$ of gps)')
                plt.xlim([0,5])
                plt.ylim([0,maxy])
                plt.xticks(np.arange(4)+1,border_name)
                plt.savefig(pr.path_out+'ssd-mean_'+add_file+'_normalized.png', bbox_inches="tight",dpi=dpi_num)

                legende = plt.legend(framealpha=1)
                frame = legende.get_frame()
                frame.set_color('white')
                pr.export_legend(legende,pr.path_out+'ssd-mean_legend_normalized.png')
                plt.legend().set_visible(False)
                plt.close("all")

            maxy=100
            # Border mean figure
            fig = plt.figure()
            ax = plt.subplot(111)
            add_file='borders_n-std'+str(n_std)+'_tot-gp'+str(tot_gp)+'_max-ssd'+str   (max_ssd)+'_'+filter_type+str(sigma)+'_'+normalise
            for isim,sim in enumerate(simulations):
                plt.plot(x1-0.1+0.05*isim,np.mean(dist_dic[isim,:,:,nnstd,ntotgp,nmaxssd],axis=0).flatten(),color=sim_colors[sim],marker=pilot_markers[sim_pilots[sim]],label=sim_names[sim], linestyle = 'None',ms=10)
            
            mm0=np.mean(dist_dic[0,:,:,nnstd,ntotgp,nmaxssd])
            mm1=np.mean(dist_dic[1,:,:,nnstd,ntotgp,nmaxssd])
            plt.ylabel('spatial spin-up distance ($\#$ of gps)')
            plt.xlim([0,5])
            plt.ylim([0,maxy])
            plt.xticks(np.arange(4)+1,border_name)
            #plt.title(add_file)
            # Shrink current axis by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            # Put a legend to the right of the current axis
            #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            plt.savefig(pr.path_out+'ssd-mean_'+add_file+'_normalized.png', bbox_inches="tight",dpi=dpi_num)

            legende = plt.legend(framealpha=1)
            frame = legende.get_frame()
            frame.set_color('white')
            pr.export_legend(legende,pr.path_out+'ssd-mean_legend_'+add_file+'_normalized.png')
            plt.legend().set_visible(False)
            plt.close("all")
            
            # Season mean figure
            fig = plt.figure()
            ax = plt.subplot(111)
            add_file='seasons_n-std'+str(n_std)+'_tot-gp'+str(tot_gp)+'_max-ssd'+str   (max_ssd)+'_'+filter_type+str(sigma)+'_'+normalise
            for isim,sim in enumerate(simulations):
                plt.plot(x1-0.1+0.05*isim,np.mean(dist_dic[isim,:,:,nnstd,ntotgp,nmaxssd],axis=1).flatten(),color=sim_colors[sim],marker=pilot_markers[sim_pilots[sim]],label=sim_names[sim], linestyle = 'None',ms=10)
            plt.ylabel('spatial spin-up distance ($\#$ of gps)')
            plt.xlim([0,5])
            plt.ylim([0,maxy])
            plt.xticks(np.arange(4)+1,seasons)
            # Shrink current axis by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            # Put a legend to the right of the current axis
            #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.savefig(pr.path_out+'ssd-mean_'+add_file+'_normalized.png', bbox_inches="tight",dpi=dpi_num)

            legende = plt.legend(loc="lower right", framealpha=1,prop={'size': 10})
            frame = legende.get_frame()
            frame.set_color('white')
            pr.export_legend(legende,pr.path_out+'ssd-mean_legend_'+add_file+'_normalized.png')
            plt.legend().set_visible(False)
            plt.close("all")

# Computational ressources figure
for ii in [0,1,2]:
    fig = plt.figure()
    ax = plt.subplot(111)
    add_file='cpu_'+str(n_std)+'_tot-gp'+str(tot_gp)+'_max-ssd'+str   (max_ssd)+'_'+filter_type+str(sigma)+'_'+normalise
    nxeff=Nx-(dist_dic[:,:,0,nnstd,ntotgp,nmaxssd]+dist_dic[:,:,1,nnstd,ntotgp,nmaxssd])
    nyeff=Ny-(dist_dic[:,:,2,nnstd,ntotgp,nmaxssd]+dist_dic[:,:,3,nnstd,ntotgp,nmaxssd])
    reduction=(Nx*Ny)/(nxeff*nyeff)
    for isim,sim in enumerate(simulations):
        if ii==0:
            add1='cpu_per_effgp'
            labely='norm. comput. cost per eff. gp'
            plt.plot(x1-0.1+0.05*isim,(1/cpu_dic[0])*cpu_dic[isim]*reduction[isim,:].flatten(),color=sim_colors[sim],marker=pilot_markers[sim_pilots[sim]],label=sim_names[sim], linestyle = 'None',ms=10)
            plt.ylim([1,1.6])
        if ii==1:
            add1='cpu_per_gp'
            labely='norm. comput. cost per gp'
            plt.plot(x1-0.1+0.05*isim,(1/cpu_dic[0])*cpu_dic[isim]*np.ones(len(seasons)),color=sim_colors[sim],marker=pilot_markers[sim_pilots[sim]], label=sim_names[sim], linestyle = 'None',ms=10)
            plt.ylim([1,1.3])

        if ii==2:
            add1='effective_domain'
            labely='norm. effective domain'
            plt.plot(x1-0.1+0.05*isim,1/reduction[isim,:].flatten(),color=sim_colors[sim],marker=pilot_markers[sim_pilots[sim]], label=sim_names[sim], linestyle = 'None',ms=10)
            plt.ylim([0.6,1.])

    plt.ylabel(labely)
    plt.xlim([0,5])
    plt.xticks(np.arange(4)+1,seasons)
    #plt.title('spatial spin up')
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(pr.path_out+add1+'_'+add_file+'_normalized.png', bbox_inches="tight",dpi=dpi_num)

    legende = plt.legend(framealpha=1)
    frame = legende.get_frame()
    frame.set_color('white')
    pr.export_legend(legende,pr.path_out+add1 + '_legend_normalized.png')
    plt.legend().set_visible(False)
    plt.close("all")
    
