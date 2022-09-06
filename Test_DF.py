# -*- coding: utf-8 -*-
#09/03/2022: created
#09/05/2022: added QC2

import os
cd=os.path.dirname(__file__)
import numpy as np
import xarray as xr
from datetime import datetime
from matplotlib import pyplot as plt
import matplotlib 
import dynamic_filter as DF
plt.close('all')
matplotlib.rcParams.update({'font.size': 12})

#%% Inputs
source='data/nwtc.lidar-halo_xrp_user_199.b0.20220421.063708.nc'

#dataset
i1=10000#initial time index
i2=10500#final time index
rmin=96#[m] blind zone of the lidar

#filter parameters
Dx=100#[m]
Dy=100#[m]
Dz=100#[m]
Dt=600#[s]
max_rws=30#[m/s]
min_p=0.25/100
min_N=10
min_ratio=0.2

#%% Initialization
C=xr.open_dataset(source)
r0=C['distance_overlapped'].values
sel_r=r0>rmin
r=r0[sel_r]

t=np.array(C['time'].values.tolist())[i1:i2]/10**9
ele=C['elevation'].values[i1:i2]
azi=C['azimuth'].values[i1:i2]
rws=np.transpose(np.array(C['wind_speed']))[sel_r,i1:i2]
snr=np.transpose(np.array(C['SNR']))[sel_r,i1:i2]

#graphics
t_plot=[datetime.utcfromtimestamp(tt) for tt in t]

#%% Main
rws_qc,rws_qc2,rws_norm,snr_norm,p_all,N_all,ratio_all,t_calc,steps=DF.dynamic_filter_v2(t, r, ele, azi, rws, snr, Dx, Dy, Dz, Dt, max_rws, min_p, min_N, min_ratio)

print('Elapsed time: ' +str(t_calc[-1]))

#%% Plots
fig=plt.figure(figsize=(16,8))
plt.subplot(1,2,1)
plt.scatter(rws.ravel(),snr.ravel(),1,p_all*100,cmap='hot',vmin=0,vmax=100,alpha=1)
plt.xlabel(r'$u_{LOS}$ [m/s]')
plt.ylabel('SNR [dB]')
plt.title('Raw')
plt.colorbar(label=r'$p$ [%]')
plt.grid()

plt.subplot(1,2,2)
plt.scatter(rws_norm.ravel(),snr_norm.ravel(),1,p_all*100,cmap='hot',vmin=0,vmax=100,alpha=1)
plt.xlabel(r'$u\prime_{LOS}$ [m/s]')
plt.ylabel('SNR\' [dB]')
plt.xlim([-5,5])
plt.ylim([-5,5])
plt.title('Normalized')
plt.colorbar(label=r'$p$ [%]')
plt.grid()

plt.figure(figsize=(18,8))
plt.subplot(1,3,1)
plt.semilogx(p_all.ravel(),rws.ravel(),'.k',alpha=0.1)
plt.semilogx([min_p,min_p],[-40,40],'--r')
plt.semilogx([np.min(p_all),np.max(p_all)],[max_rws,max_rws],'--r')
plt.semilogx([np.min(p_all),np.max(p_all)],[-max_rws,-max_rws],'--r')
plt.xlabel(r'$p$')
plt.ylabel(r'$u_{LOS}$ [m/s]')
plt.grid()
plt.subplot(1,3,2)
plt.semilogx(N_all.ravel(),rws_qc.ravel(),'.k',alpha=0.1)
plt.semilogx([min_N,min_N],[-40,40],'--r')
plt.xlabel(r'$N$')
plt.ylabel(r'$u_{LOS}$ [m/s]')
plt.grid()
plt.subplot(1,3,3)
plt.plot(ratio_all.ravel(),rws_qc.ravel(),'.k',alpha=0.1)
plt.plot([min_ratio,min_ratio],[-40,40],'--r')
plt.xlabel(r'$N_{good}/N$')
plt.ylabel(r'$u_{LOS}$ [m/s]')
plt.grid()

fig=plt.figure(figsize=(18,8))
plt.subplot(4,1,1)
plt.pcolor(t_plot,r,rws,cmap='coolwarm',vmin=-10,vmax=10,shading='nearest')
plt.ylabel(r'$r$ [m]')
plt.colorbar(label=r'$u_{LOS}$ [m/s]')
plt.subplot(4,1,2)
plt.pcolor(t_plot,r,snr,cmap='hot',vmin=-30,vmax=0,shading='nearest')
plt.colorbar(label='SNR [dB]')
plt.ylabel(r'$r$ [m]')
plt.subplot(4,1,3)
plt.pcolor(t_plot,r,rws_qc,cmap='coolwarm',vmin=-10,vmax=10,shading='nearest')
plt.ylabel(r'$r$ [m]')
plt.colorbar(label=r'$u_{LOS}$ [m/s]')
plt.xlabel('UTC time')
plt.subplot(4,1,4)
plt.pcolor(t_plot,r,rws_qc2,cmap='coolwarm',vmin=-10,vmax=10,shading='nearest')
plt.ylabel(r'$r$ [m]')
plt.colorbar(label=r'$u_{LOS}$ [m/s]')
plt.xlabel('UTC time')

plt.text(min(t_plot),200,r' $\Delta x=$'+str(Dx)+' m \n '+
                  r'$\Delta y=$'+str(Dy)+' m \n '+
                  r'$\Delta z=$'+str(Dz)+' m \n '+
                  r'$\Delta t=$'+str(Dt)+' s \n '+
                  r'$|u_{LOS}|<$'+str(max_rws)+' m/s \n '+
                  r'$p>$'+str(min_p*100)+'% \n '+
                  r'$N>$'+str(min_N)+'% \n '+
                  r'$N_{good}/N>$ = '+str(min_ratio*100)+'%')
plt.tight_layout()

plt.figure(figsize=(18,6))
plt.bar(steps,np.diff(t_calc))
plt.ylabel('Computatoinal time [s]')
plt.xticks(rotation=10)
plt.grid()