# -*- coding: utf-8 -*-
#dynamic lidar data filter based on Beck and Kuhn, Remote Sensing, 2017
import numpy as np

def lidar_xyz(r,ele,azi):
    Nr=len2(r)
    Nb=len2(ele)
    
    R=np.transpose(np.tile(r,(Nb,1)))
    A=np.tile(azi,(Nr,1))
    E=np.tile(ele,(Nr,1))
    
    X=R*cosd(E)*cosd(90-A)
    Y=R*cosd(E)*sind(90-A)
    Z=R*sind(E)
    
    return X,Y,Z

def cosd(x):
    return np.cos(x/180*np.pi)

def sind(x):
    return np.sin(x/180*np.pi)

def len2(x):
    if 'int' in str(type(x)) or 'float' in str(type(x)):
        return 1
    elif 'list' in str(type(x)) or 'array' in str(type(x))  or 'str' in str(type(x)) or 'series' in str(type(x)):
        return len(x)
    else:
        raise ValueError
        
def dynamic_filter_v2(t,r,ele,azi,rws,snr,Dx,Dy,Dz,Dt=600,max_rws=30,min_p=0.0025,min_N=10,min_ratio=0.2):
    #09/05/2022: created
    #09/06/2022: finalized
    from scipy import stats
    import time as tm
    
    #performance
    t0=tm.time()
    t_calc=0
    steps=['Stats N','Stats RWS','Index RWS','Stats SNR','Index SNR','Stats RWS/SNR','Index RWS/SNR','QC1','Stats N good','Index N','QC2']

    X,Y,Z=lidar_xyz(r, ele, azi)
    T,R=np.meshgrid(t,r)
    
    bin_x=np.arange(np.min(X),np.max(X)+Dx,Dx)
    bin_y=np.arange(np.min(Y),np.max(Y)+Dy,Dy)
    bin_z=np.arange(np.min(Z),np.max(Z)+Dz,Dz)
    bin_t=np.arange(np.min(T),np.max(T)+Dt,Dt)

    S=np.transpose(np.vstack((X.ravel(),Y.ravel(),Z.ravel(),T.ravel())))
    B=stats.binned_statistic_dd(S,[],'count',[bin_x,bin_y,bin_z,bin_t],expand_binnumbers=True)
    N=B[0]
    index=B[2]
    t_calc=np.append(t_calc,tm.time()-t0)
    
    B=stats.binned_statistic_dd(S,rws.ravel(),'mean',[bin_x,bin_y,bin_z,bin_t])
    rws_avg=B[0]
    t_calc=np.append(t_calc,tm.time()-t0)

    rws_avg_all=np.reshape(rws_avg[(index[0]-1,index[1]-1,index[2]-1,index[3]-1)],np.shape(rws))
    rws_norm=rws-rws_avg_all
    t_calc=np.append(t_calc,tm.time()-t0)

    snr_avg=stats.binned_statistic_dd(S,snr.ravel(),'mean',binned_statistic_result=B)[0]
    t_calc=np.append(t_calc,tm.time()-t0)

    snr_avg_all=np.reshape(snr_avg[(index[0]-1,index[1]-1,index[2]-1,index[3]-1)],np.shape(snr))
    snr_norm=snr-snr_avg_all
    t_calc=np.append(t_calc,tm.time()-t0)
                     
    #%probability
    sel=~np.isnan(rws_norm.ravel()+snr_norm.ravel())
    drws=3.49*np.nanstd(rws_norm.ravel()[sel])/(np.sum(sel))**(1/3)
    bin_rws=np.arange(np.min(rws_norm.ravel()[sel]),np.max(rws_norm.ravel()[sel])+drws,drws)
    dsnr=3.49*np.nanstd(snr_norm.ravel()[sel])/(np.sum(sel))**(1/3)
    bin_snr=np.arange(np.min(snr_norm.ravel()[sel]),np.max(snr_norm.ravel()[sel])+dsnr,dsnr)

    S2=np.transpose(np.vstack((rws_norm.ravel()[sel],snr_norm.ravel()[sel])))
    B2=stats.binned_statistic_dd(S2,[],'count',[bin_rws,bin_snr],expand_binnumbers=True)
    t_calc=np.append(t_calc,tm.time()-t0)
                     
    p=B2[0]/np.max(B2[0])                         
    index2=B2[2]
    p_sel=np.zeros(len(rws_norm.ravel()))
    p_sel[sel]=p[(index2[0]-1,index2[1]-1)]
    p_all=np.reshape(p_sel,np.shape(rws))   
    t_calc=np.append(t_calc,tm.time()-t0) 

    #qc1
    good=p_all>min_p
    rws_qc1=rws.copy()
    rws_qc1[np.abs(rws_qc1)>max_rws]=np.nan 
    rws_qc1[good==0]=np.nan   
    t_calc=np.append(t_calc,tm.time()-t0)

    #qc2
    N_good=stats.binned_statistic_dd(S,good.ravel(),'sum',binned_statistic_result=B)[0]
    t_calc=np.append(t_calc,tm.time()-t0)
                     
    N_all=np.reshape(N[(index[0]-1,index[1]-1,index[2]-1,index[3]-1)],np.shape(rws))
    N_good_all=np.reshape(N_good[(index[0]-1,index[1]-1,index[2]-1,index[3]-1)],np.shape(rws))
    ratio_all=N_good_all/(N_all+10**(-16))
    t_calc=np.append(t_calc,tm.time()-t0)

    rws_qc2=rws_qc1.copy()
    rws_qc2[N_all<min_N]=np.nan
    rws_qc2[ratio_all<min_ratio]=np.nan
    t_calc=np.append(t_calc,tm.time()-t0)
    
    return rws_qc1,rws_qc2,rws_norm,snr_norm,p_all,N_all,ratio_all,t_calc,steps