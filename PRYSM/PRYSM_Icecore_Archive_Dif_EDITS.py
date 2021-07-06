def densification(Tavg,bdot,rhos,z):#,model='hljohnsen'):

    '''## Calculates steady state snow/firn depth density profiles using Herron-Langway type models.
        
        # [rho,zieq]=densification(Tavg,bdot,rhos,z,model)
        
        # Herron-Langway type models. (Arthern et al. 2010 formulation).
        
        # INPUT:
        # Tavg: 10m temperature in celcius ## CELCIUS!
        # bdot: accumulation rate in mwe/yr or (kg/m2/yr)
        # rhos: surface density in kg/m3
        # z: depth in true_metres 
        
        # model can be: {'HLJohnsen' 'HerronLangway' 'LiZwally' 'Helsen' 'NabarroHerring'}
        # default is herronlangway. (The other models are tuned for non-stationary modelling (Read Arthern et al.2010 before applying in steady state).
        
        # OUTPUT: 
        # rho: density (kg/m3) for all z-values.
        # zieq: ice equivalent depth for all z-values.
        # t: age for all z-values (only taking densification into account.)
        
        # Example usage:
        # z=0:300
        # [rho,zieq,t]=densitymodel(-31.5,177,340,z,'HerronLangway')
        # plot(z,rho)
        # Aslak Grinsted, University of Copenhagen 2010
        # Adapted Sylvia Dee, Brown University, 2017
    '''

    import numpy as np
    import scipy
    from scipy import integrate
    import matplotlib.pyplot as plt
    
    Tavg = 296
    bdot = 1.28
    rho_s = 300
    z = np.arange(0, 1286, 0.12482751338928937)+0.12482751338928937

    rhoi=920. 
    rhoc=550.
    #rhoc=822.
    rhow=1000. 
    rhos = 340.
    R=8.314
    #tmodel=model

    #Tavg=248.
    #bdot=0.1
    #Herron-Langway with Johnsen et al 2000 corrections.
    #Small corrections to HL model which are not in Arthern et al. 2010 

    c0=np.float64(0.85*11*(bdot/rhow)*np.exp(-10160./(R*Tavg)))
    c1=1.15*575*np.sqrt(bdot/rhow)*np.exp(-21400./(R*Tavg))

    k0=c0/bdot ##~g4
    k1=c1/bdot

    #critical depth at which rho=rhoc
    zc=(np.log(rhoc/(rhoi-rhoc))-np.log(rhos/(rhoi-rhos)))/(k0*rhoi) #g6

    ix=z<=zc #find the z's above and below zc
    upix=np.where(ix) #indices above zc
    dnix=np.where(~ix) #indices below zc

    q=np.zeros((z.shape)) #pre-allocate some space for q, rho
    rho=np.zeros((z.shape))

    #test to ensure that this will not blow up numerically if you have a very very long core.
    # manually set all super deep layers to solid ice (rhoi=920)
    NUM=k1*rhoi*(z-zc)+np.log(rhoc/(rhoi-rhoc))

    numerical= np.where(NUM<=100.0)
    blowup=np.where(NUM>100.0)

    q[dnix]=np.exp(k1*rhoi*(z[dnix]-zc)+np.log(rhoc/(rhoi-rhoc))), dtype = 'int64' #g7
    q[upix]=np.exp(k0*rhoi*z[upix]+np.log(rhos/(rhoi-rhos))) #g7

    rho[numerical]=q[numerical]*rhoi/(1+q[numerical]) #[g8]
    rho[blowup]=rhoi

    #only calculate this if you want zieq
    tc=(np.log(rhoi-rhos)-np.log(rhoi-rhoc))/c0 #age at rho=rhoc [g17]
    t=np.zeros((z.shape)) #pre allocate a vector for age as a function of z
    t[upix]=(np.log(rhoi-rhos)-np.log(rhoi-rho[upix]))/c0 # [g16] above zc
    t[dnix]=(np.log(rhoi-rhoc)-np.log(rhoi+0.0001-rho[dnix]))/c1+tc # [g16] below zc
    tdiff=np.diff(t)

    # make sure time keeps increasing even after we reach the critical depth.
    if np.any(tdiff==0.00):
        for i in range(len(t)):
            inflection=np.where(tdiff==0.0)
            lineardepth_change=t[inflection][0]

            for i in range(len(t)):
                if t[i]<lineardepth_change:
                    t[i]=t[i]
                elif t[i]>lineardepth_change:
                    t[i]=t[i-1]+1E-5

    zieq=t*bdot/rhoi #[g15]

    rho=np.array(rho)
    return rho,zieq,t


