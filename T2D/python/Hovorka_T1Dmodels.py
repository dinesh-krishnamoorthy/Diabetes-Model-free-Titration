import numpy as np

def HovorkaModel(x,t,parm,u,D):
    """
    Hovorka model
    
    Inputs:   t: Time [min] 
              x: States
              parm: Parameters
              u: Manipulated variables [mU/min]
              D: Disturbance (meal) [g CHO/min]
           
    Outputs:  dx: Derivative wrt. time.

    Written by Dinesh Krishnamoorthy (email: dineshk@ntnu.no), Jun 2020
    based on the MATLAB code by Dimitri Boiroux (DTU)
    """
    
    dx = np.zeros((11,1))
    
    # states 
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    Q1 = x[3]   # Plasma glucose [mmol]
    Q2 = x[4]   # Glucose in peripheral tissues [mmol]
    I  = x[5]   # Plasma insulin [mU/L]
    D1 = x[6]   # Amount of glucose in gut compartment 1 [mmol]
    D2 = x[7]   # Amount of glucose in gut compartment 2 [mmol]
    S1 = x[8]   # Amount of insulin in Sc compartment 1 [mU]
    S2 = x[9]   # Amount of insulin in Sc compartment 2 [mU]
    IG = x[10]  # Measurable interstitial glucose [mmol/L]
    
    # model parameters
    EGP0 = parm.EGP0   # Endogenous glucose production [mmol/min]
    F01   = parm.F01   # Insulin independent glucose consumption [mmol/min]
    Ag    = parm.Ag    # Glucose bioavailability [-]
    k12   = parm.k12   # Transfer rate [1/min]
    ka1   = parm.ka1   # Deactivation rate [1/min]
    ka2   = parm.ka2   # Deactivation rate [1/min]
    ka3   = parm.ka3   # Deactivation rate [1/min]
    kb1   = parm.kb1   # Activation rate [(L/mU)/min]
    kb2   = parm.kb2   # Activation rate [(L/mU)/min]
    kb3   = parm.kb3   # Activation rate [(L/mU)/min]
    ke    = parm.ke    # Insulin elimination rate [1/min]
    Vi    = parm.Vi    # Insulin distribution volume [L]
    Vg    = parm.Vg    # Glucose distribution volume [L]
    tmaxI = parm.tmaxI # Insulin absorption time constant [min]
    tmaxG = parm.tmaxG # CHO absorption time constant [min]
    tauIG = 15
    
    MwG   = 180.1577e-3 # Mol. wt of glucose
    
    F01c = min(F01,F01*Q1/(Vg*4.5))
    Fr = max(0,0.003*(Q1/Vg-9)*Vg)
    Ug = D2/tmaxG
    Ui = S2/tmaxI
    G = Q1/Vg # Glucose [mmol/L]
    
    # Insulin action 
    dx[0] = -ka1*x1+kb1*I
    dx[1] = -ka2*x2+kb2*I
    dx[2] = -ka3*x3+kb3*I
    
    # Glucose compartment
    dx[3] = -F01c-Fr-x1*Q1+k12*Q2+Ug+EGP0*(1-x3)
    dx[4] = x1*Q1-(k12+x2)*Q2
    
    #Plasma insulin
    dx[5] = Ui/Vi-ke*I   #[mU/L]
    
    # Meal absorption
    dx[6] = Ag*D/MwG-D1/tmaxG
    dx[7] = (D1-D2)/tmaxG
    
    # Subcutaneous Insulin
    dx[8] = u-S1/tmaxI
    dx[9] = (S1-S2)/tmaxI
    
    # Interstitial glucose
    dx[10] = -1/tauIG*(IG-G)
    
    return dx.reshape(11,)
    


def HovorkaRandom(seed):
    """
    Function to generate the model parameters for a 
    random Hovorka virtual patient
           
    Outputs:  dx: Derivative wrt. time.

    Written by Dinesh Krishnamoorthy (email: dineshk@ntnu.no), Jun 2020
    based on the work by Dimitri Boiroux
    """
    class parameters():
        
        np.random.seed(seed)
            
        r1 = np.random.rand(2,1)
        r2 = np.random.randn(14,1)
    
        while sum(r2>1)!=0 or sum(r2<-1)!=0:
            r2 = np.random.randn(14,1)
    
        BW = 80+30*(r1[0]-0.5)  # Patient body weight [kg]
        EGP0 = (0.0161+1.6e-3*np.sqrt(6)*r2[0])*BW
        F01 = (0.0097+0.9e-3*np.sqrt(6)*r2[1])*BW
        Ag = 0.7+r1[1]*(1.2-0.7)
        k12 = 0.0649+np.sqrt(6)*0.0115*r2[2]
        ka1 = 0.0055+np.sqrt(6)*0.0023*r2[3]
        ka2 = 0.0683+np.sqrt(6)*0.0207*r2[4]
        ka3 = 0.0304+np.sqrt(6)*0.0096*r2[5]
        kb1 = ka1*(51.2+13.1*np.sqrt(6)*r2[6])*1e-4
        kb2 = ka2*(8.2+3.2*np.sqrt(6)*r2[7])*1e-4
        kb3 = ka3*(520+125*np.sqrt(6)*r2[8])*1e-4
        ke = 0.14+0.035*r2[9]
        Vi = (0.12+0.012*r2[10])*BW
        Vg = np.exp(np.log(0.15)+0.23*r2[11])*BW
        tmaxI = 1/(0.018+0.0045*r2[12])
        tmaxG = 1/(np.exp(-3.689+0.25*r2[13]))
    
        while(EGP0<0 or F01<0 or Ag<0 or k12<0 or ka1<0 or ka2<0 or ka3<0 or kb1<0 or kb2<0 or kb3<0 or ke<0 or Vi<0 or Vg<0 or tmaxI<0 or tmaxG<0):
            r1 = np.random.rand(2,1)
            r2 = np.random.randn(14,1)
    
            while sum(r2>1)!=0 or sum(r2<-1)!=0:
                r2 = np.random.randn(14,1)
                
            BW = 80+30*(r1[0]-0.5)  # Patient body weight [kg]
            EGP0 = (0.0161+1.6e-3*np.sqrt(6)*r2[0])*BW
            F01 = (0.0097+0.9e-3*np.sqrt(6)*r2[1])*BW
            Ag = 0.7+r1[1]*(1.2-0.7)
            k12 = 0.0649+np.sqrt(6)*0.0115*r2[2]
            ka1 = 0.0055+np.sqrt(6)*0.0023*r2[3]
            ka2 = 0.0683+np.sqrt(6)*0.0207*r2[4]
            ka3 = 0.0304+np.sqrt(6)*0.0096*r2[5]
            kb1 = ka1*(51.2+13.1*np.sqrt(6)*r2[6])*1e-4
            kb2 = ka2*(8.2+3.2*np.sqrt(6)*r2[7])*1e-4
            kb3 = ka3*(520+125*np.sqrt(6)*r2[8])*1e-4
            ke = 0.14+0.035*r2[9]
            Vi = (0.12+0.012*r2[10])*BW
            Vg = np.exp(np.log(0.15)+0.23*r2[11])*BW
            tmaxI = 1/(0.018+0.0045*r2[12])
            tmaxG = 1/(np.exp(-3.689+0.25*r2[13]))
    
    parm = parameters()
    return parm


def HovorkaModelWrap(x,parm,Gb):
    dx = np.zeros((11,1))
    
    # states 
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    u  = x[3]   # NB! Different from Hovorka model !
    Q2 = x[4]   
    I  = x[5]   
    D1 = x[6]   
    D2 = x[7]   
    S1 = x[8]   
    S2 = x[9]   
    IG = x[10]  
    
    # model parameters
    EGP0 = parm.EGP0   # Endogenous glucose production [mmol/min]
    F01   = parm.F01   # Insulin independent glucose consumption [mmol/min]
    Ag    = parm.Ag    # Glucose bioavailability [-]
    k12   = parm.k12   # Transfer rate [1/min]
    ka1   = parm.ka1   # Deactivation rate [1/min]
    ka2   = parm.ka2   # Deactivation rate [1/min]
    ka3   = parm.ka3   # Deactivation rate [1/min]
    kb1   = parm.kb1   # Activation rate [(L/mU)/min]
    kb2   = parm.kb2   # Activation rate [(L/mU)/min]
    kb3   = parm.kb3   # Activation rate [(L/mU)/min]
    ke    = parm.ke    # Insulin elimination rate [1/min]
    Vi    = parm.Vi    # Insulin distribution volume [L]
    Vg    = parm.Vg    # Glucose distribution volume [L]
    tmaxI = parm.tmaxI # Insulin absorption time constant [min]
    tmaxG = parm.tmaxG # CHO absorption time constant [min]
    tauIG = 15
    Q1    = Gb*parm.Vg  # NB! This is not a state anymore !
    
    MwG   = 180.1577e-3 # Mol. wt of glucose
    
    F01c = F01
    Fr = 0
    Ug = D2/tmaxG
    Ui = S2/tmaxI
    G = Q1/Vg # Glucose [mmol/L]
    
    # Insulin action 
    dx[0] = -ka1*x1+kb1*I
    dx[1] = -ka2*x2+kb2*I
    dx[2] = -ka3*x3+kb3*I
    
    # Glucose compartment
    dx[3] = -F01c-Fr-x1*Q1+k12*Q2+Ug+EGP0*(1-x3)
    dx[4] = x1*Q1-(k12+x2)*Q2
    
    #Plasma insulin
    dx[5] = Ui/Vi-ke*I   #[mU/L]
    
    # Meal absorption
    dx[6] = -D1/tmaxG
    dx[7] = (D1-D2)/tmaxG
    
    # Subcutaneous Insulin
    dx[8] = u-S1/tmaxI
    dx[9] = (S1-S2)/tmaxI
    
    # Interstitial glucose
    dx[10] = -1/tauIG*(IG-G)
    
    return dx.reshape(11,)

  
from scipy.integrate import odeint
from scipy.optimize import fsolve
import numpy.matlib

def simPatient(Bolus,Meal,patientNr):
    """
    Simulates the effect of meal and bolus insulin on a Hovorka virtual patient, 
    and returns the cumulative cost of glycemic variations.

    Written by: Dinesh Krishnamoorthy, June 2020
    
    Input Arguments:
        Bolus [mU]
        Meal [g CHO/ 5 min] (an input of 10 is equal to a 50g meal)
        Patient number
    
    Output:
        Cumulative Cost
        J = 0.5*(BG-6)**2 + 100*0.5*min(BG-4.2,0)**2

    """
    
    Ts = 5 # Sampling time
    nsimHours = 6
    nSim = int(nsimHours*60/Ts) # Simulation duration, in number of samples 
    
    parm = HovorkaRandom(patientNr)
    
    # compute steady-state
    x0 = np.zeros((11,1))
    Gb = 6  # Nominal Blood Glucose [mmol/L]
    xSS = fsolve(HovorkaModelWrap,x0,args=(parm,Gb))
    uSS = xSS[3] # Basal Insulin that maintains nominal BG 
    # initialize 4th state back to Q1 to use with HovorkaModel
    xSS[3] = Gb*parm.Vg
    
    CGMNoise = np.zeros(shape = (nSim+1,1))
    
    U = np.matlib.repmat(uSS, nSim, 1) # Basal insulin
    D = np.zeros(shape = (nSim,1))
    
    # meal and bolus taken 1 hour after the simulation start
    D[11] = Meal  # Meal [g CHO/min]
    U[11] = U[11] + Bolus
    
    x0 = xSS  # Start at steady state 

    J = 0
    y = []
    BG = []
    
    for i in range(0,nSim):
        t = np.linspace(0,Ts)+(i-1)*Ts
        
        y.append(x0[10] + CGMNoise[i])
        BG.append(x0[3]/parm.Vg)
        J = J + 0.5*(y[i]-6)**2 + 100*0.5*min(y[i]-4.2,0)**2
    
        x = odeint(HovorkaModel,x0,t,args=(parm,U[i],D[i]))
        x0 = x[-1]
        
    return J

#==========================================================================
# delete below 

def DosageFinder(PatientNr,TypicalMealSize,BolusHistory,nBudget,MaxBolus):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import sklearn

    from DiabetesModels import simPatient

    Meal = TypicalMealSize

    def obj_fun(x0,Meal = Meal,PatientNr = PatientNr):    
        f = []
        for i in range(0,len(x0)):
            Bolus = x0
            fi = simPatient(Bolus[i],Meal,PatientNr) 
            f.append(fi)
        return np.array(f)

    X_init = BolusHistory
    Y_init = obj_fun(X_init,Meal)

    Y_init[0] = 5000*Meal/10
    bounds = np.array([[0, MaxBolus]])
    noise = 0.2; # sigma 

    # ------------ Bayesian Optimization ------------- 
    import GPy
    import GPyOpt

    from GPyOpt.methods import BayesianOptimization

    kernel = GPy.kern.Matern52(input_dim=1, variance=100.0, lengthscale=10)
    bds = [{'name': 'X', 'type': 'continuous', 'domain': bounds.ravel()}]

    optimizer = BayesianOptimization(f=obj_fun, 
                                     domain=bds,
                                     model_type='GP',
                                     kernel=kernel,
                                     acquisition_type ='EI',
                                     acquisition_jitter = 0.0,
                                     X=X_init,
                                     Y=Y_init,
                                     noise_var = noise**2,
                                     exact_feval=False,
                                     normalize_Y=True,
                                     maximize=False)

    optimizer.run_optimization(max_iter=nBudget)
    
    X_init = optimizer.X
    Y_init = optimizer.Y
    
    # ------------ Gather data -------------
    Xb = X_init[0:-1]
    Y0 = Y_init[0:-1]

    X0 = np.empty([len(Xb),3])
    for i in range(0,len(Xb)):
        X0[i] = [*Xb[i],Meal,*Y0[i]]
    
    
    return X0
   
    
