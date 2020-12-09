def MVP_T2D(x,t,p):
    """
    Defines the differential equations for  the 
    4 compartment Kandarian Model modified for T2D (Aradottir et al. 2019)

    Written by: Dinesh Krishnamoorthy, May 2020
    
    Arguments:
        x :  vector of the state variables:
                  x = [I_s,I_p,I_e,G]
        t :  time
        p :  vector of the parameters:
                  p = [u,SI,pEGP,B,tau1,p2,pGEZI]

        States:
        I_s - Subcutaneous Insulin I_sc [U/day]
        I_p - Plasma Insulin I_p [U/day]
        I_e - Insulin effect on glucose I_eff [U/day]
        G   - Glucose concentration in plasma [mmol/L]
        
        Input:
        u - exogenous insulin input [U/day]
        
        Disturbances:
        SI   - Insulin sensitivity [1/U]
        pEGP - rate of endogenous glucose production [mmol/L day]
        B    - Endogenous insulin production co-eff beta [U L/mmol day]
        
        Parameters:
        tau1 - time constant [day]
        p2   - delay in insulin action [1/day]
        pGEZI-rate of glucose elimination from plasma [1/day]
        
    """
    I_s, I_p, I_e, G = x
    u, SI, pEGP, B, tau1, p2, pGEZI = p
    
    dx1 = (u - I_s)/tau1
    dx2 = (I_s - I_p)/tau1
    dx3 = p2*(I_p + B*G) - p2*I_e
    dx4 = -(pGEZI + SI*I_e)*G + pEGP
    
    f = [dx1,dx2,dx3,dx4]
    
    return f
  