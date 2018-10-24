
import numpy as np
C_to_K= 273.15
c_p_dry= 1005
c_V_dry= 717
esp= 0.622
k_dry= 0.2854


# In[2]:


def sat_vapor_pressure(T):
    sat_vapor_pressure= 6.112*np.exp((17.67*T)/(243.5+T))
    return(sat_vapor_pressure)


# In[3]:


def sat_vapor_temperature(e_s):
    sat_vapor_temperature= ((243.5*np.log(e_s)-440.8)/(19.48-np.log(e_s)))
    return(sat_vapor_temperature)


# In[4]:


def sat_mixing_ratio(p,T):
    sat_p=sat_vapor_pressure(T)
    sat_mixing_ratio= esp* (sat_p/(p-sat_p))
    return(sat_mixing_ratio)


# In[5]:


def mixing_ratio_line(p,w_s):
    e_s=(w_s*p)/(esp + w_s)
    sat_T=sat_vapor_temperature(e_s)
    mixing_ratio_line=sat_T
    return(mixing_ratio_line)


# In[6]:


def RH(T,p,w):
    sat_mr=sat_mixing_ratio(T,p)
    RH = (w/sat_mr) * 100
    return(RH)


# In[7]:


def T_LCL(T, p, w):
    rel_hum=RH(T,p,w)
    T_LCL = (1/ (1/((T + C_to_K)-55)-(np.log(rel_hum/100)/2840)))+ 55
    return(T_LCL)


# In[8]:


def theta_dry(theta, p, p_0=1000.0):
    theta_dry = (theta)* (p/p_0)**k_dry
    return(theta_dry)


# In[9]:


def pseudoeq_potential_T(T, p, w, p_0=1000.0):
    LCL_T=T_LCL(T, p, w)
    pseudoeq_potential_T = (T + C_to_K)*(p_0/p)**(0.2845*(1-(0.28*w))) * np.exp((3.373/LCL_T-0.00254) *(w*10**3) *(1+(0.81*w)))
    return(pseudoeq_potential_T)


# In[10]:


def theta_ep_field(T, p, p_0=1000.0):
    w_s= sat_mixing_ratio(T,p)
    theta_ep=pseudoeq_potential_T(T+ C_to_K, p, w_s)
    theta_ep_field= theta_ep
    return(theta_ep_field)

