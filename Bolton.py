
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
    sat_mixing_ratio= esp* (sat_vapor_pressure(T)/(p-sat_vapor_pressure(T)))
    return(sat_mixing_ratio)


# In[5]:


def mixing_ratio_line(p,w_s):
    e_s=(w_s*p)/(esp + w_s)
    mixing_ratio_line=sat_vapor_temperature(e_s)
    return(mixing_ratio_line)


# In[6]:


def RH(T,p,w):
    RH = (w/sat_mixing_ratio(T, p)) * 100
    return(RH)


# In[7]:


def T_LCL(T, p, w):
    T_LCL = (1/ (1/((T + C_to_K)-55)-(np.log(RH(T,p,w)/100)/2840)))+ 55
    return(T_LCL)


# In[8]:


def theta_dry(T, p, p_0=1000.0):
    theta_dry = (T + C_to_K)* (p_0/(p))**k_dry
    return(theta_dry)


# In[9]:


def pseudoeq_potential_T(T, p, w, p_0=1000.0):
    pseudoeq_potential_T = (T + C_to_K)*(p_0/p)**(0.2845*(1-(0.28*w))) * np.exp((3.373/T_LCL(T, p, w)-0.00254) *(w*10**3) *(1+(0.81*w)))
    return(pseudoeq_potential_T)


# In[10]:


def theta_ep_field(T, p, p_0=1000.0):
    w_s= sat_mixing_ratio(T,p)
    theta_ep_field= pseudoeq_potential_T(T + C_to_K, p, w_s)
    return(theta_ep_field)

