
import numpy as np
C_to_K= 273.15
c_p_dry= 1005
c_V_dry= 717
esp= 0.622
k_dry= 0.2854


# In[2]:


def sat_vapor_pressure(T):
    e_s= 6.112*np.exp((17.67*T)/(243.5+T))
    return(e_s)


# In[3]:


def sat_vapor_temperature(e_s):
    sat_vap_temp= (243.5*np.log(e_s)-440.8)/(19.48-np.log(e_s))
    return(sat_vap_temp)


# In[4]:


def sat_mixing_ratio(p,T):
    e_s=sat_vapor_pressure(T)
    sat_mix_rat= esp* (e_s/(p-e_s))
    return(sat_mix_rat)


# In[5]:


def mixing_ratio_line(p,w_s):
    e_s=(w_s*p)/(esp + w_s)
    mix_rat_line=sat_vapor_temperature(e_s)
    return(mix_rat_line)


mix_rat=mixing_ratio_line(900, .012)
print(mix_rat)


def RH(T,p,w):
    R_H = (w/sat_mixing_ratio(p,T)) * 100
    return(R_H)


# In[7]:


def T_LCL(T, p, w):
    R_H=RH(T, p, w)
    term_a= 1/((T+C_to_K)-55)
    term_b= np.log(R_H/100)/2840
    Temp_LCL= (1/(term_a-term_b))+55
    return(Temp_LCL)


# In[8]:


def theta_dry(theta, p, p_0=1000.0):
    theta_d = (theta)* (p/p_0)**k_dry
    return(theta_d)


# In[9]:


def pseudoeq_potential_T(T, p, w, p_0=1000.0):
    Temp_LCL=T_LCL(T+C_to_K, p, w)
    part_1=(T+C_to_K)*(p_0/p)**(0.2854*(1-0.28*(w*10**(-3))))
    part_2= (3.376/Temp_LCL)-0.00254
    part_3= (w*10**3)*(1+0.81*(w*10**-3))
    theta_ep= part_1* np.exp(part_2*part_3)
    return (theta_ep)


# In[10]:


def theta_ep_field(T, p, p_0=1000.0):
    w= sat_mixing_ratio(T,p)
    pseudo_ep_field= pseudoeq_potential_T(T, p, w, p_0)
    return(pseudo_ep_field)

