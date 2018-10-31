
import numpy as np
C_to_K= 273.15
c_p_dry= 1005
c_V_dry= 717
esp= 0.622
k_dry= 0.2854


# In[2]:


def sat_vapor_pressure(T):
    """PPK equ 10, temp is in celcius"""
    sat_vapor_pressure= 6.112 * np.exp((17.67*T)/(T+243.5))
    print("SVP=", sat_vapor_pressure)
    return(sat_vapor_pressure)


# In[3]:


def sat_vapor_temperature(e_s, T):
    """PPK equ 11, temp comes back in celcius"""
    SVP=sat_vapor_pressure(T)
    sat_vap_temp= (243.5*np.log(SVP)-440.8)/(19.48-np.log(SVP))
    return(sat_vap_temp)


# In[4]:


def sat_mixing_ratio(p,T):
    """p is in mb, temp is in celcius"""
    e_s=sat_vapor_pressure(T)
    #print("e_s=", e_s)
    sat_mixing_ratio= esp*(e_s/(p-e_s))
    print("sat_mix_rat=", sat_mixing_ratio)
    return(sat_mixing_ratio)



# In[5]:


def mixing_ratio_line(p,w_s, T):
    """equation to give the lines of mixing ratios. p in mb, w_s in kg/kg, temp is returned in celcius"""
    e_s=(w_s*p)/(esp + w_s)
    mix_rat_line=sat_vapor_temperature(e_s)
    return(mix_rat_line)



def RH(T,p,w):
    """p in mb, temp in celcius, w in kg/kg(unitless)"""
    #print("w=", w)
    #print("T=", T)
    #print("p=", p)
    sat_mix_rat=sat_mixing_ratio(p,T)
    R_H = (w/sat_mix_rat) * 100
    return(R_H)


# In[7]:


def T_LCL(T, p, w):
    """PPK equ 22, p in mb, temp in kelvin, w in kg/kg"""
    R_H=RH(T, p, w)
    term_a= 1/((T+C_to_K)-55)
    term_b= (np.log(R_H/100))/2840
    print(R_H.min())
    print(R_H.max())
    Temp_LCL= (1/(term_a-term_b))+55
    return(Temp_LCL)


# In[8]:


def theta_dry(theta, p, p_0=1000.0):
    """PPK equ 23, p in mb, temp in kelvin"""
    theta_d = (theta)* (p/p_0)**k_dry
    return(theta_d)


# In[9]:


def pseudoeq_potential_T(T, p, w, p_0=1000.0):
    """PPK equ 43, p in mb, temp in kelvin, w in kg/kg"""
    Temp_LCL=T_LCL(T, p, w)
    part_1=(T+C_to_K)*(p_0/p)**(0.2854*(1-0.28*w))
    part_2= (3.376/Temp_LCL)-0.00254
    part_3= (w*10**3)*(1+0.81*w)
    theta_ep= part_1* np.exp(part_2*part_3)
    return (theta_ep)


# In[10]:


def theta_ep_field(T, p, p_0=1000.0):
    """to make the moist adiabat lines. temp in celcius, p in mb, w in kg/kg"""
    w_s= sat_mixing_ratio(T,p)
    theta_ep_field= pseudoeq_potential_T(T, p, w_s, p_0)
    return(theta_ep_field)

