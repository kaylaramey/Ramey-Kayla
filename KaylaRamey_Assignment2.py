import math
from periodictable.constants import avogadro_number
from _testcapi import add_context_watcher

# 1g of Yb-176
mass_target_g = 1 #target mass of Yb-176 in grams
molar_mass_Yb176 = 173.05 # g/mol (atomic mass of Yb-176)
#avogadro_num = 6.022e23 #atoms/mol
avogadro_num = avogadro_number #atoms/mol 

# initial number of Yb-176 atoms at t=0
N_A0 = (mass_target_g / molar_mass_Yb176) * avogadro_num #3.422837762873707e+21
#print(N_A0)

#Production
R = 5.7e-9 

#Half-life and decay constant for Yb-177 from chart of nuclides
t_half_Yb177_s = 1.911 * 3600 #seconds
lambda_b = math.log(2)/t_half_Yb177_s

#Half-life and decay constant for Lu-177 from chart of nuclides
t_half_Lu177_s = 159.463 * 3600 # seconds
lambda_c = math.log(2)/ t_half_Lu177_s

#conversion to seconds per week
sec_per_week = 7 * 24 * 3600 

def calc_Lu177(t_sec):
    ''' Calculates the number of atoms of Lu-177 at time = t
    
    Args: t_sec: time in seconds
    
    Return: N_c: Number of Lu-177 atoms at time t_sec
    '''
    #terms for analytical solutions of numbers of atoms 
    first_term = math.exp(-R * t_sec) / ((lambda_b - R)*(lambda_c - R))
    second_term = math.exp(-lambda_b * t_sec) / ((R - lambda_b) * (lambda_c - lambda_b))
    third_term = math.exp(-lambda_c * t_sec) / ((R - lambda_c) * (lambda_b - lambda_c))
    
    N_c = R * N_A0 * lambda_b * (first_term + second_term + third_term)
    
    return N_c

max_atoms_Lu = 0
time_at_max_weeks = 0 



