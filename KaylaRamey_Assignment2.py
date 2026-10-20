'''
Created on Oct. 19, 2026
Che 6834 Assignment 2: Coding Portion
@author: kaylaramey
'''

import math
import periodictable
from periodictable.constants import avogadro_number

mass_target_g = 100 # target mass of Yb-176 in grams

molar_mass_Yb176 = periodictable.Yb[176].mass # g/mol (atomic mass of Yb-176)

# initial number of Yb-176 atoms at t=0
N_A0 = (mass_target_g / molar_mass_Yb176) * avogadro_number 

# production rate of Yb-176 in s^-1
R = 5.7e-9 

# Half-life and decay constant for Yb-177 from chart of nuclides
t_half_Yb177_s = 1.911 * 3600 #1.911 hours * sec/hour = seconds
lambda_b = math.log(2)/t_half_Yb177_s

# Half-life and decay constant for Lu-177 from chart of nuclides
t_half_Lu177_s = 159.463 * 3600 # 6.647 days * 24hrs/day = 159.463 hrs * sec/hr = secs
lambda_c = math.log(2) / t_half_Lu177_s

# Conversion to seconds per week
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

# find the maximum number of Lu-177 atoms
max_atoms_Lu = 0

for i in range(1,1000): # searching in steps of 0.02 weeks, every 3 hours for peak Lu-177
    t_weeks = i * 0.02
    t_seconds = t_weeks * sec_per_week
    atoms_current = calc_Lu177(t_seconds)
        
    if atoms_current > max_atoms_Lu:
        max_atoms_Lu = atoms_current
    else: 
        break #Lu has reached its max
        
target_atoms_90 = 0.9 * max_atoms_Lu #90% saturation of Lu 

# Guess and check loop 
atoms_at_90 = 0

for j in range(100000): # Iterate from t=0 until Lu177 90% saturation reached
    t_weeks_at_90 = j * 0.0001 #looking at every 0.0001 weeks, every 1 minute
    t_seconds = t_weeks_at_90 * sec_per_week
    atoms_calculated = calc_Lu177(t_seconds)
    
    if atoms_calculated >= target_atoms_90:
        atoms_at_90 = atoms_calculated
        break
    
print(
    f'Number of weeks: {t_weeks_at_90:.3f}, Atoms of Lu-177 at 90% saturation: {atoms_at_90:.3e}'
)


