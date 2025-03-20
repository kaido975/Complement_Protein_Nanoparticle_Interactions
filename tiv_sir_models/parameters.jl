#Kinetic Parameters
kf_C3Bsu = 4.2e8  * 1e-6
kf_pC3b = 4.2e8  * 1e-6 
kf_C3bBbC3b = 3.5e6 * 1e-6
kcat_C3bBb = 1.8 
km_C3bBb = 5.9e-6 *1e6
nfc3b = 1e-6
C3 = 6.0
k_tick = kf_C3Bsu * nfc3b
k_cat = kcat_C3bBb * C3/(km_C3bBb + C3) 
div_fac = 6
K_mod = kf_pC3b*k_cat/kf_C3bBbC3b/div_fac

Na = 6.022e23 #Avagadro Number
conc_np = 2e12/Na/1e-3 #Concentration of Nanoparticles (M)

#Diffusion Coefficient Calculation
rc3b = 3.7e-9
Kb = 1.38e-23
T = 310
μ = 4e-3
D = Kb * T / (6 * π * μ * rc3b)

#Local Concentration calculation
t_half = 2e-4
r_hem = √(6*D*t_half) #Radius of 
V_hem = 0.5*4*π*r_hem^3*1000/3 #L
A_hem = π*r_hem^2
A_c3b = 1e-16
c3b_local = A_hem / A_c3b
c_c3b = c3b_local / Na / V_hem
