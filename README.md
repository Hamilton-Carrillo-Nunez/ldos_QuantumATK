# ldos_QuantumATK
This is a python script to compute the LDOS of a 1D carbon atomic chain using the quantumATK DFT tool.
You need to perform a electronic calculation using quantumATK to obtain the *.nc file.
For this particular example the LDOS is computed as:

D(E,k) = 2*Pi sum_i Delta(E-Ei(k))

and the DOS, 

D(E) = 1/Nk * Sum_k D(E,k)

Delta function is approximated as lim_sigma->0 exp(-(x/sigma)^2) / ( sigma*sqrt(pi) ).
