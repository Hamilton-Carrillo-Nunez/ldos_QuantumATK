from NanoLanguage import *
import scipy, pylab

# Delta function approximated as lim_sigma->0 exp(-(x/sigma)^2) / ( sigma*sqrt(pi) )
def DiracDelta(x,sigma):
    eta = sigma*eV
    x = x.convertTo(eV)
    return math.exp(-(x/eta)**2) / (eta*math.sqrt(math.pi))

# oscillations in the DOS disappear when sigma > 0.02
sigma = 0.018

# energy points
Ne = 501
Emin = -1.3
Emax =  2.5
Ener = numpy.linspace(Emin,Emax,Ne)

# for extending the plot energy range
dE = 1.0;

# k points
Nk = 51
kmin = 0.0
kmax = 0.5
kvector = numpy.linspace(kmin,kmax,Nk)

# declaration of arrays for the density of states
# D(E)
dosE = numpy.zeros(len(Ener))
# D(E,k)
dos = numpy.zeros(len(Ener)*len(kvector))
# for plotting purpose
dosmax = 110

# Extracting fermi level for energy reference
bulk_configuration = nlread('bulk.nc', BulkConfiguration)[0]
FermiEnergy = ChemicalPotential(bulk_configuration)
Ef = FermiEnergy.evaluate()[0]
Ef = Ef.convertTo(eV)

print "INFO: computing density of states - begin"

# loop for the k points
for n in range(len(kvector)):

    # Extracting Hamiltonian (H) and overlap matrix (S)
    H, S = calculateHamiltonianAndOverlap(bulk_configuration,kpoint=(0.0,0.0,kvector[n]))
    H = H.inUnitsOf(eV)

    # Diagonalizing H C = E S C
    w, v = scipy.linalg.eigh(H, S)
    eigenenergies =  w*eV - Ef
  
    # loop for the energies points
    for j in range(len(Ener)):

        energy = Ener[j]*eV
        ans = 0.0

        # loop for adding the Delta functions
        for i in range(len(w)):
            ans += DiracDelta(energy-eigenenergies[i],sigma)

        # D(E,k) = 2*Pi sum_i Delta(E-Ei(k))
        dos[j + len(Ener)*n] = 2.0*math.pi*ans

# density of states is computed as D(E) = 1/Nk * Sum_k D(E,k)
for i in range(len(Ener)):
    ans = 0.0

    # loop for adding all k points 
    for n in range(len(kvector)):
        ans+= dos[i + len(Ener)*n]
    dosE[i]=ans/Nk
    
print "INFO: computing density of states - done"

#plotting the density of states
pylab.figure()
pylab.plot(Ener,dosE,label='1D carbon chain DOS')
pylab.legend()
pylab.xlim(Emin-dE, Emax+dE)
pylab.ylim(0.0, dosmax)
pylab.xlabel("Energy - Ef (eV)")
pylab.ylabel("Density of States (1/eV)")
pylab.savefig('dos.png')
pylab.show()
