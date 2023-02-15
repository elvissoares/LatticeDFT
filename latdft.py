import numpy as np
from scipy import signal
import timeit

# Author: Elvis do A. Soares
# Github: @elvissoares
# Date: 2023-01-02
# Updated: 2023-01-02


" The DFT model for Lattice fluids on 3d geometries"

class LatticeDFT():
    def __init__(self):
        print('============== The DFT 3D for Lattice fluids ==============')

    def Set_Geometry(self,L,gridsize=0.01):
        if np.isscalar(L): 
            self.L = np.array([L,L,L])
        else: 
            self.L = L
        print('Geometry properties:')
        print('Lx =', self.L[0], ' A')
        print('Ly =', self.L[1], ' A')
        print('Lz =', self.L[2], ' A')
        self.Vol = self.L[0]*self.L[1]*self.L[2]
        print('Vol =',self.Vol, ' A³')

        self.delta = gridsize
        print('The gridsize is',self.delta)
            
        self.x = np.arange(0,self.L[0],self.delta) + 0.5*self.delta
        self.y = np.arange(0,self.L[1],self.delta) + 0.5*self.delta
        self.z = np.arange(0,self.L[2],self.delta) + 0.5*self.delta

        self.N = np.array([self.x.size,self.y.size,self.z.size])

        self.rho = np.zeros((self.N[0],self.N[1],self.N[2]))
        self.rhonear = np.zeros_like(self.rho)
        self.phi = np.zeros_like(self.rho)

    def Set_FluidProperties(self,epsilon=1.0):
        self.epsilon = epsilon
        
        print('Fluid properties:')
        print('epsilon/kB =', epsilon, ' K')
    
    def Set_Temperature(self,kT):

        print('Temperature =', kT, ' K')
        self.kT = kT
        self.beta = 1/self.kT

    def Set_BulkDensity(self,rhob):

        self.rhob = rhob 
            
        self.Calculate_mu()
        self.Calculate_Pressure()
        # print('Bulk Density:',self.rhob)
        # print('mu:',self.mu.round(3))
        # print('P:',self.P.round(3))

    def Set_External_Potential(self,phi):
        self.phi[:] = phi

    def Set_InitialCondition(self):
        self.rho[:] = self.rhob
        mask = self.phi>0.0
        self.rho[mask] = 1.e-16
        self.phi[mask] = 0.0
        self.Update_System()

    def Update_System(self):
        self.Calculate_NearestNeighbor()
        self.Calculate_Omega()

    def Calculate_NearestNeighbor(self):
        self.kernel = np.array([[[0, 0, 0], [0, 1, 0], [0, 0, 0]],
                                [[0, 1, 0], [1, 0, 1], [0, 1, 0]],
                                [[0, 0, 0], [0, 1, 0], [0, 0, 0]]])

        self.rhonear[:] = signal.convolve(self.rho, self.kernel, mode='same')

    def Calculate_Free_energy(self):
        self.F = self.kT*np.sum(self.rho*np.log(self.rho) + (1-self.rho)*np.log(1-self.rho))*self.delta**3 - 0.5*self.epsilon*np.sum(self.rho*self.rhonear)*self.delta**3 

    def Calculate_Omega(self):
        self.Calculate_Free_energy()
        self.Omega = self.F + np.sum((self.phi-self.mu)*self.rho)*self.delta**3

    def Calculate_mu(self):
        self.mu = self.kT*np.log(self.rhob/(1-self.rhob)) - 6*self.epsilon*self.rhob

    def Calculate_Pressure(self):
        self.P = - self.kT*(self.rhob*np.log(self.rhob) + (1-self.rhob)*np.log(1-self.rhob)) + 3*self.epsilon*self.rhob**2 +self.mu*self.rhob

    def Calculate_Equilibrium(self,alpha0=0.2,dt=0.005,rtol=1e-3,atol=1e-5,logoutput=False):

        # print('---- Obtaining the thermodynamic equilibrium ----')

        # Fire algorithm
        Ndelay = 20
        Nmax = 10000
        finc = 1.1
        fdec = 0.5
        fa = 0.99
        Nnegmax = 2000
        dtmax = 10*dt
        dtmin = 0.02*dt
        alpha = alpha0
        Npos = 0
        Nneg = 0

        starttime = timeit.default_timer()

        lnrho = np.log(self.rho)
        V = np.zeros_like(self.rho)
        F = -self.rho*(self.kT*np.log(self.rho/(1-self.rho)) - self.epsilon*self.rhonear - self.mu + self.phi)*self.delta**3
        error0 = max(np.abs(F.min()),F.max())

        for i in range(Nmax):

            P = (F*V).sum() # dissipated power
            if (P>0):
                Npos = Npos + 1
                if Npos>Ndelay:
                    dt = min(dt*finc,dtmax)
                    alpha = alpha*fa
            else:
                Npos = 0
                Nneg = Nneg + 1
                if Nneg > Nnegmax: break
                if i> Ndelay:
                    dt = max(dt*fdec,dtmin)
                    alpha = alpha0
                lnrho[:] += - 0.5*dt*V
                V[:] = 0.0
                self.rho[:] = np.exp(lnrho)
                self.Update_System()

            V[:] += 0.5*dt*F
            V[:] = (1-alpha)*V + alpha*F*np.linalg.norm(V)/np.linalg.norm(F)
            lnrho[:] += dt*V
            self.rho[:] = np.exp(lnrho)
            self.Update_System()
            F = -self.rho*(self.kT*np.log(self.rho/(1-self.rho)) - self.epsilon*self.rhonear - self.mu + self.phi)*self.delta**3
            V[:] += 0.5*dt*F

            error = max(np.abs(F.min()),F.max())
            if error/error0 < rtol and error < atol: break
            if logoutput: print(i,self.Omega,error)
        self.Niter = i

        del V, F  

        self.Nads = self.rho.sum()*self.delta**3

        # print("Time to achieve equilibrium:", timeit.default_timer() - starttime, 'sec')
        # print('Number of iterations:', self.Niter)
        # print('error:', error)
        # print('---- Equilibrium quantities ----')
        # print('F =',self.F)
        # print('Omega =',self.Omega)
        # print('Nbulk =',self.rhob*self.Vol)
        # print('Nads =',self.Nads)
        # print('================================')