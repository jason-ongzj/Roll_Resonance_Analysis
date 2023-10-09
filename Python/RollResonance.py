import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from math import sin, radians
from scipy.interpolate import interp1d, interp2d

# Use class to avoid repeat calling of functions and variables through
# encapsulation of data within same class object
class RollResonance:
    # Initialize class with aerodynamic characteristics
    def __init__(self, cna_data, cldelta_data, clp_data, cmq_data):
        mach_cna = np.loadtxt(cna_data, skiprows=1, usecols=[0])
        cna_val = np.loadtxt(cna_data, skiprows=1, usecols=[1])*-1
        self.cna = interp1d(mach_cna, cna_val, kind="quadratic")

        mach_clp_cldelta = np.loadtxt(cldelta_data, skiprows=1, usecols=[0])
        cldelta_val = np.loadtxt(cldelta_data, skiprows=1, usecols=[1])
        self.cldelta = interp1d(mach_clp_cldelta, cldelta_val, kind="quadratic")

        clp_val = np.loadtxt(clp_data, skiprows=1, usecols=[1])
        self.clp = interp1d(mach_clp_cldelta, clp_val, kind="quadratic")

        cmq_x = np.loadtxt(cmq_data, max_rows=1)
        cmq_y = np.loadtxt(cmq_data, skiprows=1, max_rows=1)
        cmq_z = np.loadtxt(cmq_data, skiprows=2, usecols=[0,1])
        self.cmq = interp2d(cmq_x, cmq_y, cmq_z, kind='linear')

    # Feed model with Astos input data, ensure Astos output is in the order listed
    # below. Else switch the usecols number accordingly
    def InitInputData(self, astos_output):
        self.time = np.loadtxt(astos_output, skiprows=1, usecols=[0])
        self.ixx = np.loadtxt(astos_output, skiprows=1, usecols=[1])
        self.iyy = np.loadtxt(astos_output, skiprows=1, usecols=[2])
        self.Pss = np.loadtxt(astos_output, skiprows=1, usecols=[3])
        self.V = np.loadtxt(astos_output, skiprows=1, usecols=[5])*1000
        self.mach = np.loadtxt(astos_output, skiprows=1, usecols=[6])
        self.q = np.loadtxt(astos_output, skiprows=1, usecols=[7])
        self.T = np.loadtxt(astos_output, skiprows=1, usecols=[8])*1000
        self.m = np.loadtxt(astos_output, skiprows=1, usecols=[9])*1000
        self.CoM = np.loadtxt(astos_output, skiprows=1, usecols=[10])
        self.CoP = np.loadtxt(astos_output, skiprows=1, usecols=[11])
        self.Cd = np.loadtxt(astos_output, skiprows=1, usecols=[12])

        #-----------------------------------------------
        # Stable solution
        self.AstosSS = np.loadtxt(astos_output, skiprows=1, usecols=[13])
        self.AstosAoA = np.loadtxt(astos_output, skiprows=1, usecols=[15])
        self.density = np.loadtxt(astos_output, skiprows=1, usecols=[17])

    def InitializeVars(self):
        self.P = np.zeros(self.time.size)
        self.PssCalc = np.zeros(self.time.size)
        self.Prk4 = np.zeros(self.time.size)
        self.Pdot = np.zeros(self.time.size)
        self.roll = np.zeros(self.time.size)
        self.rollRK4 = np.zeros(self.time.size)
        self.phase = np.zeros(self.time.size)
        self.trimAoA = np.zeros(self.time.size)
        self.ACG = np.zeros(self.time.size)
        self.CLss = np.zeros(self.time.size)
        self.CLacg = np.zeros(self.time.size)
        self.PomegaRatio = np.zeros(self.time.size)

        # -----------------------------------------------
        # Roots of linear 2nd ODE and trimAoA
        self.RootAReal = np.zeros(self.time.size)
        self.RootAImag = np.zeros(self.time.size)
        self.RootBReal = np.zeros(self.time.size)
        self.RootBImag = np.zeros(self.time.size)
        self.trimArmReal = np.zeros(self.time.size)
        self.trimArmImag = np.zeros(self.time.size)

        # -----------------------------------------------
        # For calculation of K1, K2 and K3 arms
        self.SideslipRate = np.zeros(self.time.size)
        self.Sideslip = np.zeros(self.time.size)
        self.AoARate = np.zeros(self.time.size)
        self.AoA = np.zeros(self.time.size)
        self.K1real = np.zeros(self.time.size)
        self.K1imag = np.zeros(self.time.size)
        self.K2real = np.zeros(self.time.size)
        self.K2imag = np.zeros(self.time.size)
        self.K3real = np.zeros(self.time.size)
        self.K3imag = np.zeros(self.time.size)

    def SetRefQuantities(self, D, Sref, cant):
        self.D = D
        self.Sref = Sref
        self.cant = cant

    # mu - orientation of thrust misalignment
    def AssumeQuantities(self, eps, gamma, Toffset, CGoffset, staticAoA, mu):
        self.eps = eps
        self.gamma = gamma
        self.Toffset = Toffset
        self.CGoffset = CGoffset
        self.staticAoA = staticAoA
        self.mu = mu

    def FindNaturalFreq(self):
        self.omega = np.zeros(self.time.size)
        for i in range(0, self.omega.size):
            cna = self.CNAInterp(i);
            self.omega[i] = math.sqrt(-cna*(self.CoM[i] - self.CoP[i])*self.q[i] * self.Sref
                / self.iyy[i])

    def CMQInterpolation(self, i):
        mach = self.mach[i]
        if (mach < 0.5):
            cmq = self.cmq(self.CoM[i], 0.5)
        elif (mach > 4):
            cmq = self.cmq(self.CoM[i], 4)
        else:
            cmq = self.cmq(self.CoM[i], mach)
        return cmq

    def CNAInterp(self, i):
        mach = self.mach[i]
        if (mach < 0.5):
            cna = self.cna(0.5);
        else:
            cna = self.cna(mach)
        return cna

    def CLPInterp(self, i):
        mach = self.mach[i]
        if (mach < 0.5):
            clp = self.clp(0.5)
        else:
            clp = self.clp(mach)
        return clp

    def CLDeltaInterp(self, i):
        mach = self.mach[i]
        if (mach < 0.5):
            cldelta = self.cldelta(0.5)
        else:
            cldelta = self.cldelta(mach)
        return cldelta

    # Calculate b value based on iteration, not storing this value
    def CalcB(self, i):
        m = self.m[i]
        Iyy = self.iyy[i]
        mach = self.mach[i]
        d = self.D
        cna = self.CNAInterp(i)
        b = self.q[i] * self.Sref * (-1 * cna * (1 - self.ixx[i] / Iyy)
            - m * d ** 2 * self.CMQInterpolation(i) /Iyy)/ (m * self.V[i] * self.omega[i]);
        return b

    def CalcTrimAoA(self, i):
        P = self.P[i]/360
        w = self.omega[i]
        b = self.CalcB(i)
        denom = math.sqrt((1 - (P**2/w**2)*(1 - self.ixx[i]/self.iyy[i]))**2 + (b*P/w)**2)
        self.trimAoA[i] = self.staticAoA/denom

    def CalcPhaseAngle(self, i):
        P = self.P[i] / 360
        w = self.omega[i]
        b = self.CalcB(i)
        denom = 1 - (P**2/w**2)*(1 - self.ixx[i]/self.iyy[i])
        self.phase[i] = math.degrees(math.atan((b*P/w)/denom))
        self.PomegaRatio[i] = P/ self.omega[i]

    def CalcDynamicRoll(self, i, start):
        d = self.D
        q = self.q[i]
        Sref = self.Sref
        mach = self.mach[i]
        ixx = self.ixx[i]
        cant = self.cant
        delta_t = (self.time[i+1] - self.time[i])
        self.CalcTrimAoA(i)
        self.CalcPhaseAngle(i)

        cna = self.CNAInterp(i)
        clp = self.CLPInterp(i)
        cldelta = self.CLDeltaInterp(i)

        # Find solution to first order ODE using Laplace Transform.
        A = q*Sref*d*clp*d*0.5/(self.V[i]*ixx)
        B = q*Sref*d*(cldelta*cant- self.CGoffset*cna*self.trimAoA[i]*sin(radians(self.roll[i] - self.gamma)))/ixx

        if(self.time[i+1] <= 20.1):
            self.P[i+1] = -B/A
        else:
            self.P[i+1] = -B/A + (self.P[2099] + B/A)*math.exp(A*(self.time[i]-self.time[2099]))

        self.roll[i+1] = self.roll[i] + 0.5*(self.P[i+1] + self.P[i])*delta_t

    def CalcSteadyRoll(self, i):
        clp = self.CLPInterp(i)
        cldelta = self.CLDeltaInterp(i)
        delta_t = (self.time[i + 1] - self.time[i])

        self.PssCalc[i] = -cldelta*self.cant*2*self.V[i]/(clp*self.D)
        self.roll[i + 1] = self.roll[i] + 0.5 * (self.P[i + 1] + self.P[i]) * delta_t

    def CalcClContributions(self):
        for i in range(101, 2100):
            cna = self.CNAInterp(i)
            clp = self.CLPInterp(i)
            cldelta = self.CLDeltaInterp(i)
            mach = self.mach[i]
            d = self.D
            self.CLss[i] = self.cant * cldelta + clp*self.P[i]*d*0.5/self.V[i]
            self.CLacg[i] = -self.CGoffset * cna*self.trimAoA[i]*sin(radians(self.roll[i] - self.gamma))
            self.PomegaRatio[i] = (self.P[i]/360)/self.omega[i]

#-------------------------------------------------------------------------------------------------------------------
# 2nd part of roll resonance: finding sideslip and AoA values to determine the coning motion of the rocket in flight

    def CalcA1(self, i):
        mach = self.mach[i]
        m = self.m[i]
        d = self.D
        Sref = self.Sref
        V = self.V[i]
        I = self.iyy[i]
        g = 9.81
        Isp = 233
        cna = self.CNAInterp(i)
        A1 = self.q[i] * Sref * (cna - m* d**2 *self.cmq(self.CoM[i], mach)/I)/(m*V) \
             + self.T[i]/(m*V) + self.T[i]*self.Toffset**2 /(I*g*Isp)
        print(i, "A1:", A1)
        return A1

    def CalcA2(self, i):
        A2 = self.P[i]/360*(2- self.ixx[i]/self.iyy[i])
        print(i, "A2:", A2)
        return A2

    def CalcB1(self, i):
        Sref = self.Sref
        d = self.D
        I = self.iyy[i]
        cna = self.CNAInterp(i)
        B1 = cna*(self.CoM[i] - self.CoP[i])*self.q[i]*Sref/I \
             + (self.P[i]/360)**2 * (1 - self.ixx[i]/I)
        print(i, "B1:", B1)
        return B1

    def CalcB2(self, i):
        P = self.P[i]/360
        Ixx = self.ixx[i]
        I = self.iyy[i]
        m = self.m[i]
        V = self.V[i]
        Sref = self.Sref
        d = self.D
        g = 9.81
        Isp = 233
        T = self.T[i]
        mach = self.mach[i]
        cna = self.CNAInterp(i)
        B2 = -P*(self.q[i]*Sref/(m*V))*(cna*(1-Ixx/I) - m*d**2/I*(self.cmq(self.CoM[i], mach)))- \
             P*(T/(m*V))*(1-Ixx/I) - self.Pdot[i]/360 - P*T*self.Toffset**2/(I*g*Isp)
        return B2

    def CalcC1(self, i):
        mach = self.mach[i]
        Sref = self.Sref
        gamma = radians(self.gamma)
        mu = radians(self.mu)
        d = self.D
        I = self.iyy[i]
        lamda = radians(90 - self.staticAoA)
        cna = self.CNAInterp(i)
        Cmo = cna * (self.CoM[i] - self.CoP[i]) * self.staticAoA/d
        C1 = self.q[i]*Sref*d*(Cmo*math.cos(lamda) - self.Cd[i]*self.CGoffset*math.cos(gamma))/I + \
             self.T[i]*(self.CGoffset + self.Toffset)*self.eps*math.cos(mu)/I
        return C1

    def CalcC2(self, i):
        mach = self.mach[i]
        Sref = self.Sref
        gamma = radians(self.gamma)
        mu = radians(self.mu)
        d = self.D
        I = self.iyy[i]
        lamda = radians(90 - self.staticAoA)
        cna = self.CNAInterp(i)
        Cmo = cna * (self.CoM[i] - self.CoP[i]) * self.staticAoA/d
        C2 = self.q[i]*Sref*d*(Cmo*math.sin(lamda) - self.Cd[i]*self.CGoffset*math.sin(gamma))/I + \
             self.T[i]*(self.CGoffset + self.Toffset)*self.eps*math.sin(mu)/I
        return C2

    def CalcRoots(self, i):
        A1 = self.CalcA1(i)
        A2 = self.CalcA2(i)
        B1 = self.CalcB1(i)
        B2 = self.CalcB2(i)
        C1 = self.CalcC1(i)
        C2 = self.CalcC2(i)
        A = complex(A1, A2)
        B = complex(B1, B2)
        C = complex(C1, C2)
        a = 0.5*(-A+cmath.sqrt(A**2 + 4*B))
        b = 0.5*(-A-cmath.sqrt(A**2 + 4*B))

        print(i, "a:", a)
        print(i, "b:", b)

        self.RootAReal[i] = a.real
        self.RootAImag[i] = a.imag

        self.RootBReal[i] = b.real
        self.RootBImag[i] = b.imag

        trimArm = -C/B
        self.trimArmReal[i] = trimArm.real
        self.trimArmImag[i] = trimArm.imag

    def InitializeBC(self, i, sideslip_rate, sideslip, aoa_rate, aoa):
        self.SideslipRate[i] = sideslip_rate
        self.Sideslip[i] = sideslip
        self.AoARate[i] = aoa_rate
        self.AoA[i] = aoa

    # Refer to Harold Vaugh's paper on Tricyclic Theory.
    def CalcComplexPlane(self, i):
        # Initial BCs
        xi_rate_prev = complex(self.SideslipRate[i-1], self.AoARate[i-1])
        xi_prev = complex(self.Sideslip[i-1], self.AoA[i-1])

        rootA = complex(self.RootAReal[i], self.RootAImag[i]) # corresponds to K1 arm
        rootB = complex(self.RootBReal[i], self.RootBImag[i]) # corresponds to K2 arm
        trimArm = complex(self.trimArmReal[i], self.trimArmImag[i]) # corresponds to K3 arm

        print(i, "rootA:", rootA)
        print(i, "rootB:", rootB)

        p = self.P[i]/360
        delta_t = self.time[i] - self.time[i-1]

        # Check units for roll rate (ascertain whether its in Hz or deg/s)
        K3 = trimArm/cmath.exp(1j*p*delta_t)
        K1 = (xi_rate_prev - rootB * xi_prev - (1j*p - rootB)*K3)/(rootA - rootB)
        K2 = (xi_rate_prev - rootA * xi_prev - (1j*p - rootA)*K3)/(rootB - rootA)

        # K1 = (xi_rate_prev - rootB * xi_prev)/ (rootA - rootB)
        # K2 = (xi_rate_prev - rootA * xi_prev)/ (rootB - rootA)

        self.K1real[i] = K1.real
        self.K1imag[i] = K1.imag
        self.K2real[i] = K2.real
        self.K2imag[i] = K2.imag
        self.K3real[i] = K3.real
        self.K3imag[i] = K3.imag

        # Use epicyclic theory to derive the AoA and sideslip using only K1 and K2 arms, then superimpose with K3 to
        # get the final coning motion.
        xi_current = K1*cmath.exp(rootA*delta_t) + K2*cmath.exp(rootB*delta_t) + trimArm
        print(i, "K1:", K1)
        print(i, "K2:", K2)
        print(i, "trimArm:", trimArm)
        self.AoA[i] = xi_current.imag
        self.Sideslip[i] = xi_current.real

        # Use central difference scheme to calculate sideslip rate and AoA rate
        self.AoARate[i] = 2 * (self.AoA[i] - self.AoA[i-1])/delta_t - self.AoARate[i-1]
        self.SideslipRate[i] = 2 * (self.Sideslip[i] - self.Sideslip[i-1])/delta_t - self.SideslipRate[i-1]
