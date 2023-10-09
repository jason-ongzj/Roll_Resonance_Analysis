import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
import RollResonance
from math import sin, radians
from scipy.interpolate import interp1d, interp2d

# ------------------------------------------------------------------------------------------------------------------
x = RollResonance.RollResonance("CN_A.txt", "CL_Delta.txt", "CL_P.txt", "CMQ.txt")
x.InitInputData("Dynamic_Roll.txt")

# Diameter, Reference area, Cant angle
fin_cant = 0.35
x.SetRefQuantities(0.205, 0.03306, fin_cant)
x.FindNaturalFreq()
x.InitializeVars()

static_AoA = 0.05
cg_offset = 0.005
Toffset = 0.001
title = "Static AoA: " + str(static_AoA) + "°, CG Offset: " + str(cg_offset) + "m, Thrust Offset: " + str(Toffset) + "m"

# eps(rad), gamma(deg), Toffset(m), CGoffset(m), staticAoA(deg), mu
x.AssumeQuantities(0.00436, 22.5, Toffset, cg_offset, static_AoA, 0)

x.InitializeBC(100, -4.3138E-15, -0.003732135, 0, 2.144275746)

for i in range(101, x.time.size-1):
    x.CalcDynamicRoll(i, 101)
    x.CalcSteadyRoll(i)
    if (i <= 2100):
        x.CalcRoots(i)
        x.CalcComplexPlane(i)

x.CalcClContributions()

# ------------------------------------------------------------------------------------------------------------------
# Figure plots
plt.figure()
plt.plot(x.PomegaRatio[101:2100], x.trimAoA[101:2100]/x.staticAoA)
plt.axis([0,2,0,30])
plt.ylabel(r'$\alpha_{TRIM}/\alpha_{STATIC}$', fontsize=12)
plt.xlabel(r'$P/\omega$', fontsize=12)
plt.grid()
plt.savefig("TrimAoA",bbox_inches = 'tight', pad_inches = 0.3)

plt.clf()
plt.plot(x.time, x.PssCalc/360,label='Steady-State Equation')
plt.plot(x.time, x.P/360,label='Dynamic - Laplace Transform')
plt.plot(x.time, x.omega, label='Natural Frequency')
plt.plot(x.time, x.Pss/360, label="Astos")
plt.ylabel('Frequency (Hz)', fontsize=12)
plt.xlabel('Time (s)', fontsize=12)
plt.axis([0,100,0,8])
plt.grid()
plt.legend()
plt.savefig("Pitch_Roll_Freq",bbox_inches = 'tight', pad_inches = 0.3)

plt.clf()
plt.plot(x.time, x.P/360,label=str(fin_cant)+'\N{DEGREE SIGN} Fin Cant')
plt.plot(x.time, x.omega, label='Natural Frequency')
plt.ylabel('Frequency (Hz)', fontsize=12)
plt.xlabel('Time (s)', fontsize=12)
plt.axis([0,100,0,12])
plt.grid()
plt.legend()
plt.savefig("Design_Chart_Roll_Freq",bbox_inches = 'tight', pad_inches = 0.3)

plt.clf()
plt.plot(x.PomegaRatio[101:2100], x.CLacg[101:2100],label='Effect due to ACG')
plt.plot(x.PomegaRatio[101:2100], x.CLss[101:2100], label='Effect due to Cl_p and Cl_delta')
plt.ylabel(r'$C_L$', fontsize=12)
plt.xlabel(r'$P/\omega$', fontsize=12)
plt.axis([0,1.2,-0.04,0.12])
plt.grid()
plt.legend()
plt.savefig("CL",bbox_inches = 'tight', pad_inches = 0.3)

plt.clf()
plt.plot(x.PomegaRatio, x.phase)
plt.grid()
plt.axis([0,2,-90,90])
plt.ylabel("Phase Angle (\N{DEGREE SIGN})", fontsize=12)
plt.xlabel(r'$P/\omega$', fontsize=12)
plt.savefig("Phase_Angle",bbox_inches = 'tight', pad_inches = 0.3)

plt.clf()
fig, ax = plt.subplots(figsize=(12, 10))
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
ax.set_aspect(1)
plt.plot(x.Sideslip[101:2100], x.AoA[101:2100], label="Tricyclic")
plt.plot(x.trimArmReal[101:2100], x.trimArmImag[101:2100], label="Trim Arm")
plt.plot(x.AstosSS[101:2100], x.AstosAoA[101:2100], label="ASTOS Output")
plt.ylabel("Angle of Attack (\N{DEGREE SIGN})", fontsize=12)
plt.xlabel("Sideslip Angle (\N{DEGREE SIGN})", fontsize=12)
plt.legend(fontsize=12)
plt.title(title, fontweight="bold", fontsize=14)
plt.grid()
plt.savefig("AoA_Sideslip",bbox_inches = 'tight', pad_inches = 0.3)
