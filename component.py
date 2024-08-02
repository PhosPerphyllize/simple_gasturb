import numpy as np
import os
import copy
from gasFun import Altitude, H2OCp, N2Cp, O2Cp, CO2Cp, airk, Ma2tau, Ma2pi, pi2Ma
from log import Log
log = Log(path="./componet_log.txt")

MO2, MN2, = 32, 28.01
MCO2, MH2O = 44.0095, 18.0152
MH2, MCH4 = 2.01588, 16.043
EFS = 1e-8
FANsNUM = 10

def interpolate1D(x, point1, point2):
    # linear interpolate, point1,point2:[x1,y1], array or numpy
    # flag: if outside the points
    point1, point2 = np.array(point1), np.array(point2)
    flag = False
    ratio = (x - point1[0])/(point2[0] - point1[0])
    if ratio <0 or ratio>1:
        flag = True
    return point1[1] + ratio*(point2[1] - point1[1]), flag

def interpolate2D_charline(nl, mass, nlLine, massTable, piTable, effTable):
    # nlLine, massTable mass: from small to large
    # table not contain nlline   # nl is vertical in tables and mass in horizontal
    # flag: if outside the points
    nlLine = np.array(nlLine)
    massTable, piTable, effTable = np.array(massTable), np.array(piTable), np.array(effTable)
    if max(nlLine.shape) != massTable.shape[0] or massTable.shape != piTable.shape or massTable.shape != effTable.shape:
        raise Exception("Tables and title shape not match!")

    nlnum, massnum = massTable.shape
    # get nl left:nl_id, right value:nl_id+1, ready to interpolate
    if nl >= nlLine[-1]:
        nlid = -2
    elif nl <= nlLine[0]:
        nlid = 0
    else:
        for i in range(nlnum - 1):
            if nl >= nlLine[i] and nl <= nlLine[i+1]:
                nlid = i
                break

    flag = False
    massLine = np.zeros(massnum)
    piLine = np.zeros(massnum)
    effLine = np.zeros(massnum)
    for i in range(massnum):
        try:
            massLine[i],flagre = interpolate1D(nl, [nlLine[nlid], massTable[nlid, i]], [nlLine[nlid+1], massTable[nlid+1, i]])
        except:
            print("error")
        piLine[i],_ = interpolate1D(nl, [nlLine[nlid], piTable[nlid, i]], [nlLine[nlid+1], piTable[nlid+1, i]])
        effLine[i],_ = interpolate1D(nl, [nlLine[nlid], effTable[nlid, i]], [nlLine[nlid+1], effTable[nlid+1, i]])
        if flagre:
            flag = True

    # get mass left:massid, right value:massid+1, ready to interpolate
    if mass >= massLine[-1]:
        massid = -2
    elif mass <= massLine[0]:
        massid = 0
    else:
        for i in range(massnum - 1):
            if mass >= massLine[i] and mass <= massLine[i + 1]:
                massid = i
                break
    pi, flag1 = interpolate1D(mass, [massLine[massid], piLine[massid]], [massLine[massid+1], piLine[massid+1]])
    eff, flag2 = interpolate1D(mass, [massLine[massid], effLine[massid]], [massLine[massid+1], effLine[massid+1]])

    return pi, eff, flag or flag1 or flag2

def interpolate_2title(x, y, xtitle, ytitle, table):
    # x_title, y_title: from small to large
    # x is vertical in table and y in horizontal
    # flag: if outside the points
    xtitle, ytitle, table = np.array(xtitle), np.array(ytitle), np.array(table)
    if max(xtitle.shape) != table.shape[0] or max(ytitle.shape) != table.shape[1]:
        raise Exception("Table and title shape not match!")

    xnum,ynum = table.shape[0],table.shape[1]

    if x <= xtitle[0]:
        xid = 0
    elif x >= xtitle[-1]:
        xid = -2
    else:
        for i in range(xnum-1):
            if x >= xtitle[i] and x <= xtitle[i+1]:
                xid = i
                break

    if y <= ytitle[0]:
        yid = 0
    elif y >= ytitle[-1]:
        yid = -2
    else:
        for i in range(ynum-1):
            if y >= ytitle[i] and y <= ytitle[i+1]:
                yid = i
                break

    yleft,flag1 = interpolate1D(x, [xtitle[xid], table[xid, yid]], [xtitle[xid+1], table[xid+1, yid]])
    yright,flag2 = interpolate1D(x, [xtitle[xid], table[xid, yid+1]], [xtitle[xid+1], table[xid+1, yid+1]])

    result,flag3 = interpolate1D(y, [ytitle[yid], yleft], [ytitle[yid+1], yright])

    return result, flag1 or flag2 or flag3

class Gas():
    def __init__(self, Ts=288.15, Ps=101325.0, Ma=0.0, G=1.0):
        self.G = G   # gas mass flow kg/s

        # volume(mol) fraction
        self.N2ratio = 0.78
        self.O2ratio = 0.22
        self.CO2ratio = 0.0
        self.H2Oratio = 0.0
        self.Rg()

        # mass fraction
        self.vol2mratio()

        T = Ts
        for i in range(5):
            k = self.k(0.5*(Ts+T))
            T = Ts*Ma2tau(Ma, k)
        self.T = T                # gas total temperature K
        self.P = Ps*Ma2pi(Ma, k)  # gas total prassure    Pa

    def k(self, T):
        # T: static temprature
        Cp = self.Cp(T)
        return Cp/(Cp - self.Rg_)

    def Cp(self, T):
        # T: static temprature
        self.vol2mratio()
        if abs((self.N2mratio + self.O2mratio + self.CO2mratio + self.H2Omratio)-1) > EFS:
            raise Exception("air ratio sum != 1")
        Cp_ = self.N2mratio*N2Cp(T) + self.CO2mratio*CO2Cp(T) + self.O2mratio*O2Cp(T) + self.H2Omratio*H2OCp(T)
        return Cp_

    def Rg(self):
        # T: static temprature
        if abs((self.N2ratio + self.O2ratio + self.CO2ratio + self.H2Oratio)-1) > EFS:
            raise Exception("air ratio sum != 1")
        M = MN2*self.N2ratio + MO2*self.O2ratio + MH2O*self.H2Oratio + MCO2*self.CO2ratio
        self.Rg_ = 8314.5/M   # J/(kg K)
        return self.Rg_

    def vol2mratio(self):
        Mn2, Mo2, Mco2, Mh2o = self.N2ratio*MN2, self.O2ratio*MO2, self.CO2ratio*MCO2, self.H2Oratio*MH2O
        msum = Mn2 + Mo2 + Mco2 + Mh2o
        self.N2mratio = Mn2/msum
        self.O2mratio = Mo2/msum
        self.CO2mratio = Mco2/msum
        self.H2Omratio = 1 - self.N2mratio - self.O2mratio - self.CO2mratio
        self.Rg()
        return (self.N2mratio, self.O2mratio, self.CO2mratio, self.H2Omratio)

    def m2volratio(self):
        Vn2, Vo2, Vco2, Vh2o = self.N2mratio/MN2, self.O2mratio/MO2, self.CO2mratio/MCO2, self.H2Omratio/MH2O
        Vsum = Vn2 + Vo2 + Vco2 + Vh2o
        self.N2ratio = Vn2/Vsum
        self.O2ratio = Vo2/Vsum
        self.CO2ratio = Vco2/Vsum
        self.H2Oratio = 1 - self.N2ratio - self.O2ratio - self.CO2ratio
        self.Rg()
        return (self.N2ratio, self.O2ratio, self.CO2ratio, self.H2Oratio)

    def expand(self, pi, iteration=5):
        # pi > 1 : self.P/expand P
        if pi < 1.0:
            raise Exception("pi must > 1 pi:%f"%(pi))
        Ts = self.T  # complete expand
        Ma,k = 0.0, 0.0
        for i in range(iteration):
            k = airk(0.5 * (self.T + Ts))
            Ma = pi2Ma(pi, k)
            Ts = self.T / Ma2tau(Ma, k)
        Ps = self.P / pi
        V = Ma*np.sqrt(k*self.Rg_*Ts)   # Velocity
        rhos = Ps/(self.Rg_*Ts)
        A = self.G/(rhos*V)  # Aera
        return Ts, Ps, Ma, V, rhos, A

    def static(self, Ma, iteration=5):
        if Ma <= 0.0:
            return self.T, self.P, 0.0, self.P/(self.Rg_*self.T), -1
        Ts = self.T
        k = 1.4
        for i in range(iteration):
            k = airk(0.5 * (self.T + Ts))
            Ts = self.T / Ma2tau(Ma, k)
        Ps = self.P / Ma2pi(Ma, k)
        V = Ma*np.sqrt(k*self.Rg_*Ts)
        rhos = Ps / (self.Rg_ * Ts)
        A = self.G / (rhos * V)  # Aera
        return Ts, Ps, V, rhos, A

class Fan05():
    def __init__(self):
        self.Aout = 1.0
        self.nozzleCf = 0.995

        self.Gd = 1.0   # design para without table
        self.pi = 1.6
        self.eff = 0.78

    def workDesign(self, Height, Ma, G, pi=1.6, eff=0.78, Pamb=101325, iteration=5):
        # design Aout, para will insert in para without table
        gas  = Gas(Height, Ma, G)
        _, _, Vin, _, _ = gas.static(Ma)

        Pout = pi*gas.P
        Tout = gas.T
        for i in range(iteration):
            Cp = gas.Cp(0.5*(Tout+gas.T))
            k = Cp/(Cp - gas.Rg_)
            Tout_ad = gas.T * np.power(pi, (k-1)/k)
            W_ad = gas.G*Cp*(Tout_ad - gas.T)
            W = W_ad/eff
            Tout = gas.T + W/(gas.G*Cp)

        gas.P, gas.T = Pout, Tout
        Ts, Ps, Ma, V, rhos, A = gas.expand(Pout/Pamb)
        self.Aout = A
        self.Gd, self.pi, self.eff = gas.G, pi, eff
        thrust = gas.G*V*self.nozzleCf + A*(Ps-Pamb) - gas.G*Vin
        return gas, thrust

    def getTable(self, path=os.path.join(os.getcwd(), "CharData2"), Gd=1.0, pid=1.0):
        # return the total corrected table
        mass = np.loadtxt(os.path.join(path, "FanMassFlow.txt"))
        eff = np.loadtxt(os.path.join(path, "FanEfficiency.txt"))
        pi = np.loadtxt(os.path.join(path, "FanPressureRatio.txt"))

        if mass.shape != eff.shape or mass.shape != pi.shape or pi.shape != eff.shape:
            raise Exception("MassFlow.txt and Efficiency.txt and PressureRatio.txt shape not match!")

        nl_num = mass.shape[0]
        nl = np.zeros(nl_num)
        for i in range(nl_num):
            if mass[i, 0] != eff[i, 0] or mass[i, 0] != pi[i, 0] or pi[i, 0] != eff[i, 0]:
                raise Exception("MassFlow.txt and Efficiency.txt and PressureRatio.txt nl not match!")
            nl[i] = mass[i, 0]

        Rg, P, T = 287.06, 100818, 288.15  # table calculate condition

        self.nlLine, self.massTable = nl/np.sqrt(Rg*T), mass[:, 1:]*np.sqrt(Rg*T)/P*Gd
        self.effTable, self.piTable = eff[:, 1:], pi[:, 1:]*pid
        return self.nlLine, self.massTable, self.piTable, self.effTable

    def workTable(self, gas:Gas, nl, Gd=1.0, iteration=5):
        nlcor = nl/np.sqrt(gas.Rg_*gas.T)
        nldot = nlcor / (1/np.sqrt(287.06*288.15))
        Gcor = gas.G*np.sqrt(gas.Rg_*gas.T)/gas.P
        Gdot = Gcor / (Gd*np.sqrt(287.06*288.15)/101325)
        pi, eff, flag = interpolate2D_charline(nlcor, Gcor, self.nlLine, self.massTable, self.piTable, self.effTable)

        Pout = pi * gas.P
        Tout = gas.T
        for i in range(iteration):
            Cp = gas.Cp(0.5 * (Tout + gas.T))
            k = Cp / (Cp - gas.Rg_)
            Tout_ad = gas.T * np.power(pi, (k - 1) / k)
            W_ad = gas.G * Cp * (Tout_ad - gas.T)
            W = W_ad / eff
            Tout = gas.T + W / (gas.G * Cp)

        gas.P, gas.T = Pout, Tout
        return gas, W, Gdot, nldot, pi, eff

class Fan():
    def __init__(self,leakG_per=0.00,
                 ifTable=False, path=os.path.join(os.getcwd(), "CharData2"), ifShape=False, pi=4.0, mass=1.0, eff=0.8):
        self.leakG_per = leakG_per  # leakage gas massflow/total gas massflow, G will take outlet the fan
        if ifTable:
            # use characteristic table to calculate mass, pi, eff
            self.getTable(path)
        if ifShape:
            # re determined the pi, mass, eff
            massd = 1.0
            gas = Gas(G=massd)
            _, _, _, _, pid, effd = self.workTable(gas, nl=1.0)
            self.massTable, self.piTable, self.effTable = self.massTable*mass/massd, self.piTable*pi/pid, self.effTable*eff/effd

    def workDesign(self, gas:Gas, pi, eff, iteration=5):
        Pout = pi * gas.P
        Tout = gas.T
        for i in range(iteration):
            Cp = gas.Cp(0.5 * (Tout + gas.T))
            k = Cp / (Cp - gas.Rg_)
            Tout_ad = gas.T * np.power(pi, (k - 1) / k)
            W_ad = gas.G * Cp * (Tout_ad - gas.T)
            W = W_ad / eff
            Tout = gas.T + W / (gas.G * Cp)

        gas.G = gas.G * (1 - self.leakG_per)
        gas.P, gas.T = Pout, Tout
        return gas, W,

    def getTable(self, path=os.path.join(os.getcwd(), "CharData2")):
        # return the total corrected table
        mass = np.loadtxt(os.path.join(path, "FanMassFlow.txt"))
        eff = np.loadtxt(os.path.join(path, "FanEfficiency.txt"))
        pi = np.loadtxt(os.path.join(path, "FanPressureRatio.txt"))

        if mass.shape != eff.shape or mass.shape != pi.shape or pi.shape != eff.shape:
            raise Exception("MassFlow.txt and Efficiency.txt and PressureRatio.txt shape not match!")

        nl_num = mass.shape[0]
        nl = np.zeros(nl_num)
        for i in range(nl_num):
            if mass[i, 0] != eff[i, 0] or mass[i, 0] != pi[i, 0] or pi[i, 0] != eff[i, 0]:
                raise Exception("MassFlow.txt and Efficiency.txt and PressureRatio.txt nl not match!")
            nl[i] = mass[i, 0]

        self.Rg, self.P, self.T = 287.06, 101325.0, 288.15  # table calculate condition

        self.nlLine, self.massTable = nl/np.sqrt(Rg*T), mass[:, 1:]*np.sqrt(Rg*T)/P
        self.effTable, self.piTable = eff[:, 1:], pi[:, 1:]
        return self.nlLine, self.massTable, self.piTable, self.effTable

    def workTable(self, gas:Gas, nl, Gd=1.0, iteration=5):
        nlcor = nl/np.sqrt(gas.Rg_*gas.T)
        nldot = nlcor / (1/np.sqrt(self.Rg*self.T))
        Gcor = gas.G*np.sqrt(gas.Rg_*gas.T)/gas.P
        Gdot = Gcor / (Gd*np.sqrt(self.Rg*self.T)/self.P)
        pi, eff, flag = interpolate2D_charline(nlcor, Gcor, self.nlLine, self.massTable, self.piTable, self.effTable)
        if pi<1.00001:
            pi=1.00001
        if eff<0.1:
            eff=0.1

        Pout = pi * gas.P
        Tout = gas.T
        for i in range(iteration):
            Cp = gas.Cp(0.5 * (Tout + gas.T))
            k = Cp / (Cp - gas.Rg_)
            Tout_ad = gas.T * np.power(pi, (k - 1) / k)
            W_ad = gas.G * Cp * (Tout_ad - gas.T)
            W = W_ad / eff
            Tout = gas.T + W / (gas.G * Cp)

        gas.G = gas.G * (1-self.leakG_per)
        gas.P, gas.T = Pout, Tout
        return gas, W, Gdot, nldot, pi, eff

class Intake():
    def __init__(self, sigma=0.995, Ain=0.1):
        self.sigma = sigma
        self.Ain = Ain

    def work(self, gas:Gas, Ma):
        if Ma < 0:
            raise Exception("Ma must > 0")
        elif Ma <= 1:
            sigma = self.sigma
        else:
            sigma = self.sigma * (1 - 0.075 * np.power(Ma - 1, 1.35))

        gas.P = gas.P*sigma
        if Ma<=0.0:
            inlet_drag = 0.0
            self.Gin = 0.0
        else:
            Ts, Ps, V, rhos, A = gas.static(Ma)
            inlet_drag = gas.G * V
            self.Gin = self.Ain*rhos*V
        return gas, inlet_drag

class Duct():
    def __init__(self, sigma=0.995):
        self.sigma = sigma
    def work(self, gas:Gas):
        gas.P = gas.P * self.sigma
        return gas

class Compressor():
    def __init__(self,leakG_per=0.01, ifTable=False, path=os.path.join(os.getcwd(), "CharData2"), ifShape=False, pi=4.0, mass=1.0, eff=0.8):
        self.leakG_per = leakG_per  # leakage gas massflow/total gas massflow, G will take outlet the comp
        if ifTable:
            # use characteristic table to calculate mass, pi, eff
            self.getTable(path)
        if ifShape:
            # re determined the pi, mass, eff
            massd = 0.49
            gas = Gas(G=massd)
            _, _, _, _, pid, effd = self.workTable(gas, nl=1.0)
            self.massTable, self.piTable, self.effTable = self.massTable*mass/massd, self.piTable*pi/pid, self.effTable*eff/effd

    def workDesign(self, gas:Gas, pi, eff, iteration=5):
        Pout = pi * gas.P
        Tout = gas.T
        for i in range(iteration):
            Cp = gas.Cp(0.5 * (Tout + gas.T))
            k = Cp / (Cp - gas.Rg_)
            Tout_ad = gas.T * np.power(pi, (k - 1) / k)
            W_ad = gas.G * Cp * (Tout_ad - gas.T)
            W = W_ad / eff
            Tout = gas.T + W / (gas.G * Cp)

        gas.G = gas.G * (1 - self.leakG_per)
        gas.P, gas.T = Pout, Tout
        return gas, W,

    def getTable(self, path=os.path.join(os.getcwd(), "CharData2")):
        # return the total corrected table
        mass = np.loadtxt(os.path.join(path, "CompMassFlow.txt"))
        eff = np.loadtxt(os.path.join(path, "CompEfficiency.txt"))
        pi = np.loadtxt(os.path.join(path, "CompPressureRatio.txt"))

        if mass.shape != eff.shape or mass.shape != pi.shape or pi.shape != eff.shape:
            raise Exception("MassFlow.txt and Efficiency.txt and PressureRatio.txt shape not match!")

        nl_num = mass.shape[0]
        nl = np.zeros(nl_num)
        for i in range(nl_num):
            if mass[i, 0] != eff[i, 0] or mass[i, 0] != pi[i, 0] or pi[i, 0] != eff[i, 0]:
                raise Exception("MassFlow.txt and Efficiency.txt and PressureRatio.txt nl not match!")
            nl[i] = mass[i, 0]

        self.Rg, self.P, self.T = 287.06, 100818, 288.15  # table calculate condition

        self.nlLine, self.massTable = nl/np.sqrt(self.Rg*self.T), mass[:, 1:]*np.sqrt(self.Rg*self.T)/self.P
        self.effTable, self.piTable = eff[:, 1:], pi[:, 1:]
        return self.nlLine, self.massTable, self.piTable, self.effTable

    def workTable(self, gas:Gas, nl, Gd=0.49, iteration=5):
        nlcor = nl/np.sqrt(gas.Rg_*gas.T)
        nldot = nlcor / (1/np.sqrt(self.Rg*self.T))
        Gcor = gas.G*np.sqrt(gas.Rg_*gas.T)/gas.P
        Gdot = Gcor / (Gd*np.sqrt(self.Rg*self.T)/self.P)
        pi, eff, flag = interpolate2D_charline(nlcor, Gcor, self.nlLine, self.massTable, self.piTable, self.effTable)
        if pi<1.00001:
            pi=1.00001
        if eff<0.1:
            eff=0.1

        Pout = pi * gas.P
        Tout = gas.T
        for i in range(iteration):
            Cp = gas.Cp(0.5 * (Tout + gas.T))
            k = Cp / (Cp - gas.Rg_)
            Tout_ad = gas.T * np.power(pi, (k - 1) / k)
            W_ad = gas.G * Cp * (Tout_ad - gas.T)
            W = W_ad / eff
            Tout = gas.T + W / (gas.G * Cp)

        gas.G = gas.G * (1-self.leakG_per)
        gas.P, gas.T = Pout, Tout
        return gas, W, Gdot, nldot, pi, eff

class Combustor():
    def __init__(self, sigma=0.95, eff=0.98, LHV=119.988*1e6, fuelType="H2"):
        self.sigma = sigma
        self.eff = eff
        self.fuelLHV = LHV   # fuel low heat value
        self.fuelType = fuelType  # "H2" or "CH4"
    def workT4(self, gas:Gas, Tout=1130.0, iteration=8):
        Cpin = gas.Cp(gas.T)

        # N2 mass of inlet gas
        moln2in = gas.N2mratio * gas.G / MN2
        molo2in = gas.O2mratio * gas.G / MO2
        molh2oin = gas.H2Omratio * gas.G / MH2O
        molco2in = gas.CO2mratio * gas.G / MCO2

        Gfuel=0.0
        for i in range(iteration):
            G=gas.G+Gfuel
            Cp = 0.5*(Cpin + gas.Cp(Tout))
            Qad = G *Cp*(Tout - gas.T)
            Gfuel = Qad/self.eff/self.fuelLHV  # kg/s

            if self.fuelType == "H2":
                molFuel = Gfuel/MH2
                molh2o = molh2oin + molFuel
                molo2 = molo2in - 0.5*molFuel
                moln2, molco2 = moln2in, molco2in
                if molo2 < 0:
                    log("Warning: O2 is no enough to burn.")
                    molo2 = 0.0

                molsum = moln2 + molo2 + molh2o + molco2
                gas.N2ratio, gas.O2ratio, gas.H2Oratio = moln2/molsum, molo2/molsum, molh2o/molsum
                gas.CO2ratio = 1 - gas.N2ratio - gas.O2ratio - gas.H2Oratio
            else:
                # fuel tpye == "CH4"
                molFuel = Gfuel / MCH4
                molh2o = molh2oin + 2*molFuel
                molo2 = molo2in - 2*molFuel
                molco2 = molco2in + molFuel
                moln2 = moln2in
                if molo2 < 0:
                    log("Warning: O2 is no enough to burn.")
                    molo2 = 0.0

                molsum = moln2 + molo2 + molh2o + molco2
                gas.N2ratio, gas.O2ratio, gas.H2Oratio = moln2/molsum, molo2/molsum, molh2o/molsum
                gas.CO2ratio = 1 - gas.N2ratio - gas.O2ratio - gas.H2Oratio

            gas.vol2mratio()

        gas.G += Gfuel
        gas.P, gas.T = gas.P*self.sigma, Tout
        return gas, Gfuel

class Turbine():
    def __init__(self, leakG_per=0.0, ifTable=False, path=os.path.join(os.getcwd(), "CharData2"), ifShape=False, pi=4.0, mass=1.0, eff=0.8):
        self.leakG_per = leakG_per  # leakage gas massflow/total gas massflow, G will take out before turb
        if ifTable:
            # use characteristic table to calculate mass, pi, eff
            self.getTable(path)
        if ifShape:
            # re determined the pi, mass, eff
            massd = 0.49
            gas = Gas(G=massd)
            _, _, _, _, pid, effd = self.workTable(gas, nl=1.0)
            self.massTable, self.piTable, self.effTable = self.massTable*mass/massd, self.piTable*pi/pid, self.effTable*eff/effd

    def workDesign_pi(self, gas:Gas, pi, eff, iteration=5):
        gas.G = gas.G * (1 - self.leakG_per)
        Pout = gas.P / pi
        Tout = gas.T
        for i in range(iteration):
            Cp = gas.Cp(0.5 * (Tout + gas.T))
            k = Cp / (Cp - gas.Rg_)
            Tout_ad = gas.T * np.power(1 / pi, (k - 1) / k)
            W_ad = gas.G * Cp * (gas.T - Tout_ad)
            W = W_ad * eff
            Tout = gas.T - W / (gas.G * Cp)

        gas.P, gas.T = Pout, Tout
        return gas, W,

    def workDesign_work(self, gas:Gas, work, eff, iteration=5):
        # work, unit: W
        gas.G = gas.G * (1 - self.leakG_per)
        Tout = gas.T
        for i in range(iteration):
            Cp = gas.Cp(0.5 * (Tout + gas.T))
            k = Cp / (Cp - gas.Rg_)
            W_ad = work/eff
            Tout_ad = gas.T - W_ad / (gas.G * Cp)
            pi = np.power(Tout_ad/gas.T, k/(k-1))
            Tout = gas.T - work / (gas.G * Cp)

        gas.P, gas.T = gas.P*pi, Tout
        return gas, pi

    def getTable(self, path=os.path.join(os.getcwd(), "CharData2")):
        # return the total corrected table
        mass = np.loadtxt(os.path.join(path, "TurbMassFlow.txt"))
        eff = np.loadtxt(os.path.join(path, "TurbEfficiency.txt"))
        pi = np.loadtxt(os.path.join(path, "TurbPressureRatio.txt"))

        if mass.shape != eff.shape or mass.shape != pi.shape or pi.shape != eff.shape:
            raise Exception("MassFlow.txt and Efficiency.txt and PressureRatio.txt shape not match!")

        nl_num = mass.shape[0]
        nl = np.zeros(nl_num)
        for i in range(nl_num):
            if mass[i, 0] != eff[i, 0] or mass[i, 0] != pi[i, 0] or pi[i, 0] != eff[i, 0]:
                raise Exception("MassFlow.txt and Efficiency.txt and PressureRatio.txt nl not match!")
            nl[i] = mass[i, 0]

        self.Rg, self.P, self.T = 287.06, 391079, 1017.83  # table calculate condition

        self.nlLine, self.massTable = nl/np.sqrt(self.Rg*self.T), mass[:, 1:]*np.sqrt(self.Rg*self.T)/self.P
        self.effTable, self.piTable = eff[:, 1:], pi[:, 1:]
        return self.nlLine, self.massTable, self.piTable, self.effTable

    def workTable(self, gas:Gas, nl, Gd=0.49, iteration=5):
        gas.G = gas.G * (1 - self.leakG_per)

        nlcor = nl/np.sqrt(gas.Rg_*gas.T)
        nldot = nlcor / (1 / np.sqrt(self.Rg * self.T))
        Gcor = gas.G*np.sqrt(gas.Rg_*gas.T)/gas.P
        Gdot = Gcor / (Gd * np.sqrt(self.Rg * self.T) / self.P)
        pi, eff, flag = interpolate2D_charline(nlcor, Gcor, self.nlLine, self.massTable, self.piTable, self.effTable)
        if pi<1.00001:
            pi=1.00001
        if eff<0.1:
            eff=0.1

        Pout = gas.P/pi
        Tout = gas.T
        for i in range(iteration):
            Cp = gas.Cp(0.5 * (Tout + gas.T))
            k = Cp / (Cp - gas.Rg_)
            Tout_ad = gas.T * np.power(1/pi, (k - 1) / k)
            W_ad = gas.G * Cp * (gas.T - Tout_ad)
            W = W_ad * eff
            Tout = gas.T - W / (gas.G * Cp)

        gas.P, gas.T = Pout, Tout
        return gas, W, Gdot, nldot, pi, eff

class Nozzle():
    def __init__(self, type, Pamb=101325.0, Cd=0.98, Cf=0.98, sigma=0.99):
        if type not in ["con", "con/div"]:
            raise Exception("type not in [con, con/div]")
        self.type_ = type
        self.A8 = 1.0
        self.A9 = 1.1
        self.Pamb = Pamb
        self.sigma = sigma  # pressure loss, might overloap Cf
        self.Cf = Cf  # Thrust coef
        self.Cd = Cd  # massflow coef, deu to boundary layer

    def workDesign(self, gas:Gas):
        # input: inlet Gas, output: outlet Gas, grass Thrust (without -W2*V0)
        gas.P = gas.P*self.sigma
        P8, T8 = gas.P, gas.T   # throat

        Ts_ad, Ps_ad, Ma_ad, V_ad, _, A_ad = gas.expand(P8/self.Pamb)  # complete expend

        if Ma_ad >= 1.0:
            # throat para
            T8s, P8s, V8, _, A8 = gas.static(1.0)

            if self.type_ == "con":
                Ts_out, Ps_out = T8s, P8s
                Ma_out, V_out = 1.0, V8
                A_out = A8
            else:
                Ts_out, Ps_out = Ts_ad, Ps_ad
                Ma_out, V_out =  Ma_ad, V_ad
                A_out = A_ad
        elif Ma_ad < 1.0:
            Ts_out, Ps_out, Ma_out, V_out = Ts_ad, Ps_ad, Ma_ad, V_ad
            A_out = A_ad
            A8 = A_out
        else:
            raise Exception("Ma_ad must > 0.0, sth wrong with P8/self.Pamb")

        gro_thrust = gas.G*V_out*self.Cf + A_out/self.Cd*(Ps_out - self.Pamb)
        self.A8, self.A9 = A8/self.Cd, A_out/self.Cd
        list = (Ts_out, Ps_out, Ma_out, V_out)
        return gas, gro_thrust, self.A8, self.A9, list

    def work(self, gas:Gas):
        # work with fix nozzle geom, only
        gas.P = gas.P * self.sigma
        P8, T8 = gas.P, gas.T  # throat

        if P8>=self.Pamb:
            T8s, P8s, V8, rhos, A8 = gas.static(1.0)
            if self.type_ == "con" and self.Pamb<P8s:
                # not complete expand
                Ts_ad, Ps_ad, Ma_ad, V_ad = T8s, P8s, 1.0, V8
                rhos_ad, A_ad = rhos, A8
            else:
                Ts_ad, Ps_ad, Ma_ad, V_ad, rhos_ad, A_ad = gas.expand(P8/self.Pamb)  # complete expend
        else:
            # back flow from amb
            _, _, Ma_ad, V_ad, rhos_ad, A_ad = gas.expand(self.Pamb/P8)  # complete expend
            Ts_ad, Ps_ad = T8, self.Pamb
            Ma_ad, V_ad = -Ma_ad, -V_ad

        G = rhos_ad*V_ad*self.A9*self.Cd
        gas.G = G
        gro_thrust = G * V_ad * self.Cf + self.A9 * (Ps_ad - self.Pamb)
        # gas, grass Thrust (without -W2*V0),
        # outlet static T, ..P, outlet Ma, velocity, density, best nozzle area
        return gas, gro_thrust, G, (Ts_ad, Ps_ad, Ma_ad, V_ad, rhos_ad, A_ad)

class Detector():
    def __init__(self):
        pass
    def get(self, gas:Gas, name="detector", isPrint=False):
        if isPrint:
            log(name+":\t%.5f kg/s \t%.2f K \t%.5f kPa"%(gas.G, gas.T, gas.P/1000))

class Mixer():
    def __init__(self, sigma=0.985):
        self.sigma = sigma
    def work(self, gas1:Gas, gas2:Gas):
        gas = Gas()
        gas.G = gas1.G + gas2.G
        gas.P = (gas1.G*gas1.P + gas2.G*gas2.P)/gas.G * self.sigma
        gas.T = (gas1.G*gas1.T + gas2.G*gas2.T)/gas.G

        Mn2, Mo2  = gas1.G*gas1.N2mratio + gas2.G*gas2.N2mratio, gas1.G*gas1.O2mratio + gas2.G*gas2.O2mratio
        Mco2 = gas1.G*gas1.CO2mratio + gas2.G*gas2.CO2mratio
        gas.N2mratio = Mn2 / gas.G
        gas.O2mratio = Mo2 / gas.G
        gas.CO2mratio = Mco2 / gas.G
        gas.H2Omratio = 1 - gas.N2mratio - gas.O2mratio - gas.CO2mratio

        gas.m2volratio()

        return gas

class Seprator():
    def __init__(self, sigma=0.99):
        self.sigma = sigma
    def work(self, gas:Gas, ratio):
        # ratio: outlet gas.G/main gas.G
        gas.P *= self.sigma
        Gtotal = gas.G
        gas2 = copy.deepcopy(gas)
        gas.G = Gtotal/(ratio+1)
        gas2.G = Gtotal * ratio/(ratio+1)
        return gas, gas2

class HeatExchenger():
    def __init__(self, heatEff=1.0, coe=0.85, sigma=0.95):
        self.heatEff = heatEff
        self.coe = coe
        self.sigma = sigma

    def workDesign(self, gas1:Gas, gas2:Gas, iteration=8):
        T1out, T2out = gas2.T, gas1.T
        for i in range(iteration):
            Cp1 = gas1.Cp(0.5*(gas1.T+T1out))
            Cp2 = gas2.Cp(0.5*(gas2.T+T2out))
            Qexc = min(gas1.G*Cp1*abs(gas1.T-gas2.T), gas2.G*Cp2*abs(gas1.T-gas2.T))
            Qexc = Qexc*self.coe*self.heatEff

            if gas1.T > gas2.T:
                T1out = gas1.T - Qexc/(gas1.G*Cp1)
                T2out = gas2.T + Qexc/(gas2.G*Cp2)
            else:
                T1out = gas1.T + Qexc/(gas1.G*Cp1)
                T2out = gas2.T - Qexc/(gas2.G*Cp2)
        gas1.T, gas2.T = T1out, T2out
        gas1.P, gas2.P = self.sigma*gas1.P, self.sigma*gas2.P
        return gas1, gas2, Qexc

class Gennerator():
    def __init__(self, root='./GenData'):
        Table = np.loadtxt(os.path.join(root, "gennerator_eff.csv"), delimiter=',')
        self.effTor, self.effN, self.effTable = Table[1:,0], Table[0,1:], Table[1:,1:]

        Table = np.loadtxt(os.path.join(root, "gennerator_volt.csv"), delimiter=',')
        self.voltTor, self.voltN, self.voltTable = Table[1:, 0], Table[0, 1:], Table[1:, 1:]

    def model(self, n, Tor):
        # n:prm, Torque:N*m
        eff,_ = interpolate_2title(Tor, n, self.effTor, self.effN, self.effTable)
        volt,_ = interpolate_2title(Tor, n, self.voltTor, self.voltN, self.voltTable)
        return volt, eff

    def loadR(self, U, R=1.0):
        return U/R

    def workR(self, n, R, lr=1.0, maxStep=3000):
        # n:prm  R: omi
        Tor1 = 3.0
        detaTor = 0.0001
        Tor2 = Tor1 + detaTor

        step = 1
        while True:
            U1, eff1 = self.model(n, Tor1)
            U2, eff2 = self.model(n, Tor2)

            Wgen1, Wgen2 = Tor1*n* np.pi/30 * eff1, Tor2*n* np.pi/30 * eff2
            Wload1, Wload2 = U1 * self.loadR(U1, R), U2 * self.loadR(U2, R)
            err1, err2 = Wgen1 - Wload1, Wgen2 - Wload2
            if (abs(err1) < 0.000001):
                break

            grad = (err2 - err1) / detaTor
            if (grad > 1e9):
                log("\nWarning: gradient > 1e9!\n")

            Tor1 = Tor1 - lr * err1 / grad
            Tor2 = Tor1 + detaTor

            if step > maxStep:
                log("\nstep %d grad: %f\n" % (step, grad))
                break
            step = step + 1

        return Tor1, U1, Wload1, Tor1*n*np.pi/30

    def loadFans(self, U, Radd=np.zeros(FANsNUM)):
        fans = ElecFans(Uratio=5.0, Ueff=0.99)
        return fans.work(U, Radd)

    def workFans(self, n, Radd=np.zeros(FANsNUM), lr=1.0, maxStep=3000):
        # input gen n, return work state params
        # n:prm  R: omi
        Tor1 = 3.0
        detaTor = 0.0001
        Tor2 = Tor1 + detaTor

        step = 1
        while True:
            U1, eff1 = self.model(n, Tor1)
            U2, eff2 = self.model(n, Tor2)

            Wgen1, Wgen2 = Tor1*n* np.pi/30 * eff1, Tor2*n* np.pi/30 * eff2
            Wload1, Wload2 = U1 * self.loadFans(U1, Radd), U2 * self.loadFans(U2, Radd)
            err1, err2 = Wgen1 - Wload1, Wgen2 - Wload2
            if (abs(err1) < 0.000001):
                break

            grad = (err2 - err1) / detaTor
            if (grad > 1e9):
                log("\nWarning: gradient > 1e9!\n")

            Tor1 = Tor1 - lr * err1 / grad
            Tor2 = Tor1 + detaTor

            if step > maxStep:
                log("\nstep %d grad: %f\n" % (step, grad))
                break
            step = step + 1

        fans = ElecFans(Uratio=5.0, Ueff=0.99)
        Tlist = fans.work(U1, Radd, ifThrust=True)
        # return: gen input Tor, gen output U, gen ouput work, gen input work, fans(10) output thrust list
        return Tor1, U1, Wload1, Tor1*n*np.pi/30, Tlist

    def designR(self, n, W):
        # design rpm, design input gennerator power, unit:W
        Tor = W/n/(np.pi/30)
        U, eff = self.model(n, Tor)
        R = U*U/(W*eff)
        return R, U

    def designFans(self, n, W, lr=1.0, maxStep=3000):
        # get the designe Radd of fans
        # design rpm, design input gennerator power, unit:W

        R1 = 0
        detaR = 0.0001
        R2 = R1 + detaR

        step = 1
        while True:
            _, U1, _, W1, Tlist = self.workFans(n, R1*np.ones(FANsNUM))
            _, _, _, W2, _ = self.workFans(n, R2*np.ones(FANsNUM))

            err1, err2 = W - W1, W - W2
            if (abs(err1) < 0.0001):
                break

            grad = (err2 - err1) / detaR
            if (grad > 1e9):
                log("\nWarning: gradient > 1e9!\n")

            R1 = R1 - lr * err1 / grad
            R2 = R1 + detaR

            if step > maxStep:
                log("\nstep %d grad: %f\n" % (step, grad))
                break
            step = step + 1

        return R1, U1, W1, Tlist.sum()

class ElecFans():
    def __init__(self, Uratio=5.0, Ueff=0.99):
        self.Uratio = Uratio   # Uin=290V, Uout=58V
        self.Ueff = Ueff

    def fanI(self, U):
        # Ulist = [51.8, 56, 58.8]  V  U=23.022,I=0
        # Ilist = [107, 123, 133]   A
        # U < 23, I<0
        return 3.7218 * U - 85.684

    def fanT(self, U):
        # Ulist = [51.8, 56, 58.8]  V  U=20.838,T=0
        # Thrust list = [78.48, 90.252, 96.138]   N
        return 2.5447 * U - 53.026

    def fanElec(self, U, Radd=0.1):
        U1 = 0
        U2 = U
        Usum1 = U1 + self.fanI(U1) * Radd   # U small
        Usum2 = U2 + self.fanI(U2) * Radd
        while True:
            Umid = 0.5*(U1+U2)
            Usummid = Umid + self.fanI(Umid)*Radd
            if abs(U - Usummid)<1e-6:
                break

            if U - Usummid < 0:
                U2 = Umid
            else:
                U1 = Umid

        return self.fanI(Umid), self.fanT(Umid)

    def work(self, U, Rlist=np.zeros(FANsNUM), ifThrust=False):
        Uin = U/self.Uratio
        Iin = 0
        for i in range(FANsNUM):
            Iin += (self.fanElec(Uin, Radd=Rlist[i]))[0]

        Win = Iin*Uin
        Wout = Win/self.Ueff
        if ifThrust:
            Tlist = np.zeros_like(Rlist)
            for i in range(FANsNUM):
                Tlist[i] = (self.fanElec(Uin, Radd=Rlist[i]))[1]
            return Tlist
        else:
            return Wout/U

class Model():
    def __init__(self):
        pass

    def model(self, Ts=288.15, Ps=101325.0, Ma=0.0, G=0.49, nl=96464.0, T4d=1017.0, isPrint:bool=False):
        Gd = 0.495
        nld = 94636
        leak_per=0.00
        nl = nl/nld
        detector = Detector()

        gas = Gas(Ts, Ps, Ma, G)
        detector.get(gas, "engine inlet", isPrint)

        intake = Intake(sigma=0.995)
        gas, inlet_drag = intake.work(gas, Ma)

        detector.get(gas, "Compressor inlet", isPrint)
        comp = Compressor(ifTable=True, path=os.path.join(os.getcwd(), "CharData2",), leakG_per=leak_per)
        gas, Wcomp, Gcdot, nlcdot, pic, effc = comp.workTable(gas, nl, Gd=Gd)
        gas_leak = copy.deepcopy(gas)
        gas_leak.G = gas.G/(1-leak_per)*leak_per

        duct1 = Duct(sigma=0.995)
        gas = duct1.work(gas)

        detector.get(gas, "Combustor inlet", isPrint)
        H2LHV = 119.988 * 1e6  # H2 lower heat value ，J/kg
        combu = Combustor(sigma=0.95, eff=0.98, LHV=H2LHV)
        if T4d<0:
            gas, Gfuel = combu.workT4(gas, Tout=gas.T)
        else:
            gas, Gfuel = combu.workT4(gas, Tout=T4d)

        detector.get(gas, "Mixer inlet", isPrint)
        mixer = Mixer(sigma=0.999)
        gas_leak.P = gas.P
        gas = mixer.work(gas, gas_leak)

        detector.get(gas, "Turbine inlet", isPrint)
        turb = Turbine(ifTable=True, path=os.path.join(os.getcwd(), "CharData2",), leakG_per=0.0)
        gas, Wturb, Gtdot, nltdot, pit, efft = turb.workTable(gas, nl, Gd=Gd)

        duct2 = Duct(sigma=0.998)
        gas = duct2.work(gas)

        G0, Pnoz = gas.G, gas.P
        detector.get(gas, "Nozzle inlet", isPrint)
        nozzle = Nozzle(type="con", Pamb=Ps, Cd=0.98, Cf=0.98, sigma=0.995)
        nozzle.A9 = 6.939778172e-3 *1  # outlet Area m2
        gas, gro_thrust, G1, list = nozzle.work(gas)
        detector.get(gas, "Engine outlet", isPrint)
        Ts_ad, Ps_ad, Ma_ad, V_ad, rhos_ad, A_ad = list

        err = G0-G1

        work = (0.99*Wturb-Wcomp)

        if isPrint:
            log(" ")
            log("Engine massflow:   %.5f kg/s"%(G))
            log("Engine Work:       %.5f kw"%(work/1000))
            log("Engine Efficiency: %.5f "%(work/(Gfuel*H2LHV)))
            log("Engine Fuleflow:   %.5f kg/s"%(Gfuel))
            log("Engine SFC:        %.5f kg/(kwh)"%(Gfuel*3600*1000/work))
            log("Engine nl:         %.5f rpm"%(nl*nld))
            log("Engine Thrust:     %.5f N"%(gro_thrust-inlet_drag))
            log("Engine outlet Ps:  %.5f Pa"%(Ps_ad))
            log("Engine outlet V:   %.5f m/s"%(V_ad))
            log(" ")
            log("Compressor Gdot: %.5f  nldot: %.5f"%(Gcdot, nlcdot))
            log("Compressor pi:   %.5f   eff:   %.5f"%(pic, effc))
            log("Compressor work: %.5f kw"%(Wcomp/1000))
            log("Turbine    Gdot: %.10f  nldot: %.10f"%(Gtdot, nltdot))
            log("Turbine    pi:   %.10f  eff:   %.10f"%(pit, efft))
            log("Turbine    work: %.5f kw"%(Wturb/1000))
            log(" ")

        return err, (work, Gfuel), (pic, pit, Pnoz)

    def solve(self, T0s=288.15, P0s=101325.0, Ma=0.0, nlratio=1.0, Gratio=1.0, T4d=1017.0, lr=0.1, maxStep=1000, isPrint:bool=False):
        # make the engine mass equal
        gas = Gas(T0s, P0s, Ma)
        intake = Intake(sigma=0.995)
        gas, _ = intake.work(gas, Ma)

        P0d = 101325  # design inlet tatal temperature/pressure
        T0d = 288.15

        Gd = 0.5006  # design mass flow rate, kg/s
        nld = 94636  # design speed spool init, rpm/min
        nl = nld * np.sqrt(gas.T / T0d) * nlratio  # speed spool init, rpm/min

        G1 = Gd * gas.P / P0d * np.sqrt(T0d / gas.T) * Gratio
        detaR = 0.0001
        G2 = (1-detaR)*G1  # to get the derivative of G
        lr = lr  # learning rate

        errLog, picLog, pitLog = np.zeros(maxStep), np.zeros(maxStep), np.zeros(maxStep)
        GLog = np.zeros(maxStep)
        PnozLog = np.zeros(maxStep)
        step = 0
        while True:
            err1, _, (pic, pit, Pnoz) = self.model(T0s, P0s, Ma, G1, nl, T4d, isPrint=False)
            if (abs(err1) < 0.000001):
                break

            err2, _, _ = self.model(T0s, P0s, Ma, G2, nl, T4d, isPrint=False)

            grad = (err1 - err2) / (G1 - G2)
            if (grad > 1e9):
                log("\nWarning: gradient > 1e9!\n")

            errLog[step],picLog[step],pitLog[step] = err1, pic, pit
            GLog[step], PnozLog[step] = G1, Pnoz
            # Narrow the learning rate to prevent fluctuations
            if step>3 and abs(err1)>abs(errLog[step-1]):
                lr *= 0.5
                G1 = GLog[step-2]

            G1 = G1 - lr * err1 / grad
            G2 = (1 - detaR) * G1
            step = step + 1

            if step >= maxStep:
                log("\nstep %d grad: %f err1: %f min err: %f\n lr: %f" % (step, grad, err1, errLog.min(), lr))
                import matplotlib.pyplot as plt
                fig = plt.figure("state")
                line_color = 'b-'
                axis_font = {'size': '14'}
                fig.set_size_inches(18, 9)

                axes = plt.subplot(2, 3, 1)
                axes.plot(errLog, line_color)
                axes.set_xlabel("step", axis_font)
                axes.set_ylabel("errLog", axis_font)
                set_axes(axes)

                axes = plt.subplot(2, 3, 2)
                axes.plot(picLog, line_color)
                axes.set_xlabel("step", axis_font)
                axes.set_ylabel("picLog", axis_font)
                set_axes(axes)

                axes = plt.subplot(2, 3, 3)
                axes.plot(pitLog, line_color)
                axes.set_xlabel("step", axis_font)
                axes.set_ylabel("pitLog", axis_font)
                set_axes(axes)

                axes = plt.subplot(2, 3, 4)
                axes.plot(GLog, line_color)
                axes.set_xlabel("step", axis_font)
                axes.set_ylabel("GLog", axis_font)
                set_axes(axes)

                axes = plt.subplot(2, 3, 5)
                axes.plot(PnozLog, line_color)
                axes.set_xlabel("step", axis_font)
                axes.set_ylabel("Pout", axis_font)
                set_axes(axes)

                plt.tight_layout()
                plt.show()
                break

        _, (W, Gfuel), _ = self.model(T0s, P0s, Ma, G1, nl, T4d, isPrint=isPrint)
        Gratio = G1/(Gd * gas.P / P0d * np.sqrt(T0d / gas.T))
        if isPrint:
            log("Engine Gratio: %.5f  nlratio: %.5f" % (Gratio, nlratio))
            log(" ")
        return W, Gfuel, (G1, nl, Gratio)

    def workElec(self, nlratio, T0s=288.15, P0s=101325.0, Ma=0.0, R=4.5, lr=1.0, maxStep=1000, clip=50.0, isPrint:bool=False):
        gen = Gennerator()
        T4d1 = 1200.0
        detaT4 = 0.001
        T4d2 = T4d1 + detaT4
        Gratio = 1.5 * nlratio - 0.7  # T4d=1130  nl:0.47 - 1.01
        # adjust lr, when input lr=1.0, nlratio:1--0.1, 0.7--0.3, clip
        if nlratio<=1.0:
            lr *= -2/3 * nlratio + 0.77
        elif nlratio<=0.7:
            lr *= 0.3
        else:
            lr = 0.1

        W1, _, res = self.solve(T0s=T0s, P0s=P0s, Ma=Ma, nlratio=nlratio, Gratio=Gratio, T4d=T4d1, lr=lr, isPrint=False)
        _, nl, Gratio = res
        _, U, _, Wgen = gen.workR(nl, R)

        step = 1
        while True:
            try:
                W1, _, res = self.solve(T0s=T0s, P0s=Ps, Ma=Ma, nlratio=nlratio, Gratio=0.9*Gratio, T4d=T4d1, lr=lr, isPrint=False)
                _, _, Gratio = res
            except:
                print("error")
            err1 = W1 - Wgen
            if (abs(err1) < 0.000001):
                break

            W2, _, _ = self.solve(T0s=T0s, P0s=Ps, Ma=Ma, nlratio=nlratio, Gratio=0.9*Gratio, T4d=T4d2, lr=lr, isPrint=False)
            err2 = W2 - Wgen

            grad = (err2 - err1) / detaT4
            if (grad > 1e9):
                log("\nWarning: gradient > 1e9!\n")

            plus = - 1.0 * err1 / grad
            if plus>clip:
                plus = clip
            elif plus <-clip:
                plus = -clip
            T4d1 = T4d1 + plus
            T4d2 = T4d1 + detaT4

            if step > maxStep:
                log("\nstep %d grad: %f\n" % (step, grad))
                break
            step = step + 1
        Weng, Gfuel, (G, nl, Gratio) \
            = self.solve(T0s=T0s, P0s=Ps, Ma=Ma, nlratio=nlratio, Gratio=Gratio, T4d=T4d1, lr=lr, isPrint=isPrint)
        return (G, Gfuel, nl, T4d1), (Weng, U), Gratio

    def stablizeT4(self, nlstart=0.5, nlend=0.9, R=4.5, maxStep=300, isPlot=True):
        Ts, Ps, _ = Altitude(0)
        model = Model()
        log("########## start state ##########")
        (G, Gfuel, nl, T4), (Weng, U), Gratio \
            = model.workElec(nlstart, T0s=Ts, P0s=Ps, Ma=0.0, R=R, lr=0.2, maxStep=300, clip=200, isPrint=True)
        # nlratio=0.5, R=4.60 get the good condition

        gen = Gennerator()
        _, U, _, Wgen = gen.workR(nl, R)

        timeArr = np.zeros(maxStep)
        GArr, WArr, GfuelArr = np.zeros(maxStep), np.zeros(maxStep), np.zeros(maxStep)
        nlArr, UArr, WgenArr = np.zeros(maxStep), np.zeros(maxStep), np.zeros(maxStep)
        T4Arr = np.zeros(maxStep)
        err = np.zeros(maxStep)

        GArr[0], WArr[0], GfuelArr[0] = G, Weng, Gfuel
        nlArr[0], UArr[0], WgenArr[0],T4Arr[0] = nl, U, Wgen, T4
        err[0] = nlend-nlstart

        kp, ki, kd = 200, 0, 0.0

        detaTime = 0.01  # deta time, s
        J = 0.00119002  # rotational inertia kg/(m2)
        step, time = 1, 0.0
        while True:
            nlratio = nl /(94636*np.sqrt(Ts/288.15))
            err[step] = nlend - nlratio
            T4 = T4 + kp*(err[step] + ki*detaTime*err.sum() + kd*(err[step]-err[step-1])/detaTime)
            if T4>2500:
                T4 = 2500
            elif T4 < 500:
                T4 = 500

            Weng, Gfuel, (G, _, Gratio) = model.solve(Ts, Ps, 0.0, nlratio, 0.94*Gratio, T4d=T4, lr=0.1, isPrint=False)

            _, U, _, Wgen = gen.workR(nl, R)  # gennerator intake work in current nl

            log("step:%d  time:%.3fs  nl:%.1frpm  nlratio:%.3f  Gratio:%.3f  G:%.3fkg/s  Weng:%.3fw  Wgen:%.3fw  Gfuel:%.3fg/s  T4:%.2fK"
                  %(step, time, nl, nlratio, Gratio, G, Weng, Wgen, Gfuel*1000, T4))

            nl += detaTime * ((Weng - Wgen) / (J * nl) / ((np.pi / 30) ** 2))
            GArr[step], WArr[step], GfuelArr[step] = G, Weng, Gfuel
            nlArr[step], UArr[step], WgenArr[step], T4Arr[step] = nl, U, Wgen, T4

            time+=detaTime
            timeArr[step] = time
            step+=1
            if step>=maxStep:
                break

        if isPlot:
            import matplotlib.pyplot as plt
            line_color = 'b-'
            line_color2 = 'r-'
            axis_font = {'size': '14'}

            fig = plt.figure("state")
            fig.set_size_inches(15, 9)

            axes = plt.subplot(2, 3, 1)
            axes.plot(timeArr, nlArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("nl prm", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 2)
            axes.plot(timeArr, GArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("massflow kg/s", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 3)
            axes.plot(timeArr, WArr/1000, line_color)
            axes.plot(timeArr, WgenArr/1000, line_color2)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("W kw", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 4)
            axes.plot(timeArr, T4Arr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("T4 K", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 5)
            axes.plot(timeArr, GfuelArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("Gfuel kg/s", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 6)
            axes.plot(timeArr, UArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("voltage V", axis_font)
            set_axes(axes)

            plt.tight_layout()

            plt.figure("Work")
            plt.plot(timeArr, WArr/1000, lw="1.2", color="m", label="Engine work")
            plt.plot(timeArr, WgenArr/1000, lw="1.2", color="b", label="Gennerator work")
            plt.xlabel("Time/s", fontdict=axis_font)
            plt.ylabel("Work kw", fontdict=axis_font)
            plt.grid("on")
            plt.legend()
            plt.tight_layout()

            plt.show()

        return timeArr, GArr, WArr, GfuelArr, nlArr, UArr, WgenArr

    def stablizeT4constW(self, nlstart=0.5, nlend=0.9, R=4.5, T4d=1130, maxStep=300, isPlot=True):
        Ts, Ps, _ = Altitude(0)
        model = Model()
        log("\n########## start state ##########")
        Gratio = 1.5 * nlstart - 0.7  # T4d=1130  nl:0.47 - 1.01
        Weng, Gfuel, (G, nl, Gratio) = model.solve(T0s=Ts, P0s=Ps, Ma=0.0, nlratio=nlstart, Gratio=Gratio, T4d=T4d,
                                                   isPrint=True, lr=0.1, maxStep=1000)
        WgenEnd, _, _ = model.solve(T0s=Ts, P0s=Ps, Ma=0.0, nlratio=nlend, Gratio=1.5*nlend-0.7, T4d=T4d,
                                                   isPrint=True, lr=0.1, maxStep=1000)

        timeArr = np.zeros(maxStep)
        GArr, WArr, GfuelArr = np.zeros(maxStep), np.zeros(maxStep), np.zeros(maxStep)
        nlArr, WgenArr = np.zeros(maxStep), np.zeros(maxStep)
        T4Arr = np.zeros(maxStep)
        err = np.zeros(maxStep)

        GArr[0], WArr[0], GfuelArr[0] = G, Weng, Gfuel
        nlArr[0], WgenArr[0],T4Arr[0] = nl, Weng, T4d
        err[0] = nlend-nlstart

        kp, ki, kd = 1000, 0.0, 0.5

        detaTime = 0.01  # deta time, s
        J = 0.00119002  # rotational inertia kg/(m2)
        step, time = 1, 0.0
        T4 = T4d
        startStep = 50
        endStep = 0.7 * maxStep - startStep
        WgenStart = Weng
        while True:
            nlratio = nl /(94636*np.sqrt(Ts/288.15))
            err[step] = nlend - nlratio
            if step > startStep and step<endStep:
                # Wgen = WgenEnd
                Wgen = (WgenEnd-WgenStart)/(endStep-startStep) * (step-startStep) + WgenStart
                T4 = T4 + kp*(err[step] + ki*detaTime*err.sum() + kd*(err[step]-err[step-1])/detaTime)
                if T4>1700:
                    T4 = 1700
                elif T4 < 800:
                    T4 = 800
            elif step > endStep:
                Wgen = WgenEnd
                T4 = T4 + kp * (err[step] + ki * detaTime * err.sum() + kd * (err[step] - err[step - 1]) / detaTime)
                if T4 > 1700:
                    T4 = 1700
                elif T4 < 800:
                    T4 = 800
            else:
                Wgen = Weng

            Weng, Gfuel, (G, _, Gratio) = model.solve(Ts, Ps, 0.0, nlratio, 0.94*Gratio, T4d=T4, lr=0.1, isPrint=False)

            log("step:%d  time:%.3fs  nl:%.1frpm  nlratio:%.3f  Gratio:%.3f  G:%.3fkg/s  Weng:%.3fw  Wgen:%.3fw  Gfuel:%.3fg/s  T4:%.2fK"
                  %(step, time, nl, nlratio, Gratio, G, Weng, Wgen, Gfuel*1000, T4))

            nl += detaTime * ((Weng - Wgen) / (J * nl) / ((np.pi / 30) ** 2))
            GArr[step], WArr[step], GfuelArr[step] = G, Weng, Gfuel
            nlArr[step], T4Arr[step] = nl, T4
            WgenArr[step] = Wgen

            time+=detaTime
            timeArr[step] = time
            step+=1
            if step>=maxStep:
                break

        if isPlot:
            import matplotlib.pyplot as plt
            line_color = 'b-'
            line_color2 = 'r-'
            axis_font = {'size': '14'}

            fig = plt.figure("state")
            fig.set_size_inches(15, 9)

            axes = plt.subplot(2, 3, 1)
            axes.plot(timeArr, nlArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("nl prm", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 2)
            axes.plot(timeArr, GArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("massflow kg/s", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 3)
            axes.plot(timeArr, WArr/1000, line_color)
            axes.plot(timeArr, WgenArr/1000, line_color2)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("W kw", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 4)
            axes.plot(timeArr, T4Arr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("T4 K", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 5)
            axes.plot(timeArr, GfuelArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("Gfuel kg/s", axis_font)
            set_axes(axes)

            plt.tight_layout()

            plt.figure("Work")
            plt.plot(timeArr, WArr/1000, lw="1.2", color="m", label="Engine work")
            plt.plot(timeArr, WgenArr/1000, lw="1.2", color="b", label="Gennerator work")
            plt.xlabel("Time/s", fontdict=axis_font)
            plt.ylabel("Work kw", fontdict=axis_font)
            plt.grid("on")
            plt.legend()
            plt.tight_layout()

            plt.show()

        return timeArr, GArr, WArr, GfuelArr, nlArr, WgenArr

    def stablizeR(self, nlstart=0.5, nlend=0.9, T4d=1130, maxStep=300, isPlot=True):
        Ts, Ps, _ = Altitude(0)
        model = Model()
        log("\n########## start state ##########")
        Gratio = 1.5 * nlstart - 0.7  # T4d=1130  nl:0.47 - 1.01
        Weng, Gfuel, (G, nl, Gratio) = model.solve(T0s=Ts, P0s=Ps, Ma=0.0, nlratio=nlstart, Gratio=Gratio, T4d=T4d,
                                                 isPrint=True, lr=0.1,maxStep=1000)

        gen = Gennerator()
        R, U = gen.designR(nl, Weng)
        _, U, _, Wgen = gen.workR(nl, R)

        timeArr = np.zeros(maxStep)
        GArr, WArr, GfuelArr = np.zeros(maxStep), np.zeros(maxStep), np.zeros(maxStep)
        nlArr, UArr, WgenArr = np.zeros(maxStep), np.zeros(maxStep), np.zeros(maxStep)
        RArr = np.zeros(maxStep)
        err = np.zeros(maxStep)

        GArr[0], WArr[0], GfuelArr[0] = G, Weng, Gfuel
        nlArr[0], UArr[0], WgenArr[0],RArr[0] = nlstart, U, Wgen, R
        err[0] = nlend-nlstart

        kp, ki, kd = 300, 0.1, 0

        detaTime = 0.01  # deta time, s
        J = 0.00119002  # rotational inertia kg/(m2)
        step, time = 1, 0.0
        while True:
            nlratio = nl /(94636*np.sqrt(Ts/288.15))
            err[step] = nlend - nlratio
            R = kp*(err[step] + ki*detaTime*err.sum() + kd*(err[step]-err[step-1])/detaTime)
            if R>5000:
                R = 5000
            elif R < 0.2:
                R = 0.2

            Weng, Gfuel, (G, _, Gratio) = model.solve(Ts, Ps, 0.0, nlratio, 0.94*Gratio, T4d=T4d, lr=0.1, isPrint=False)

            _, U, _, Wgen = gen.workR(nl, R)  # gennerator intake work in current nl

            log("step:%d  time:%.3fs  nl:%.1frpm  nlratio:%.3f  Gratio:%.3f  G:%.3fkg/s  Weng:%.3fw  Wgen:%.3fw  Gfuel:%.3fg/s  R:%.2fomi"
                  %(step, time, nl, nlratio, Gratio, G, Weng, Wgen, Gfuel*1000, R))

            nl += detaTime * ((Weng - Wgen) / (J * nl) / ((np.pi / 30) ** 2))
            GArr[step], WArr[step], GfuelArr[step] = G, Weng, Gfuel
            nlArr[step], UArr[step], WgenArr[step], RArr[step] = nlratio, U, Wgen, R

            time+=detaTime
            timeArr[step] = time
            step+=1
            if step>=maxStep:
                break

        log("\n########## end state ##########")
        model.solve(Ts, Ps, 0.0, nlratio, 0.94 * Gratio, T4d=T4d, lr=0.1, isPrint=True)
        if isPlot:
            import matplotlib.pyplot as plt
            line_color = 'b-'
            line_color2 = 'r-'
            axis_font = {'size': '14'}

            fig = plt.figure("state")
            fig.set_size_inches(15, 9)

            axes = plt.subplot(2, 3, 1)
            axes.plot(timeArr, nlArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("nlratio prm", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 2)
            axes.plot(timeArr, GArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("massflow kg/s", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 3)
            axes.plot(timeArr, WArr/1000, line_color)
            axes.plot(timeArr, WgenArr/1000, line_color2)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("W kw", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 4)
            axes.plot(timeArr, RArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("R omi", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 5)
            axes.plot(timeArr, GfuelArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("Gfuel kg/s", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 6)
            axes.plot(timeArr, UArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("voltage V", axis_font)
            set_axes(axes)

            plt.tight_layout()

            plt.figure("Work")
            plt.plot(timeArr, WArr/1000, lw="1.2", color="m", label="Engine work")
            plt.plot(timeArr, WgenArr/1000, lw="1.2", color="b", label="Gennerator work")
            plt.xlabel("Time/s", fontdict=axis_font)
            plt.ylabel("Work kw", fontdict=axis_font)
            plt.grid("on")
            plt.legend()
            plt.tight_layout()

            plt.show()

        return timeArr, GArr, WArr, GfuelArr, nlArr, UArr, WgenArr

    def stablizeFans(self, nlstart=0.5, nlend=0.9, T4d=1130, maxStep=300, isPlot=True):
        Ts, Ps, _ = Altitude(0)
        model = Model()
        log("\n########## start state ##########")
        Gratio = 1.5 * nlstart - 0.7  # T4d=1130  nl:0.47 - 1.01
        Weng, Gfuel, (G, nl, Gratio) = model.solve(T0s=Ts, P0s=Ps, Ma=0.0, nlratio=nlstart, Gratio=Gratio, T4d=T4d,
                                                 isPrint=True, lr=0.1,maxStep=1000)

        gen = Gennerator()
        R, U, Wgen, thrust = gen.designFans(nl, Weng)

        timeArr = np.zeros(maxStep)
        GArr, WArr, GfuelArr = np.zeros(maxStep), np.zeros(maxStep), np.zeros(maxStep)
        nlArr, UArr, WgenArr = np.zeros(maxStep), np.zeros(maxStep), np.zeros(maxStep)
        RArr = np.zeros(maxStep)
        Thrust = np.zeros(maxStep)
        err = np.zeros(maxStep)

        GArr[0], WArr[0], GfuelArr[0] = G, Weng, Gfuel
        nlArr[0], UArr[0], WgenArr[0],RArr[0] = nlstart, U, Wgen, R
        Thrust[0] = thrust
        err[0] = nlend-nlstart

        kp, ki, kd = 10, 0.0, 0

        detaTime = 0.01  # deta time, s
        J = 0.00119002  # rotational inertia kg/(m2)
        step, time = 1, 0.0
        while True:
            nlratio = nl /(94636*np.sqrt(Ts/288.15))
            err[step] = nlend - nlratio
            if step>100:
                R = kp*(err[step] + ki*detaTime*err.sum() + kd*(err[step]-err[step-1])/detaTime)
                if R>5000:
                    R = 5000
                elif R < 0.0:
                    R = 0.0

            Weng, Gfuel, (G, _, Gratio) = model.solve(Ts, Ps, 0.0, nlratio, 0.94*Gratio, T4d=T4d, lr=0.1, isPrint=False)

            _, U, _, Wgen, Tlist = gen.workFans(nl, R*np.ones(FANsNUM))  # gennerator intake work in current nl

            log("step:%d  time:%.3fs  nl:%.1frpm  nlratio:%.3f  Gratio:%.3f  G:%.3fkg/s  Weng:%.3fw  Wgen:%.3fw  Gfuel:%.3fg/s  R:%.2fomi"
                  %(step, time, nl, nlratio, Gratio, G, Weng, Wgen, Gfuel*1000, R))

            nl += detaTime * ((Weng - Wgen) / (J * nl) / ((np.pi / 30) ** 2))
            GArr[step], WArr[step], GfuelArr[step] = G, Weng, Gfuel
            nlArr[step], UArr[step], WgenArr[step], RArr[step] = nlratio, U, Wgen, R
            Thrust[step] = Tlist.sum()

            time+=detaTime
            timeArr[step] = time
            step+=1
            if step>=maxStep:
                break

        log("\n########## end state ##########")
        model.solve(Ts, Ps, 0.0, nlratio, 0.94 * Gratio, T4d=T4d, lr=0.1, isPrint=True)
        if isPlot:
            import matplotlib.pyplot as plt
            line_color = 'b-'
            line_color2 = 'r-'
            axis_font = {'size': '14'}

            fig = plt.figure("state")
            fig.set_size_inches(15, 9)

            axes = plt.subplot(2, 3, 1)
            axes.plot(timeArr, nlArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("nlratio prm", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 2)
            axes.plot(timeArr, GArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("massflow kg/s", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 3)
            axes.plot(timeArr, WArr/1000, line_color)
            axes.plot(timeArr, WgenArr/1000, line_color2)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("W kw", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 4)
            axes.plot(timeArr, RArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("R omi", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 5)
            axes.plot(timeArr, GfuelArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("Gfuel kg/s", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 6)
            axes.plot(timeArr, UArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("voltage V", axis_font)
            set_axes(axes)

            plt.tight_layout()

            Signal =WArr[-1]* np.ones_like(WArr)
            Signal[:100] = WArr[0]
            plt.figure("Work")
            plt.plot(timeArr, WArr/1000, lw="1.2", color="m", label="Engine work")
            plt.plot(timeArr, WgenArr/1000, lw="1.2", color="b", label="Gennerator work")
            plt.plot(timeArr, Signal/1000, lw="1.2", color="k", label="Signal")
            plt.xlabel("Time/s", fontdict=axis_font)
            plt.ylabel("Work kw", fontdict=axis_font)
            plt.grid("on")
            plt.legend()
            plt.tight_layout()

            plt.figure("Thrust")
            plt.plot(timeArr, Thrust/9.8, lw="1.2", color="b", label="Fans thrust")
            plt.xlabel("Time/s", fontdict=axis_font)
            plt.ylabel("Thrust kg", fontdict=axis_font)
            plt.grid("on")
            plt.tight_layout()

            plt.show()

        return timeArr, GArr, WArr, GfuelArr, nlArr, UArr, WgenArr, Thrust

    def stablizeFansConst(self, nlstart=0.5, nlend=0.9, T4d=1130, maxStep=300, isPlot=True):
        Ts, Ps, _ = Altitude(0)
        model = Model()

        Gratio = 1.5 * nlend - 0.7  # T4d=1130  nl:0.47 - 1.01
        Weng, _, (_, nl, _) = model.solve(T0s=Ts, P0s=Ps, Ma=0.0, nlratio=nlend, Gratio=Gratio, T4d=T4d,
                                                   isPrint=True, lr=0.1, maxStep=1000)

        gen = Gennerator()
        Rend, _, _, _ = gen.designFans(nl, Weng)

        log("\n########## start state ##########")
        Gratio = 1.5 * nlstart - 0.7  # T4d=1130  nl:0.47 - 1.01
        Weng, Gfuel, (G, nl, Gratio) = model.solve(T0s=Ts, P0s=Ps, Ma=0.0, nlratio=nlstart, Gratio=Gratio, T4d=T4d,
                                                 isPrint=True, lr=0.1,maxStep=1000)

        gen = Gennerator()
        R, U, Wgen, thrust = gen.designFans(nl, Weng)

        timeArr = np.zeros(maxStep)
        GArr, WArr, GfuelArr = np.zeros(maxStep), np.zeros(maxStep), np.zeros(maxStep)
        nlArr, UArr, WgenArr = np.zeros(maxStep), np.zeros(maxStep), np.zeros(maxStep)
        RArr = np.zeros(maxStep)
        Thrust = np.zeros(maxStep)
        err = np.zeros(maxStep)

        GArr[0], WArr[0], GfuelArr[0] = G, Weng, Gfuel
        nlArr[0], UArr[0], WgenArr[0],RArr[0] = nlstart, U, Wgen, R
        Thrust[0] = thrust
        err[0] = nlend-nlstart

        detaTime = 0.01  # deta time, s
        J = 0.00119002  # rotational inertia kg/(m2)
        step, time = 1, 0.0
        while True:
            nlratio = nl /(94636*np.sqrt(Ts/288.15))
            err[step] = nlend - nlratio

            if step>=100:
                R = Rend

            Weng, Gfuel, (G, _, Gratio) = model.solve(Ts, Ps, 0.0, nlratio, 0.94*Gratio, T4d=T4d, lr=0.1, isPrint=False)

            _, U, _, Wgen, Tlist = gen.workFans(nl, R*np.ones(FANsNUM))  # gennerator intake work in current nl

            log("step:%d  time:%.3fs  nl:%.1frpm  nlratio:%.3f  Gratio:%.3f  G:%.3fkg/s  Weng:%.3fw  Wgen:%.3fw  Gfuel:%.3fg/s  R:%.2fomi"
                  %(step, time, nl, nlratio, Gratio, G, Weng, Wgen, Gfuel*1000, R))

            nl += detaTime * ((Weng - Wgen) / (J * nl) / ((np.pi / 30) ** 2))
            GArr[step], WArr[step], GfuelArr[step] = G, Weng, Gfuel
            nlArr[step], UArr[step], WgenArr[step], RArr[step] = nlratio, U, Wgen, R
            Thrust[step] = Tlist.sum()

            time+=detaTime
            timeArr[step] = time
            step+=1
            if step>=maxStep:
                break

        log("\n########## end state ##########")
        model.solve(Ts, Ps, 0.0, nlratio, 0.94 * Gratio, T4d=T4d, lr=0.1, isPrint=True)
        if isPlot:
            import matplotlib.pyplot as plt
            line_color = 'b-'
            line_color2 = 'r-'
            axis_font = {'size': '14'}

            fig = plt.figure("state")
            fig.set_size_inches(15, 9)

            axes = plt.subplot(2, 3, 1)
            axes.plot(timeArr, nlArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("nlratio prm", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 2)
            axes.plot(timeArr, GArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("massflow kg/s", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 3)
            axes.plot(timeArr, WArr/1000, line_color)
            axes.plot(timeArr, WgenArr/1000, line_color2)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("W kw", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 4)
            axes.plot(timeArr, RArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("R omi", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 5)
            axes.plot(timeArr, GfuelArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("Gfuel kg/s", axis_font)
            set_axes(axes)

            axes = plt.subplot(2, 3, 6)
            axes.plot(timeArr, UArr, line_color)
            axes.set_xlabel("time s", axis_font)
            axes.set_ylabel("voltage V", axis_font)
            set_axes(axes)

            plt.tight_layout()

            plt.figure("Work")
            plt.plot(timeArr, WArr/1000, lw="1.2", color="m", label="Engine work")
            plt.plot(timeArr, WgenArr/1000, lw="1.2", color="b", label="Gennerator work")
            plt.xlabel("Time/s", fontdict=axis_font)
            plt.ylabel("Work kw", fontdict=axis_font)
            plt.grid("on")
            plt.legend()
            plt.tight_layout()

            plt.figure("Thrust")
            plt.plot(timeArr, Thrust/9.8, lw="1.2", color="b", label="Fans thrust")
            plt.xlabel("Time/s", fontdict=axis_font)
            plt.ylabel("Thrust kg", fontdict=axis_font)
            plt.grid("on")
            plt.tight_layout()

            plt.show()

        return timeArr, GArr, WArr, GfuelArr, nlArr, UArr, WgenArr, Thrust

def set_axes(axes):
    # This sets the axis parameters for all plots
    axes.minorticks_on()
    axes.grid(which='major', linestyle='-', linewidth=0.5, color='grey')
    axes.grid(which='minor', linestyle=':', linewidth=0.5, color='grey')
    axes.grid(True)
    axes.get_yaxis().get_major_formatter().set_scientific(False)
    axes.get_yaxis().get_major_formatter().set_useOffset(False)
    return

def checkTable(path="./CHarData2", isFan=True, isComp=True, isTurb=True ):
    if isFan:
        fan = Fan05()


if __name__ == '__main__':
    # Nozzle
    if False:
        gas = Gas()
        gas.G, gas.T, gas.P, = 13.650,    998.14,   220477
        nozzle = Nozzle(type="con", Pamb=101325, Cf=0.98, Cd=1.0, sigma=1)
        gas, gro_thrust, A8, A9, list = nozzle.workDesign(gas)
        Ts_out, Ps_out, Ma_out, V_out = list

        gas25 = Gas(Height=0, Ma=0.2)
        Ts, Ps, V, _, _ = gas25.static(Ma=0.2)
        print(gro_thrust-5.121*V, A8, A9)

    # Nozzle with Intake
    if False:
        Pamb = 22632
        gas0 = Gas()
        gas0.G, gas0.T, gas0.P, = 3.794,    223.60,    25272

        intake = Intake(sigma=1)
        gas1, intake_drag = intake.work(gas0, Ma=0.4)
        print(gas1.T, gas1.P)

        gas = Gas()
        gas.G, gas.T, gas.P, = 3.876,   1038.25,    66912
        nozzle = Nozzle(type="con", Pamb=Pamb, Cf=0.98, Cd=0.985, sigma=1)
        gas, gro_thrust, A8, A9, list = nozzle.workDesign(gas)
        Ts_out, Ps_out, Ma_out, V_out = list
        thrust = (gro_thrust - intake_drag)/10

        print("Engine Thrust:       %.5f daN" % (thrust))
        print("Engine Thrust:       %.5f kN" % (thrust / 100))
        print("Engine Thrust:       %.5f kgf" % (thrust * 10 / 9.8))
        print("Engine gro Thrust:   %.5f kN" % (gro_thrust / 1000))
        print("Engine inlet Thrust: %.5f kN" % (intake_drag / 1000))
        print("Engine outlet Pi:    %.5f " % (gas.P / Pamb))
        print("Engine outlet Ma:    %.5f " % (Ma_out))
        print("Engine outlet Vel:   %.5f m/s" % (V_out))
        print("Engine outlet Ps:    %.5f kPa" % (Ps_out / 1000))
        print("Engine A8:           %.5f m2" % (A8))
        print("Engine A9:           %.5f m2" % (A9))
        print("Engine A9 diam:      %.5f m" % (np.sqrt(A9 / np.pi)))

    # gas chemical
    if False:
        T = 300
        gas1 = Gas()
        print(gas1.Cp(T=T))
        gas2 = copy.deepcopy(gas1)
        gas2.O2ratio, gas2.N2ratio = 0.5, 0.5
        gas2.vol2mratio()
        print(gas2.Cp(T=T))
        print(gas2.Rg())

    # Fans simulation
    if False:
        T,P,den = Altitude(25000)
        fan = Fan05()
        gas, thrust = fan.workDesign(Height=20000, Ma=0, G=1.0, pi=1.6, eff=0.78, Pamb=101325)
        print(thrust)

    # Fans simulation
    if False:
        Ts, Ps, _ = Altitude(25000)
        Ma = 0.6

        gas = Gas(Ts, Ps, Ma, G=1.06)

        intake = Intake(sigma=0.99)
        gas, inlet_drag = intake.work(gas, Ma)

        fan = Fan()
        gas2, W = fan.workDesign(gas, pi=2.3, eff=0.75)

        noz = Nozzle(type="con", Pamb=Ps, Cd=0.98, Cf=0.98, sigma=1.0)
        gas, gro_thrust, A8, A9, _ = noz.workDesign(gas2)
        # 0.3395919805107334 0.3395919805107334

        print((gro_thrust - inlet_drag) / 10, W / 1000, A8, A9, (gro_thrust - inlet_drag) / W * 1000)
        print(gro_thrust, inlet_drag)

        noz.A8 = 0.3395919805107334
        _, gro_thrust, G, _ = noz.work(gas2)

        print(" ")
        print((gro_thrust-inlet_drag)/10, W/1000, (gro_thrust-inlet_drag)/W*1000)
        print(gro_thrust, inlet_drag)

    if False:
        turb = Turbine()
        nlLine, massTable, piTable, effTable = turb.getTable()
        nl, mass = 0.95, 0.41645
        Rg, P, T = 287.06,  391079,  1017.89
        pi, eff, flag = interpolate2D_charline(nl/np.sqrt(Rg*T), mass*np.sqrt(Rg*T)/P, nlLine, massTable, piTable, effTable)
        print(pi, eff, flag)

        import matplotlib.pyplot as plt
        num = 100
        mass = np.linspace(0.25, 0.51, num)
        pi50 = np.zeros_like(mass)
        pi63 = np.zeros_like(mass)
        pi66 = np.zeros_like(mass)
        pi75 = np.zeros_like(mass)
        pi825 = np.zeros_like(mass)
        pi975 = np.zeros_like(mass)
        pi110 = np.zeros_like(mass)

        eff50 = np.zeros_like(mass)
        eff63 = np.zeros_like(mass)
        eff66 = np.zeros_like(mass)
        eff75 = np.zeros_like(mass)
        eff825 = np.zeros_like(mass)
        eff975 = np.zeros_like(mass)
        eff110 = np.zeros_like(mass)

        for i in range(num):
            pi50[i], eff50[i], _ = interpolate2D_charline(0.5/np.sqrt(Rg*T), mass[i]*np.sqrt(Rg*T)/P, nlLine, massTable, piTable, effTable)
            pi63[i], eff63[i], _ = interpolate2D_charline(0.63/np.sqrt(Rg*T), mass[i]*np.sqrt(Rg*T)/P, nlLine, massTable, piTable, effTable)
            pi66[i], eff66[i], _ = interpolate2D_charline(0.66/np.sqrt(Rg*T), mass[i]*np.sqrt(Rg*T)/P, nlLine, massTable, piTable, effTable)
            pi75[i], eff75[i], _ = interpolate2D_charline(0.75/np.sqrt(Rg*T), mass[i]*np.sqrt(Rg*T)/P, nlLine, massTable, piTable, effTable)
            pi825[i], eff825[i], _ = interpolate2D_charline(0.825/np.sqrt(Rg*T), mass[i]*np.sqrt(Rg*T)/P, nlLine, massTable, piTable, effTable)
            pi975[i], eff975[i], _ = interpolate2D_charline(0.975/np.sqrt(Rg*T), mass[i]*np.sqrt(Rg*T)/P, nlLine, massTable, piTable, effTable)
            pi110[i], eff110[i], _ = interpolate2D_charline(1.1/np.sqrt(Rg*T), mass[i]*np.sqrt(Rg*T)/P, nlLine, massTable, piTable, effTable)

        plt.figure("pi")
        plt.plot(massTable[0, :]/(np.sqrt(Rg*T)/P), piTable[0, :], lw="1.6", color="m", ls="-", label="n=0.6")
        plt.plot(massTable[1, :]/(np.sqrt(Rg*T)/P), piTable[1, :], lw="1.6", color="b", ls="-", label="n=0.7")
        plt.plot(massTable[2, :]/(np.sqrt(Rg*T)/P), piTable[2, :], lw="1.6", color="g", ls="-", label="n=0.8")
        plt.plot(massTable[3, :]/(np.sqrt(Rg*T)/P), piTable[3, :], lw="1.6", color="k", ls="-", label="n=0.9")
        plt.plot(massTable[4, :]/(np.sqrt(Rg*T)/P), piTable[4, :], lw="1.6", color="r", ls="-", label="n=1.0")

        plt.plot(mass, pi50, lw="1.2", color="m", ls="--", label="n=0.5")
        plt.plot(mass, pi63, lw="1.2", color="b", ls="--", label="n=0.63")
        plt.plot(mass, pi66, lw="1.2", color="g", ls="--", label="n=0.66")
        plt.plot(mass, pi75, lw="1.2", color="k", ls="--", label="n=0.75")
        plt.plot(mass, pi825, lw="1.2", color="r", ls="--", label="n=0.825")
        plt.plot(mass, pi975, lw="1.2", color="c", ls="--", label="n=0.975")
        plt.plot(mass, pi110, lw="1.2", color="y", ls="--", label="n=1.1")

        plt.ylim([1,4.7])
        plt.legend()

        plt.figure("eff")
        plt.plot(massTable[0, :]/(np.sqrt(Rg*T)/P), effTable[0, :], lw="1.6", color="m", ls="-", label="n=0.6")
        plt.plot(massTable[1, :]/(np.sqrt(Rg*T)/P), effTable[1, :], lw="1.6", color="b", ls="-", label="n=0.7")
        plt.plot(massTable[2, :]/(np.sqrt(Rg*T)/P), effTable[2, :], lw="1.6", color="g", ls="-", label="n=0.8")
        plt.plot(massTable[3, :]/(np.sqrt(Rg*T)/P), effTable[3, :], lw="1.6", color="k", ls="-", label="n=0.9")
        plt.plot(massTable[4, :]/(np.sqrt(Rg*T)/P), effTable[4, :], lw="1.6", color="r", ls="-", label="n=1.0")

        plt.plot(mass, eff50, lw="1.2", color="m", ls="--", label="n=0.5")
        plt.plot(mass, eff63, lw="1.2", color="b", ls="--", label="n=0.63")
        plt.plot(mass, eff66, lw="1.2", color="g", ls="--", label="n=0.66")
        plt.plot(mass, eff75, lw="1.2", color="k", ls="--", label="n=0.75")
        plt.plot(mass, eff825, lw="1.2", color="r", ls="--", label="n=0.825")
        plt.plot(mass, eff975, lw="1.2", color="c", ls="--", label="n=0.975")
        plt.plot(mass, eff110, lw="1.2", color="y", ls="--", label="n=1.1")

        plt.ylim([0.5, 1.02])
        plt.legend()

        plt.show()

    if False:
        gas1 = Gas(Ts=300, Ps=1e5, G=1.0)
        gas1.O2ratio, gas1.N2ratio = 1.0, 0.0
        gas1.vol2mratio()
        gas2 = Gas(Ts=600, Ps=2e5, G=2.0)
        gas2.O2ratio, gas2.N2ratio = 0.0, 1.0
        gas2.vol2mratio()

        mixer = Mixer(sigma=1.0)
        gasout = mixer.workT4(gas1, gas2)
        print(gasout.G, gasout.T, gasout.P)
        print(gasout.N2mratio, gasout.O2mratio)

    if False:
        Ts, Ps, _ = Altitude(0)
        nlratio = 1.0
        Gratio = 1.5 * nlratio - 0.7  # T4d=1130  nl:0.47 - 1.01
        model = Model()
        W,_,res, = model.solve(T0s=Ts, P0s=Ps, Ma=0.0, nlratio=nlratio, Gratio=Gratio, T4d=1130, isPrint=True, lr=0.1, maxStep=1000)
        _,nl,_ = res

        gen = Gennerator()
        R, U = gen.designR(nl, W)
        print(R, U, U*U/R)

        Radd = 0*np.ones(FANsNUM)
        print(Radd)
        print(gen.workFans(nl, Radd))
        print(gen.designFans(nl, W))

    if False:
        Ts, Ps, _ = Altitude(0)
        nlratio = 0.6
        model = Model()
        _, (Weng, Wcomp, Wturb, U), _ = model.workElec(nlratio, T0s=Ts, P0s=Ps, Ma=0.0, R=1.15, lr=0.1, maxStep=1000, clip=200, isPrint=True)
        print(U, Wgen, U*U/1.15)

    if False:
        Ts, Ps, _ = Altitude(0)
        model = Model()
        model.stablizeT4(nlstart=0.75, nlend=0.8, R=1.6, maxStep=2000, isPlot=True)

    if True:
        Ts, Ps, _ = Altitude(0)
        model = Model()
        model.stablizeT4constW(nlstart=0.78, nlend=0.8, R=1.6, maxStep=1800, isPlot=True)

    if False:
        Ts, Ps, _ = Altitude(0)
        model = Model()
        model.stablizeR(nlstart=0.6, nlend=0.8, T4d=1130, maxStep=2000, isPlot=True)

    if False:
        Ts, Ps, _ = Altitude(0)
        model = Model()
        nllist = [0.52, 0.56, 0.61, 0.665, 0.728, 0.785, 0.83, 0.869, 0.908, 0.95, 1.0]
        for i in range(len(nllist)-1):
            model.stablizeR(nlstart=nllist[i], nlend=nllist[i+1], T4d=1130, maxStep=2000, isPlot=False)

    if False:
        Ts, Ps, _ = Altitude(0)
        model = Model()
        model.stablizeFans(nlstart=0.6, nlend=0.8, T4d=1130, maxStep=400, isPlot=True)

    if False:
        Ts, Ps, _ = Altitude(0)
        model = Model()
        model.stablizeFansConst(nlstart=0.6, nlend=0.8, T4d=1130, maxStep=1000, isPlot=True)


    print("code end~")
