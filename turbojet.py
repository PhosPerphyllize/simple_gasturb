import numpy as np
from component import *
log.changeDir("./turbojet_log.txt")

def effReCor(Height, Ma=0.8):
    # cor the eff due to height
    # height: unit:m
    compEff, turbEff, combuEff = 0.78, 0.86, 0.98
    compEff25, turbEff25, combuEff25 = 0.68, 0.75, 0.90

    # if Height > 22000:
    #     compEff = (compEff25 - compEff)/(25000-10000)*(Height-10000) + compEff
    #     turbEff = (turbEff25 - turbEff)/(25000-10000)*(Height-10000) + turbEff
    #     combuEff = (combuEff25 - combuEff)/(25000-10000)*(Height-10000) + combuEff

    if Height==11000:
        return 0.7774, 0.8540, 0.9626, 5.917

    if Height==21500:
        return 0.7346, 0.7565, 0.8631, 3.841

    return compEff, turbEff, combuEff, 5.47

class Turbojet():
    def __init__(self):
        pass

    def design(self,state=[0.78, 0.86, 0.98, 5.47],
               Ts=288.15, Ps=101325.0, Ma=0.0, Height=0.0, Gd=1.0, Tout=1200, isPrint: bool = True):
        detector = Detector()
        leakG_per = 0.013
        compEff, turbEff, combuEff,compPi = state

        gas = Gas(Ts, Ps, 0, Gd)
        detector.get(gas, "Ambient", isPrint)

        gas = Gas(Ts, Ps, Ma, Gd)
        detector.get(gas, "engine inlet", isPrint)

        intake = Intake(sigma=1)
        gas, inlet_drag = intake.work(gas, Ma)

        detector.get(gas, "Compressor inlet", isPrint)
        comp = Compressor(ifTable=False, leakG_per=leakG_per)
        gas, Wcomp, = comp.workDesign(gas, pi=compPi, eff=compEff)
        gas_leak = copy.deepcopy(gas)
        gas_leak.G = gas.G / (1-leakG_per) * leakG_per

        duct1 = Duct(sigma=1)
        gas = duct1.work(gas)

        detector.get(gas, "Combustor inlet", isPrint)
        H2LHV = 43.325 * 1e6  # lower heat value ，J/kg
        combu = Combustor(sigma=0.96, eff=combuEff, LHV=H2LHV, fuelType="CH4")
        gas, Gfuel = combu.workT4(gas, Tout=Tout)

        detector.get(gas, "Mixer inlet", isPrint)
        mixer = Mixer(sigma=0.99)
        gas_leak.P = gas.P
        gas = mixer.work(gas, gas_leak)

        detector.get(gas, "Turbine inlet", isPrint)
        nleff = 0.99
        Wele = 20 *1000
        turb = Turbine(ifTable=False, leakG_per=0.0)
        gas, pit, = turb.workDesign_work(gas, work=1.002*(Wcomp/nleff+Wele), eff=turbEff)

        detector.get(gas, "Turbine outlet", isPrint)
        duct2 = Duct(sigma=0.99)
        gas = duct2.work(gas)

        detector.get(gas, "Nozzle inlet", isPrint)
        nozzle = Nozzle(type="con", Pamb=Ps, Cd=0.985, Cf=0.985, sigma=1)
        gas, gro_thrust, A8, A9, list = nozzle.workDesign(gas)
        detector.get(gas, "Engine outlet", isPrint)
        Ts_out, Ps_out, Ma_out, V_out = list

        thrust = (gro_thrust-inlet_drag)/10
        if isPrint:
            log(" ")
            log("Engine Thrust:       %.5f daN"%(thrust))
            log("Engine Thrust:       %.5f kN"%(thrust/100))
            log("Engine Thrust:       %.5f kgf"%(thrust*10/9.8))
            log("Engine SFC:          %.5f kg/(daN*h)"%(Gfuel*3600/thrust))
            log("Engine SFC:          %.5f g/(kN*s)"%(Gfuel*1000*100/thrust))
            log("Engine SFC:          %.5f kg/(kgf*s)"%(Gfuel*3600*0.98/thrust))
            log("Engine Fuleflow:     %.5f kg/s"%(Gfuel))
            log("Engine outlet Pi:    %.5f "%(gas.P/Ps))
            log("Engine outlet Ma:    %.5f "%(Ma_out))
            log("Engine outlet Vel:   %.5f m/s"%(V_out))
            log("Engine outlet Ps:    %.5f kPa"%(Ps_out/1000))
            log("Engine A8:           %.5f m2"%(A8))
            log("Engine A9:           %.5f m2"%(A9))
            log("Engine A9 diam:      %.5f m"%(np.sqrt(A9/np.pi)))
            log(" ")
            log("Compressor pi:   %.5f   eff:   %.5f"%(compPi, compEff))
            log("Compressor work: %.5f kw"%(Wcomp/1000))
            log("Turbine    pi:   %.5f  eff:   %.5f"%(1/pit, turbEff))
            log("Turbine    work: %.5f kw"%(Wcomp/nleff/1000+Wele/1000))
            log("compEff: %f turbEff: %f combuEff: %f" % (compEff, turbEff, combuEff))
            log(" ")
            log("Engine gro Thrust:   %.5f kN" % (gro_thrust / 1000))
            log("Engine inlet Thrust: %.5f kN" % (inlet_drag / 1000))
            log(" ")

        return 0

    def work(self,state=[0.78, 0.86, 0.98, 5.47],
             Ts=288.15, Ps=101325.0, Ma=0.0, Gd=1.0, Tout=1200, A9=1.0, isPrint: bool = True):
        detector = Detector()
        leakG_per = 0.013
        compEff, turbEff, combuEff, compPi = state

        gas = Gas(Ts, Ps, 0, Gd)
        detector.get(gas, "Ambient", isPrint)

        gas = Gas(Ts, Ps, Ma, Gd)
        detector.get(gas, "engine inlet", isPrint)

        intake = Intake(sigma=0.99, Ain=0.10406)
        gas, inlet_drag = intake.work(gas, Ma)

        detector.get(gas, "Compressor inlet", isPrint)
        comp = Compressor(ifTable=False, leakG_per=leakG_per)
        gas, Wcomp, = comp.workDesign(gas, pi=compPi, eff=compEff)
        gas_leak = copy.deepcopy(gas)
        gas_leak.G = gas.G / (1-leakG_per) * leakG_per

        duct1 = Duct(sigma=0.99)
        gas = duct1.work(gas)

        detector.get(gas, "Combustor inlet", isPrint)
        H2LHV = 43.325 * 1e6  # lower heat value ，J/kg
        combu = Combustor(sigma=0.96, eff=combuEff, LHV=H2LHV, fuelType="CH4")
        gas, Gfuel = combu.workT4(gas, Tout=Tout)

        detector.get(gas, "Mixer inlet", isPrint)
        mixer = Mixer(sigma=0.99)
        gas_leak.P = gas.P
        gas = mixer.work(gas, gas_leak)

        detector.get(gas, "Turbine inlet", isPrint)
        nleff = 0.99
        Wele = 20 *1000
        turb = Turbine(ifTable=False, leakG_per=0.0)
        gas, pit, = turb.workDesign_work(gas, work=1.002*(Wcomp/nleff+Wele), eff=turbEff)

        detector.get(gas, "Turbine outlet", isPrint)
        duct2 = Duct(sigma=0.99)
        gas = duct2.work(gas)
        Gin = gas.G

        detector.get(gas, "Nozzle inlet", isPrint)
        nozzle = Nozzle(type="con", Pamb=Ps, Cd=0.985, Cf=0.985, sigma=0.985)
        nozzle.A8, nozzle.A9 = A9, A9
        gas, gro_thrust, Gout, list = nozzle.work(gas)
        detector.get(gas, "Engine outlet", isPrint)
        Ts_out, Ps_out, Ma_out, V_out,_,_ = list

        thrust = (gro_thrust-inlet_drag)/10*0.98
        Gfuel /=0.98
        if isPrint:
            log(" ")
            log("Engine Thrust:       %.5f daN"%(thrust))
            log("Engine Thrust:       %.5f kN"%(thrust/100))
            log("Engine Thrust:       %.5f kgf"%(thrust*10/9.8))
            log("Engine SFC:          %.5f kg/(daN*h)"%(Gfuel*3600/thrust))
            log("Engine SFC:          %.5f g/(kN*s)"%(Gfuel*1000*100/thrust))
            log("Engine SFC:          %.5f kg/(kgf*h)"%(Gfuel*3600*0.98/thrust))
            log("Engine Fuleflow:     %.5f kg/s"%(Gfuel))
            log("Engine outlet Pi:    %.5f "%(gas.P/Ps))
            log("Engine outlet Ma:    %.5f "%(Ma_out))
            log("Engine outlet Vel:   %.5f m/s"%(V_out))
            log("Engine outlet Ps:    %.5f kPa"%(Ps_out/1000))
            log("Engine A8:           %.5f m2"%(nozzle.A8))
            log("Engine A9:           %.5f m2"%(nozzle.A9))
            log("Engine A9 diam:      %.5f m"%(np.sqrt(A9/np.pi)))
            log(" ")
            log("Compressor pi:   %.5f   eff:   %.5f"%(compPi, compEff))
            log("Compressor work: %.5f kw"%(Wcomp/1000))
            log("Turbine    pi:   %.5f  eff:   %.5f"%(1/pit, turbEff))
            log("Turbine    work: %.5f kw"%(Wcomp/nleff/1000+Wele/1000))
            log("compEff: %f turbEff: %f combuEff: %f" % (compEff, turbEff, combuEff))
            log(" ")
            log("Engine intake G:     %.5f kg/s" % (intake.Gin))
            log("Engine gro Thrust:   %.5f kN" % (gro_thrust / 1000))
            log("Engine inlet Thrust: %.5f kN" % (inlet_drag / 1000))
            log("Engine Thrust:       %.5f kgf" % (thrust/0.98 * 10 / 9.8))
            log("Engine SFC:          %.5f kg/(kgf*h)" % (Gfuel * 3600 * 0.98 / thrust*0.98*0.98))
            log(" ")

        return Gin-Gout



if __name__ == '__main__':
    if True:
        state = [0.8, 0.86, 0.98, 5.78]
        #       compEff, turbEff, combuEff, compPi
        Ma = 0.7
        Ts, Ps, _ = Altitude(5000)
        Gd = 13.4

        gas = Gas(Ts, Ps, Ma, Gd)
        T0d, P0d, _ = Altitude(0)
        G1 = Gd * gas.P / P0d * np.sqrt(T0d / gas.T) * 1

        jet = Turbojet()
        jet.design(state, Ts, Ps, Ma, Gd=G1, Tout=1180, isPrint=True)


    if True:
        print("########################\n######################")
        state = [0.75, 0.83, 0.98, 6.4]
        #       compEff, turbEff, combuEff, compPi
        Ma = 0.8
        height = 20000
        A9 = 0.0520929
        # A9 = 0.05641043
        # A9 = 0.054
        T4 = 1200
        Ts, Ps, _ = Altitude(height)

        Gd = 13.4
        gas = Gas(Ts, Ps, Ma, Gd)
        T0d, P0d, _ = Altitude(0)
        G1 = Gd * gas.P / P0d * np.sqrt(T0d / gas.T) * 1
        G0 = G1
        T1,P1 = gas.T, gas.P

        jet = Turbojet()
        detaR = 0.00001
        G2 = (1 - detaR) * G1  # to get the derivative of G
        lr = 0.2  # learning rate
        step = 1
        while True:
            err1 = jet.work(state, Ts, Ps, Ma, Gd=G1, Tout=T4, A9=A9, isPrint=False)
            if (abs(err1) < 0.000001):
                break

            err2 = jet.work(state, Ts, Ps, Ma, Gd=G2, Tout=T4, A9=A9, isPrint=False)

            grad = (err1 - err2) / (G1 - G2)
            if (grad > 1e9):
                print("\nWarning: gradient > 1e9!\n")

            G1 = G1 - lr * err1 / grad
            G2 = (1 - detaR) * G1

            if step > 2000:
                print("\nstep %d grad: %f\n" % (step, grad))
                break
            step = step + 1
        jet.work(state, Ts, Ps, Ma, Gd=G1, Tout=T4, A9=A9, isPrint=True)
        print("G1/Gcor: %f"%(G1/G0))
        print("nl/nl0: %f"%(np.sqrt(T1/T0d)))
        print("nl0/nl: %f"%(np.sqrt(T0d/T1)))
        print("best T4: %f"%(T1/T0d*1200))