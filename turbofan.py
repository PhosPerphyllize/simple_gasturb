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

class Turbofan():
    def __init__(self):
        pass

    def design(self, state=[1.5, 0.86,   0.57, 17,   0.85,   0.95,     0.985,   1400,  0.90,    0.90],
               Ts=288.15, Ps=101325.0, Ma=0.0, Gd=1.0, isPrint: bool = True):
        pif, fanEff, ratio, pic, compEff, picombu, combuEff, Tout, turbEff, turbEff2  = state
        leakFan_per, leakComp_per =0.0, 0.013
        detector = Detector()

        gas = Gas(Ts, Ps, Ma, Gd)
        detector.get(gas, "engine inlet", isPrint)

        intake = Intake(sigma=1)
        gas, inlet_drag = intake.work(gas, Ma)

        detector.get(gas, "Fan inlet", isPrint)
        fan = Compressor(ifTable=False, leakG_per=leakFan_per)
        gas, Wfan, = fan.workDesign(gas, pi=pif, eff=fanEff)
        fan_leak = copy.deepcopy(gas)
        fan_leak.G = gas.G / (1 - leakFan_per) * leakFan_per

        seprator = Seprator(sigma=1.0)
        gas, gasOut = seprator.work(gas, ratio=ratio)

        detector.get(gasOut, "Bypass inlet", isPrint)
        bypass = Nozzle(type="con", Pamb=Ps, Cd=0.985, Cf=0.985, sigma=1)
        gasOut, gro_thrust2, bypassA8, bypassA9, list = bypass.workDesign(gasOut)
        detector.get(gasOut, "Bypass outlet", isPrint)
        byTs_out, byPs_out, byMa_out, byV_out = list

        duct0 = Duct(sigma=1)
        gas = duct0.work(gas)

        detector.get(gas, "Compressor inlet", isPrint)
        comp = Compressor(ifTable=False, leakG_per=leakComp_per)
        gas, Wcomp, = comp.workDesign(gas, pi=pic, eff=compEff)
        comp_leak = copy.deepcopy(gas)
        comp_leak.G = gas.G / (1-leakComp_per) * leakComp_per

        duct1 = Duct(sigma=1)
        gas = duct1.work(gas)

        detector.get(gas, "Combustor inlet", isPrint)
        H2LHV = 43.325 * 1e6  # lower heat value ，J/kg
        combu = Combustor(sigma=picombu, eff=combuEff, LHV=H2LHV, fuelType="CH4")
        gas, Gfuel = combu.workT4(gas, Tout=Tout)

        detector.get(gas, "Mixer inlet", isPrint)
        mixer = Mixer(sigma=1)
        comp_leak.P = gas.P
        gas = mixer.work(gas, comp_leak)

        Wfan, Wcomp = 1.00*Wfan, 1.008*Wcomp
        detector.get(gas, "Turbine1 inlet", isPrint)
        nleff = 0.99
        turb = Turbine(ifTable=False, leakG_per=0.0)
        gas, pit, = turb.workDesign_work(gas, work=1.005*(Wcomp/nleff), eff=turbEff)

        detector.get(gas, "Turbine2 inlet", isPrint)
        nl2eff = 0.99
        Wele = 0*1000
        turb2 = Turbine(ifTable=False, leakG_per=0.0)
        gas, pit2, = turb2.workDesign_work(gas, work=1.005*(Wfan/nl2eff + Wele), eff=turbEff2)

        detector.get(gas, "Turbine outlet", isPrint)
        duct2 = Duct(sigma=0.99)
        gas = duct2.work(gas)

        detector.get(gas, "Nozzle inlet", isPrint)
        nozzle = Nozzle(type="con", Pamb=Ps, Cd=0.985, Cf=0.985, sigma=1)
        gas, gro_thrust, A8, A9, list = nozzle.workDesign(gas)
        detector.get(gas, "Engine outlet", isPrint)
        Ts_out, Ps_out, Ma_out, V_out = list

        thrust = (gro_thrust+gro_thrust2-inlet_drag)/10
        if isPrint:
            log(" ")
            log("Engine Thrust:       %.5f daN"%(thrust))
            log("Engine Thrust:       %.5f kN"%(thrust/100))
            log("Engine Thrust:       %.5f kgf"%(thrust*10/9.8))
            log("Engine SFC:          %.5f kg/(daN*h)"%(Gfuel*3600/thrust))
            log("Engine SFC:          %.5f g/(kN*s)"%(Gfuel*1000*100/thrust))
            log("Engine SFC:          %.5f kg/(kgf*s)"%(Gfuel*3600*0.98/thrust))
            log("Engine Fuleflow:     %.5f kg/s"%(Gfuel))
            log(" ")
            log("Engine outlet Pi:    %.5f " % (gas.P / Ps))
            log("Engine outlet Ma:    %.5f " % (Ma_out))
            log("Engine outlet Vel:   %.5f m/s" % (V_out))
            log("Engine outlet Ps:    %.5f kPa" % (Ps_out / 1000))
            log("Engine A8:           %.5f m2" % (A8))
            log("Engine A9:           %.5f m2" % (A9))
            log("Engine A9 diam:      %.5f m" % (np.sqrt(A9 / np.pi)))
            log(" ")
            log("Engine bypass outlet Pi:    %.5f " % (gasOut.P / Ps))
            log("Engine bypass outlet Ma:    %.5f " % (byMa_out))
            log("Engine bypass outlet Vel:   %.5f m/s" % (byV_out))
            log("Engine bypass outlet Ps:    %.5f kPa" % (byPs_out / 1000))
            log("Engine bypass A8:           %.5f m2" % (bypassA8))
            log("Engine bypass A9:           %.5f m2" % (bypassA9))
            log("Engine bypass A9 diam:      %.5f m" % (np.sqrt(bypassA9 / np.pi)))
            log(" ")
            log("Fan        pi:   %.5f   eff:   %.5f" % (pif, fanEff))
            log("Fan        work: %.5f kw" % (Wfan/1000))
            log("Compressor pi:   %.5f   eff:   %.5f"%(pic, compEff))
            log("Compressor work: %.5f kw"%(Wcomp/1000))
            log("Turbine1   pi:   %.5f  eff:   %.5f"%(1/pit, turbEff))
            log("Turbine1   work: %.5f kw"%(Wcomp/nleff/1000+Wele/1000))
            log("Turbine2   pi:   %.5f  eff:   %.5f" % (1 / pit2, turbEff2))
            log("Turbine2   work: %.5f kw" % (Wfan/nleff/1000 + Wele/1000))
            log("compEff: %f turbEff: %f combuEff: %f" % (compEff, turbEff, combuEff))
            log(" ")
            log("Engine gro Thrust:        %.5f kN" % ((gro_thrust+gro_thrust2) / 1000))
            log("Engine core Thrust:       %.5f kN" % (gro_thrust / 1000))
            log("Engine bypass Thrust:     %.5f kN" % (gro_thrust2 / 1000))
            log("Engine inlet Thrust:      %.5f kN" % (inlet_drag / 1000))
            log("Engine total outlet Area: %.5f m2" % (bypassA9+A9))
            log("bypass V8out/core V8out:  %.5f" % (byV_out/V_out))
            log(" ")


        return 0

    def work(self, state=[1.5, 0.86,   0.57, 17,   0.85,   0.95,     0.985,   1400,  0.90,    0.90],
                design=[1.0, 1.0],
                Ts=288.15, Ps=101325.0, Ma=0.0, Gd=1.0, isPrint: bool = True):
        pif, fanEff, ratio, pic, compEff, picombu, combuEff, Tout, turbEff, turbEff2  = state
        A9, bypassA9 = design
        leakFan_per, leakComp_per =0.0, 0.013
        detector = Detector()

        gas = Gas(Ts, Ps, Ma, Gd)
        detector.get(gas, "engine inlet", isPrint)

        intake = Intake(sigma=1)
        gas, inlet_drag = intake.work(gas, Ma)

        detector.get(gas, "Fan inlet", isPrint)
        fan = Compressor(ifTable=False, leakG_per=leakFan_per)
        gas, Wfan, = fan.workDesign(gas, pi=pif, eff=fanEff)
        fan_leak = copy.deepcopy(gas)
        fan_leak.G = gas.G / (1 - leakFan_per) * leakFan_per

        seprator = Seprator(sigma=1.0)
        gas, gasOut = seprator.work(gas, ratio=ratio)

        detector.get(gasOut, "Bypass inlet", isPrint)
        bypass = Nozzle(type="con", Pamb=Ps, Cd=0.985, Cf=0.985, sigma=1)
        # bypass.A8, bypass.A9 = bypassA9, bypassA9
        gasOut, gro_thrust2, bypassA8, bypassA9, list = bypass.workDesign(gasOut)
        detector.get(gasOut, "Bypass outlet", isPrint)
        byTs_out, byPs_out, byMa_out, byV_out = list

        duct0 = Duct(sigma=1)
        gas = duct0.work(gas)

        detector.get(gas, "Compressor inlet", isPrint)
        comp = Compressor(ifTable=False, leakG_per=leakComp_per)
        gas, Wcomp, = comp.workDesign(gas, pi=pic, eff=compEff)
        comp_leak = copy.deepcopy(gas)
        comp_leak.G = gas.G / (1-leakComp_per) * leakComp_per

        duct1 = Duct(sigma=1)
        gas = duct1.work(gas)

        detector.get(gas, "Combustor inlet", isPrint)
        H2LHV = 43.325 * 1e6  # lower heat value ，J/kg
        combu = Combustor(sigma=picombu, eff=combuEff, LHV=H2LHV, fuelType="CH4")
        gas, Gfuel = combu.workT4(gas, Tout=Tout)

        detector.get(gas, "Mixer inlet", isPrint)
        mixer = Mixer(sigma=1)
        comp_leak.P = gas.P
        gas = mixer.work(gas, comp_leak)

        Wfan, Wcomp = 1.00*Wfan, 1.008*Wcomp
        detector.get(gas, "Turbine1 inlet", isPrint)
        nleff = 0.99
        turb = Turbine(ifTable=False, leakG_per=0.0)
        gas, pit, = turb.workDesign_work(gas, work=1.005*(Wcomp/nleff), eff=turbEff)

        detector.get(gas, "Turbine2 inlet", isPrint)
        nl2eff = 0.99
        Wele = 0*1000
        turb2 = Turbine(ifTable=False, leakG_per=0.0)
        gas, pit2, = turb2.workDesign_work(gas, work=1.005*(Wfan/nl2eff + Wele), eff=turbEff2)

        detector.get(gas, "Turbine outlet", isPrint)
        duct2 = Duct(sigma=0.99)
        gas = duct2.work(gas)

        Gin = gas.G
        detector.get(gas, "Nozzle inlet", isPrint)
        nozzle = Nozzle(type="con", Pamb=Ps, Cd=0.985, Cf=0.985, sigma=1)
        nozzle.A8, nozzle.A9 = A9, A9
        gas, gro_thrust, Gout, list = nozzle.work(gas)
        detector.get(gas, "Engine outlet", isPrint)
        Ts_out, Ps_out, Ma_out, V_out, _, _ = list

        thrust = (gro_thrust+gro_thrust2-inlet_drag)/10
        if isPrint:
            log(" ")
            log("Engine Thrust:       %.5f daN"%(thrust))
            log("Engine Thrust:       %.5f kN"%(thrust/100))
            log("Engine Thrust:       %.5f kgf"%(thrust*10/9.8))
            log("Engine SFC:          %.5f kg/(daN*h)"%(Gfuel*3600/thrust))
            log("Engine SFC:          %.5f g/(kN*s)"%(Gfuel*1000*100/thrust))
            log("Engine SFC:          %.5f kg/(kgf*s)"%(Gfuel*3600*0.98/thrust))
            log("Engine Fuleflow:     %.5f kg/s"%(Gfuel))
            log(" ")
            log("Engine outlet Pi:    %.5f " % (gas.P / Ps))
            log("Engine outlet Ma:    %.5f " % (Ma_out))
            log("Engine outlet Vel:   %.5f m/s" % (V_out))
            log("Engine outlet Ps:    %.5f kPa" % (Ps_out / 1000))
            log("Engine A8:           %.5f m2" % (nozzle.A8))
            log("Engine A9:           %.5f m2" % (nozzle.A9))
            log("Engine A9 diam:      %.5f m" % (np.sqrt(A9 / np.pi)))
            log(" ")
            log("Engine bypass outlet Pi:    %.5f " % (gasOut.P / Ps))
            log("Engine bypass outlet Ma:    %.5f " % (byMa_out))
            log("Engine bypass outlet Vel:   %.5f m/s" % (byV_out))
            log("Engine bypass outlet Ps:    %.5f kPa" % (byPs_out / 1000))
            log("Engine bypass A8:           %.5f m2" % (bypass.A8))
            log("Engine bypass A9:           %.5f m2" % (bypass.A9))
            log("Engine bypass A9 diam:      %.5f m" % (np.sqrt(bypass.A9 / np.pi)))
            log(" ")
            log("Fan        pi:   %.5f   eff:   %.5f" % (pif, fanEff))
            log("Fan        work: %.5f kw" % (Wfan/1000))
            log("Compressor pi:   %.5f   eff:   %.5f"%(pic, compEff))
            log("Compressor work: %.5f kw"%(Wcomp/1000))
            log("Turbine1   pi:   %.5f  eff:   %.5f"%(1/pit, turbEff))
            log("Turbine1   work: %.5f kw"%(Wcomp/nleff/1000+Wele/1000))
            log("Turbine2   pi:   %.5f  eff:   %.5f" % (1 / pit2, turbEff2))
            log("Turbine2   work: %.5f kw" % (Wfan/nleff/1000 + Wele/1000))
            log("compEff: %f turbEff: %f combuEff: %f" % (compEff, turbEff, combuEff))
            log(" ")
            log("Engine gro Thrust:        %.5f kN" % ((gro_thrust+gro_thrust2) / 1000))
            log("Engine core Thrust:       %.5f kN" % (gro_thrust / 1000))
            log("Engine bypass Thrust:     %.5f kN" % (gro_thrust2 / 1000))
            log("Engine inlet Thrust:      %.5f kN" % (inlet_drag / 1000))
            log("Engine total outlet Area: %.5f m2" % (bypassA9+A9))
            log("bypass V8out/core V8out:  %.5f" % (byV_out/V_out))
            log(" ")


        return Gin-Gout

    def work_ori(self,state, Ts=288.15, Ps=101325.0, Ma=0.0, Gd=1.0, A9=1.0, isPrint: bool = True):
        detector = Detector()
        leakG_per = 0.013
        pic, compEff, picombu, combuEff, Tout, turbEff,  = state

        gas = Gas(Ts, Ps, Ma, Gd)
        detector.get(gas, "engine inlet", isPrint)

        intake = Intake(sigma=1)
        gas, inlet_drag = intake.work(gas, Ma)

        detector.get(gas, "Compressor inlet", isPrint)
        comp = Compressor(ifTable=False, leakG_per=leakG_per)
        gas, Wcomp, = comp.workDesign(gas, pi=pic, eff=compEff)
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
        efft = turbEff
        turb = Turbine(ifTable=False, leakG_per=0.0)
        gas, pit, = turb.workDesign_work(gas, work=1.00*(Wcomp/nleff+Wele), eff=efft)

        detector.get(gas, "Turbine outlet", isPrint)
        duct2 = Duct(sigma=0.99)
        gas = duct2.work(gas)
        Gin = gas.G

        detector.get(gas, "Nozzle inlet", isPrint)
        nozzle = Nozzle(type="con", Pamb=Ps, Cd=0.985, Cf=0.985, sigma=1)
        nozzle.A8, nozzle.A9 = A9, A9
        gas, gro_thrust, Gout, list = nozzle.work(gas)
        detector.get(gas, "Engine outlet", isPrint)
        Ts_out, Ps_out, Ma_out, V_out,_,_ = list

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
            log("Engine A8:           %.5f m2"%(nozzle.A8))
            log("Engine A9:           %.5f m2"%(nozzle.A9))
            log("Engine A9 diam:      %.5f m"%(np.sqrt(A9/np.pi)))
            log(" ")
            log("Compressor pi:   %.5f   eff:   %.5f"%(pic, compEff))
            log("Compressor work: %.5f kw"%(Wcomp/1000))
            log("Turbine    pi:   %.5f  eff:   %.5f"%(1/pit, efft))
            log("Turbine    work: %.5f kw"%(Wcomp/nleff/1000+Wele/1000))
            log("compEff: %f turbEff: %f combuEff: %f" % (compEff, turbEff, combuEff))
            log(" ")
            log("Engine gro Thrust:   %.5f kN" % (gro_thrust / 1000))
            log("Engine inlet Thrust: %.5f kN" % (inlet_drag / 1000))
            log(" ")

        return Gin-Gout



if __name__ == '__main__':
    if True:
        state = [1.5, 0.86,   0.57, 17,   0.85,   0.95,     0.985,   1400,  0.90,    0.90]
        # state=pif, fanEff, ratio, pic, compEff, picombu, combuEff, Tout, turbEff, turbEff2
        Ma = 0.0
        Ts, Ps, _ = Altitude(0)
        Gd = 120

        gas = Gas(Ts, Ps, Ma, Gd)
        T0d, P0d, _ = Altitude(0)
        G1 = Gd * gas.P / P0d * np.sqrt(T0d / gas.T) * 1

        jet = Turbofan()
        jet.design(state, Ts, Ps, Ma, Gd=G1, isPrint=True)


    if True:
        state = [1.5, 0.86,   0.57, 17,   0.85,   0.95,     0.985,   1400,  0.90,    0.90]
        # state=pif, fanEff, ratio, pic, compEff, picombu, combuEff, Tout, turbEff, turbEff2
        design = [0.15827, 0.13676]
        #          A9,     bypassA9
        Ma = 0
        Ts, Ps, _ = Altitude(0)

        Gd = 120
        gas = Gas(Ts, Ps, Ma, Gd)
        T0d, P0d, _ = Altitude(0)
        G1 = Gd * gas.P / P0d * np.sqrt(T0d / gas.T) * 1
        G0 = G1

        jet = Turbofan()
        detaR = 0.00001
        G2 = (1 - detaR) * G1  # to get the derivative of G
        lr = 0.2  # learning rate
        step = 1
        while True:
            err1 = jet.work(state, design, Ts, Ps, Ma, Gd=G1, isPrint=False)
            if (abs(err1) < 0.000001):
                break

            err2 = jet.work(state, design, Ts, Ps, Ma, Gd=G2, isPrint=False)

            grad = (err1 - err2) / (G1 - G2)
            if (grad > 1e9):
                print("\nWarning: gradient > 1e9!\n")

            G1 = G1 - lr * err1 / grad
            G2 = (1 - detaR) * G1

            if step > 2000:
                print("\nstep %d grad: %f\n" % (step, grad))
                break
            step = step + 1
        jet.work(state, design, Ts, Ps, Ma, Gd=G1, isPrint=True)
        print("G1/Gcor: %f"%(G1/G0))