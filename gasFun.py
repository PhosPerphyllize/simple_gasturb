import numpy as np

def Altitude(Height):
    # input: unit, m   0-30km
    Height = Height/1000
    if Height < 0:
        print("Warning: Height must >= 0, Height: %f km"%(Height))
        Height = 0

    if Height < 11:
        T = 288.15 - 6.5*Height
        P = 101325*np.power(1-0.0225577*Height, 5.25588)
    elif Height<25:
        T = 216.65
        P = 22632*np.power(np.e, (11-Height)/6.34162)
    else:
        T = 216.65 + 3*(Height-25)
        P = 2488.6 * np.power(216.15/T, 11.8)

    den = P/(287.06*T)
    return T,P,den

def airCp(T):
    # T unit:K, suggest: 250-1200K, return unit:J/(kg.K)
    fai = T / 1000
    Cp = 1.05 - 0.365 * fai + 0.85 * fai * fai - 0.39 * fai * fai * fai

    # T suggest: 250-2500K, return unit:J/(kg.K)
    # fai = T-273
    # Cp = 1.00261 - 3.48511e-5*fai + 1.05471e-7*fai**2 - 6.43638e-11*fai**3 + 1.13284e-14*fai**4
    return  Cp*1000

def H2OCp(T):
    # T unit:K, suggest: 250-1200K, return unit:J/(kg.K)
    fai = T/1000
    Cp = 1.79 - 0.107*fai + 0.586*fai*fai - 0.2*fai*fai*fai
    return  Cp*1000

def N2Cp(T):
    # T unit:K, suggest: 250-1200K, return unit:J/(kg.K)
    fai = T/1000
    Cp = 1.11 - 0.48*fai + 0.96*fai*fai - 0.42*fai*fai*fai
    return  Cp*1000

def O2Cp(T):
    # T unit:K, suggest: 250-1200K, return unit:J/(kg.K)
    fai = T/1000
    Cp = 0.88 - 0.0001*fai + 0.54*fai*fai - 0.33*fai*fai*fai
    return  Cp*1000

def CO2Cp(T):
    # T unit:K, suggest: 250-1200K, return unit:J/(kg.K)
    fai = T / 1000
    Cp = 0.45 + 1.67 * fai - 1.27 * fai * fai + 0.39 * fai * fai * fai
    return Cp * 1000

def airk(T, R=287.06):
    Cp = airCp(T)
    return Cp/(Cp-R)

def Ma2tau(Ma, k=1.4):
    # total temperature/static temperature
    return 1 + 0.5*(k-1)*Ma*Ma

def Ma2pi(Ma, k=1.4):
    # total pressure/static pressure
    return np.power(1+0.5*(k-1)*Ma*Ma, k/(k-1))

def pi2Ma(pi, k=1.4):
    # total pressure/static pressure
    return np.sqrt(2/(k-1)*(np.power(pi, (k-1)/k) - 1))

if __name__ == '__main__':
    if False:
        import matplotlib.pyplot as plt
        num = 200
        height = np.linspace(0, 30000, num)
        T, P = np.zeros_like(height), np.zeros_like(height)
        for i in range(num):
            T[i],P[i],_ = Altitude(height[i])

        fig, ax1 = plt.subplots()
        ax1.plot(height, T, 'g-', label='T')
        ax1.set_xlabel('X data')
        ax1.set_ylabel('Y1 data', color='g')

        # 创建共享x轴的第二个坐标轴
        ax2 = ax1.twinx()
        ax2.plot(height, P, 'b-', label='P')
        ax2.set_ylabel('Y2 data', color='b')

        # 添加图例
        ax1.legend()
        ax2.legend()
        plt.tight_layout()

        plt.show()

    if True:
        T = 500
        Cp = 0.78*N2Cp(T)+0.22*O2Cp(T)
        M = 0.78*28+0.22*32
        Rg = 8314.5/M
        k = Cp/(Cp-Rg)
        print(Cp, k)
        print(Rg)


    if False:
        T = 300
        Cp = N2Cp(T)
        M = 28
        Rg = 8314.5 / M
        k = Cp / (Cp - Rg)
        print(Cp, k)

    if False:
        print(pi2Ma(0.8))