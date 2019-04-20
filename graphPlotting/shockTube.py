import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import fsolve

xS       = []
density  = []
velocity = []
pressure = []
intEnergy= []

def riemann(num):
    # --------o-- Input Section --o-------
    # Sod Shock Tube Problem
    PrL = 1.0; PrR = 0.1; rhoL = 1.0; rhoR = 0.125;
    uL = 0.0; uR = 0.0; t_end = num; mu = 0.35; # mu = dt/dx
    # --------o--END Input Section --o-----------------
    
    gamma = 1.4
    gammaA = gamma - 1.0
    gammaB = 1/ gammaA
    gammaC = gamma + 1.0
    PRL = PrR/PrL
    cR = np. sqrt ( gamma * PrR/ rhoR )
    cL = np. sqrt ( gamma * PrL/ rhoL )
    CRL = cR/cL
    machL = (uL - uR )/ cL

    def func (p34 ):
        wortel = np.sqrt (2 * gamma * ( gammaA + gammaC * p34 ))
        yy = ( gammaA * CRL * (p34 - 1)) / wortel
        yy = (1 + machL * gammaA /2 - yy )**(2 * gamma / gammaA )
        y = yy/p34 - PRL
        return y

    p34 = fsolve (func ,3.0) # p34 = p3/p4
    print (p34)
    p3 = p34 * PrR
    alpha = gammaC / gammaA
    rho3 = rhoR * (1 + alpha * p34 )/( alpha + p34 )
    rho2 = rhoL * (p34 * PrR /PrL )**(1/ gamma )
    u2 = uL -uR +(2/ gammaA )* cL *(1 - (p34 * PrR /PrL )**( gammaA /(2 * gamma )))
    c2 = np. sqrt ( gamma * p3/ rho2 )
    spos = (0.5 + t_end * cR *np. sqrt ( gammaA /(2 * gamma ) + gammaC /(2 * gamma ) * p34 ) + t_end * uR) # Shock position

    conpos = 0.5 + u2 * t_end + t_end * uR # Position of contact discontinuity
    pos1 = 0.5 + (uL - cL) * t_end # Start of expansion fan
    pos2 = 0.5 + (u2 + uR - c2) * t_end # End of expansion fan

    xgrid = np. linspace (0 ,1 ,100)
    PrE = np. zeros ((1 , len( xgrid )))
    uE= np. zeros ((1 , len ( xgrid )))

    rhoE = np. zeros ((1 , len ( xgrid )))
    machE = np. zeros ((1 , len( xgrid )))
    cE = np. zeros ((1 , len( xgrid )))
    xgrid = np. matrix ( xgrid )

    for i in range (0, xgrid . size ):
        if xgrid [0,i] <= pos1 :
            PrE [0,i] = PrL
            rhoE [0,i] = rhoL
            uE [0,i] = uL
            cE [0,i]= np. sqrt ( gamma *PrE [0,i]/ rhoE [0,i]);
            machE [0,i]= uE [0,i]/ cE [0,i];
        elif xgrid [0,i] <= pos2 :
            PrE [0,i] = (PrL *(1 + ( pos1 - xgrid [0,i]) / ( cL * alpha * t_end ))**(2 * gamma / gammaA ))
            rhoE [0,i] = ( rhoL *(1+( pos1 - xgrid [0,i]) / ( cL * alpha * t_end ))**(2/ gammaA ))
            uE [0,i] = uL + (2/ gammaC )*( xgrid [0,i] - pos1 )/ t_end
            cE [0,i] = np. sqrt ( gamma * PrE [0,i]/ rhoE [0,i])
            machE [0,i] = uE [0,i]/ cE [0,i]
        elif xgrid [0,i] <= conpos :
            PrE [0,i] = p3
            rhoE [0,i] = rho2
            uE [0,i] = u2 + uR;
            cE [0,i]= np. sqrt ( gamma * PrE [0,i]/ rhoE [0,i]);
            machE [0,i] = uE [0,i]/ cE [0,i]
        elif xgrid [0,i] <= spos :
            PrE [0,i] = p3
            rhoE [0,i] = rho3
            uE [0,i] = u2+uR
            cE [0,i] = np. sqrt ( gamma * PrE [0,i]/ rhoE [0,i])
            machE [0,i] = uE [0,i]/ cE [0,i]
        else :
            PrE [0,i] = PrR
            rhoE [0,i] = rhoR ;
            uE [0,i] = uR
            cE [0,i] = np. sqrt ( gamma *PrE [0,i]/ rhoE [0,i]);
            machE [0,i] = uE [0,i]/ cE [0,i];

    entropy_E = np.log(PrE/ rhoE ** gamma );

    xS.clear()
    density.clear()
    velocity.clear()
    pressure.clear()
    intEnergy.clear()

    for i in range (0, xgrid.size ):
        xS.append(xgrid[0,i])
        density.append(rhoE [0,i])
        velocity.append(machE [0,i])
        pressure.append(PrE [0,i])
        intEnergy.append(uE [0,i])

    
fig = plt.figure()

def animate(i):
    
    riemann(i/100)

    plt.subplot(2,2,1)
    plt.subplot(2,2,1).clear()
    plt.plot(xS, density)

    #plt.subplot(2,2,2)
    #plt.subplot(2,2,2).clear()
    #plt.plot(xS, velocity)

    plt.subplot(2,2,3)
    plt.subplot(2,2,3).clear()
    plt.plot(xS, pressure)

    plt.subplot(2,2,2)
    plt.subplot(2,2,2).clear()
    plt.plot(xS, intEnergy)


ani = animation.FuncAnimation(fig, animate, frames=25, interval=500)

plt.show()