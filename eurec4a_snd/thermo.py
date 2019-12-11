#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thermodynamics functions

From Chris Bretherton, Bolton 1980, etc.

Simon de Szoeke
Created on Mon Aug 27 21:30:59 2018

@author: sdeszoek
"""

import numpy as np

#  Important global constants (Emanuel, Appendix 2)
g     = 9.8
stefan= 5.67e-8
r_earth = 6370*1000
Omega   = 2*np.pi/86400

Rd    = 287.04
Rv    = 461.50
Cp    = 1005.7
Cpv   = 1870
Cw    = 4190
L     = 2.501e6
delta = 0.608 # Rv/Rd - 1
CtoK  = 273.15

# Blackbody radiative constants
kb=1.381e-23     # J/K   Boltzmann
hplanck=6.6261e-34   # J s   Planck
clight=2.998e8   # m/s   speed of light
hc=1.98645e-25   # J m   Planck*speed of light

# Usually functions describing water and water vapor depend on T in Celsius
# and functions of generalized gases, thermodynmics, and physics depend on 
# absolute T in Kelvin.

def Lv(T=0):
    #  Latent heat of vaporization as a function of T
    #  Lv(Temp[C]) = 2.501e6 + (Cpv-Cw)*(Temp)  [J/kg]
    #  From Bolton, 1980, MWR, 108, 1046-1053.
    Lvap = L + (Cpv-Cw) * T 
    return Lvap

def es(T,p=1.0e5):
    #   es(T[C],p[Pa]) [Pa] is saturation vapor pressure based on Wexler's formula,
    #   with enhancement factor for moist air rather than water vapor.
    #   The enhancement factor requires a pressure.
    #   T [degrees C], p [Pa] (note the reversed input order), es [Pa]
    #   From A. L. Buck 1981: JAM, 20, 1527-1532.
    #   SPdeS 7 July 2004
    esat = 1e2*6.1121*(1.0007 + 3.46e-8*p)*np.exp((17.502*T)/(240.97 + T))
    return esat

def qs(T,p=1.0e5):
    #  qs(p,T) is saturation specific humidity based on Wexler's formula for es
    #   with enhancement factor (see es.m).
    #   T [degrees C], p [Pa], qs [kg/kg]
    #   From A. L. Buck 1981: JAM, 20, 1527-1532.
    #   SPdeS 7 July 2004
    esat = es(T,p)
    qsat = (Rd/Rv)*esat/(p+(Rd/Rv-1)*esat)
    return qsat

# qsea = lambda T,p: 0.98 * qs(T,p)
def qsea(T,p):
    # seawater specific humidity [kg/kg]
    return 0.98*qs(T,p)

def Tdew(e,p):
    #   Tdew(e,p)
    #   e [Pa] vapor pressure, p pressure [Pa]
    #   Dewpoint Temperature [Celsius] from inversion of es formula of
    #   A. L. Buck 1981: JAM, 20, 1527-1532.
    #   SPdeS 2 August 2010
    a  = 17.502
    b  = 240.97
    k  = np.log( e/6.1121 / (1.0007+3.46e-6*p) ) # depends only on e,p
    Td = b*k / (a-k)
    return Td

def dqsdT(T,p):
    # dqsdT(T[C],p[Pa]) [(kg/kg) / K]
    # saturation specific humidity change with T from Clausius-Clapeyron
    dqsatdT = qs(T,p)*Lv(T) / ( Rv * (T+CtoK)**2 )
    return dqsatdT
    
def dqsdp(Temp,p):
    #  dqsdp(p,Temp) = d(qs)/d(Temp) = -qs(p,Temp)/(p - es(Temp,p))   [(kg/kg) Pa^{-1}]
    #  p[Pa], Temp[degree C]
    dqsatdp = -qs(Temp,p) / (p - es(Temp,p));
    return dqsatdp

def rs(T,p=1.0e5):
    # rs(p,T) is saturation mixing ratio based on Wexler's formula for es
    # with enhancement factor (see es.m).
    # p [Pa], T [degrees C], rs [kg/kg]
    # From A. L. Buck 1981: JAM, 20, 1527-1532.
    # SPdeS 7 July 2004
    esat = es(T,p)
    rsat = (Rd/Rv)*esat/(p-esat)
    return rsat

def Tlcl(Temp,p,qv):
    #   Tl=Tlcl(p[Pa], Temp[K], qv[kg/kg])
    #   Temperature at the LCL
    #   From Bolton, 1980, MWR, 108, 1046-1053.
    ev = p * qv / (Rd/Rv + qv);
    Tl = 2840/(3.5*np.log(Temp) - np.log(0.01*ev) - 4.805) + 55;
    return Tl


def Exner(p,qv=0):
    # Exner function Pi = (p/1.e5)**((Rd/Cp)*(1-0.28*qv))
    Pi = (p/1.e5)**((Rd/Cp)*(1-0.28*qv))
    return Pi

def theta(Temp,p,qv=0):
    #   Potential temperature Temp*(1e5/p)^(Rd/Cp)
    th = Temp / Exner(p,qv)
    return th

def temp(th,p,qv=0):
    # temp(theta[K],p[Pa],qv[kg/kg]=0)
    return th * Exner(p,qv)

def theta_e(Temp,p,qv):
    # theta_e(Temp[K],p[Pa],qv[kg/kg]) from eqn. 43 of Bolton, 1980, MWR, 108, 1046-1053.
    Tl = Tlcl(Temp,p,qv)
    thetae = Temp/Exner(p,qv) * np.exp((3376/Tl - 2.54) * qv * (1 + 0.81*qv))
    return thetae

def theta_es(Temp,p):
    # theta_es(Temp[K],p[Pa],qv[kg/kg]) from eqn. 43 of Bolton, 1980, MWR, 108, 1046-1053.
    return theta_e(Temp,p,qs(Temp-CtoK,p))

def Twet(T,qv,p=1.0e5,niter=5):
    # Tw=Twet(Tdry[C], qv[kg/kg], p[Pa], niter=5)
    # Computes the wet bulb temperature from an isobaric process
    # by using heat internal energy -Cp*dT=L*dq to evaporate
    # and raise qv to qs(Twb).
    #
    # Simon de Szoeke 2016-04-07

    # initialize Tw, qw
    Tw=T
    # problem for temperatures greater than 100 C!
    qw=min(qs(Tw,p),qv)

    # Steps by 1/4 of saturation deficit and enforces -Cp*dT=L*dq by successive approximations
    for iter in range(1,niter+1):
        dq = 0.25 * (qs(p,Tw) - qw)
        if dq < 0:
            break
        dT = -Lv(Tw) / Cp*dq
        qw = qw + dq
        Tw = Tw + dT
        
    return Tw

def theta_w(Temp,p,qv):
    # theta_w(Temp[K],p[Pa],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008 and 
    # eqn. 43 of Bolton, 1980, MWR, 108, 1046-1053.
    
    # Davies-Jones rational function coefficients
    a0=   7.101574
    a1= -20.68208
    a2=  16.11182
    a3=   2.574631
    a4=  -5.205688
    b1=  -3.552497
    b2=   3.781782
    b3=  -0.6899655
    b4=  -0.5929340
    
    # compute Theta_e from Bolton as first approx for thw
    th = Temp / Exner(p,qv)
    Tl = Tlcl(Temp,p,qv)
    thw = th * np.exp( (3376/Tl - 2.54) *qv * (1 + 0.81*qv) );
    
    # refine thw
    # evaluate Davies-Jones 2008 rational function fit for thw [here in K]
    if thw >= 173.15:
        X = thw / 273.15
        X2 = X  * X
        X3 = X2 * X
        X4 = X2 * X2
        thw = thw - np.exp( (a0+a1*X+a2*X2+a3*X3+a4*X4) / (1+b1*X+b2*X2+b3*X3+b4*X4) );
    
    return thw

def Tv(T,qv):
    # Tv(T[K],qv[kg/kg]) = T * (1 + delta*qv); delta = 0.608
    Tvirt = T * (1 + delta*qv)
    return Tvirt
    
def density(T,p,qv=0):
    # rho = p./(Rd*Tv(T,qv))
    rho = p / (Rd * Tv(T,qv))
    return rho


# lapse rate and sounding methods
    
def gamma(Temp,p):
    # dT/dz = gamma(T[C],p[Pa])
    # moist adiabatic lapse rate [degree C/m]
    gam = ( Lv(Temp) / Cp ) * dqsdT(Temp,p)
    return gam

def Hscale(Temp,qv=0):
    # Hscale(Temp[K],qv[kh/kh]=0) = Rd*Temp/g [m] is pressure scale height
    return Rd * Tv(Temp,qv) / g;
    
def dqsdzs(Temp,p):
    # dqsdzs(p,Temp) = dqs/dz [K/m] on a moist adiabat.
    # Temp [C], p [Pa]
    #
    # dqsatdzs = dqsdzu(Temp,p) / (1 + gamma(Temp,p))
    return dqsdzu(Temp,p) / (1 + gamma(Temp,p))

def dqsdzu(p,Temp):
    #  dqsdzu = dqs/dz [K/m] on a dry adiabat.
    # Temp [C], p [Pa]
    dqsatdzu = ((Rd*(Temp+CtoK)/Cp) * dqsdT(Temp,p) + p * dqsdp(Temp,p)) / Hscale(Temp+CtoK)
    
    return dqsatdzu


def buoyancy_flux(shf, lhf, T=25.0, p=1.0e5, q=0):
    # [bflx, bflx_s]=buoyancy_flux(shf, lhf, T, q, p)  [shf, lhf W/m^2; T °C; q kg/kg; p hPa]
    # Buoyancy flux from sensible and latent heat flux. Optionally also returns the sensible
    # part of the buoyancy flux. bflx = bflx_virt + bflx_s
    #
    # Simon de Szoeke  2018-09-01
    
    virtfac = 1.0 + delta*q
    Tvirt   = (T+CtoK) * virtfac
    rho     = density(Tvirt,p)
    wt      = shf / (rho*Cp)
    wtv     = wt * virtfac + lhf/(rho * Lv(T)) * (T+CtoK)*delta
    bflx    = g * wtv / Tvirt
    bflx_s  = g * wt  / Tvirt
    
    return bflx, bflx_s


# Planck functions
def Planck_la(la,T):
    # B_la=Planck_la(lambda[m],T[K]) [W/sr/m^3]
    # Computes Planck's blackbody spectral radiance as a function of wavelength lambda [m].
    # (c) Simon de Szoeke 2016-04-18
    
    B_la = 2 * hc * clight / ( (la**5) * ( np.exp(hc / (la*kb*T)) - 1 ) ) # W/sr/m^3
    return B_la

def Planck_nu(nu,T):
    # B_nu=Planck_nu(nu[1/s],T[K]) [W/sr/m^3]
    # Computes Planck's blackbody spectral radiance as a function of wavelength lambda [m].
    # (c) Simon de Szoeke 2016-04-18
    B_nu = 2 * hplanck * nu*nu*nu / ( clight*clight * ( np.exp( hplanck*nu / (kb*T) ) - 1 ) ) # W/sr/m^2 s
    return B_nu
