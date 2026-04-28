#!/usr/bin/env python3

import math

def eta(theta):
    # Added to prevent 0 crashes
    if theta < 1e-10:
        theta = 1e-10
    elif (math.pi - theta) < 1e-10:
        theta = math.pi - 1e-10

    return -math.log(math.tan(theta / 2.0))

def theta_region(theta):
    regs = []
    if 0.70 < theta < 2.44:
        regs.append('barrel')
    if 0.99 < theta < 2.15:
        regs.append('centbarrel')
    if (0.7 < theta < 0.99) or (2.15 < theta < 2.44):
        regs.append('transition')
    if (0.175 < theta < 0.7) or (2.44 < theta < 2.96):
        regs.append('endcap')
    if len(regs) > 0: return regs
    return None

def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2*math.pi
    while dphi < -math.pi:
        dphi += 2*math.pi
    return dphi
