import numpy as np
from scipy.special import roots_legendre, eval_legendre
def FuncaoForma(nint,nint_1d,elemtype="Q1"):
    if elemtype=="P1":
        nsf = 3 #número de funções de forma (shape function)
    elif elemtype=="Q1":
        nsf = 4
    elif elemtype=="tri":
        nsf = 4
    elif elemtype=="Q2":
        nsf = 9
        
    
    phi = np.zeros([nsf,nint*nint])
    dphidxi = np.zeros([nsf,nint*nint])
    dphideta = np.zeros([nsf,nint*nint])
    wk = np.zeros(nint*nint)
    roots, weights = roots_legendre(nint)
    
    phixi = np.zeros([3,1])
    phieta = np.zeros([3,1])
    dphixi = np.zeros([3,1])
    dphieta = np.zeros([3,1])
    
    for j in range(0,nint):
        xi = roots[j]
        for k in range(0,nint):
            wk[j*nint+k] = weights[j]*weights[k]
            eta = roots[k]
            if elemtype=="Q1":
                phi[0,j*nint + k] = (1-xi)*(1-eta)/4
                dphidxi[0,j*nint + k] = -(1-eta)/4
                dphideta[0,j*nint + k] = -(1-xi)/4

                phi[1,j*nint + k] = (1+xi)*(1-eta)/4
                dphidxi[1,j*nint + k] = (1-eta)/4
                dphideta[1,j*nint + k] = -(1-xi)/4

                phi[2,j*nint + k] = (1+xi)*(1+eta)/4
                dphidxi[2,j*nint + k] = (1-eta)/4
                dphideta[2,j*nint + k] = (1-xi)/4

                phi[3,j*nint + k] = (1-xi)*(1+eta)/4
                dphidxi[3,j*nint + k] = -(1-eta)/4
                dphideta[3,j*nint + k] = (1-xi)/4
            
            if elemtype=="Q2":
                
                phixi[0] = 0.5*xi*(xi-1)
                phixi[1] = 0.5*xi*(xi+1)
                phixi[2] = (1-xi*xi)
                phieta[0] = 0.5*eta*(eta-1)
                phieta[1] = 0.5*eta*(eta+1)
                phieta[2] = (1-eta*eta)
                dphixi[0] = xi-0.5
                dphixi[1] = xi+0.5
                dphixi[2] = -2*xi
                dphieta[0] = eta-0.5
                dphieta[1] = eta+0.5
                dphieta[2] = -2*eta
                
                
                phi[0,j*nint + k] = phixi[0]*phieta[0]
                phi[1,j*nint + k] = phixi[1]*phieta[0]
                phi[2,j*nint + k] = phixi[1]*phieta[1]
                phi[3,j*nint + k] = phixi[0]*phieta[1]
                phi[4,j*nint + k] = phixi[2]*phieta[0]
                phi[5,j*nint + k] = phixi[1]*phieta[2]
                phi[6,j*nint + k] = phixi[2]*phieta[1]
                phi[7,j*nint + k] = phixi[0]*phieta[2]
                phi[8,j*nint + k] = phixi[2]*phieta[2]

                dphidxi[0,j*nint + k] = dphixi[0]*phieta[0]
                dphidxi[1,j*nint + k] = dphixi[1]*phieta[0]
                dphidxi[2,j*nint + k] = dphixi[1]*phieta[1]
                dphidxi[3,j*nint + k] = dphixi[0]*phieta[1]
                dphidxi[4,j*nint + k] = dphixi[2]*phieta[0]
                dphidxi[5,j*nint + k] = dphixi[1]*phieta[2]
                dphidxi[6,j*nint + k] = dphixi[2]*phieta[1]
                dphidxi[7,j*nint + k] = dphixi[0]*phieta[2]
                dphidxi[8,j*nint + k] = dphixi[2]*phieta[2]

                dphideta[0,j*nint + k] = phixi[0]*dphieta[0]
                dphideta[1,j*nint + k] = phixi[1]*dphieta[0]
                dphideta[2,j*nint + k] = phixi[1]*dphieta[1]
                dphideta[3,j*nint + k] = phixi[0]*dphieta[1]
                dphideta[4,j*nint + k] = phixi[2]*dphieta[0]
                dphideta[5,j*nint + k] = phixi[1]*dphieta[2]
                dphideta[6,j*nint + k] = phixi[2]*dphieta[1]
                dphideta[7,j*nint + k] = phixi[0]*dphieta[2]
                dphideta[8,j*nint + k] = phixi[2]*dphieta[2]

    
    roots, weights = roots_legendre(nint_1d)
    phi_1d,dphidxi_1d,wk_1d = np.zeros([nint_1d,nint_1d]),np.zeros([nint_1d,nint_1d]),np.zeros(nint_1d)
    for j in range(nint_1d):
        xi = roots[j]
        if elemtype=="Q1":
            phi_1d[0,j] = (1-xi)/2
            phi_1d[1,j] = (1+xi)/2
            wk_1d[j] = weights[j]
        if elemtype=="Q2":
            phi_1d[0,j] = 0.5*(xi-1)*xi
            phi_1d[1,j] = 0.5*(xi+1)*xi
            phi_1d[2,j] = (1-xi*xi)

            dphidxi_1d[0,j] = (xi-0.5)
            dphidxi_1d[1,j] = (xi+0.5)
            dphidxi_1d[2,j] = -2*xi
            wk_1d[j] = weights[j]

    return phi,dphidxi,dphideta,wk,phi_1d,dphidxi_1d,wk_1d

