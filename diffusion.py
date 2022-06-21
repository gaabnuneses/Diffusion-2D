import numpy as np
from scipy.sparse import lil_matrix


def matriz_elementar(dphidxi,dphideta,wk,nint,elem,nodes,gamma,elemType):
    # Inicializações
    if elemType =="Q1":
        nl = 4
    elif elemType=="Q2":
        nl = 9
    
    mate = np.zeros([nl,nl])
    vete = np.zeros([nl])
    # Retomada dos nós do elemento
    
    if elemType=="Q1":
        no1 = elem[2]
        no2 = elem[3]
        no3 = elem[4]
        no4 = elem[5]
        # Obtenção das coordenadas dos pontos
        X = np.array([nodes[no1][1],nodes[no2][1],nodes[no3][1],nodes[no4][1]])    
        Y = np.array([nodes[no1][2],nodes[no2][2],nodes[no3][2],nodes[no4][2]])
    elif elemType=="Q2":
        no1 = elem[2]
        no2 = elem[3]
        no3 = elem[4]
        no4 = elem[5]
        no5 = elem[6]
        no6 = elem[7]
        no7 = elem[8]
        no8 = elem[9]
        no9 = elem[10]
        # Obtenção das coordenadas dos pontos
        X = np.array([nodes[no1][1],nodes[no2][1],nodes[no3][1],nodes[no4][1],nodes[no5][1],nodes[no6][1],nodes[no7][1],nodes[no8][1],nodes[no9][1]])    
        Y = np.array([nodes[no1][2],nodes[no2][2],nodes[no3][2],nodes[no4][2],nodes[no5][2],nodes[no6][2],nodes[no7][2],nodes[no8][2],nodes[no9][2]])
    
    
    # Cálculo das derivadas de x e y em xi e eta
    dxdxi = np.matmul(np.transpose(X),dphidxi)
    dxdeta = np.matmul(np.transpose(X),dphideta)
    dydxi = np.matmul(np.transpose(Y),dphidxi)
    dydeta = np.matmul(np.transpose(Y),dphideta)
    # Cálculo da matriz elementar
    for k in range(0,nint*nint):             # Laço sobre os pontos de integração
        J = np.matrix([[dxdxi[k],dydxi[k]],[dxdeta[k],dydeta[k]]])
        Jac = J[0,0]*J[1,1] - J[0,1]*J[1,0]
        Jinv = np.matrix([[J[1,1]/Jac,-J[0,1]/Jac],[-J[1,0]/Jac,J[0,0]/Jac]])

        for i in range(0,nl):           # Laço sobre os nós (linha)
            for j in range(0,nl):       # Laço sobre os nós (coluna)
                dphidx2 = (Jinv[0,0]*dphidxi[i,k] + Jinv[0,1]*dphideta[i,k])*(Jinv[0,0]*dphidxi[j,k] + Jinv[0,1]*dphideta[j,k])
                dphidy2 = (Jinv[1,0]*dphidxi[i,k] + Jinv[1,1]*dphideta[i,k])*(Jinv[1,0]*dphidxi[j,k] + Jinv[1,1]*dphideta[j,k])
                mate[i,j] += gamma*(dphidx2+dphidy2)*wk[k]*Jac
    return mate,vete

def matriz_elementar_1d(phi,dphidxi,wk,nint,elem,nodes,elemType):
    # Inicializações
    if elemType=="Q1":
        nl = 2
    if elemType=="Q2":
        nl = 3
    mate = np.zeros([nl,nl])
    vete = np.zeros([nl])
    # Retomada dos nós do elemento
    L = np.zeros(nl)
    if elemType=="Q1":
        no1 = elem[3]
        no2 = elem[4]
        X = np.array([nodes[no1][1],nodes[no2][1]])    
        Y = np.array([nodes[no1][2],nodes[no2][2]])
        L[0] = np.sqrt((X[1]-X[0])**2+(Y[1]-Y[0])**2)
        L[1] = np.sqrt((X[1]-X[0])**2+(Y[1]-Y[0])**2)
    elif elemType=="Q2":
        no1 = elem[3]
        no2 = elem[4]
        no3 = elem[5]
        X = np.array([nodes[no1][1],nodes[no2][1],nodes[no3][1]])    
        Y = np.array([nodes[no1][2],nodes[no2][2],nodes[no3][2]])
        dxdxi = np.transpose(X)@dphidxi
        dydxi = np.transpose(Y)@dphidxi
        for k in range(3):
            L[k] = 2*np.sqrt((dxdxi[k])**2+(dydxi[k])**2)
    
    phy = elem[2]
    h=0
    Tinf=0
    if phy == 10:
        h = 1000
        Tinf = 1700
    elif phy == 9:
        h = 200
        Tinf = 400
    else:
        return mate,vete


    for k in range(nint):             # Laço sobre os pontos de integração
        for i in range(nl):           # Laço sobre os nós (linha)
            for j in range(nl):       # Laço sobre os nós (coluna)
                mate[i,j] += h*phi[i,k]*phi[j,k]*wk[k]*L[k]/2
            vete[i] += h*phi[i,k]*wk[k]*L[k]/2*Tinf
    return mate,vete

def montar_sistema(phi,dphidxi,dphideta,wk,phi_1d,dphidxi_1d,wk_1d,nint,nint_1d,elem,bound_elem,elemType,nodes,gamma):
    ngnode = len(nodes) # Número de nós no problema
    if elemType =="Q1":
        nldof = 4 # Número de nós por elemento
    if elemType =="Q2":
        nldof = 9 # Número de nós por elemento
    # Vetor b
    b = np.zeros([ngnode])
    A = lil_matrix((ngnode,ngnode))
    
    # CONSTRUÇÃO DA PARTE DIFUSIVA DO DOMÍNIO
    for nel in range(len(elem)):
        # Parte ∇T∇w
        mate,vete = matriz_elementar(dphidxi,dphideta,wk,nint,elem[nel],nodes,gamma,elemType)
        iq = []
        for i in range(nldof):
            iq.append(elem[nel][2+i])
        
        m = 0
        for i in iq:
            n = 0
            for j in iq:
                A[i,j] = A[i,j] + mate[m,n]
                n+=1
            b[i] = b[i] + vete[m]
            m+=1

    # CONSTRUÇÃO DAS CONDIÇÕES DE CONTORNO - NEUMANN
    ngnode = len(nodes) # Número de nós no problema
    if elemType =="Q1":
        nldof = 2 # Número de nós por elemento
    if elemType =="Q2":
        nldof = 3 # Número de nós por elemento
    for nel in range(len(bound_elem)):
        # Parte w∇Tn
        mate,vete = matriz_elementar_1d(phi_1d,dphidxi_1d,wk_1d,nint_1d,bound_elem[nel],nodes,elemType)
        iq = []
        for i in range(nldof):
            iq.append(bound_elem[nel][3+i])
        m = 0
        for i in iq:
            n = 0
            for j in iq:
                A[i,j] = A[i,j] + mate[m,n]
                n+=1
            b[i] = b[i] + vete[m]
            m+=1
    
    # CONSTRUÇÃO DAS CONDIÇÕES DE CONTORNO - DIRICHLET
    for nel in range(len(bound_elem)):
        phys = bound_elem[nel][2]   # Condição de contorno
        #if phys == 5 or phys == 6 or phys == 7 or phys == 8:
        if phys == 70 or phys == 80 or phys == 90 or phys == 100:
            iq = []
            for i in range(nldof):
                iq.append(bound_elem[nel][3+i])
            for i in iq:
                A[i,:] = 0
                A[i,i] = 1 
        #        if phys == 5 or phys == 6 or phys == 7: # Adiabáticas
                if phys == 80 or phys == 90 or phys == 100: # Adiabáticas
                    b[i] = 1
        #        if phys == 8:   # T(x,0) = 1
                if phys == 70:   # T(x,0) = 1
                    b[i] = 0
                

    return A,b
