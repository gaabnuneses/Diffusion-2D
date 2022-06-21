##------------------------------- Bibliotecas utilizadas -------------------------------#
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.special import roots_legendre, eval_legendre
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
from scipy.interpolate import griddata
##--------------------------- Funções em outros Arquivos -------------------------------#
from loadMeshes import *
from diffusion import *
from plotmalha import *
from funcaoforma import *

##------------------------------- Execução do programa -------------------------------#
# Arquivo .txt com dados da malha gerada no gmsh Versão 2 gmsh
# numel can vary in this list: [12,48,300,1200,2700,4800,10800]
numel = 1800
elemType = "Q1"

#nmalhas = [50,200,450,800,1800,3200,5000,8450,12800]
#for numel in nmalhas:
if elemType=="Q1":
    nint = 2
    nint_1d = 2
    file = "Malhas/PaTurbina"+str(numel)+".txt"
elif elemType=="Q2":
    nint = 3
    nint_1d = 3
    file = "Malhas/PaTurbina"+str(numel)+"-Q2.txt"

# Carregamento da Malha
nodes, bound_elem, elem = LoadMesh(file,1e-3)

# Mostrando a malha
plotMalha(nodes,elem,False)

# Cálculo das funções de forma
phi,dphidxi,dphideta,wk,phi_1d,dphidxi_1d,wk_1d = FuncaoForma(nint,nint_1d,elemType)

# Propriedades Físicas
gamma = 25

# Cálculo da matriz global
t = time.time()
A,b = montar_sistema(phi,dphidxi,dphideta,wk,phi_1d,dphidxi_1d,wk_1d,nint,nint_1d,elem,bound_elem,elemType,nodes,gamma)
# Solução do sistema
A = A.tocsr()
T = spsolve(A, b)
elapsed = time.time() - t
print("Tempo decorrido da montagem à solução: "+str(round(elapsed,ndigits=3))+"s")
    
def valMaxMin():
    print(val_max)
    print(val_min)
    erro5000_max = []
    erro5000_min = []
    for i in range(len(val_max)):
        erro5000_max.append(100*np.abs(val_max[i]-val_max[len(val_max)-1])/val_max[len(val_max)-1])
        erro5000_min.append(100*np.abs(val_min[i]-val_min[len(val_min)-1])/val_min[len(val_min)-1])
    f = plt.figure()
    plt.plot(nmalhas,erro5000_max,marker='o',label="Valor Máximo")
    plt.plot(nmalhas,erro5000_min,marker='o',label="Valor Mínimo")
    plt.xlabel("Número de Elementos")
    plt.ylabel("Diferença absoluta relativa com o mais refinado [%]")
    plt.title("Refinamento de Malha para Elementos do tipo Q2")
    plt.legend()
    plt.show()

def plot_campos():
    # Mostrando o Resultado na tela
    x = []
    y = []
    for i in range(len(nodes)):
        x.append(nodes[i][1])
        y.append(nodes[i][2])

    f = plt.figure()
    f.set_figwidth(10)
    f.set_figheight(6)    
    plt.tricontourf(x,y,T,15,cmap='jet')
    plt.colorbar()

    xq = np.linspace(0,.005,30)
    yq = np.linspace(0,.003,30)
    Xq,Yq = np.meshgrid(xq,yq)

    Tq = griddata((x, y), T, (Xq, Yq))

    qx,qy=np.gradient(Tq,xq,yq)

    # Ordem inversa de qx e qy por conta da alocação na matriz gerada no meshgrid

    plt.fill([.002,.002,.005,.005],[0,.001,.001,0],color='white')
    plt.quiver(xq,yq,-qy,-qx)
    
    plt.show()

#valMaxMin()
plot_campos()