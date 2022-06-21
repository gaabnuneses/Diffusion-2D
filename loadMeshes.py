import numpy as np

# Função que carrega malha e separa informações necessárias
def LoadMesh(file,factor):
    # Inicialização da Leitura e Armazenamento de Dados
    lines = open(file).readlines()
    nodes = []
    elem = []
    num_nodes = 0
    num_elem = 0
    phys = []
    # Laço sobre as linhas do arquivo
    for i in range(0,len(lines)):
        # Armazenamento dos nós
        if "$PhysicalNames" in lines[i]:
            num_phys = int(lines[i+1])
            m=0
            for j in range(i+2,i+2+num_phys):
                phys.append(lines[j])

        if "$Nodes" in lines[i]:
            num_nodes = int(lines[i+1])
            m=0
            for j in range(i+2,i+2+num_nodes):
                aux = lines[j].strip().split(' ')
                aux2 = [m]
                for k in aux[1:]:
                    aux2.append(float(k)*factor)
                nodes.append(aux2)
                m+=1

        # Armazenamento dos elementos
        if "$Elements" in lines[i]:
            num_elem = int(lines[i+1])
            for j in range(i+2,i+2+num_elem):
                aux = lines[j].strip().split(' ')
                aux2 = []
                for k in aux:
                    aux2.append(int(k))
                elem.append(aux2)
    # Separação dos elementos por tipo - contorno/domínio
    BOUND = []          # Lista de Elementos de Contorno
    ELEM = []           # Lista de Elementos de Domínio
    # contadores:
    counter_elem = 0    
    counter_bound = 0
    for i in range(0,len(elem)):
        if elem[i][1] == 1: # Tipo: Contorno/Linear
            BOUND.append([counter_bound,1,elem[i][3],elem[i][5]-1,elem[i][6]-1])
            counter_bound += 1
        if elem[i][1] == 8: # Tipo: Contorno/Quadrático
            BOUND.append([counter_bound,8,elem[i][3],elem[i][5]-1,elem[i][6]-1,elem[i][7]-1])
            counter_bound += 1
        if elem[i][1] == 2: # Tipo: Domínio/Triangular P1
            ELEM.append([counter_elem,2,elem[i][5]-1,elem[i][6]-1,elem[i][7]-1] )
            counter_elem += 1 
        if elem[i][1] == 3: # Tipo: Domínio/Quadrilateral Q1
            ELEM.append([counter_elem,3,elem[i][5]-1,elem[i][6]-1,elem[i][7]-1,elem[i][8]-1] )
            counter_elem += 1 
        if elem[i][1] == 10: # Tipo: Domínio/Quadrilateral Q2
            aux = [counter_elem,10]
            for k in range(5,14):
                aux.append(elem[i][k]-1)
            ELEM.append(aux)
            counter_elem += 1 
    
    # Exibição de dados no Terminal
    print("Phyiscal Groups:")
    for i in phys:
        print(i[0:(len(i)-1)])
    print("Número de nós = "+str(len(nodes)))
    print("Número de elementos de domínio = "+str(len(ELEM)))
    print("Número de elementos de contorno = "+str(len(BOUND)))
    return nodes,   BOUND,  ELEM
