import matplotlib.pyplot as plt
import numpy as np
from random import random

def plotMalha(nodes, elem,showIndex=False):
    f = plt.figure()
    #f.set_figwidth(10)
    #f.set_figheight(6)
    
    for i in range(0,len(elem)):
        x = []
        y = []
        for k in range(2,len(elem[i])):
            x.append(nodes[elem[i][k]][1])
            y.append(nodes[elem[i][k]][2])
        xc = np.mean(x)
        yc = np.mean(y)
        x.append(nodes[elem[i][2]][1])
        y.append(nodes[elem[i][2]][2])
        #(0.1, 0.2, 0.5)
        plt.fill(x,y,color='k',alpha=random()*.6+.1, linewidth=1)
        if showIndex == True:
            plt.annotate(str(i),(xc,yc),color='k',weight='bold', ha='center', va='center',size=8)
    if showIndex == True:
        for i in range(0,len(nodes)):
            plt.annotate(str(i),(nodes[i][1],nodes[i][2]),size=6)
    #plt.xlabel("x")
    #plt.ylabel("y")
    plt.axis('off')
    p = plt.title("Malha com "+str(len(elem))+" elementos")
    #plt.show()
    return p

def plotContorno(nodes, bound_elem):
    f = plt.figure()
    f.set_figwidth(10)
    f.set_figheight(6)
    for i in range(len(bound_elem)):
        x = []
        y = []
        for k in range(3,len(bound_elem[i])):
            x.append(nodes[bound_elem[i][k]][1])
            y.append(nodes[bound_elem[i][k]][2])
        x.append(nodes[bound_elem[i][3]][1])
        y.append(nodes[bound_elem[i][3]][2])

        plt.plot(x,y,color=(0.1, 0.2, 0.5),marker='o', linewidth=1)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()