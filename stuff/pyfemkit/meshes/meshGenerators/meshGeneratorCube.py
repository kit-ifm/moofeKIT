import numpy as np

def deleteNodes(EDOF,NODES,delNodes): #--> noch zu erledigen
    #deletes all nodes from given in delNodes
    delNodes = np.unique(delNodes)

    #delete all the node --> nargout ist bei python sichergestellt

    for i in range(len(delNodes)):
        #delete nodes
        actDupl = delNodes[i]
        #EDOF(EDOF>=actDul) = EDOF(EDOF>=actDul)-1
        delNodes[i:] = delNodes[i:]-1

    EDOF = EDOF.astype(int)
    return EDOF,Nodes

def meshGeneratorCube(lengthX, lengthY, lengthZ, nelX,nelY,nelZ,order,serendipity,*args):
    
    assert lengthX > 0, 'lengthX must be positive'
    assert lengthY > 0, 'lengthY must be positive'
    assert lengthZ > 0, 'lengthZ must be positive'
    assert (isinstance(nelX,int) and  nelX> 0), 'nelx must be an integer and >0'
    assert (isinstance(nelY,int) and  nelX> 0), 'nely must be an integer and >0'
    assert (isinstance(nelZ,int) and  nelX> 0), 'nelz must be an integer and >0'

    #set standard arguments --> unnötig?
    if (order==1 and serendipity ==True):
        print('setting serendipity=true has no effect if order is 1')

    #NODES
    #number of nodes
    nel = nelX*nelY*nelZ
    nnoX = nelX*order+1
    nnoY = nelY*order+1
    nnoZ = nelZ*order+1

    #geometry
    geometry = [[-lengthX/2,lengthX/2],[-lengthY/2,lengthY/2],[-lengthZ/2,lengthZ/2]]
    nodesXDir = np.array([np.linspace(geometry[0][0],geometry[0][1],int(nnoX))]).T
    nodesYDir = np.array([np.linspace(geometry[1][0],geometry[1][1],nnoY)]).T
    nodesZDir = np.array([np.linspace(geometry[2][0],geometry[2][1],nnoZ)]).T
    zw = np.concatenate((np.kron(np.ones((nnoY,1)),nodesXDir),np.kron(nodesYDir,np.ones((nnoX,1)))),axis=1) #nodes2D
    NODES = np.concatenate((np.kron(np.ones((nnoZ,1)),zw),np.kron(nodesZDir,np.ones((zw.shape[0],1)))),axis=1)

    #EDOF
    #e.dof 2D
    if order == 1:
        edof2D=np.zeros((nelX*nelY,4))
        a=np.kron(np.ones((nelY,1)),np.arange(1,nnoX,1).reshape((nnoX-1,1)))+np.kron(np.arange(0,nelX**2,nnoX).reshape((nelX,1))+0*nnoX,np.ones((nelY,1)))
        b=np.kron(np.ones((nelY,1)),np.arange(2,nnoX+1,1).reshape((nnoX-1,1)))+np.kron(np.arange(0,nelX**2,nnoX).reshape((nelX,1))+0*nnoX,np.ones((nelY,1)))
        c=np.kron(np.ones((nelY,1)),np.arange(2,nnoX+1,1).reshape((nnoX-1,1)))+np.kron(np.arange(0,nelX**2,nnoX).reshape((nelX,1))+1*nnoX,np.ones((nelY,1)))
        d=np.kron(np.ones((nelY,1)),np.arange(1,nnoX,1).reshape((nnoX-1,1)))+np.kron(np.arange(0,nelX**2,nnoX).reshape((nelX,1))+1*nnoX,np.ones((nelY,1)))

        i=0
        while i < len(a):
            edof2D[i,0]=a[i]
            edof2D[i,1]=b[i]
            edof2D[i,2]=c[i]
            edof2D[i,3]=d[i]
            i+=1
        EDOF=np.zeros((nel,8),dtype=int)
        EDOF=np.concatenate((edof2D,edof2D+nnoX*nnoY),axis=1)
        EDOF=np.kron(np.ones((nelZ,1)),EDOF)+np.kron(np.arange(0,nnoX-1,1).reshape((nnoX-1,1)),1*nnoX*nnoY*np.ones((nelX*nelY,8)))  
      
    elif order==2: #funktioniert
        edof2D=np.zeros((nelY*nelY,9))    
        a = np.kron(np.ones(nelY),np.arange(1,nnoX,2)).reshape((nelX*nelY,1))+ np.kron((np.arange(0,(nnoX-1)*(nnoY-1),nnoX*2) + 0*nnoX),np.ones(nelX)).reshape((nelY*nelY,1))
        b = np.kron(np.ones(nelY),np.arange(3,nnoX+1,2)).reshape((nelX*nelY,1))+ np.kron((np.arange(0,(nnoX-1)*(nnoY-1),nnoX*2) + 0*nnoX),np.ones(nelX)).reshape((nelY*nelY,1))
        c = np.kron(np.ones(nelY),np.arange(3,nnoX+1,2)).reshape((nelX*nelY,1))+ np.kron((np.arange(0,(nnoX-1)*(nnoY-1),nnoX*2) + 2*nnoX),np.ones(nelX)).reshape((nelY*nelY,1))
        d = np.kron(np.ones(nelY),np.arange(1,nnoX,2)).reshape((nelX*nelY,1))+ np.kron((np.arange(0,(nnoX-1)*(nnoY-1),nnoX*2) + 2*nnoX),np.ones(nelX)).reshape((nelY*nelY,1))
        e = np.kron(np.ones(nelY),np.arange(2,nnoX,2)).reshape((nelX*nelY,1))+ np.kron((np.arange(0,(nnoX-1)*(nnoY-1),nnoX*2) + 0*nnoX),np.ones(nelX)).reshape((nelY*nelY,1))
        f = np.kron(np.ones(nelY),np.arange(3,nnoX+1,2)).reshape((nelX*nelY,1))+ np.kron((np.arange(0,(nnoX-1)*(nnoY-1),nnoX*2) + 1*nnoX),np.ones(nelX)).reshape((nelY*nelY,1))
        g = np.kron(np.ones(nelY),np.arange(2,nnoX,2)).reshape((nelX*nelY,1))+ np.kron((np.arange(0,(nnoX-1)*(nnoY-1),nnoX*2) + 2*nnoX),np.ones(nelX)).reshape((nelY*nelY,1))
        h = np.kron(np.ones(nelY),np.arange(1,nnoX,2)).reshape((nelX*nelY,1))+ np.kron((np.arange(0,(nnoX-1)*(nnoY-1),nnoX*2) + 1*nnoX),np.ones(nelX)).reshape((nelY*nelY,1))
        j = np.kron(np.ones(nelY),np.arange(2,nnoX,2)).reshape((nelX*nelY,1))+ np.kron((np.arange(0,(nnoX-1)*(nnoY-1),nnoX*2) + 1*nnoX),np.ones(nelX)).reshape((nelY*nelY,1))
        
        i=0
        while i < len(a):
            edof2D[i,0]=a[i]
            edof2D[i,1]=b[i]
            edof2D[i,2]=c[i]
            edof2D[i,3]=d[i]
            edof2D[i,4]=e[i]
            edof2D[i,5]=f[i]
            edof2D[i,6]=g[i]
            edof2D[i,7]=h[i]
            edof2D[i,8]=j[i]           
            i+=1
        #edof 3D unstructured for one layer (only 2d structuring applied)
        edofZw = np.zeros((nelX*nelY,27),dtype=int)
        edofZw = np.concatenate((edof2D,edof2D+nnoX*nnoY,edof2D+2*nnoX*nnoY),axis=1) # passt auch
        
        #structure first layer with node numbering of esra
        #1-8 corner nodes
        #9-16 midside nodes
        #17 midplane node (irregular)
        #18-21 midside nodes
        #22-26 midplane nodes
        #27 center nonlocal
        #EDOF = np.zeros((nel,27),dtype=int)
        EDOF=np.concatenate((edofZw[:,[0,1,2,3]],edofZw[:,[18,19,20,21]],edofZw[:,[9,10,11,12,4,5,6,7,8,22,23,24,25,26,13,14,15,16,17]]),axis=1)#,[18]]],axis=1)
        EDOF=np.concatenate((EDOF,np.zeros((nel-nelX*nelY,27))),axis=0) #funktioniert
        #complete EDOF --> funktioniert
        EDOF=np.concatenate((EDOF[0:nelX*nelY,:],EDOF[0:nelX*nelY,:],EDOF[0:nelX*nelY,:]),axis=0)+np.kron(np.arange(0,nelZ,1).reshape((nelZ,1)),order*nnoX*nnoY*np.ones((nelX*nelY,27)))
    else:
        print('order not implemented')

    #get the elements that touch the border, then select the nodes on the boundary in the correct order
    if order==1:
        
        elementSetSX = np.kron(np.ones(nelZ),np.arange(1,nelX*nelY-nelX+2,nelX)).reshape((nelZ*nelY,1)) + np.kron(np.arange(0,nelZ,1),np.ones(nelY)).reshape((nelZ*nelY,1))*nelX*nelY-1
        elementSetSY = np.kron(np.ones(nelZ),np.arange(1,nelX+1,1)).reshape((nelX*nelY,1)) + np.kron(np.arange(0,nelZ,1),np.ones(nelX)).reshape((nelZ*nelY,1))*nelX*nelY-1
        elementSetSZ = np.arange(1,nelX*nelY+1,1).reshape((nelX*nelY,1))-1
        bounEDOFsSX1 = []
        bounEDOFsSX2 = []
        bounEDOFsSY1 = []
        bounEDOFsSY2 = []
        bounEDOFsSZ1 = []
        bounEDOFsSZ2 = []
        xcoordSX1=np.array([3,0,4,7])
        xcoordSX2=np.array([1,2,6,5])
        xcoordSY1=np.array([0,1,5,4])
        xcoordSY2=np.array([2,3,7,6])
        xcoordSZ1=np.array([3,2,1,0])
        xcoordSZ2=np.array([4,5,6,7])
        i=0
        while (i<nelX*nelY):
            j=0
            zeileSX1=[]
            zeileSX2=[]
            zeileSY1=[]
            zeileSY2=[]
            zeileSZ1=[]
            zeileSZ2=[]
            while(j<len(xcoordSX1)):
                zeileSX1.append(EDOF[int(elementSetSX[i,0])][int(xcoordSX1[j])])
                zeileSX2.append(EDOF[int(elementSetSX[i,0])+nelX-1][int(xcoordSX2[j])])
                zeileSY1.append(EDOF[int(elementSetSY[i,0])][int(xcoordSY1[j])])
                zeileSY2.append(EDOF[int(elementSetSY[i,0])+nelX*nelY-nelX][int(xcoordSY2[j])])
                zeileSZ1.append(EDOF[int(elementSetSZ[i,0])][int(xcoordSZ1[j])])
                zeileSZ2.append(EDOF[int(elementSetSZ[i,0])+nelX*nelY*nelZ-nelX*nelY][int(xcoordSZ2[j])])
                j+=1
            bounEDOFsSX1.extend([zeileSX1])
            bounEDOFsSX2.extend([zeileSX2])
            bounEDOFsSY1.extend([zeileSY1])
            bounEDOFsSY2.extend([zeileSY2])
            bounEDOFsSZ1.extend([zeileSZ1])
            bounEDOFsSZ2.extend([zeileSZ2])
            i+=1
        bounEDOFsSX1=np.array(bounEDOFsSX1)
        bounEDOFsSX2=np.array(bounEDOFsSX2)
        bounEDOFsSY1=np.array(bounEDOFsSY1)
        bounEDOFsSY2=np.array(bounEDOFsSY2)
        bounEDOFsSZ1=np.array(bounEDOFsSZ1)
        bounEDOFsSZ2=np.array(bounEDOFsSZ2)
        
        bounEDOFs={'bounEDOFsSX1':bounEDOFsSX1,
                   'bounEDOFsSX2':bounEDOFsSX2,
                   'bounEDOFsSY1':bounEDOFsSY1,
                   'bounEDOFsSY2':bounEDOFsSY2,
                   'bounEDOFsSZ1':bounEDOFsSZ1,
                   'bounEDOFsSZ2':bounEDOFsSZ2,}
        #print(bounEDOFs['bounEDOFsSX1'])

    elif order==2: #funktioniert
        elementSetSX = np.kron(np.ones(nelZ),np.arange(1,nelX*nelY-nelX+2,nelX)) + np.kron(np.arange(0,nelZ,1),np.ones(nelY))*nelX*nelY-1
        elementSetSY = np.kron(np.ones(nelZ),np.arange(1,nelX+1,1)) + np.kron(np.arange(0,nelZ,1),np.ones(nelX))*nelX*nelY-1
        elementSetSZ = np.arange(1,nelX*nelY+1,1)-1
        bounEDOFsSX1 = []
        bounEDOFsSX2 = []
        bounEDOFsSY1 = []
        bounEDOFsSY2 = []
        bounEDOFsSZ1 = []
        bounEDOFsSZ2 = []
        xcoordSX1=np.array([3,0,4,7,15,8,20,11,25])
        xcoordSX2=np.array([1,2,6,5,13,10,18,9,23])
        xcoordSY1=np.array([0,1,5,4,12,9,17,8,22])
        xcoordSY2=np.array([2,3,7,6,14,11,19,10,24])
        xcoordSZ1=np.array([3,2,1,0,14,13,12,15,16])
        xcoordSZ2=np.array([4,5,6,7,17,18,19,20,21])
        i=0
        while (i<nelX*nelY):
            j=0
            zeileSX1=[]
            zeileSX2=[]
            zeileSY1=[]
            zeileSY2=[]
            zeileSZ1=[]
            zeileSZ2=[]
            while(j<len(xcoordSX1)):
                zeileSX1.append(EDOF[int(elementSetSX[i])][int(xcoordSX1[j])])
                zeileSX2.append(EDOF[int(elementSetSX[i])+nelX-1][int(xcoordSX2[j])]) #
                zeileSY1.append(EDOF[int(elementSetSY[i])][int(xcoordSY1[j])])
                zeileSY2.append(EDOF[int(elementSetSY[i])+nelX*nelY-nelX][int(xcoordSY2[j])]) #
                zeileSZ1.append(EDOF[int(elementSetSZ[i])][int(xcoordSZ1[j])])
                zeileSZ2.append(EDOF[int(elementSetSZ[i])+nelX*nelY*nelZ-nelX*nelY][int(xcoordSZ2[j])]) #
                j+=1
            bounEDOFsSX1.extend([zeileSX1])
            bounEDOFsSX2.extend([zeileSX2])
            bounEDOFsSY1.extend([zeileSY1])
            bounEDOFsSY2.extend([zeileSY2])
            bounEDOFsSZ1.extend([zeileSZ1])
            bounEDOFsSZ2.extend([zeileSZ2])
            i+=1
        bounEDOFsSX1=np.array(bounEDOFsSX1)
        bounEDOFsSX2=np.array(bounEDOFsSX2)
        bounEDOFsSY1=np.array(bounEDOFsSY1)
        bounEDOFsSY2=np.array(bounEDOFsSY2)
        bounEDOFsSZ1=np.array(bounEDOFsSZ1)
        bounEDOFsSZ2=np.array(bounEDOFsSZ2)

        bounEDOFs={'SX1':bounEDOFsSX1,
                   'SX2':bounEDOFsSX2,
                   'SY1':bounEDOFsSY1,
                   'SY2':bounEDOFsSY2,
                   'SZ1':bounEDOFsSZ1,
                   'SZ2':bounEDOFsSZ2,}
    else:
        print('order not implemented')

    #change mesh for serendipity 
    if (order==2 and serendipity):
        #reduce edofs--> funktioniert
        EDOF= np.concatenate((np.array([i[:8] for i in EDOF]),np.array([i[12:16] for i in EDOF]),np.array([i[17:21] for i in EDOF]),np.array([i[8:12] for i in EDOF])),axis=1)
        bounEDOFs['SX1'] = np.array([i[:1] for i in bounEDOFs['SX1']]).reshape((1,nelX*nelY))
        bounEDOFs['SX2'] = np.array([i[:1] for i in bounEDOFs['SX2']]).reshape((1,nelX*nelY))
        bounEDOFs['SY1'] = np.array([i[:1] for i in bounEDOFs['SY1']]).reshape((1,nelX*nelY))
        bounEDOFs['SY2'] = np.array([i[:1] for i in bounEDOFs['SY2']]).reshape((1,nelX*nelY))
        bounEDOFs['SZ1'] = np.array([i[:1] for i in bounEDOFs['SZ1']]).reshape((1,nelX*nelY))
        bounEDOFs['SZ2'] = np.array([i[:1] for i in bounEDOFs['SZ2']]).reshape((1,nelX*nelY))
        
        #find surplus nodes --> funktioniert
        delNodes = np.arange(1,len(NODES)+1).reshape((len(NODES),1))
        delNodes= np.delete(delNodes,list(map(int,np.unique(EDOF)-1)))
        
        #delete surplus nodes --> noch zu erledigen
        EDOF,NODES = deleteNodes(EDOF,NODES,delNodes)
        """
        bounEDOFs.SX1 = deleteNodes(bounEDOFs.SX1,[],delNodes)
        bounEDOFs.SX2 = deleteNodes(bounEDOFs.SX2,[],delNodes)
        bounEDOFs.SY1 = deleteNodes(bounEDOFs.SY1,[],delNodes)
        bounEDOFs.SY2 = deleteNodes(bounEDOFs.SY2,[],delNodes)
        bounEDOFs.SZ1 = deleteNodes(bounEDOFs.SZ1,[],delNodes)
        bounEDOFs.SZ2 = deleteNodes(bounEDOFs.SZ2,[],delNodes)
        """
    EDOF = EDOF.astype(int)
    return NODES,EDOF,bounEDOFs   

        
    
