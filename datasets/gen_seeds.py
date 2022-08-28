import sys
import random
from scipy import sparse
from sklearn import preprocessing
import numpy as np
import math
import networkx as nx

if __name__ == '__main__':
    edge_file = sys.argv[1]
    edge_seed_flag = int(sys.argv[2])

    # graph = nx.read_edgelist(edge_file, create_using=nx.Graph(), nodetype=int)
    # adj = nx.adjacency_matrix(graph, nodelist=graph.nodes())

    print("processing %s"%(edge_file))

    row=[]
    col=[]
    n=0
    with open(edge_file,'r') as fin:
        for line in fin:
            u,v = line.strip().split()
            u,v=int(u),int(v)

            if u>n:
                n=u

            if v>n:
                n=v

            row.append(u)
            col.append(v)

            row.append(v)
            col.append(u)

    n=n+1

    adj = sparse.csr_matrix(([1]*len(row), (row, col)),shape=(n,n),dtype=np.float)

    P = preprocessing.normalize(adj, norm='l1', axis=1)
    P = sparse.csr_matrix(P)
    P.eliminate_zeros()
    del adj
    #print(P)

    # node_map={}
    # i=0
    # for node in graph.nodes():
        # node_map[i]=node
        # i+=1

    n = P.shape[0]
    m = P.nnz

    # if edge_seed_flag:
        # (row,col)=P.nonzero()

    S=[]
    for i in range(5):
        if edge_seed_flag:
            x=random.randint(0,m)
            u = row[x]
            v = col[x]
        else:
            u=random.randint(0,n)
            v=random.randint(0,n)

        e = np.zeros(n)
        e[u] = 1.0
        e[v] = -1.0
        p = e.copy()
        p[u] = p[u]/P[u].nnz #graph.degree(u)
        p[v] = p[v]/P[v].nnz #graph.degree(v)
        pr = p.copy()
        for j in range(1000):
            p = P.dot(p)
            pr = pr + p

        er = e.T.dot(pr)

        #print(er)

        ui = u #node_map[u]
        vi = v #node_map[v]
        sys.stdout.write(str(u)+" "+str(v)+";")
        sys.stdout.flush()
        s = str(ui)+" "+str(vi)+" "+str(er)+"\n"
        sys.stdout.write(s)
        sys.stdout.flush()
        S.append(s)

    folder = edge_file.split('/')[0]
    if edge_seed_flag:
        seed_file = folder + "/edge_seeds.txt"
    else:
        seed_file = folder + "/rand_seeds.txt"

    with open(seed_file, 'w') as f:
        for s in S:
            f.write(s)
