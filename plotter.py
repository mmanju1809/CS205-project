import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
%matplotlib inline
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import seaborn as sns
sns.set_style('ticks',{'xtick.direction': u'in','ytick.direction': u'in'})

files = [['BFS_3.txt']]
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
current_palette = sns.color_palette()
col_map = plt.get_cmap('Oranges')

P = np.zeros((60248,5))##19993
for t in range(1):
    sys=open(files[t][0]);
    holder = sys.readline()
    holder = sys.readline()
    for j in range(60248):
        sys_line = sys.readline()
        sys_line = [float(k) for k in sys_line.split()]
        for q in range(5):
            P[j,q] = sys_line[q]
    sys.close()
I = np.argsort(P[:,3])
B = P[I,:]
Q = np.zeros((60248,4))
Q2 = np.zeros((60248,4))
done1 = []
done2 = []
s = 1
s2 = 1
Q[0,:] = np.array([B[0,0],B[0,1],B[0,3],B[0,4]])
Q2[0,:] = np.array([B[0,0],B[0,2],B[0,3],B[0,4]])

for kk in range(1,len(B[:,0])):
    t = B[kk,3]
    BB = B[np.where(abs(B[:,3] - t)<0.00001),:][0]
    CC = BB[np.where(abs(BB[:,0] - B[kk,0])<0.0001),:][0]
    DD = CC[np.where(abs(CC[:,1] - B[kk,1])<0.0001),:][0]
    EE = CC[np.where(abs(CC[:,2] - B[kk,2])<0.0001),:][0]
    if [B[kk,0],B[kk,1],B[kk,3]] not in done1:
        Q[s,:] = np.array([B[kk,0],B[kk,1],B[kk,3],np.sum(DD[:,4])])
        done1 += [[B[kk,0],B[kk,1],B[kk,3]]]   
        s = s+1
    if [B[kk,0],B[kk,2],B[kk,3]] not in done2:
        Q2[s2,:] = np.array([B[kk,0],B[kk,2],B[kk,3],np.sum(EE[:,4])])
        done2 += [[B[kk,0],B[kk,2],B[kk,3]]]   
        s2 = s2+1
Q = Q[np.where(Q[:,3] > 0.0),:][0]
Q2 = Q2[np.where(Q2[:,3] > 0.0),:][0] 

U = np.zeros((len(Q[:,0]),4))
U2 = np.zeros((len(Q2[:,0]),4))
done1 = []
done2 = []
s = 1
s2 = 1
U[0,:] = np.array([Q[0,0],Q[0,1],Q[0,2],Q[0,3]])
U2[0,:] = np.array([Q2[0,0],Q2[0,1],Q2[0,2],Q[0,3]])

for ll in range(1,7):
    BB = Q[np.where((Q[:,2] > (ll-1)*0.001905) * (Q[:,2] <=ll*0.001905)),:][0]
    for kk in range(len(BB[:,0])):
        CC = BB[np.where(abs(BB[:,0] - BB[kk,0])<0.00001),:][0]
        DD = CC[np.where(abs(CC[:,1] - BB[kk,1])<0.00001),:][0]
        if [BB[kk,0],B[kk,1],np.mean(BB[:,2])] not in done1:
            U[s,:] = np.array([BB[kk,0],BB[kk,1],np.mean(BB[:,2]),np.mean(DD[:,3])])
            done1 += [[BB[kk,0],BB[kk,1],np.mean(BB[:,2])]]   
            s = s+1
    BB = Q2[np.where((Q2[:,2] > (ll-1)*0.001905) * (Q2[:,2] <=ll*0.001905)),:][0]
    for kk in range(len(BB[:,0])):
        CC = BB[np.where(abs(BB[:,0] - BB[kk,0])<0.00001),:][0]
        DD = CC[np.where(abs(CC[:,1] - BB[kk,1])<0.00001),:][0]
        if [BB[kk,0],B[kk,1],np.mean(BB[:,2])] not in done1:
            U2[s2,:] = np.array([BB[kk,0],BB[kk,1],np.mean(BB[:,2]),np.mean(DD[:,3])])
            done2 += [[BB[kk,0],BB[kk,1],np.mean(BB[:,2])]]   
            s2 = s2+1
  
U = U[np.where(U[:,3] > 0.0),:][0]
U2 = U2[np.where(U2[:,3] > 0.0),:][0] 


nd = 0.308473;
a = 1.0104;
c = 1.0074;
spos_Si = [[0.33333333,  0.66666667,  0.93750000-1], [0.00000, 0.00000, 0.1875], [0.66666667,  0.33333333,  0.43750000], [0.000000,  0.000000,  0.68750000], [0.33333333,  0.66666667,  0.93750000], [0.00000, 0.00000, 0.1875+1]]
cell2 = [[nd, 0*a, 0*a], [-nd/2,nd/2*np.sqrt(3), 0*a], [0*1.0, 0*1.0, 10.086*c*nd/3.078*a/2.57218587467527*2.51866888630220]]
t = -1;
fig = plt.figure()
ax = fig.add_subplot(111)
qq = 1
maxer = U[0,3];
done = []
for jj in range(len(U)):
    if U[jj,2] not in done:
        maxer = max(U[jj:,3]);
        U3 = Q[np.where(Q[:,2]<=U[jj,2]),:][0]
        cax = ax.scatter(U3[:,0],U3[:,1],s=100, c=U3[:,3], vmin = 0, vmax = maxer, edgecolors='none', cmap = 'Oranges')
        
        cbar = fig.colorbar(cax, ticks=[0, maxer/4,maxer/2,3*maxer/4, maxer])
        plt.axis('equal')
        plt.title('Probability map at: %.4f s (1300 K)' %U[jj,2],fontsize=18)
        plt.xlabel('$\\Delta$ x (nm)',fontsize=16);
        plt.ylabel('$\\Delta$ y (nm)',fontsize=16);
        ax.set_xbound([-1.8, 1.8]);
        ax.set_ybound([-1.8, 1.8]);
        ax.set_xlim([-1.8, 1.8]);
        ax.set_ylim([-1.8, 1.8]);
        plt.savefig('PXY%02d.png'%qq)
        plt.close()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        qq = qq +1
        done += [U[jj,2]]
        print(U[jj,2])
t = -1;
fig = plt.figure()
ax = fig.add_subplot(111)
qq = 1
maxer = U2[0,3];
done = []
for jj in range(len(U2)):
    if U2[jj,2] not in done:
        maxer = max(U2[jj:,3]);
        U3 = Q2[np.where(Q2[:,2]<=U2[jj,2]),:][0]
        cax = ax.scatter(U3[:,0],U3[:,1],s=100, c=U3[:,3], vmin = 0, vmax = maxer, edgecolors='none', cmap = 'Oranges')
        
        cbar = fig.colorbar(cax, ticks=[0, maxer/4,maxer/2,3*maxer/4, maxer])
        plt.axis('equal')
        plt.title('Probability map at: %.4f s (1300 K)' %Q[jj,2],fontsize=18)
        plt.xlabel('$\\Delta$ x (nm)',fontsize=16);
        plt.ylabel('$\\Delta$ z (nm)',fontsize=16);
        ax.set_xbound([-1.8, 1.8]);
        ax.set_ybound([-1.8, 1.8]);
        ax.set_xlim([-1.8, 1.8]);
        ax.set_ylim([-1.8, 1.8]);
        plt.savefig('PXZ%02d.png'%qq)
        plt.close()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        qq = qq +1
        done += [U2[jj,2]]
        print(U2[jj,2])

f, (ax0) = plt.subplots(1, 1)
T0_0 = 256
T0 = 226310.000000/1000
ax0.plot([1,4,8,16,36],[1,T0_0/103/4,T0_0/83/8,T0_0/59/16,T0_0/62/36],label="$(12^6-1)/11$ paths - recursive")
ax0.plot([1,4,8,16,36],[1,T0/104/4,T0/83/8,T0/66/16,T0/67/36],label="$(12^6-1)/11$ paths - BFS")
ax0.set_xlabel("Number of cores",fontsize=16)
ax0.set_ylabel("Efficiency",fontsize=16)
ax0.legend(loc='best',fontsize=16)
plt.savefig('efficiency.png')

f, (ax0) = plt.subplots(1, 1)

ax0.plot([1,4,8,16,36],[1,T0_0/103,T0_0/83,T0_0/59,T0_0/62],label="$(12^6-1)/11$ paths - recursive")
ax0.plot([1,4,8,16,36],[1,T0/104,T0/83,T0/66,T0/67],label="$(12^6-1)/11$ paths - BFS")
ax0.set_xlabel("Number of cores",fontsize=16)
ax0.set_ylabel("Speedup",fontsize=16)
ax0.legend(loc='best',fontsize=16)
plt.savefig('speedup.png')
