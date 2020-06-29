import sys
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["text.usetex"] = "True"
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16
plt.rcParams["axes.labelsize"] = 18
plt.rcParams["figure.figsize"] = [5, 7]

def parse_out_eig(lines):
    n = 0
    for line in lines:
        if len(line.split()) == 0:
            break
        else:
            n = n + 1
            continue

    ks = get_kpts(lines[0:n])

    nband = int( len(lines) / n )
    bands = []
    for ib in range(nband):
        n1 = ib * (n + 1)
        n2 = n1 + n

        bands.append( get_bands(lines[n1:n2]) )
        
    return ks, bands 

def get_bands(lines):
    es = []

    for line in lines:
        sp = line.split()
        es.append( float(sp[3]) )

    return np.array( es )

def get_kpts(lines):
    ks = []
    
    sp = lines[0].split()
    k0 = np.array( [float(sp[0]), float(sp[1]), float(sp[2])] )

    kdist = 0.0
    for line in lines:
        sp = line.split()
        k1 = np.array( [float(sp[0]), float(sp[1]), float(sp[2])] )
        kdist = kdist + np.linalg.norm( k1-k0 )
        ks.append( kdist )
        k0 = k1

    return np.array( ks )

def get_nkpt(line):
    sp = line.split()

    return int( sp[1] )

def get_labels(lines):
    labels = []
    for line in lines:
        sp = line.split()
        labels.append( sp[0] )

    sp = lines[-1].split()
    labels.append( sp[4] )

    return labels

######

out_eig_file = sys.argv[1]
in_k_band_file = sys.argv[2]

ymin = float(sys.argv[3])
ymax = float(sys.argv[4])

f = open( out_eig_file )
lines = f.readlines()
f.close()

ks, es = parse_out_eig( lines )

for v in es:
    plt.plot( ks, v, 'black', lw = 1.6 )

plt.xlim( 0, ks[-1] )
plt.ylabel( "Energy (ev)" )

plt.xticks([])

######

f = open( in_k_band_file )
lines = f.readlines()
f.close()

nkpt = get_nkpt( lines[0] )
labels = get_labels( lines[1:] )

plt.ylim(ymin,ymax)

for ilab in range(len(labels)):
    if ilab == 0:
        ik = 0
    else:
        ik = ilab * nkpt - 1
        
    if labels[ilab] == 'G':
        plt.text(ks[ik], ymin - (ymax - ymin) * 0.01, \
                 '$\Gamma$', fontsize = 18, \
                 horizontalalignment= 'center', verticalalignment = 'top')
    else:
        plt.text(ks[ik], ymin - (ymax - ymin) * 0.01, \
                 labels[ilab], fontsize = 18, \
                 horizontalalignment= 'center', verticalalignment = 'top')

    plt.axvline( ks[ik], color = 'black' ) 


        
plt.tight_layout()

##plt.show()
plt.savefig( "out.eig.pdf" )
