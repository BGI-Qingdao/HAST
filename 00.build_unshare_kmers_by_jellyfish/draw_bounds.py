#!/usr/bin/python3

import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend("agg")

maternal_xy = np.loadtxt("maternal.kmercount.count.txt",dtype=int)
paternal_xy = np.loadtxt("paternal.kmercount.count.txt",dtype=int)

maternal_min=0
maternal_man=0
maternal_lower=0
maternal_upper=0

paternal_min=0
paternal_man=0
paternal_lower=0
paternal_upper=0

with open('maternal.bounds.txt', 'r') as f:
    for c in f.readlines():
        c_array=c.split("=")
        if c_array[0] == 'MIN_INDEX':
            maternal_min = int(c_array[1])
        elif c_array[0] == 'MAX_INDEX':
            maternal_max = int(c_array[1])
        elif c_array[0] == 'LOWER_INDEX':
            maternal_lower = int(c_array[1])
        elif c_array[0] == 'UPPER_INDEX':
            maternal_upper = int(c_array[1])
        else :
            print("ERROR : unknow line %s"%c)

with open('paternal.bounds.txt', 'r') as f:
    for c in f.readlines():
        c_array=c.split("=")
        if c_array[0] == 'MIN_INDEX':
            paternal_min = int(c_array[1])
        elif c_array[0] == 'MAX_INDEX':
            paternal_max = int(c_array[1])
        elif c_array[0] == 'LOWER_INDEX':
            paternal_lower = int(c_array[1])
        elif c_array[0] == 'UPPER_INDEX':
            paternal_upper = int(c_array[1])
        else :
            print("ERROR : unknow line %s"%c)

plt.figure()
plt.subplot(2,1,1)
plt.plot(maternal_xy[:,0],maternal_xy[:,1])
plt.axvline(x=maternal_min,ls='--',c='r',label="MIN INDEX %d "%maternal_min)
plt.axvline(x=maternal_max,ls='--',c='g',label="MAX INDEX %d "%maternal_max)
plt.axvline(x=maternal_lower,ls='-.',c='r',label="LOWER INDEX %d "%maternal_lower)
plt.axvline(x=maternal_upper,ls='-.',c='g',label="UPPER INDEX %d "%maternal_upper)
plt.legend(loc='best')
plt.xlim(1,150)
plt.xlabel('kmer depth')
plt.ylabel('count')
#plt.yscale('log')
plt.title('maternal kmer-depth count')
plt.subplot(2,1,2)
plt.plot(paternal_xy[:,0],paternal_xy[:,1])
plt.axvline(x=paternal_min,ls='--',c='r',label="MIN INDEX %d "%paternal_min)
plt.axvline(x=paternal_max,ls='--',c='g',label="MAX INDEX %d "%paternal_max)
plt.axvline(x=paternal_lower,ls='-.',c='r',label="LOWER INDEX %d "%paternal_lower)
plt.axvline(x=paternal_upper,ls='-.',c='g',label="UPPER INDEX %d "%paternal_upper)
plt.xlim(1,150)
#plt.yscale('log')
plt.xlabel('kmer depth')
plt.ylabel('count')
plt.legend(loc='best')
plt.title('paternal kmer-depth count')
plt.subplots_adjust(hspace=0.4)
plt.savefig('test.png')

