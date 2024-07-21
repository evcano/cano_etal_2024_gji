import numpy as np

stacode = np.loadtxt('./STATIONS', usecols=(0),dtype='str')
netcode = np.loadtxt('./STATIONS', usecols=(1),dtype='str')

for i in range(0, len(stacode)):
    sta = stacode[i]
    net = netcode[i]
    _f = open(f'irec_main_noise_{net}.{sta}', 'w')
    _f.write(str(i+1))
    _f.close()
