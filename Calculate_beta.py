import numpy
import math
import re
import matplotlib.pyplot as plt


# Read frequencies from fildyn files
def dmat_getfreqs(filepath):
    regex = re.compile('freq.')
    freqs = []
    with open(filepath) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            match = re.match(regex, line)
            if match:
                freqs.append(float(line[40:52].strip()) * 3.0E10)  # Convert cm^-1 to s^-1
    return freqs


def dm_getfreqs(filepath):
    freqs = []
    regex1 = re.compile('^(?!0)([0-9]{1,2})')
    regex2 = re.compile('((?!\.).)*$')
    with open(filepath) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if re.match(regex1, line):
                if re.match(regex2, line[:4]):
                    freqs.append(float(line[4:12].strip()) * 3.0E10)
    return freqs


# Calculate 1000ln_beta from fildyn file
# file_path1: the directory containing heavy isotopes
def dmat_calcbeta(t_low, t_high, t_step, q_point, isotope_sites, file_path1, file_path2, file_name):
    h = -6.62607004E-34
    k = 1.38064852E-23
    beta_ts = []
    betas = [1] * int(math.ceil((t_high - t_low) / t_step))
    for t, T in enumerate(numpy.arange(t_low, t_high, t_step)):
        for q in range(1, q_point+1):
            dmat1 = file_path1 + file_name + str(q)
            dmat2 = file_path2 + file_name + str(q)

            freqs1 = dmat_getfreqs(dmat1)
            freqs2 = dmat_getfreqs(dmat2)

            # Product over all frequencies. The first three translational frequencies are neglected.
            for i in range(3, len(freqs1)):
                betas[t] *= \
                    freqs1[i] * math.exp(h * freqs1[i] / (2 * k * T)) * (1 - math.exp(h * freqs2[i] / (k * T))) / \
                    (freqs2[i] * (1 - math.exp(h * freqs1[i] / (k * T))) * math.exp(h * freqs2[i] / (2 * k * T)))

    for beta in betas:
        beta_ts.append(1000.0 * math.log(math.pow(beta, 1 / (q_point * isotope_sites))))

    return beta_ts


# Calculate 1000ln_beta from dm.x output
# dm1 & dm2 require two files, not file directory
# material parameter ('crystal', 'linear', 'nonlinear') depends on the material
def dm_calcbeta(t_low, t_high, t_step, isotope_sites, dm1, dm2, material='linear'):
    h = -6.62607004E-34
    k = 1.38064852E-23
    beta_ts = []
    betas = [1] * int(math.ceil((t_high - t_low) / t_step))
    for t, T in enumerate(numpy.arange(t_low, t_high, t_step)):
        if material == 'linear':
            cut = 5
        elif material == 'nonlinear':
            cut = 6
        elif material == 'crystal':
            cut = 3

        freqs1 = dm_getfreqs(dm1)
        freqs2 = dm_getfreqs(dm2)

        # Product over all frequencies. The first three translational frequencies are neglected.
        for i in range(cut, len(freqs1)):
            betas[t] *= \
                freqs1[i] * math.exp(h * freqs1[i] / (2 * k * T)) * (1 - math.exp(h * freqs2[i] / (k * T))) / \
                (freqs2[i] * (1 - math.exp(h * freqs1[i] / (k * T))) * math.exp(h * freqs2[i] / (2 * k * T)))

    for beta in betas:
        beta_ts.append(1000.0 * math.log(math.pow(beta, 1 / isotope_sites)))

    return beta_ts


# Generate 1000ln_beta versus 10e6/T^2 data output for Excel/OriginLab
def outbeta(out_temperatures, out_betas, outdir):
    with open(outdir, 'w') as f:
        for out_t, out_beta in zip(out_temperatures, out_betas):
            f.write(str(out_t)+', '+str(out_beta)+'\n')


# Generate temperature coordinates (X axis)
T_low = 50.0
T_high = 5250.0
T_step = 50.0
temperatures = []

for T in numpy.arange(T_low, T_high, T_step):
    temperatures.append(1.0E6 / T ** 2)

# Computation information
q_point1 = 1    # No. of q points in the calculation (same as # of fildyn files)
q_point2 = 1
q_point3 = 1
q_point4 = 1

isotope1 = 2    # No. of isotopic sites in the unit cell
isotope2 = 2
isotope3 = 4
isotope4 = 1

file_name1 = 'dm.K41.out'        # used for dmat_calcbeta
file_name2 = 'dm.K39.out'        # Not used in dm_calcbeta
file_name3 = 'dynmat.out'         # Not used in dm_calcbeta
file_name4 = 'dynmat.out'         # Not used in dm_calcbeta

file1 = 'D:/Research/K-Bearing/microcline/Na-ONCV/dm.K41.out'   # Heavy isotope
file2 = 'D:/Research/K-Bearing/microcline/Na-ONCV/dm.K39.out'   # Light isotope
out1 =  'D:/Research/K-Bearing/microcline/microcline-ONCV_Na-K.txt'

file3 = 'D:/Research/K-Bearing/microcline/Na-GBRV/dm.K41.out'   # Heavy isotope
file4 = 'D:/Research/K-Bearing/microcline/Na-GBRV/dm.K39.out'   # Light isotope
out2 =  'D:/Research/K-Bearing/microcline/microcline-GBRV_Na-K.txt'

file5 = 'D:/Research/K-Bearing/microcline/GBRV-LDA/dm.K41.out'   # Heavy isotope
file6 = 'D:/Research/K-Bearing/microcline/GBRV-LDA/dm.K39.out'   # Light isotope
out3 =  'D:/Research/K-Bearing/microcline/microcline-GBRV_K.txt'

file7 = 'D:/Research/solvation/config6/K41/'   # Heavy isotope
file8 = 'D:/Research/solvation/config6/K39/'   # Light isotope
out4 =  'D:/Research/solvation/config6/'

#beta_ts1 = dmat_calcbeta(T_low, T_high, T_step, q_point1, isotope1, file1, file2, file_name3)
beta_ts1 = dm_calcbeta(T_low, T_high, T_step, isotope1, file1, file2, 'crystal')
#beta_ts2 = dmat_calcbeta(T_low, T_high, T_step, q_point2, isotope2, file3, file4, file_name2)
beta_ts2 = dm_calcbeta(T_low, T_high, T_step, isotope2, file3, file4, 'crystal')
# beta_ts3 = dmat_calcbeta(T_low, T_high, T_step, q_point3, isotope3, file5, file6, file_name3)
beta_ts3 = dm_calcbeta(T_low, T_high, T_step, isotope3, file5, file6, 'crystal')
# beta_ts4 = dmat_calcbeta(T_low, T_high, T_step, q_point4, isotope4, file7, file8, file_name4)

outbeta(temperatures, beta_ts1, out1)
outbeta(temperatures, beta_ts2, out2)
# outbeta(temperatures, beta_ts3, out3)
# outbeta(temperatures, beta_ts4, out4)

plt.plot(
   temperatures, beta_ts1, 'b--',
   temperatures, beta_ts2, 'g--',
   temperatures, beta_ts3, 'y--',
#   temperatures, beta_ts4, 'r'
)

plt.xlabel('10^6/T^2 (K^-2)')
plt.ylabel('1000ln(beta)')
plt.axis([0, 12, 0, 5])
plt.show()
