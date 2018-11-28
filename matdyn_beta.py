import re
import math
import numpy
import matplotlib.pyplot as plt

def get_freqs(file):
    regex = re.compile('^        .')
    freqs = []
    remove = []
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if re.search('&.', line):
                freq_num = int(re.sub(',', '', line.split()[2]))
                nq = int(line.split()[4])
                n = nq - 1
                continue
            if re.search(regex, line) is None:
                for freq in line.split():
                    freqs.append(float(freq) * 3.0E10)
        while n >= 0:
            remove.append(n * freq_num + 2)
            remove.append(n * freq_num + 1)
            remove.append(n * freq_num)
            n -= 1
        for i in remove:
            del freqs[i]

        return freqs


def get_qnum(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if re.search('&.', line):
                q_point = int(line.split()[4])
    return q_point


def matdyn_calcbeta(t_low, t_high, t_step, isotope_sites, file_path1, file_path2):
    h = -6.62607004E-34
    k = 1.38064852E-23
    beta_ts = []
    q_point = get_qnum(file_path1)
    betas = [1] * int(math.ceil((t_high - t_low) / t_step))
    freqs1 = get_freqs(file_path1)
    freqs2 = get_freqs(file_path2)
    for t, T in enumerate(numpy.arange(t_low, t_high, t_step)):
        # Product over all frequencies. The first three translational frequencies are neglected.
        for i in range(1, len(freqs1)):
            betas[t] *= \
                freqs1[i] * math.exp(h * freqs1[i] / (2 * k * T)) * (1 - math.exp(h * freqs2[i] / (k * T))) / \
                (freqs2[i] * (1 - math.exp(h * freqs1[i] / (k * T))) * math.exp(h * freqs2[i] / (2 * k * T)))

    for beta in betas:
        beta_ts.append(1000.0 * math.log(math.pow(beta, 1 / (q_point * isotope_sites))))

    return beta_ts


def outbeta(out_temperatures, out_betas, outdir):
    with open(outdir, 'w') as f:
        for out_t, out_beta in zip(out_temperatures, out_betas):
            f.write(str(out_t)+', '+str(out_beta)+'\n')

# Generate temperature coordinates (X axis)
T_low = 250.0
T_high = 5250.0
T_step = 250.0
temperatures = []

for T in numpy.arange(T_low, T_high, T_step):
    temperatures.append(1.0E6 / T ** 2)

isotope1 = 2    # No. of isotopic sites in the unit cell
isotope2 = 1
isotope3 = 1


file_name1 = 'D:/Research/K-Bearing Minerals/muscovite/muscovite.freq-K41'  # Heavy isotope
file_name2 = 'D:/Research/K-Bearing Minerals/muscovite/muscovite.freq-K39'  # Light isotope
outdir1 =    'D:/Research/K-Bearing Minerals/muscovite/muscovite_K.txt'

file_name3 = 'D:/Research/K-Bearing Minerals/illite/illite.freq-Rb87'   # Heavy isotope
file_name4 = 'D:/Research/K-Bearing Minerals/illite/illite.freq-Rb85'   # Light isotope
outdir2 =    'D:/Research/K-Bearing Minerals/illite/illite_Rb.txt'

file_name5 = 'D:/Research/Diamond/N1V0/diamond.freq-N15'   # Heavy isotope
file_name6 = 'D:/Research/Diamond/N1V0/diamond.freq-N14'   # Light isotope
outdir3 =    'D:/Research/Diamond/N1V0/diamond_N.txt'

beta_ts1 = matdyn_calcbeta(T_low, T_high, T_step, isotope1, file_name1, file_name2)
outbeta(temperatures, beta_ts1, outdir1)

beta_ts2 = matdyn_calcbeta(T_low, T_high, T_step, isotope2, file_name3, file_name4)
outbeta(temperatures, beta_ts2, outdir2)

#beta_ts3 = matdyn_calcbeta(T_low, T_high, T_step, isotope3, file_name5, file_name6)
#outbeta(temperatures, beta_ts3, outdir3)


plt.plot(
    temperatures, beta_ts1, 'b',
#    temperatures, beta_ts2, 'g--',
#    temperatures, beta_ts3, 'r--',
#    temperatures, beta_ts4, 'r'
)

plt.xlabel('10^6/T^2 (K^-2)')
plt.ylabel('1000ln(beta)')
plt.axis([0, 14, 0, 2])
plt.show()