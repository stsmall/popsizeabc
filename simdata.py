# -*- coding: utf-8 -*-
"""
This script simulates genomic samples, and the corresponding AFS and LD summary
statistics, as described in Boitard 2016. Results will be written in the 'res'
directory, in two files with suffices: FOO.params and FOO.stat. Each of this
file includes one line per simulated dataset.

Each column in the params file corresponds to a different parameter.
    mutation_rate recombination_rate pop_size_1 pop_size_2 .... pop_size_I

    ,where I is the number of intervals

Each column in the FOO.stat file corresponds to a different statistic.
    AFS0 AFS1 ... AFSn LD1 LD2 ... LDp

    ,where n is the number of diploid individuals
    ,p is the number of physical distance bins

If option save_ms is set to True, a compressed tar file example.ms.tar.bz2
including the raw ms output for all simulated datasets and all segments is
also created.

Simulation settings are similar to those used in our study, but this can
easily be changed by modifying variables in the parameters section below.
In practice, it is strongly recommended to excute this script using many
processors in parallel (with different output names!) and to merge all output
files at the end.

"""
import numpy as np
import popgen_abc
import tarfile
import os
import argparse
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

parser = argparse.ArgumentParser()
parser.add_argument('-o', "--outfile", type=str, required=True,
                    help='outfile name')
parser.add_argument('-c', "--configfile", type=str, required=True,
                    help='config file')
args = parser.parse_args()


def timewindows(nb_times, Tmax, a):
    """The length of time windows increases when time increases at a speed that
    is determined by this coefficient computation of time windows based on the
    above parameters

    Parameters:
        nb_times: int, number of time windows
        Tmax: int, the oldest time window will start at Tmax
        a: float, the length of time windows increases when time increases,
             at a speed that is determined by this coefficient
    Returns:
        times: numpy array

    """
    times = -np.ones(shape=nb_times, dtype='float')
    for i in range(nb_times):
        times[i] = (np.exp(np.log(1 + a * Tmax) *
                           i/(nb_times - 1)) - 1)/a
    print("Population size changes at the following times (in generations): {}"
          .format(times))
    return(times)


def ldstats(nb_times, time, r, L, per_err, Tmax):
    """Creation of the bins of physical distance for which the average LD will
    be computed, based on the time windows defined above.

   Parameters:
        nb_times: int, number of time windows
        times: array,
        r: float,recomb rate per generation per bp
        L: int, size of each segment, in bp.
        per_err: int, the length of each interval, as a per of the target
            distance
        Tmax: int, the oldest time window will start at Tmax

    Returns:
        intervals_list: list

    """
    interval_list = []
    for i in range(nb_times - 1):
        t = (times[i + 1] + times[i])/2
        d = 1/(2*r*t)
        if d <= L:
            interval_list.append([d - per_err * d/100,
                                  d + per_err * d/100])
    t = Tmax + times[nb_times - 1] - times[nb_times - 2]
    d = 10.0**8/(2 * t)
    # d = 1/(2*r*t)
    interval_list.append([d-per_err * d/100, d + per_err * d/100])
    print("Average LD will be computed for the following distance bins (in bp)"
          ": {}".format(interval_list))
    return(interval_list)


def simdata(interval_list, nb_rep, nb_times, haps, mmu, r_min, r_max, Nmin,
            Nmax, outfile, nb_seg, L, times, mac, mac_ld, save_msp):
    """Simulates data under model of popsizeabc

    Parameters:
        interval_list: list
        nb_rep: int, number of simulated datasets. This must be the same as
            that used in simul_data.py.
        nb_times: int, number of time windows
        haps: int, haploid sample size. This must be the same as that used in
            simul_data.py.
        mmu: float, per generation per bp mutation rate
        r_min: float, lower bound for the per generation per bp recomb rate
        r_max: float, upper bound for the per generation per bp recomb rate
        Nmin: int, lower bound for the pop size in each time window log10 scale
        Nmax: int, upper bound for the pop size in each time window log10 scale
        outfile: file. root of the output files. This must be the same as
            that used in simul_data.py.
        nb_seg: int, number of independent segments in each dataset. This must
            be the same as that used in simul_data.py.
        L: int, length of one segment, in bp. This must be the same as that
            used in simul_data.py.
        times: array,
        mac: int, minor allele count threshold for AFS statistics computation
        mac_ld: int, minor allele count threshold for LD statistics computation
        save_msp: bool, if True, snp positions and haplotypes corresponding to
            the same dataset will be stored in a compressed tar file. This
            allows to keep and potentially re-use the exact genomic samples,
            rather than just the summary statistics, but this require high
            memory ressources (approx 1 Mo per simulated dataset, on average,
            with current parameter values).
    """

    # create the matrices where results (parameters and statistics) are stored
    nb_dist = len(interval_list)
    params = -np.ones(shape=[nb_rep, nb_times + 2], dtype='float')
    results = -np.ones(shape=[nb_rep, nb_dist + haps/2 + 1], dtype='float')
    # missing statistic values will be set to -1. LD statistics, because it is
    # sometimes impossible to find SNP pairs in a given distance bin.
    print("Total number of statistics:{}\n".format(nb_dist + haps/2 + 1))

    # simulate parameters from the prior
    mu = mmu * np.ones([nb_rep])
    params[:, 0] = mu
    r = np.random.uniform(low=r_min, high=r_max, size=nb_rep)
    params[:, 1] = r
    for i in range(nb_times):
        if i == 0:
            pop_size = 10**(np.random.uniform(low=Nmin, high=Nmax,
                                              size=nb_rep))
        else:
            multiplication_factor = 10**(np.random.uniform(-1, 1, size=nb_rep))
            pop_size = params[:, 1 + i] * multiplication_factor
            pop_size = np.maximum(10**Nmin, pop_size)
            pop_size = np.minimum(pop_size, 10**Nmax)
        params[:, 2 + i] = pop_size

    # simulate summary statistics
    directory = "./res"
    if not os.path.exists(directory):
        os.makedirs(directory)
    out_name_2 = "./res/{}_n{}_s{}".format(outfile, haps, nb_seg)
    print('Started the simualtions')
    for i in range(nb_rep):
        print("Simulating replicate, {}".format(i + 1))
        try:
            results[i, :] = popgen_abc.simul_stats_one_rep_macld(
                out_name_2, i, nb_seg, L, haps, times, params[i, :],
                interval_list, mac=mac, mac_ld=mac_ld, save_msp=save_msp)
        except:  # this should really be a specific error
            print("Problem with replicate, {}".format(i + 1))
            continue

    # print the results
    fname = "{}_mac{}_macld{}.stat".format(out_name_2, mac, mac_ld)
    np.savetxt(fname, results[0:nb_rep, :], fmt='%.3e')
    np.savetxt("{}.params".format(out_name_2), params[0:nb_rep, :], fmt='%.3e')
    print("Printed the results")

    # tar ms files (if they exist)
    if save_msp:
        mytar = tarfile.open("{}.msp.tar.bz2".format(out_name_2), 'w:bz2')
        for i in range(nb_rep):
            try:
                for j in range(nb_seg):
                    filename = "{}_rep{}_seg{}.msp".format(out_name_2, i + 1,
                                                           j + 1)
                    mytar.add(filename)
                    os.remove(filename)
            except OSError:
                continue
        mytar.close()
        print("Created tar for msp files")


if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read(args.configfile)
    sh = "simulation"
    nb_times = config.getint(sh, "nbtimes")
    Tmax = config.getint(sh, "tmax")
    a = config.getfloat(sh, "a")
    times = timewindows(nb_times, Tmax, a)
    per_err = config.getint(sh, "pererr")
    r = config.getfloat(sh, "recomb")
    L = config.getint(sh, "seg_size")
    interval_list = ldstats(nb_times, times, r, L, per_err, Tmax)
    outfile_name = args.outfile
    nb_rep = config.getint(sh, "nb_rep")
    nb_seg = config.getint(sh, "nb_seg")
    haps = config.getint(sh, "haps")
    mac = config.getint(sh, "minallelcount_sfs")
    mac_ld = config.getint(sh, "minallelcount_LD")
    save_msp = config.getboolean(sh, "savemsp")
    # prior distributions
    r_min = config.getfloat(sh, "recombmin")
    r_max = config.getfloat(sh, "recombmax")
    mmu = config.getfloat(sh, "mmu")
    Nmin = config.getfloat(sh, "Nemin")
    Nmax = config.getfloat(sh, "Nemax")
    simdata(interval_list, nb_rep, nb_times, haps, mmu, r_min, r_max, Nmin,
            Nmax, outfile_name, nb_seg, L, times, mac, mac_ld, save_msp)
