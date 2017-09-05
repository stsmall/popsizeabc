# -*- coding: utf-8 -*-
"""
This script computes the AFS and LD summary statistics described in
Boitard 2016, from a set of genomic samples that have already been produced by
the script simul_data.py (with option save_ms=True). It can be used for
instance to compute summary statistics using different MAF thresholds, or a
different sample size, than in the first analysis.
With the current parameter values, the only difference with simul_data.py
is that we use a MAF threshol of 20%, similar to what is done in the manuscript
for the cattle data analysis.
This can easily be changed by modifying variables in the parameters section
below. In practice, it is strongly recommended to excute this script using many
processors in parallel (with different output names!) and to merge all output
files at the end.
"""
import numpy as np
import popgen_abc
import argparse
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--infile", type=str, required=True,
                    help='infile name')
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
        times[i] = (np.exp(np.log(1.0 + a * Tmax)
                           * i/(nb_times - 1.0)) - 1.0)/a
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
        d = 1/(2 * r * t)
        if d <= L:
            interval_list.append([d - per_err * d/100.0,
                                  d + per_err * d/100.0])
    t = Tmax + times[nb_times - 1.0] - times[nb_times - 2.0]
    # d = 10.0**8/(2.0 * t)
    d = 1/(2*r*t)
    interval_list.append([d-per_err * d/100.0, d + per_err * d/100.0])
    print("Average LD will be computed for the following distance bins (in bp)"
          ": {}".format(interval_list))
    return(interval_list)


def calcsimstats(interval_list, nb_rep, haps2, infile, nb_seg, nb_seg2, mac,
                 mac_ld, L, haps):
    """

    Parameters:
        interval_list: list
        nb_rep: int, number of simulated datasets. This must be the same as
            that used in simul_data.py.
        haps: int, haploid sample size. This must be the same as that used in
            simul_data.py.
        haps2: int, number of haploid genomes used to compute summary
            statistics. must be less than parameter haps and even.
        infile: file,  root of the output files. This must be the same as
            that used in simul_data.py.
        nb_seg: int, number of independent segments in each dataset. This must
            be the same as that used in simul_data.py.
        nb_seg2: int, number of independent segments used to compute summary
            statistics. must be lower than nb_seg.
        mac: int, minor allele count threshold for AFS statistics computation
        mac_ld: int, minor allele count threshold for LD statistics computation
        L: int, length of one segment, in bp. This must be the same as that
            used in simul_data.py.

    """
    nb_dist = len(interval_list)
    results = -np.ones(shape=[nb_rep, nb_dist + haps2/2 + 1], dtype='float')
    # missing statistic values will be set to 1.
    # most common in LD statistics, because it is sometimes impossible to find\
    # SNP pairs in a given distance bin.
    print("Total number of statistics:{}".format(nb_dist + haps2/2 + 1))
    # compute summary statistics
    out_name_2_read = "{}_n{}_s{}".format(infile, haps, nb_seg)
    out_name_2_write = "{}_n{}_s{}_mac{}_macld{}".format(infile,
                                                         haps2, nb_seg2,
                                                         mac, mac_ld)
    for i in range(nb_rep):
        try:
            results[i, :] = popgen_abc.comp_stats_one_rep_macld(
                   out_name_2_read, i, nb_seg2, L, haps2, interval_list,
                   mac=mac, mac_ld=mac_ld)
        except:  # this should really be a specific error
            raise Exception("Problem with replicate, {}".format(i + 1))
            pass
    # print the result
    np.savetxt("{}.stat".format(out_name_2_write), results[0:nb_rep, :],
               fmt='%.3e')


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
    nb_rep = config.getint(sh, "nb_rep")
    nb_seg = config.getint(sh, "nb_seg")
    haps = config.getint(sh, "haps")
    haps2 = config.getint(sh, "haps2")
    nb_seg2 = config.getint(sh, "nb_seg2")
    mac = config.getint(sh, "minallelcount_sfs")
    mac_ld = config.getint(sh, "minallelcount_LD")
    infile_name = args.infile
    calcsimstats(interval_list, nb_rep, haps2, infile_name, nb_seg, nb_seg2,
                 mac, mac_ld, L, haps)
