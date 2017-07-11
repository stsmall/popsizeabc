# -*- coding: utf-8 -*-
"""
This script computes the AFS, LD and IBS summary statistics described in the
manuscript of Boitard 2016.

input: vcf format, inds.txt, pops.txt

inds.txt:
Angus_6
Angus_7
Angus_8

pops.txt as ped file format
Angus Angus_6 0 0 1 -999
Angus Angus_7 0 0 1 -999
Angus Angus_8 0 0 1 -999

output: Summary statistics are stored in a file with suffix .stat, which has a
single column with one statistic per line.
"""
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import os
from gwas import IO
from gwas import data
from gwas import summary_stat as ss
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-pop', "--population", type=str,
                    help='population to calculate')
parser.add_argument('-ind', "--individuals", type=str, required=True,
                    help="path to individuals file")
parser.add_argument('-vcf', "--vcffile", type=str, required=True,
                    help="path to vcf directory, vcfs should be named as"
                    "chr.vcf.gz")
parser.add_argument('-ped', "--pedfile", type=str, required=True,
                    help="path to pedfile")
parser.add_argument('-chr', "--chrlist", type=str, required=True,
                    help="path to chromfile")
parser.add_argument('-o', "--out", type=str, required=True,
                    help="outfile prefix")
parser.add_arguement('-c', "--config", type=str, required=True)
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
        times[i] = (np.exp(np.log(1 + a * Tmax)
                            * i/(nb_times - 1)) - 1)/a
    print("Population size changes at the following times (in generations): {}"
          .format(times))
    return(times)


def ldstats(nb_times, time, r, L, per_err, Tmax):
    """Creation of the bins of physical distance for which the average LD will
    be computed, based on the time windows defined above.

   Parameters:
        nb_times: int,
        times: array,
        r: float,recomb rate per generation per bp
        L: int, size of each segment, in bp.
        per_err: int,
        Tmax: int,

    Returns:
        intervals_list: list

    """
    interval_list = []
    for i in range(nb_times - 1):
        t = (times[i + 1] + times[i])/2
        d = 1/(2 * r * t)
        if d <= L:
            interval_list.append([d - per_err * d/100,
                                  d + per_err * d/100])
    t = Tmax + times[nb_times - 1] - times[nb_times - 2]
    d = 10**8/(2 * t)
    # d = 1/(2*r*t)
    interval_list.append([d-per_err * d/100, d + per_err * d/100])
    print("Average LD will be computed for the following distance bins (in bp)"
          ": {}".format(interval_list))
    return(interval_list)


def obsstats(haps, interval_list, chromlist, vcf, list_ani, popsped,
             pop, mac, mac_ld, L, outfile):
    """program to calc summary stats from vcf

    Parameters:
        haps: int, number of haps, dips * 2
        interval_list: list, list of time intervals
        chromlist: list, list of chromosomes; name chrom.vcf.gz
        vcf: file, vcf file
        list_ani: list, individuals from the vcf
        popsped: file, ped file
        pop: str, population outfile id
        mac: int, minor allele count
        mac_ld: int, minor allele count LD
        L: int, length of regions

    Returns:
        file: file with string of stats

    """
    print("Number of diploid individuals:{}\n".format(haps/2))
    nb_dist = len(interval_list)
    print("Total number of statistics:{}".format(nb_dist + haps/2 + 1))
    # summarize each chromosome
    count_list = []
    pos_list = []  # for ld
    geno_list = []  # for ld
    nb_snp = 0
    Lchr = 0
    for chrom in chromlist:
        # store data and pre-analyse it
        print("Processing chromosome,{}".format(chrom))
        infile_vcf = "{}/{}.vcf.gz".format(vcf, chrom)
        [mydata, mymap] = IO.parseVcfFile(infile_vcf, includeInd=list_ani)
        pedfile = popsped
        IO.parsePedFile_nogeno(pedfile, mydata)
        u = data.Genos_and_counts(mydata, mymap, pop, mac=mac)
        count_list.append(u[1][0])
        p = len(u[0][0])
        nb_snp += p
        Lchr += u[0][0][p - 1]
        u = data.Genos_and_counts(mydata, mymap, pop, mac=mac_ld)
        pos_list.append(u[0][0])
        geno_list.append(u[2][0])
    # compute summary statistics
    res_afs = ss.histo(count_list, int(haps/2))
    u = ss.break_chr(pos_list, geno_list, L)
    pos_list = u[0]
    geno_list = u[1]
    res_ld_zyg = ss.distrib_zyg_r2(pos_list, geno_list, interval_list)
    # print the result
    directory = "{}/res".format(vcf)
    if not os.path.exists(directory):
        os.makedirs(directory)
    fname = "{}/res/{}_{}_n{}_mac{}_macld{}.stat".format(vcf, outfile, pop,
                                                         haps, mac, mac_ld)
    fnp = np.array([np.float(nb_snp)/np.float(Lchr)], dtype='float')
    np.savetxt(fname, np.concatenate((fnp, res_afs, res_ld_zyg[0])),
               fmt='%.3e')

if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read(config.args)
    sh = "simulation"
    nb_times = config.getint(sh, "nbtimes")
    Tmax = config.getint(sh, "tmax")
    a = config.getfloat(sh, "a")
    times = timewindows(nb_times, Tmax, a)
    per_err = config.getint(sh, "per_err")
    r = config.getfloat(sh, "recomb")
    L = config.getint(sh, "seg_size")
    interval_list = ldstats(nb_times, times, r, L, per_err, Tmax)
    chromlist = []
    with open(args.chrlist, 'r') as chrm:
        for line in chrm:
            chromlist.append(line.strip())
    popsped = args.pedfile
    vcf = args.vcffile
    pop = args.population  # pop id
    list_ani = IO.read_list(args.individuals)  # list of inds for obs stats
    mac = config.getint(sh, "minallelcount_sfs")
    mac_ld = config.getint(sh, "minallelcount_LD")
    haps = len(list_ani) * 2  # haploid sample size
    obsstats(haps, interval_list, chromlist, vcf, list_ani, popsped,
             pop, mac, mac_ld, L, args.out)
