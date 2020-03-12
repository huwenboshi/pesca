import numpy as np
import pandas as pd
import argparse, math, sys

# main function
def main():
    
    args = get_command_line()
 
    # parse log file
    chrom_cau_cnt_est = []
    for i in range(1,23):
        chrom_cau_cnt_est.append([])
        log_file = open('{}{}.log'.format(args.prefix, i), 'r')
        for line in log_file:
            cols = line.strip().split()
            if cols[0] == 'sigma_sq' or cols[0] == 'iter':
                continue
            chrom_cau_cnt_est[i-1].append([float(cols[i]) for i in range(2,6)])
        log_file.close()

    # extract last 50 entries
    null_cvg_est = []
    eas_cvg_est = []
    eur_cvg_est = []
    both_cvg_est = []
    for j in range(22):
        sum_null = []
        sum_eas = []
        sum_eur = []
        sum_both =[]
        nentry = len(chrom_cau_cnt_est[j])
        for i in range(50):
            sum_null.append(chrom_cau_cnt_est[j][nentry-1-i][0])
            sum_eas.append(chrom_cau_cnt_est[j][nentry-1-i][1])
            sum_eur.append(chrom_cau_cnt_est[j][nentry-1-i][2])
            sum_both.append(chrom_cau_cnt_est[j][nentry-1-i][3])
        null_cvg_est.append(sum_null)
        eas_cvg_est.append(sum_eas)
        eur_cvg_est.append(sum_eur)
        both_cvg_est.append(sum_both)

    # covert python array to numpy array
    null_cvg_est = np.array(null_cvg_est)
    eas_cvg_est = np.array(eas_cvg_est)
    eur_cvg_est = np.array(eur_cvg_est)
    both_cvg_est = np.array(both_cvg_est)

    # compute the sum of each chromosome
    null_est = np.sum(null_cvg_est, axis=0)
    eas_est = np.sum(eas_cvg_est, axis=0)
    eur_est = np.sum(eur_cvg_est, axis=0)
    both_est = np.sum(both_cvg_est, axis=0)

    # compute the mean
    null_est_mean = np.mean(null_est)
    eas_est_mean = np.mean(eas_est)
    eur_est_mean = np.mean(eur_est)
    both_est_mean = np.mean(both_est)

    # get the parameter
    f00 = 0.0
    f01 = np.log(eas_est_mean / null_est_mean)
    f10 = np.log(eur_est_mean / null_est_mean)
    f11 = np.log(both_est_mean / null_est_mean) - f01 - f10

    # print to screen
    print("number of SNPs")
    print("q00 {:.2f}".format(null_est_mean))
    print("q01 {:.2f}".format(eas_est_mean))
    print("q10 {:.2f}".format(eur_est_mean))
    print("q11 {:.2f}".format(both_est_mean))
    print()
    print("MVB parameters")
    print("f00 {:.4f}".format(f00))
    print("f01 {:.4f}".format(f01))
    print("f10 {:.4f}".format(f10))
    print("f11 {:.4f}".format(f11))

# get command line
def get_command_line():
    
    parser = argparse.ArgumentParser(description='plot convergence')

    parser.add_argument('--prefix', dest='prefix', type=str,
        help='prefix file', required=True)

    args = parser.parse_args()
    
    return args

if(__name__ == '__main__'):
    main()
