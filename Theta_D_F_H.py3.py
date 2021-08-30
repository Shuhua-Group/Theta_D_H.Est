# -*- coding: utf-8 -*-

#####################################################################################################
## By:             Pan Yuwen, 05/2021
## Contact:        panyuwen.x@gmail.com
#####################################################################################################

import numpy as np
import pandas
import gzip
import re
from functools import reduce
import math
from math import sqrt
import argparse
import sys
import socket
import os
import time
import gc
from scipy import integrate
from scipy.special import gamma
#from rpy2.robjects.packages import importr
#from rpy2.robjects.vectors import FloatVector
#stats = importr('stats')

## 1993-Statistical tests of neutrality of mutations
## under the neutral model
## E(s) = an * theta; 
## E(pi) = theta
## E(η) = n/(n-1) * theta

## FU & Li's D: K vs. singleton
## FU & Li's F: pi vs. singleton
## Fay & Wu's H: pi vs. sum(#mutant^2)
## Tajima's D: pi vs. K

## Fu and Li's D, with or without outgroup
## #singletons may overestimate the #mutations in the external branches without outgroup
def calculate_Dfuli_outgroup(s,n,m):
    ## s: num of segregating sites
    ## n: num of DNA sequences
    ## m: num of derived singletons

    an = reduce(lambda x,y: x+1.0/y, range(1,n))
    bn = reduce(lambda x,y: x+1.0/pow(y,2), range(1,n))
    if n==2:
        cn = 1
    else:
        cn = 2.0 * (n*an - 2.0*(n-1.0)) / ((n-1.0)*(n-2.0))

    v = 1.0 + (pow(an,2) / (bn+pow(an,2))) * (cn - (n+1.0)/(n-1.0))
    u = an - 1.0 - v
    D = (s - m*an) / sqrt(u*s + v*pow(s,2))

    return D

def calculate_Dfuli_no_outgroup(s,n,m):
    ## s: num of segregating sites
    ## n: num of DNA sequences
    ## m: num of singletons

    an = reduce(lambda x,y: x+1.0/y, range(1,n))
    an1 = an + 1.0/n
    bn = reduce(lambda x,y: x+1.0/pow(y,2), range(1,n))
    if n==2:
        cn = 1
    else:
        cn = 2.0 * (n*an - 2.0*(n-1.0)) / ((n-1.0)*(n-2.0))
    dn = cn + (n-2.0)/pow(n-1.0,2) + 2.0/(n-1.0)*(3.0/2-(2.0*an1-3.0)/(n-2.0)-1.0/n)

    v = (pow(n*1.0/(n-1.0),2)*bn + pow(an,2)*dn - 2.0*n*an*(an+1.0)/pow(n-1,2)) / (pow(an,2)+bn)
    u = n*1.0/(n-1.0) * (an-n*1.0/(n-1.0)) - v

    D = (s*n*1.0 / (n-1.0) - an*m) / sqrt(u*s + v*pow(s,2))

    return D

def calculate_Ffuli_outgroup(s,n,m,pi):
    ## s: num of segregating sites
    ## n: num of DNA sequences
    ## m: num of derived singletons

    an = reduce(lambda x,y: x+1.0/y, range(1,n))
    an1 = an + 1.0/n
    bn = reduce(lambda x,y: x+1.0/pow(y,2), range(1,n))
    if n==2:
        cn = 1
    else:
        cn = 2.0 * (n*an - 2.0*(n-1.0)) / ((n-1.0)*(n-2.0))

    v = (cn+2.0*(pow(n,2)+n+3.0)/(9.0*n*(n-1.0))-2.0/(n-1.0)) / (pow(an,2)+bn)
    u = (1.0 + (n+1.0)/(3.0*(n-1.0)) - 4.0*(n+1.0)/pow(n-1,2)*(an1-2.0*n/(n+1.0))) / an - v
    F = (pi - m) / sqrt(u*s + v*pow(s,2))

    return F

## 1995-Properties of Statistical Tests of Neutrality for DNA Polymorphism Data
def calculate_Ffuli_no_outgroup(s,n,m,pi):
    ## s: num of segregating sites
    ## n: num of DNA sequences
    ## m: num of singletons

    an = reduce(lambda x,y: x+1.0/y, range(1,n))
    an1 = an + 1.0/n
    bn = reduce(lambda x,y: x+1.0/pow(y,2), range(1,n))

    #v = (dn + 2.0*(pow(n,2)+n+3.0)/(9.0*n*(n-1.0)) - 2.0/(n-1.0)*(4.0*bn-6.0+8.0/n)) / (pow(an,2) + bn)
    #u = (n*1.0/(n-1.0)+(n+1.0)/(3.0*(n-1.0))-4.0/(n*(n-1.0))+2.0*(n+1.0)/pow(n-1,2)*(an1-2.0*n/(n+1.0))) / an -v
    v = ((2.0*pow(n,3)+110.0*pow(n,2)-255.0*n+153.0) / (9.0*pow(n,2)*(n-1.0)) + 2.0*(n-1.0)*an/pow(n,2) - 8.0*bn/n) / (pow(an,2) + bn)
    u = (4.0*pow(n,2)+19.0*n+3.0-12.0*(n+1.0)*an1) / (3.0*n*(n-1)) / an - v

    F = (pi - m*1.0*(n-1.0)/n) / sqrt(u*s + v*pow(s,2))

    return F

## ancestral stat required
## 2000-Hitchhiking Under Positive Darwinian Selection
## 2006-Statistical Tests for Detecting Positive Selection by Utilizing High-Frequency Variants, (8, 11, 12)
def calculate_Hfaywu(s,n,pi,hapcount):
    ## s: num of segregating sites
    ## n: num of DNA sequences

    ## original Fay and Wu's H
    #count = pandas.value_counts(hapcount['1'])
    #thetaH = 2.0*(count.values * np.power(np.array(count.index),2)).sum() / (n*(n-1))
    thetaH = 2.0 * np.power(hapcount['1'].values,2).sum() / (n*(n-1.0))
    H = pi - thetaH

    ## normalized Fay and Wu's H
    thetaL = hapcount['1'].sum()*1.0 / (n-1.0)
    an = reduce(lambda x,y: x+1.0/y, range(1,n))
    bn = reduce(lambda x,y: x+1.0/pow(y,2), range(1,n))
    bn1 = bn + 1.0/pow(n,2)
    thetaW = s*1.0 / an
    theta_squre = s*1.0*(s-1.0) / (pow(an,2)+bn)
    var = thetaW*(n-2.0)/(6.0*(n-1.0)) + theta_squre * (18.0*pow(n,2)*(3.0*n+2.0)*bn1 - (88.0*pow(n,3)+9.0*pow(n,2)-13.0*n+6.0)) / (9.0*n*pow(n-1,2))
    normH = (pi - thetaL)*1.0 / sqrt(var)
    
    return H, normH

## 1989-Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism
def calculate_Dtajima(s,n,pi):
    ## s: num of segregating sites
    ## n: num of DNA sequences

    ## Tajima's D
    a1 = reduce(lambda x,y: x+1.0/y, range(1,n)); a2 = reduce(lambda x,y: x+1.0/pow(y,2), range(1,n))
    b1 = (n+1.0)/(3.0*n-3.0); b2 = 2.0*(pow(n,2)+n+3.0)/(9.0*n*(n-1.0))
    c1 = b1-1.0/a1; c2 = b2-(n+2.0)/(a1*n)+a2/pow(a1,2)
    e1 = c1/a1; e2 = c2/(pow(a1,2)+a2)
    D = (pi-s/a1)/sqrt(e1*s+e2*s*(s-1.0))
    ## P value for Tajima's D, assuming that D follows the beta distribution
    if n%2 == 0:
        Dmax = (n/(2.0*(n-1))-1.0/a1)/sqrt(e2)
    else:
        Dmax = ((n+1.0)/(2.0*n)-1.0/a1)/sqrt(e2)
    Dmin = (2.0/n-1.0/a1)/sqrt(e2)
    a = Dmin; b = Dmax
    alpha = -(1.0+a*b)*b/(b-a); beta = (1.0+a*b)*a/(b-a)
    func = lambda d: gamma(alpha+beta)*pow(b-d,alpha-1.0)*pow(d-a,beta-1.0)/(gamma(alpha)*gamma(beta)*pow(b-a,alpha+beta-1.0))
    pvalue = 2*min(integrate.quad(func, Dmin, D)[0], integrate.quad(func, D, Dmax)[0])
    #pvalue = np.nan

    return D, pvalue

## Nei, M., and Li, W.H. (1979). MATHEMATICAL-MODEL FOR STUDYING GENETIC-VARIATION IN TERMS OF RESTRICTION ENDONUCLEASES. Proc Natl Acad Sci U S A 76, 5269-5273
def theta_pi_k(hapcount,s,n):
    pi = (hapcount['0'].values * hapcount['1'].values).sum() * 1.0/(n*(n-1.0)/2.0)
    k = s * 1.0 / reduce(lambda x,y: x+1.0/y, range(1,n))

    return pi, k

## Nei, M., and Tajima, F. (1981). DNA POLYMORPHISM DETECTABLE BY RESTRICTION ENDONUCLEASES. Genetics 97, 145-163
def haplotype_diversity(haps):    
    ## Haplotype Diversity (H), H = N/(N-1) * (1-sigma(x^2))
    ## x is the haplotype frequency of each haplotype
    ## N is the sample size (haplotypes)
    ## This measure of gene diversity is analogous to the heterozygosity at a single locus
    haplist = haps.apply(lambda x: "".join(list(x)),axis=0) # assemble each hap to string
    nsample = haps.shape[1]
    sigmax2 = reduce(lambda x,y: x+pow(y,2), [0]+[z*1.0/nsample for z in list(haplist.value_counts())])
    nhap = len(set(haplist))
    H = nsample*1.0/(nsample-1.0)*(1.0-sigmax2)

    return nhap, H

def calculate_one_region_stat(hap_df,hapcount_df,nseq,chromid,start,end,outgroup):
    ## info in the given region
    haps = hap_df[(hap_df['#CHROM']==str(chromid)) & (hap_df['POS']>=int(start)) & (hap_df['POS']<=int(end))].copy()
    hapcount = hapcount_df[(hapcount_df['#CHROM']==str(chromid)) & (hapcount_df['POS']>=int(start)) & (hapcount_df['POS']<=int(end))].copy()

    if haps.empty:
        nmarker = 0; sigtn = 0; thetaPI = 0; thetaK = 0; seg = 0; nhap = 0; H = 0; 
        Dtajima = np.nan; DtajimaP = np.nan; 
        Hfaywu = np.nan; normHfaywu = np.nan;
        Ffuli = np.nan; Dfuli = np.nan
        return nmarker, sigtn, thetaPI, thetaK, seg, nhap, H, Hfaywu, normHfaywu, Ffuli, Dfuli, Dtajima, DtajimaP
    else:
        nmarker = hapcount.shape[0]
        haps.drop(['#CHROM','POS'],axis=1,inplace=True); hapcount.drop(['#CHROM','POS'],axis=1,inplace=True)

        ## check non-biallelic (or missing) sites + homozygotes
        site2rm = list(hapcount[((hapcount['0']+hapcount['1']) <nseq) | (hapcount['1']==0) | (hapcount['0']==0)].index)
        if len(site2rm) == nmarker:
            sigtn = 0; thetaPI = 0; thetaK = 0; seg = 0; nhap = 0; H = 0; 
            Dtajima = np.nan; DtajimaP = np.nan; 
            Hfaywu = np.nan; normHfaywu = np.nan;
            Ffuli = np.nan; Dfuli = np.nan
            return nmarker, sigtn, thetaPI, thetaK, seg, nhap, H, Hfaywu, normHfaywu, Ffuli, Dfuli, Dtajima, DtajimaP
        else:
            if len(site2rm) >0:
                haps.drop(site2rm,inplace=True)
                hapcount.drop(site2rm,inplace=True)
            else:
                pass

            seg = hapcount.shape[0] ## num of segregating site
            thetaPI, thetaK = theta_pi_k(hapcount,seg,nseq)
            nhap, H = haplotype_diversity(haps)

            Dtajima, DtajimaP = calculate_Dtajima(seg,nseq,thetaPI)
            
            if outgroup == 'Y':
                sigtn = hapcount[(hapcount['1']==1)].shape[0] ## num of derived singleton
                Hfaywu, normHfaywu = calculate_Hfaywu(seg,nseq,thetaPI,hapcount)
                Ffuli = calculate_Ffuli_outgroup(seg,nseq,sigtn,thetaPI)
                Dfuli = calculate_Dfuli_outgroup(seg,nseq,sigtn)
            else:
                sigtn = hapcount[(hapcount['0']==1) | (hapcount['1']==1)].shape[0] ## num of singleton
                Hfaywu, normHfaywu = np.nan, np.nan
                Ffuli = calculate_Ffuli_no_outgroup(seg,nseq,sigtn,thetaPI)
                Dfuli = calculate_Dfuli_no_outgroup(seg,nseq,sigtn)

            return nmarker, sigtn, thetaPI, thetaK, seg, nhap, H, Hfaywu, normHfaywu, Ffuli, Dfuli, Dtajima, DtajimaP

## remain required geno data
## convert to ped format
## count alleles
def convert_vcf(vcf,regionfile,window_shift,info,haplist):
    if window_shift == 'target_region':
        windowsize = 5000
    else:
        windowsize = int(window_shift.split('@')[0])
    
    region = pandas.read_csv(regionfile,sep='\s+',header=None,usecols=[0,1,2,3],names=['regionID','chr','start','end'])
    region['chr'] = region['chr'].astype(str)
    region['start'] = region['start'] - windowsize
    region['end'] = region['end'] + windowsize
    region.sort_values(by=['chr','start','end'],ascending=True,inplace=True)
    region.reset_index(inplace=True,drop=True)

    ## merge regions
    if region.shape[0] == 1:
        pass
    else:
        for index in list(region.index)[:-1]:
            chrom1, start1, end1 = list(region.loc[index])[1:]
            chrom2, start2, end2 = list(region.loc[index+1])[1:]

            if ((chrom2 == chrom1) & (start2 <= end1+1)):
                new_end = max(end1,end2)
                region.loc[index+1,'start'] = start1
                region.loc[index+1,'end'] = new_end
                region.drop(index,inplace=True)
            else:
                pass

    ## extract vcf, and convert format
    geno = pandas.concat(list(region.apply(lambda x: vcf[(vcf['#CHROM']==x['chr']) & (vcf['POS']>=x['start']) & (vcf['POS']<=x['end'])].copy(), axis=1)),ignore_index=True)
    
    if geno.empty:
        haps = pandas.DataFrame()
        hapcount = pandas.DataFrame()
        nseq = 0
        return haps, hapcount, nseq
    else:        
        slist = list(geno.columns)[2:]
        mlist = [s for s in slist if info[s]==1]
        ## convert format, for male individuals
        if ((len(mlist) >0) & ('X' in list(geno['#CHROM'].unique()))):
            geno[mlist] = geno[mlist].applymap(lambda x: x[0]).astype('category')
            hapnames = []
            for s in slist:
                if info[s] == 1:
                    hapnames += [s+'_1']
                else:
                    hapnames += [s+'_1',s+'_2']
        else:
            hapnames = [s+'_'+i for s in slist for i in ['1','2']]

        ## to ped format, as type of category
        haps = pandas.DataFrame(geno.apply(lambda x: '|'.join(x[2:]).split('|'),axis=1).tolist(),columns=hapnames).astype('category')
        if (len(set(hapnames) & set(haplist)) < len(hapnames)):
            hapnames = [h for h in hapnames if h in haplist]
            haps = haps[hapnames]
        else:
            pass
        nseq = len(hapnames)
        haps.columns = ['hap'+str(x) for x in range(1,nseq+1)]
        #haps = pandas.DataFrame(geno.apply(lambda x: list(''.join(x[2:]).replace('|','')),axis=1).tolist(),columns=['hap'+str(x) for x in range(1,nseq+1)]).astype('category')
        ## allele count, using bool types to speed up (the same as summing up numbers)
        #hapcount = haps.apply(pandas.value_counts,axis=1) # Matrix[#site number, 2]
        count1 = (haps=='1').apply(sum,axis=1); count0 = (haps=='0').apply(sum,axis=1)
        hapcount = pandas.concat([count0, count1],axis=1).astype('int64')
        hapcount.rename(columns=lambda x: str(x),inplace=True)

        haps['#CHROM'] = hapcount['#CHROM'] = geno['#CHROM'].values
        haps['POS'] = hapcount['POS'] = geno['POS'].values

        ## compress
        haps = haps.astype({'#CHROM':'category','POS':'int32'})
        hapcount = hapcount.astype({'#CHROM':'category','POS':'int32'})

        return haps, hapcount, nseq


def split_window(regionID,chromID,start,end,window_shift):
    windowsize = int(window_shift.split('@')[0])
    stepsize = int(window_shift.split('@')[1])
    overlapsize = windowsize - stepsize

    length = end - start + 1
    bin_num = max(int(math.ceil((length - overlapsize)*1.0 / stepsize)),1)
    ex_len = bin_num * stepsize + overlapsize
    ex_start = int(max(start-(ex_len-length)/2.0, 1.0))
    ex_end = int(end + (ex_len-length)/2.0)

    region = pandas.DataFrame(columns=['regionID','chr','start','end'])
    region['regionID'] = [regionID] * bin_num
    region['chr'] = chromID
    region['start'] = [ex_start + num*stepsize for num in range(bin_num)]
    region['end'] = region['start'] + windowsize - 1
    
    return region

def make_regions(regionfile,window_shift):
    region = pandas.read_csv(regionfile,sep='\s+',header=None,usecols=[0,1,2,3],names=['regionID','chr','start','end'])
    region['chr'] = region['chr'].astype(str)

    if window_shift == 'target_region':
        pass
    else:
        region['tmp'] = region.apply(lambda x: split_window(x['regionID'],x['chr'],x['start'],x['end'],window_shift),axis=1)
        region = pandas.concat(list(region['tmp']),ignore_index=True)
    region.sort_values(by=['chr','start','end'],ascending=True,inplace=True)
    
    return region

def fdr(pvaluelist):
    ## numpy.array format
    ## should be sorted, decreasing order (ascending Pvalues)
    n = len(pvaluelist)
    pvalues = pvaluelist[~np.isnan(pvaluelist)]
    if len(pvalues) <= 1:
        return list(pvaluelist)
    else:
        num = len(pvalues)
        adj_pvalues = pvalues * num / range(1,num+1)
        if adj_pvalues[-1] > 1.0:
            adj_pvalues[-1] = 1.0
        for i in range(num-2,-1,-1):
            adj_pvalues[i] = min(adj_pvalues[i+1],adj_pvalues[i])
        adj_pvalues = list(adj_pvalues) + [np.nan] * (n-num)
        return adj_pvalues

def make_sample_hap(samplefile, hapfile, allsamplelist):
    if samplefile == 'all':
        sampleinfo = pandas.DataFrame(columns=[0,1])
        sampleinfo[0] = allsamplelist
        sampleinfo[1] = 2
    else:
        sampleinfo = pandas.read_csv(samplefile,header=None,sep='\s+')
        if sampleinfo.shape[1] == 1:
            sampleinfo[1] = 2
        else:
            pass
        
        sampleinfo[1] = sampleinfo[1].astype(int)
        if len(set([1,2]) | set(sampleinfo[1])) == 2:
            pass
        else:
            print('something wrong with the genders. only accept 1 and 2.')
            exit()
    
    if hapfile == 'all':
        hapinfo = pandas.DataFrame(columns=[0,1])
        hapinfo[0] = allsamplelist + allsamplelist
        hapinfo[1] = [1]*len(allsamplelist) + [2]*len(allsamplelist)
    else:
        hapinfo = pandas.read_csv(hapfile,header=None,sep='\s+')
        
        hapinfo[1] = hapinfo[1].astype(int)
        if len(set([1,2]) | set(hapinfo[1])) == 2:
            pass
        else:
            print('something wrong with the hap index. only accept 1 and 2.')
            exit()

    s2r = list(set(sampleinfo[0]) & set(hapinfo[0]) & set(allsamplelist))
    if len(s2r) == 0:
        print('NO sample included.')
        exit()
    else:
        pass
    samplelist = [s for s in allsamplelist if s in s2r]
    
    sampleinfo = sampleinfo[sampleinfo[0].isin(s2r)]
    sampleinfo = dict(zip(list(sampleinfo[0]), list(sampleinfo[1])))
    
    hapinfo = hapinfo[hapinfo[0].isin(s2r)]
    haplist = list(hapinfo.apply(lambda x: '{}_{}'.format(x[0],x[1]), axis=1))
    haplist = [s+'_'+i for s in samplelist for i in ['1','2'] if s+'_'+i in haplist]

    return samplelist, sampleinfo, haplist

def main():
    parser = argparse.ArgumentParser(description='Theta_D_H.Est, https://github.com/Shuhua-Group/Theta_D_H.Est for more details')
    parser.add_argument("--gzvcf", type=str, required = True, \
                        help="phased.vcf.gz, format:GT (i.e., 0|1). able to deal with diploids and haploids, seperately.")
    parser.add_argument("--samples", type=str, required = False, default='all', \
                        help="sample_ID.list, samples to be remained. 1 or 2 column: <sample ID> <gender 1/2, optional, useful for X chromosome>, no header line")
    parser.add_argument("--haps", type=str, required = False, default='all', \
                        help="haplotype info, haplotypes to be used. 2 column: <sample ID> <haplotype index, 1/2 denotes the first/second hap>, no header line")
    parser.add_argument("--region", type=str, required = True, \
                        help="region.bed, variants in regions to be used, 4 columns: <region ID> <chrom ID> <start pos> <end pos>, no header line, tab or space sperated")
    parser.add_argument("--window_shift", type=str, required = False, default='target_region', \
                        help="windowsize@increment, for example, 50000@10000.")
    parser.add_argument("--outgroup",type=str, required = False, choices=['Y','N'], default='N', \
                        help="whether the state of ancestral/derived allele is determined, required for FU&Li's and Fay&Wu's tests")
    parser.add_argument("--out", type=str, required = False, default='out.txt', \
                        help="output file name. default: out.txt, will be automatically gzipped.")
    args = parser.parse_args()

    ## log
    with open(args.out+'.logfile','w') as log:
        log.write('{}\n'.format(sys.argv[0]))
        log.write('{}--gzvcf        {}\n'.format(' '*8, args.gzvcf))
        log.write('{}--samples      {}\n'.format(' '*8, args.samples))
        log.write('{}--haps         {}\n'.format(' '*8, args.haps))
        log.write('{}--region       {}\n'.format(' '*8, args.region))
        log.write('{}--window_shift {}\n'.format(' '*8, args.window_shift))
        log.write('{}--outgroup     {}\n'.format(' '*8, args.outgroup))
        log.write('{}--out          {}\n\n'.format(' '*8, args.out))
        
        log.write('Hostname: '+socket.gethostname()+'\n')
        log.write('Working directory: '+os.getcwd()+'\n')
        log.write('Start time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

    ## sample info
    with gzip.open(args.gzvcf) as f:
        headerline = 0
        line = f.readline().decode("utf-8")
        while line[:2] == "##":
            headerline += 1
            line = f.readline().decode("utf-8")

    allsamplelist = line.strip().split('\t')[9:]
    samplelist, sampleinfo, haplist = make_sample_hap(args.samples, args.haps, allsamplelist)  ## gender not consider for haplist
    if ((len(samplelist) <=1) | (len(haplist) <=3)):
        print('No enough sequences, at least 4 sequences are required.')
        exit()
    else:
        pass        

    ## input
    ## read str using Categorical dtypes, to save memory
    datatype = dict(zip(['#CHROM','POS']+samplelist, ['category','int32']+['category']*len(samplelist)))
    vcfdata = pandas.read_csv(args.gzvcf,sep='\t',skiprows=range(headerline),usecols=['#CHROM','POS']+samplelist,dtype=datatype)
    if ((vcfdata.shape[0] <=1) | (vcfdata.shape[1] <=4)):
        print('There may be something wrong with the input data.')
        print('Plz check the #sites and #samples.')
    else:
        pass

    ## convert
    hapdata, hapcountdata, nseq = convert_vcf(vcfdata, args.region, args.window_shift, sampleinfo, haplist)
    if nseq <=3:
        print('No enough sequences, at least 4 sequences are required.')
        exit()
    else:
        pass

    ## get output
    result = make_regions(args.region, args.window_shift)
    result['#sequence'] = nseq
    result['tmp'] = result.apply(lambda x: calculate_one_region_stat(hapdata,hapcountdata,nseq, x['chr'],x['start'],x['end'],args.outgroup),axis=1)

    result = pandas.concat([result[['regionID','chr','start','end','#sequence']], pandas.DataFrame(result['tmp'].tolist(),columns=['#marker','#singleton','ThetaPI','ThetaK','#segregating','#haplotype','Hap_diversity',"Hfaywu","norm_Hfaywu","Ffuli","Dfuli","Dtajima",'Dtajima_P'], index=list(result.index))],axis=1)

    result.sort_values(by='Dtajima_P',ascending=True,inplace=True)
    #result['Dtajima_adj.P'] = stats.p_adjust(FloatVector(result['Dtajima_P']), method='BH')
    result['Dtajima_adj.P'] = fdr(result['Dtajima_P'].values)
    result.sort_values(by=['chr','start','end'],ascending=True,inplace=True)

    result.to_csv(args.out+'.gz',sep='\t',index=None,compression='gzip',na_rep='NA')

    with open(args.out+'.logfile','a') as log:
        log.write("Done.\n")
        log.write("Output "+args.out+'.gz\n')
        log.write('End time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

    print('Done.')
    print('Have a Nice Day!')

if __name__ == '__main__':
    main()

