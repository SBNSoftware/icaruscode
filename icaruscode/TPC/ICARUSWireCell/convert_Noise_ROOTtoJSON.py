#!/usr/bin/env python
# coding: utf-8
import ROOT
from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile
import sys
import json
import os
import numpy

#Run Script in the same directory with 4 Fillipo's noise files : "dataFFTHistosEE.root","dataFFTHistosEW.root","dataFFTHistosWE.root","dataFFTHistosWW.root"
#Script will produce 8 .json.bz2 files (4 for int noise and 4 for coh). These files can be copied to WireCellData 

def get_freq_amps(hist):
        freq=[]
        amps=[]
        for i in range(0,hist.GetXaxis().GetNbins()):
            freq.append(hist.GetBinLowEdge(i+1)*pow(10,-6))
            amps.append(hist.GetBinContent(i+1)*pow(10,-9)*3.3*pow(10,3)/4095.0)#*pow(10,-9))#*pow(10,-9))
        for i in range(hist.GetXaxis().GetNbins(),0,-1):
            amps.append(hist.GetBinContent(i)*pow(10,-9)*3.3*pow(10,3)/4095.0)#*pow(10,-9))
        for i in range(0,hist.GetXaxis().GetNbins()):
            freq.append(hist.GetBinLowEdge(hist.GetXaxis().GetNbins()+1)*pow(10,-6)+hist.GetBinLowEdge(i+1)*pow(10,-6))
        return freq,amps

def load_noise_spectra(filenames):
    noise_type = "int"
    for filename in filenames:
        
        noises = list()
        modelFile = TFile.Open(filename,"READ")
        tpcName=filename.replace('.root','')
        tpcName=tpcName[-2:]
        powerHist_ind1 = modelFile.Get(noise_type+"powerI1")
        fr_ind1,amps_ind1 = get_freq_amps(powerHist_ind1)
        powerHist_ind2 = modelFile.Get(noise_type+"powerI2")
        fr_ind2,amps_ind2 = get_freq_amps(powerHist_ind2)
        powerHist_coll = modelFile.Get(noise_type+"powerC")
        fr_coll,amps_coll = get_freq_amps(powerHist_coll)
        freq = [fr_ind1,fr_ind2,fr_coll]
        amps = [amps_ind1,amps_ind2,amps_coll]
        gain = [8.811970678500002e-10,8.811970678500002e-10,8.811970678500002e-10]
        period = [400.0,400.0,400.0]
        nsamples = [powerHist_ind1.GetXaxis().GetNbins()*2,powerHist_ind2.GetXaxis().GetNbins()*2,powerHist_coll.GetXaxis().GetNbins()*2]
        shaping = [1.3,1.3,1.3]
        wirelen = [8949.51,3658.0939799999996,3658.0939799999996]
        const = [0.0,0.0,0.0]
        for n in range(0,3):
            ns = { 'period' : period[n],
              'nsamples' : nsamples[n],
              'gain' : gain[n],
              'shaping' : shaping[n],
              'wirelen': wirelen[n],
              'const' : const[n],
              'plane' : n,
              'tpcname' : tpcName,
              'freqs' : list(freq[n]),
              'amps' : list(amps[n])
            }
            noises.append(ns)
        json_format = json.dumps(noises,indent=4)
        int_noise_file = open('icarus_noise_model_'+noise_type+'_TPC'+tpcName+'.json','w')
        int_noise_file.write(json_format)
        int_noise_file.close()
    return 1

def load_coherent_noise_spectra(filenames):
    noise_type = "coh"
    for filename in filenames:
        noises = list()
        modelFile = TFile.Open(filename,"READ")
        tpcName=filename.replace('.root','')
        tpcName=tpcName[-2:]
        powerHist_ind1 = modelFile.Get(noise_type+"powerI1")
        fr_ind1,amps_ind1 = get_freq_amps(powerHist_ind1)
        powerHist_ind2 = modelFile.Get(noise_type+"powerI2")
        fr_ind2,amps_ind2 = get_freq_amps(powerHist_ind2)
        powerHist_coll = modelFile.Get(noise_type+"powerC")
        fr_coll,amps_coll = get_freq_amps(powerHist_coll)
        freq = [fr_ind1,fr_ind2,fr_coll]
        amps = [amps_ind1,amps_ind2,amps_coll]
        gain = [8.811970678500002e-10,8.811970678500002e-10,8.811970678500002e-10]
        period = [400.0,400.0,400.0]
        nsamples = [powerHist_ind1.GetXaxis().GetNbins()*2,powerHist_ind2.GetXaxis().GetNbins()*2,powerHist_coll.GetXaxis().GetNbins()*2]
        shaping = [1.3,1.3,1.3]
        const = [0.0,0.0,0.0]
        
        #Due to the bug in WC code 1st set of channels has only channel 0 (it will have NO coherent noise)
        #the next group would be 31 channel and the rest are all 32 channel group
        ns_0 = { 'period' : period[0],
              'nsamples' : nsamples[0],
              'gain' : gain[0],
              'shaping' : shaping[0],
              'wire-delta': 1.0,
              'const' : const[0],
              'tpcname' : tpcName,
              'freqs' : list(freq[0]),
              'amps' : list(amps[0])
            }
        noises.append(ns_0)
        ns_1 = { 'period' : period[0],
              'nsamples' : nsamples[0],
              'gain' : gain[0],
              'shaping' : shaping[0],
              'wire-delta': 31.0,
              'const' : const[0],
              'tpcname' : tpcName,
              'freqs' : list(freq[0]),
              'amps' : list(amps[0])
            }
        noises.append(ns_1)
        #Group for induction 1 tarts from 1 to accomodate for the bug
        for group in range(1,66):#72
            ns = { 'period' : period[0],
              'nsamples' : nsamples[0],
              'gain' : gain[0],
              'shaping' : shaping[0],
              'wire-delta': 32.0,
              'const' : const[0],
              'tpcname' : tpcName,
              'freqs' : list(freq[0]),
              'amps' : list(amps[0])
            }
            noises.append(ns)
        for group in range(0,175):#180
            ns = { 'period' : period[0],
              'nsamples' : nsamples[0],
              'gain' : gain[0],
              'shaping' : shaping[0],
              'wire-delta': 32.0,
              'const' : const[0],
              'tpcname' : tpcName,
              'freqs' : list(freq[1]),
              'amps' : list(amps[1])
            }
            noises.append(ns)
        for group in range(0,175):#180
            ns = { 'period' : period[0],
              'nsamples' : nsamples[0],
              'gain' : gain[0],
              'shaping' : shaping[0],
              'wire-delta': 32.0,
              'const' : const[0],
              'tpcname' : tpcName,
              'freqs' : list(freq[2]),
              'amps' : list(amps[2])
            }
            noises.append(ns)
        json_format = json.dumps(noises,indent=4)
        int_noise_file = open('icarus_noise_model_'+noise_type+'_TPC'+tpcName+'.json','w')
        int_noise_file.write(json_format)
        int_noise_file.close()
    return 1

print("Start Job")
json_out_int = load_noise_spectra(["dataFFTHistosEE.root","dataFFTHistosEW.root","dataFFTHistosWE.root","dataFFTHistosWW.root"])
json_out_coh = load_coherent_noise_spectra(["dataFFTHistosEE.root","dataFFTHistosEW.root","dataFFTHistosWE.root","dataFFTHistosWW.root"])
#compress files to make tham wirecell compatable
bashCommand = "bzip2 -zf icarus_noise_model_*"
os.system(bashCommand)
print("Finished Job")
