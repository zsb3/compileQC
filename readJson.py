#!/usr/bin/env python3

import pandas as pd
import json
import sys
import argparse
import warnings
warnings.filterwarnings("ignore")

#Settle command line business
parser=argparse.ArgumentParser()
parser.add_argument("--jsonFile")
parser.add_argument("--coverageFile")
parser.add_argument("--minoralleleFile")
parser.add_argument("--nonACGTNsThreshold")
parser.add_argument("--nonACGTNsSpikeThreshold")
parser.add_argument("--medDepthThreshold")
parser.add_argument("--minorFrequencyThreshold")
parser.add_argument("--minorCountThreshold")
parser.add_argument("--missingSpikeThreshold")

args=parser.parse_args()




#make function to read and parse jsons
def readParse(jsonFile):
    with open(jsonFile) as json_data:
        data = json.load(json_data)
        df = pd.json_normalize(data['results']).T

    query = ['seqName','refName','customNodeAttributes.Nextclade_pango','clade', #'qc.overallScore','qc.overallStatus',
             'totalNonACGTNs','nonACGTNs','qc.frameShifts.frameShifts','missing', #'totalSubstitutions','totalDeletions','totalInsertions','totalMissing',
             'deletions','alignmentScore', 'privateNucMutations.totalPrivateSubstitutions', #'alignmentRange.begin','alignmentRange.end',
             'privateNucMutations.totalPrivateDeletions']
    exportdf = df.loc[query]

    #parse the deletions to be a list of strings? of intervals (e.g. start-stop, start-stop)
    intervList = list()
    for x in range(len(exportdf.loc['deletions'][0])):
        strt = exportdf.loc['deletions'][0][x]['range']['begin']
        stp = exportdf.loc['deletions'][0][x]['range']['end']
        intervList.append(str(strt) + '-' + str(stp))

    #frameShiftList = list()
    #for x in range(len(exportdf.loc['qc.frameShifts.frameShifts'][0])):
    #    strt = exportdf.loc['qc.frameShifts.frameShifts'][0][x]['range']['begin']
    #    stp = exportdf.loc['qc.frameShifts.frameShifts'][0][x]['range']['end']
    #    intervList.append(str(strt) + '-' + str(stp))

    #Repeat for nonACGTNs
    nonACList = list()
    nonACspikecheck = list()
    for x in range(len(exportdf.loc['nonACGTNs'][0])):
        strt = exportdf.loc['nonACGTNs'][0][x]['range']['begin']
        stp = exportdf.loc['nonACGTNs'][0][x]['range']['end']
        nonACList.append(str(strt) + '-' + str(stp))
        nonACspikecheck.append(list(range(strt,stp)))

    #Repeat for missing
    missingList = list()
    missspikecheck = list()
    for x in range(len(exportdf.loc['missing'][0])):
        strt = exportdf.loc['missing'][0][x]['range']['begin']
        stp = exportdf.loc['missing'][0][x]['range']['end']
        missingList.append(str(strt) + '-' + str(stp))
        missspikecheck.append(list(range(strt,stp)))

    exportdf.loc['deletions'][0] = intervList
    exportdf.loc['nonACGTNs'][0] = nonACList
    exportdf.loc['missing'][0] = missingList

    #Find missing in spike
    spike = list(range(21563,25384))
    missingCount = len(set(sum(missspikecheck,[])) & set(spike))
    exportdf = pd.concat([exportdf,pd.DataFrame([missingCount],index = ['missing_in_spike'])])
    #Find nonACGTN in spike
    nonACspikeCount = len(set(sum(nonACspikecheck,[])) & set(spike))
    exportdf = pd.concat([exportdf, pd.DataFrame([nonACspikeCount], index=['nonACGTN_in_spike'])])
    return exportdf[0]


#Make a matching parser for coverage csvs
def covParse(csvFile):
    covInpt = pd.read_csv(csvFile, sep = "\t")

    cov = int(covInpt['Coverage Depth'].median())

    expCovDf = pd.DataFrame([cov], index = ['MedianDepth'])

    return expCovDf

#Pull minor alleles csv
def minParse(csvFile, minThresh=float(args.minorFrequencyThreshold) ):
    minInpt = pd.read_csv(csvFile, sep = "\t")

    minorv = minInpt['Minority_Frequency']
    minorloc = minInpt['Position']

    #count minor alleles above threshold
    minAbvThresh = sum(minorv > minThresh)

    #return a complex list of position and frequency
    posFreqList = list()
    for x in range(len(minorv)) :
        posFreqList.append([str(minorloc[x]) + ':' + str(minorv[x])])

    mindf = pd.DataFrame([[posFreqList]],index=['MinorVariantsFreqs'])
    mindf.loc['minorAboveThreshold'] = minAbvThresh

    return mindf


#json table
jsonTabl = readParse(args.jsonFile)

#depth table
covTabl = covParse(args.coverageFile)

#minor alleles table
minTabl = minParse(args.minoralleleFile)

merge_df = pd.concat([jsonTabl, covTabl, minTabl])


#now do booleans
#Need to make these cutoffs command line addable
NonACTGNs_pf = (merge_df.loc['totalNonACGTNs'] < int(args.nonACGTNsThreshold) ).map({False:'Fail',True:'Pass'})
NonACTGNS_spike_pf = (merge_df.loc['nonACGTN_in_spike'] < int(args.nonACGTNsSpikeThreshold) ).map({False:'Fail',True:'Pass'})
medDepth_pf = (merge_df.loc['MedianDepth'] > int(args.medDepthThreshold) ).map({False:'Fail',True:'Pass'})
minCount_pf = (merge_df.loc['minorAboveThreshold'] < int(args.minorCountThreshold) ).map({False:'Fail',True:'Pass'})
miss_spike_pf = (merge_df.loc['missing_in_spike'] < int(args.missingSpikeThreshold) ).map({False:'Fail',True:'Pass'})

bool_tabl = pd.DataFrame([NonACTGNs_pf,NonACTGNS_spike_pf,medDepth_pf,miss_spike_pf,merge_df.loc['seqName']])

bool_merge = pd.merge(merge_df.T,bool_tabl.T,on='seqName').T


#Now build a report sheet
#Going to attack this by first adding a Description column
#Do this with a mapping dict
descript_dict = {'seqName':'Sequence Name', 'refName':'Reference Name',
                 'customNodeAttributes.Nextclade_pango':'Pangolin Lineage', 'clade':'Nextstrain Clade',
                 'totalNonACGTNs_x':'Total Non-ACGTNs', 'nonACGTNs':'Non-ACGTN Locations',
                 'qc.frameShifts.frameShifts':'Frameshift Locations', 'missing':'Missing Locations',
                 'deletions':'Deletion Locations', 'alignmentScore':'Alignment Score',
                 'privateNucMutations.totalPrivateSubstitutions':'Total Private Mutations',
                 'privateNucMutations.totalPrivateDeletions':'Total Private Deletions',
                 'missing_in_spike_x':'Total Missing Data in Spike','nonACGTN_in_spike_x':'Total Non-ACGTNs in Spike',
                 'MedianDepth_x':'Median Depth', 'totalNonACGTNs_y':'Total Non-ACGTNs (Pass/Fail)',
                 'nonACGTN_in_spike_y':'Total Non-ACGTNs in Spike (Pass/Fail)',
                 'MedianDepth_y': 'Median Depth (Pass/Fail)',
                 'missing_in_spike_y':'Total Missing Data in Spike (Pass/Fail)',
                 'minorAboveThreshold':'Minor Alleles Below Threshold (Pass/Fail)',
                 'MinorVariantsFreqs': 'Minor Alleles Reported From IRMA (Position:Frequency)'}

descr = (bool_merge.index).map(descript_dict)
desc_table = pd.DataFrame(descr,index=bool_merge.index)
format_table = bool_merge.merge(desc_table,left_index=True,right_index=True)
format_table = format_table[['0_y','0_x']]
format_table.rename(columns={'0_x':'Value','0_y':'Description'})


format_table.to_csv('qc_table.csv',index=False,header=False)