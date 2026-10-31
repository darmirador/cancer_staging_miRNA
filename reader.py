#!/usr/bin/env python3
"""
cancer_staging_reader_miRNA.py

Goal
----
For a single miRNA, read per-cancer expression files and
stratify samples by TCGA clinical tumor stage using TCGA_clinical_data.tsv.

Output groups
-------------
normal
stage i
stage ii
stage iii
stage iv
"""

import os
import sys
import csv
import numpy
import pandas
from scipy import stats

# User parameters
pathname = '/Users/darmirador/Desktop/Python Scripts/cancer_staging_miRNA'

if "-i" not in sys.argv:
    print("Usage: python cancer_staging_reader_miRNA.py -i <input_miRNA>")
    sys.exit(1)

miRNAInput = sys.argv[sys.argv.index("-i") + 1]


# Discover input files
directoryFiles = os.listdir(pathname)
filenameInput = []
filenameOutput = []  # unused; kept
cancerTypes = []

for cancerFilename in directoryFiles:
    if cancerFilename.split("_")[0] == miRNAInput + "-output":
        filenameInput.append(cancerFilename)
        cancerTypes.append(cancerFilename.split("_")[1].split(".")[0])


class stage:
    def __init__(self, normal, stage_i, stage_ii, stage_iii, stage_iv):
        self.normal = normal
        self.stage_i = stage_i
        self.stage_ii = stage_ii
        self.stage_iii = stage_iii
        self.stage_iv = stage_iv


class output:
    def __init__(self, barcode, staging, rpm):
        self.barcode = barcode
        self.staging = staging
        self.rpm = rpm


# Load clinical annotations
clinicalPath = os.path.join(pathname, "TCGA_clinical_data.tsv")
tsvExtract = pandas.read_csv(clinicalPath, sep='\t')
patient = tsvExtract["submitter_id"].tolist()
cancerStage = tsvExtract["tumor_stage"].fillna("").astype(str).str.lower().tolist()


# Initialize (create output files + headers)
for cancerType in cancerTypes:
    csvOutput = open(os.path.join(pathname, "staging-" + miRNAInput + "_" + cancerType + ".csv"), "w")
    csvOutput.write("barcode,stage,rpm\n")
    csvOutput.close()

filenameInput.sort()
cancerTypes.sort()


# Main processing: stage stratification
totalCancerTypes = len(filenameInput)
processCounter = 0

for filename in filenameInput:
    currentCancerType = cancerTypes[processCounter]
    print(f"[{processCounter + 1}/{totalCancerTypes}] Processing {currentCancerType} ({filename})...")

    barcodelist = stage([], [], [], [], [])
    transcript_rpm = stage([], [], [], [], [])
    finalList = output([], [], [])

    csvData = pandas.read_csv(os.path.join(pathname, filename))
    cancerState = csvData["tissue"].astype(str).str.lower().tolist()
    patientMatched = csvData["barcode"].tolist()
    transcriptMatched = csvData["rpm"].tolist()

    for counter in range(0, len(patientMatched)):
        barcode = patientMatched[counter]
        sampleType = cancerState[counter]
        expressionValue = transcriptMatched[counter]

        matchedPatient = barcode[:12]

        if sampleType == "normal":
            barcodelist.normal.append(barcode)
            transcript_rpm.normal.append(expressionValue)

        elif sampleType == "cancer":
            if matchedPatient in patient:
                position_patientId = patient.index(matchedPatient)
                currentStage = cancerStage[position_patientId]

                if "stage iv" in currentStage:
                    barcodelist.stage_iv.append(barcode)
                    transcript_rpm.stage_iv.append(expressionValue)
                elif "stage iii" in currentStage:
                    barcodelist.stage_iii.append(barcode)
                    transcript_rpm.stage_iii.append(expressionValue)
                elif "stage ii" in currentStage:
                    barcodelist.stage_ii.append(barcode)
                    transcript_rpm.stage_ii.append(expressionValue)
                elif "stage i" in currentStage:
                    barcodelist.stage_i.append(barcode)
                    transcript_rpm.stage_i.append(expressionValue)

    finalList.barcode.extend(barcodelist.normal)
    finalList.barcode.extend(barcodelist.stage_i)
    finalList.barcode.extend(barcodelist.stage_ii)
    finalList.barcode.extend(barcodelist.stage_iii)
    finalList.barcode.extend(barcodelist.stage_iv)

    for _ in barcodelist.normal:
        finalList.staging.append("normal")
    for _ in barcodelist.stage_i:
        finalList.staging.append("stage i")
    for _ in barcodelist.stage_ii:
        finalList.staging.append("stage ii")
    for _ in barcodelist.stage_iii:
        finalList.staging.append("stage iii")
    for _ in barcodelist.stage_iv:
        finalList.staging.append("stage iv")

    finalList.rpm.extend(transcript_rpm.normal)
    finalList.rpm.extend(transcript_rpm.stage_i)
    finalList.rpm.extend(transcript_rpm.stage_ii)
    finalList.rpm.extend(transcript_rpm.stage_iii)
    finalList.rpm.extend(transcript_rpm.stage_iv)

    for counter in range(0, len(finalList.barcode)):
        csvOutput = open(os.path.join(pathname, "staging-" + miRNAInput + "_" + currentCancerType + ".csv"), "a")
        csvOutput.write(
            finalList.barcode[counter] + "," +
            finalList.staging[counter] + "," +
            str(finalList.rpm[counter]) + "\n"
        )
        csvOutput.close()

    print(f"    Done {currentCancerType}: wrote {len(finalList.barcode)} rows -> staging-{miRNAInput}_{currentCancerType}.csv")
    processCounter += 1

print("All cancer types processed.")