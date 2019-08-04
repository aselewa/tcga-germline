#!/bin/python

import sys
import json
import numpy as np

# Steps
# 1. Get JSON file from GDC-legacy, extract file ID and TCGA barcode (extract_from_json.py)
# 2. Get file ID and TCGA barcode from GDC for somatic mutations
# 3. Use this script to merge the two (each individual's genotype and SSM file IDs)

outName = sys.argv[2]

jfile = open(sys.argv[1])
json_file = json.load(jfile)
fnames = []
ids = []
project_id = []
gender = []
age = []
race = []
for i in json_file:
    fnames.append(i['file_id'])

    ids.append(i['associated_entities'][0]['entity_submitter_id'])
    project_id.append(i['cases'][0]['project']['project_id'])

    if 'demographic' in i['cases'][0].keys():
        gender.append(i['cases'][0]['demographic']['gender'])
        race.append(i['cases'][0]['demographic']['race'])
    else:
        gender.append('NA')
        race.append('NA')
        
    if 'diagnoses' in i['cases'][0].keys():
        age.append(i['cases'][0]['diagnoses'][0]['age_at_diagnosis'])
    else:
        age.append('NA')

keep = []
for i in ids:
    types = i.split('-')[3]
    a = int(list(types)[0]+list(types)[1])
    if a <= 9:
        keep.append(False)
    else:
        keep.append(True)

keep = np.array(keep)

ids = np.array(ids)[keep]
gender = np.array(gender)[keep]
age = np.array(age)[keep]
fnames = np.array(fnames)[keep]
project_id = np.array(project_id)[keep]
race = np.array(race)[keep]

age_mat = np.column_stack((ids, fnames, project_id, age, gender, race)) 
np.savetxt(outName,age_mat,delimiter="\t", fmt="%s")

