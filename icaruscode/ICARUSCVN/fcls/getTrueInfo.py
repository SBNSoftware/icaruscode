import re
import glob
import os

def normalize(value):
    if value > 2:
        return 3
    return value 

def normalize2(value):
    if value < 0:
        return 1
    return 0 

directory = '/pnfs/icarus/scratch/users/fwieler/pixelMaps_score2/W_set'  # current directory
pattern = re.compile(r'e(\d+)')
files = list(glob.iglob(directory + "/**/*.gz", recursive=True))

output_lines = []
extracted_e_values = []
i = 0

for filename in files:
    match = pattern.search(filename)
    
    if not match:
        continue
    i += 1
    print(f'Processing [{i}/{len(files)}]')
    e_value = f"e{match.group(1)}"
    extracted_e_values.append(e_value)

    file_ID = filename.split("/")[-1][:-3]
    infofile = file_ID + '.info'
    dir_path = os.path.dirname(filename)
     
    try: 
        with open(dir_path + "/" + infofile, 'r') as f:
            info = f.readlines()
    except FileNotFoundError:
        continue
    
    if not len(info):
        continue
    
    fInt = int(info[0].strip())
    flavour = fInt // 4
    if flavour == 3:
        flavour = 2;
    interaction = fInt % 4

    fNuEnergy = float(info[1].strip())
    fLepEnergy = float(info[2].strip())
    fRecoNueEnergy = float(info[3].strip())
    fRecoNumuEnergy = float(info[4].strip())
    fEventWeight = float(info[6].strip())

    fNuPDG = normalize2(int(info[7].strip()))
    fNProton = normalize(int(info[8].strip()))
    fNPion = normalize(int(info[9].strip()))
    fNPizero = normalize(int(info[10].strip()))
    fNNeutron = normalize(int(info[11].strip()))
    
    # Create a CSV line: eXX, flavour, interaction, energies, weights, normalized particles
    line = f"{file_ID},{flavour},{interaction},{fNuEnergy},{fLepEnergy},{fRecoNueEnergy},{fRecoNumuEnergy},{fEventWeight},{fNuPDG},{fNProton},{fNPion},{fNPizero},{fNNeutron}"
    output_lines.append(line)
    
# Write all results to a file
with open("event_info_wset_score2.csv", "w") as outfile:
    outfile.write("fid,flavour,interaction,fNuEnergy,fLepEnergy,fRecoNueEnergy,fRecoNumuEnergy,fEventWeight,fNuPDG,fNProton,fNPion,fNPizero,fNNeutron\n")
    for line in output_lines:
        outfile.write(line + "\n")

print("\nSaved event info to 'event_info_wset_score2.csv'.")
