# %%
import os
import shutil
import glob

folder_path = "tmp/FATCAT_aligments_twisted/"

# Check if the folder exists
if os.path.exists(folder_path):
    # Remove the folder and its contents
    shutil.rmtree(folder_path)
    
# Create the folder
os.makedirs(folder_path)




def process_file(file_path):

    output_name, output_name2  = file_path.split('/')[-1][:-14].split('_')
    atom_a_file = 'tmp/FATCAT_aligments_twisted/' + output_name + '_' + output_name + '_' + output_name2 +"_twisted.pdb"
    atom_b_file = 'tmp/FATCAT_aligments_twisted/' + output_name2 + '_' + output_name + '_' + output_name2 + "_twisted.pdb"
    
    with open(file_path, "r") as file:
        with open(atom_a_file, "w") as atom_a:
            with open(atom_b_file, "w") as atom_b:
                for line in file:
                    if line.startswith("ATOM"):
                        elements = line.strip().split()
                        if len(elements) >= 5:
                            if elements[4] == "A":
                                atom_a.write(line)
                            elif elements[4] == "B":
                                atom_b.write(line)
                atom_a.write("END")
                atom_b.write("END")


# %%


file_paths = glob.glob("tmp/FATCAT_aligments/*.opt.twist.pdb")

for file_path in file_paths:
    process_file(file_path)



