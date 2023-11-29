# %%
def process_file(file_path):

    output_name, output_name2  = file_path.split('/')[-1][:-14].split('_')
    atom_a_file = output_name + "_twisted.pdb"
    atom_b_file = output_name2 + "_twisted.pdb"
    
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


file_paths = snakemake.input #glob.glob("../git_repo_cloned/FATCAT/Examples_FATCAT/*.opt.twist.pdb")

for file_path in file_paths:
    process_file(file_path)



