# %%
import multiprocessing
import gemmi
import os
import glob

# %%
# Specify the number of cores to use
num_cores = snakemake.threads


cif_files = glob.glob('tmp/FATCAT_pdb_files/*.cif')
pdb_files = [file.replace('.cif', '.pdb') for file in cif_files]



# %%




def convert_cif_to_pdb(cif_file, pdb_file):
    structure = gemmi.read_structure(cif_file)
    structure.write_pdb(pdb_file)

def parallel_convert_cif_to_pdb(cif_files, pdb_files, num_cores):
    # Create a pool of worker processes
    pool = multiprocessing.Pool(processes=num_cores)
    
    # Map the convert_cif_to_pdb function to the pool of worker processes
    results = pool.starmap(convert_cif_to_pdb, zip(cif_files, pdb_files))
    
    # Close the pool of worker processes
    pool.close()
    
    # Wait for all the worker processes to finish
    pool.join()
    
    # Return the results
    return results


# Call the parallel_convert_cif_to_pdb function
results = parallel_convert_cif_to_pdb(cif_files, pdb_files, num_cores)






# %%
# Delete the CIF files
for cif_file in cif_files:
    os.remove(cif_file)


