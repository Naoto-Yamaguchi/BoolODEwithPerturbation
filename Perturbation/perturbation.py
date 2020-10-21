import BoolODE as bo
import numpy as np
import pandas as pd
import re

def perturbation(    
    model_definition, 
    simulation_time,
    num_cells,
    sampling_time,
    perturbed_transcription, 
    perturbed_translation,
    output_dir,
    output_filename,
    modeltype
):
    # parameter settings
    model_dir = "data"
    gs = bo.GlobalSettings(model_dir, output_dir, True, False, modeltype)
    model_df = pd.read_csv(model_dir + "/" + model_definition, sep='\t', engine='python')
    genelist = list(model_df['Gene'].values)

    # parameter validation
    if max(sampling_time) > simulation_time * 100:
        print("Sampling_time is larger than simulation_time")
        return

    if type(perturbed_transcription) is float:
        rate = perturbed_transcription
        perturbed_transcription = []
        for gene in genelist:
            perturbed_transcription.append({ gene: rate })
    
    if type(perturbed_translation) is float:
        rate = perturbed_translation
        perturbed_translation = []
        for gene in genelist:
            perturbed_translation.append({ gene: rate })

    if len(perturbed_transcription) != 0 and len(perturbed_translation) != 0:
        print("Either of transcription or translation is allowed")
        return

    # Normal simulation
    js_normal = bo.JobSettings(
    [{
        "name": model_definition[:-4],
        "model_definition": model_definition,
        "model_initial_conditions": model_definition + "_ics.txt",
        "modeltype": modeltype,
        "simulation_time": simulation_time,
        "num_cells": num_cells,
        "do_parallel": False,
        "sample_cells": False,
        "perturbation": False,
        "perturbation_control": False,
        "perturbed_transcription": {},
        "perturbed_translation": {},
        "perturbation_input": "", 
        "perturbation_sampling_time": [],
        "perturbation_sampling_filename": ""
    }]
    )
    boolodejobs = bo.BoolODE(job_settings=js_normal, global_settings=gs, postproc_settings="")
    boolodejobs.execute_jobs()
    print("\n \n \n")
    print("Normal simulation ... Done")

    print("==========================================")
    print("Starting perturbation control simulation...")
    print("==========================================")

    # Perturbation control simulation
    js_ctrl = bo.JobSettings(
    [{
        "name": model_definition[:-4] + "-perturbation-control",
        "model_definition": model_definition,
        "model_initial_conditions": model_definition + "_ics.txt",
        "modeltype": modeltype,
        "simulation_time": simulation_time,
        "num_cells": num_cells,
        "do_parallel": False,
        "sample_cells": False,
        "perturbation": True,
        "perturbation_control": True,
        "perturbation_input": model_definition[:-4] + "/simulations/",

        "perturbation_sampling_time": sampling_time,
        "perturbation_sampling_filename": output_dir + "/" + model_definition[:-4] + "-" + output_filename
    }]
    )
    boolodejobs = bo.BoolODE(job_settings=js_ctrl, global_settings=gs, postproc_settings="")
    boolodejobs.execute_jobs()


    
    print("Starting perturbation simulation...")

    # Perturbation simulation
    if len(perturbed_transcription) != 0:
        print("\n \n \n")
        print("==========================================")
        print("Starting transcription perturbation simulation...")
        print("==========================================")
        for trans in perturbed_transcription:
            if list(trans.keys())[0] in genelist:
                print("\n")
                print("---------")
                print("{} perturbation".format(list(trans.keys())[0]))
                print("---------")

                js_transcription = bo.JobSettings(
                [{
                    "name": model_definition[:-4] + "-perturbation-transcription-" + list(trans.keys())[0],
                    "model_definition": model_definition,
                    "model_initial_conditions": model_definition + "_ics.txt",
                    "modeltype": modeltype,
                    "simulation_time": simulation_time,
                    "num_cells": num_cells,
                    "do_parallel": False,
                    "sample_cells": False,
                    "perturbation": True,
                    "perturbation_control": False,
                    "perturbed_transcription": trans,
                    #"perturbed_translation": ,
                    "perturbation_input": model_definition[:-4] + "/simulations/",
                    "perturbation_sampling_time": sampling_time,
                    "perturbation_sampling_filename": output_dir + "/" + model_definition[:-4] + "-" + output_filename
                }]
                )
                boolodejobs = bo.BoolODE(job_settings=js_transcription, global_settings=gs, postproc_settings="")
                boolodejobs.execute_jobs()
            else:
                print("Please specify a genename in the genes in your model definition")
        
    elif len(perturbed_translation) != 0:
        print("\n \n \n")
        print("==========================================")
        print("Starting translation perturbation simulation...")
        print("==========================================")
        for trans in perturbed_translation:
            if list(trans.keys())[0] in genelist:
                print("\n")
                print("---------")
                print("{} perturbation".format(list(trans.keys())[0]))
                print("---------")

                # job実行
                js_translation = bo.JobSettings(
                [{
                    "name": model_definition[:-4] + "-perturbation-translation-" + list(trans.keys())[0],
                    "model_definition": model_definition,
                    "model_initial_conditions": model_definition + "_ics.txt",
                    "modeltype": modeltype,
                    "simulation_time": simulation_time,
                    "num_cells": num_cells,
                    "do_parallel": False,
                    "sample_cells": False,
                    "perturbation": True,
                    "perturbation_control": False,
                    #"perturbed_transcription": ,
                    "perturbed_translation": trans,
                    "perturbation_input": model_definition[:-4] + "/simulations/",

                    "perturbation_sampling_time": sampling_time,
                    "perturbation_sampling_filename": output_dir + "/" + model_definition[:-4] + "-" + output_filename
                }]
                )
                boolodejobs = bo.BoolODE(job_settings=js_translation, global_settings=gs, postproc_settings="")
                boolodejobs.execute_jobs()
            else:
                print("Please specify a genename in the genes in your model definition")
        
            
    E = np.load(output_dir + "/" + model_definition[:-4] + "-" + output_filename)

    refNetwork_filepath = output_dir + "/" + model_definition[:-4] +  "/refNetwork.csv" # perturbationの有無はgene_names, grn_mtxに影響を与えない。
    with open(refNetwork_filepath, 'r') as f:
        df = pd.read_csv(refNetwork_filepath)
        gene_names = sorted(list(set(df["Gene2"]) | set(df["Gene1"])))
        grn_mtx = np.zeros((len(gene_names), len(gene_names)))
        
        for row in df.iterrows():
            i = gene_names.index(row[1]["Gene2"])
            j = gene_names.index(row[1]["Gene1"])
            if row[1]["Type"] == "+":
                grn_mtx[i][j] = 1
            elif row[1]["Type"] == "-":
                grn_mtx[i][j] = -1

        grn_mtx_df = pd.DataFrame(grn_mtx, dtype="int32")
        grn_mtx_df.index = gene_names
        grn_mtx_df.columns = gene_names
        grn_mtx_df.to_csv(output_dir + "/" + "GRN_mtx_columns_regulate_rows.csv") 
        print("\n \n \n ")
        print("Writing GRN matrix to csv file...") 
        
    print("It's all Done!!")
    return E, gene_names, grn_mtx
    
# Run like this

# import BoolODE as bo
# import numpy as np
# import pandas as pd
# from Perturbation import perturbation as pt
# # user setting
# model_definition = "HSC.txt" #"dyn-cycle.txt"

# # perturbation parameter validation
# # either of perturbation_transcription and perturbation_translation should be []
# # either of perturbation_transcription and perturbation_translation should be float or list of dict
# perturbed_transcription = 0.1 #[] #[{ 'g1': 10.0 }, { 'g2': 10.0 }]
# perturbed_translation = [] #[{ 'Gata2': 0.5 }, { 'Cebpa': 0.3 }] # [{ 'g3': 5.0 }, { 'g4': 5.0 }] # 0.5

# simulation_time = 1 # * 100 step
# num_cells = 10 # = # of sampling_cells
# sampling_time = [10, 30, 50] # 1~simulation_time*100-1?
# output_dir = "Run1736" # change every simulation
# output_filename = "PerturbationSampling.npy" # npy file
# modeltype = "hill" # or heaviside

# # default
# do_parallel = False
# sample_cells = False

# # execution from module
# E, gene_names, grn_mtx = pt.perturbation(
#     model_definition, 
#     simulation_time,
#     num_cells,
#     sampling_time,
#     perturbed_transcription, 
#     perturbed_translation,
#     output_dir,
#     output_filename,
#     modeltype
# )