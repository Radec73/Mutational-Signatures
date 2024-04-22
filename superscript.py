import csv
import os
import subprocess
import signatureanalyzer as sa
import pandas as pd
import sigProfilerPlotting as sigPlt
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze

input_raw_data_folder = 'RawData'
output_processed_data_folder = 'ProcessedInputData'
output_directory_SignatureAnalyzer = "SignatureAnalyzerOutput"
output_directory_sigProfilerExtractor = "SigProfilerExtractorOutput"
output_directory_sigProfilerAssignment = "SigProfilerAssignmentOutput"


def merge_files(input_folder, output_file):
    first_file = True
    with open(output_file, 'w') as merged_file:
        for file_name in os.listdir(input_folder):
            if file_name.endswith('.maf'):
                file_path = os.path.join(input_folder, file_name)
                with open(file_path, 'r') as input_file:
                    if not first_file:
                        next(input_file)  # Skip header except for the first file
                    merged_file.write(input_file.read())
                first_file = False


def form_input(parameters_string, software):
    default_dict_parameters = {}
    parameters_dictionary = {}
    parameters = parameters_string.strip().split(',')
    if software == "Extractor":
        default_dict_parameters = set_default_parameters_Extractor()
    if software == "Assignment":
        default_dict_parameters = set_default_parameters_Assignment()
    if software == "Analyzer":
        default_dict_parameters = set_default_parameters_Analyzer()
    if software == "SigMiner":
        default_dict_parameters = set_default_parameters_SigMiner()
    for pair in parameters:
        key, value = pair.split('=')
        parameters_dictionary[key] = value
    for key, value in parameters_dictionary.items():
        default_dict_parameters[key] = value
    return default_dict_parameters


def delete_files_in_folder(folder):
    # Check if the folder exists
    if os.path.exists(folder):
        # Iterate over files in the folder
        for file_name in os.listdir(folder):
            file_path = os.path.join(folder, file_name)
            # Check if it's a file and not a directory
            if os.path.isfile(file_path):
                # Delete the file
                os.remove(file_path)


def set_default_parameters_Assignment():
    return {
        'input_type': 'matrix',
        'context_type': '96',
        'collapse_to_SBS96': True,
        'cosmic_version': 3.4,
        'exome': False,
        'genome_build': 'GRCh37',
        'signature_database': None,
        'exclude_signature_subgroups': None,
        'export_probabilities': False,
        'export_probabilities_per_mutation': False,
        'make_plots': False,
        'sample_reconstruction_plots': False,
        'verbose': False
    }


def set_default_parameters_Extractor():
    return {"reference_genome": "GRCh38",
            "opportunity_genome": "GRCh38",
            "cosmic_version": 3.4,
            "context_type": "default",
            "exome": False,
            "minimum_signatures": 1,
            "maximum_signatures": 25,
            "nmf_replicates": 100,
            "resample": True,
            "batch_size": 1,
            "cpu": -1,
            "gpu": False,
            "nmf_init": "random",
            "precision": "single",
            "matrix_normalization": "gmm",
            "seeds": "random",
            "min_nmf_iterations": 10000,
            "max_nmf_iterations": 1000000,
            "nmf_test_conv": 10000,
            "nmf_tolerance": 1e-15,
            "nnls_add_penalty": 0.05,
            "nnls_remove_penalty": 0.01,
            "initial_remove_penalty": 0.05,
            "collapse_to_SBS96": True,
            "clustering_distance": "cosine",
            "export_probabilities": True,
            "make_decomposition_plots": True,
            "stability": 0.8,
            "min_stability": 0.2,
            "combined_stability": 1.0,
            "allow_stability_drop": False,
            "get_all_signature_matrices": False
            }


def set_default_parameters_Analyzer():
    return {
        "maf": "input",
        "outdir": "output",
        "cosmic": 'cosmic3',
        "hg_build": 'hg38.2bit',
        "nruns": 10,
        "verbose": False,
        "plot_results": True,
        "K0": None,
        "objective": 'poisson',
        "max_iter": 10000,
        "del_": 1,
        "tolerance": 1e-6,
        "phi": 1.0,
        "a": 10.0,
        "b": None,
        "prior_on_W": 'L1',
        "prior_on_H": 'L1',
        "report_freq": 100,
        "active_thresh": 1e-2,
        "cut_norm": 0.5,
        "cut_diff": 1.0,
        "cuda_int": 0,
        "tag": ""
    }


def set_default_parameters_SigMiner():
    return {
    "nruns" : "30",
    "nsigs" : "5",
    "data_folder" : "RawData",
    "method" : "brunet",
    "threshold" : 0.05
    }


def SigAnalyzerMafInput(input_folder, output_folder):
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    else:
        delete_files_in_folder(output_folder)

    for file_name in os.listdir(input_folder):
        if file_name.endswith('.vcf'):
            input_file_path = os.path.join(input_folder, file_name)
            output_file_path = os.path.join(output_folder,
                                            file_name[11:-4].replace("filtered", "").replace(".bam", "") + '.maf')

            with open(input_file_path, 'r') as vcf_file, open(output_file_path, 'w') as output_file:
                output_file.write(
                    "Chromosome\tStart_Position\tTumor_Sample_Barcode\tReference_Allele\tTumor_Seq_Allele2\tVariant_Type\n")
                for line in vcf_file:
                    if not line.startswith('##') and not line.startswith('#'):
                        columns = line.strip().split('\t')
                        selected_columns = [columns[i] for i in [0, 1, 2, 3, 4]]
                        selected_columns.append('SNP')
                        selected_columns[0] = columns[0].replace('chr', '')
                        selected_columns[2] = file_name[11:-4].replace("filtered", "").replace(".bam", "")
                        output_file.write('\t'.join(selected_columns) + '\n')


def SigProfilerVcfInput(input_folder, output_folder):
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    else:
        delete_files_in_folder(output_folder)

    for file_name in os.listdir(input_folder):
        if file_name.endswith('.vcf'):
            input_file_path = os.path.join(input_folder, file_name)
            output_file_path = os.path.join(output_folder, file_name[11:].replace("filtered", "").replace(".bam", ""))
            with open(input_file_path, 'r') as vcf_file, open(output_file_path, 'w') as output_file:
                for line in vcf_file:
                    if not line.startswith('##') and not line.startswith('#'):
                        columns = line.strip().split('\t')
                        selected_columns = [columns[i] for i in [0, 1, 2, 3, 4]]
                        selected_columns[0] = columns[0].replace('chr', '')
                        output_file.write('\t'.join(selected_columns) + '\n')


def option_1():
    print("*** SigProfiler Tools ***")
    SigProfilerVcfInput(input_raw_data_folder, output_processed_data_folder)
    print("\nRAW DATA SUCCESSFULLY PROCESSED!\n")

    while True:
        print("Choose your option:")
        print("a - SigProfilerAssignment - Signature Exposure Fitting and Optimization")
        print("e - SigProfilerExtractor - De Novo Signature Identification")
        print("b - Back to menu")
        choice = input()
        if choice == 'a':
            while True:
                print("Options to run SigProfilerAssignment")
                print("h - help and instructions for running SigProfilerAssignment")
                print("r - run SigProfilerAssignment with your own setup of parameters")
                print("t - run SigProfilerAssignment with our default setup of parameters")
                print("b - Back to menu\n")
                choice1 = input()
                if choice1 == 'h':
                    print(f"""INPUT PARAMETERS AND HYPERPARAMETERS WITH INSTRUCTIONS\n
samples                  String      Path to the input somatic mutations file (if using segmentation file/mutational matrix) or input folder (mutation calling file/s).
output                   String      Path to the output folder.
input_type               String      Three accepted input types:

                                     "vcf": if using mutation calling file/s (VCF, MAF, simple text file) as input
                                     "seg:TYPE": if using a segmentation file as input. Please check the required format at https://github.com/AlexandrovLab/SigProfilerMatrixGenerator#copy-number-matrix-generation. The accepted callers for TYPE are the following {"ASCAT", "ASCAT_NGS", "SEQUENZA", "ABSOLUTE", "BATTENBERG", "FACETS", "PURPLE", "TCGA"}. For example:"seg:BATTENBERG"
                                     "matrix": if using a mutational matrix as input

                                     The default value is "matrix".
context_type             String      Required context type if input_type is "vcf". context_type takes which context type of the input data is considered for assignment. Valid options include "96", "288", "1536", "DINUC", and "ID". The default value is "96".
cosmic_version           Float       Defines the version of the COSMIC reference signatures. Takes a positive float among 1, 2, 3, 3.1, 3.2, 3.3, and 3.4. The default value is 3.4.
exome                    Boolean     Defines if the exome renormalized COSMIC signatures will be used. The default value is False.
genome_build             String      The reference genome build, used for select the appropriate version of the COSMIC reference signatures, as well as processing the mutation calling file/s. Supported genomes include "GRCh37", "GRCh38", "mm9", "mm10" and "rn6". The default value is "GRCh37". If the selected genome is not in the supported list, the default genome will be used.
signature_database       String      Path to the input set of known mutational signatures (only in case that COSMIC reference signatures are not used), a tab delimited file that contains the signature matrix where the rows are mutation types and columns are signature IDs.
exclude_signature_subgroups List    Removes the signatures corresponding to specific subtypes to improve refitting (only available when using default COSMIC reference signatures). The usage is explained below. The default value is None, which corresponds to use all COSMIC signatures.
export_probabilities     Boolean     Defines if the probability matrix per mutational context for all samples is created. The default value is True.
export_probabilities_per_mutation Boolean Defines if the probability matrices per mutation for all samples are created. Only available when input_type is "vcf". The default value is False.
make_plots               Boolean     Toggle on and off for making and saving plots. The default value is True.
sample_reconstruction_plots String  Select the output format for sample reconstruction plots. Valid inputs are {'pdf', 'png', 'both', None}. The default value is None.
verbose                  Boolean     Prints detailed statements. The default value is False.\n
*** DEFAULT SETUP *** for Assignment of known mutational signatures to individual samples\n
Function -> cosmic_fit(samples, output, input_type="matrix", context_type="96",
                   collapse_to_SBS96=True, cosmic_version=3.4, exome=False,
                   genome_build="GRCh37", signature_database=None,
                   exclude_signature_subgroups=None, export_probabilities=False,
                   export_probabilities_per_mutation=False, make_plots=False,
                   sample_reconstruction_plots=False, verbose=False)\n""")

                elif choice1 == 'r':
                    print(f"""*** Enter parameters for *** SigProfilerAssignment *** separated with comma (,) without whitespaces \nExample: "
input_type="vcf",
context_type="96",
collapse_to_SBS96=True,
cosmic_version=3.4,
exome=False,
genome_build="GRCh38",
signature_database=None,
exclude_signature_subgroups=[],
export_probabilities=True,
export_probabilities_per_mutation=True,
make_plots=True,
sample_reconstruction_plots=True,
verbose=False\n""")
                    print("ENTER PARAMETERS:\n")
                    parameters_string = input()
                    sigProfAssignment_parameters = form_input(parameters_string, "Assignment")
                    print("\nRunning *** SigProfilerAssignment *** with these parameters:\n" + parameters_string)
                    # try catch blok a vypisat zly input
                    # Analyze.cosmic_fit(samples=output_processed_data_folder,
                    #                    output=output_directory_sigProfilerAssignment,
                    #                    input_type=sigProfAssignment_parameters["input_type"],
                    #                    context_type=sigProfAssignment_parameters["context_type"],
                    #                    collapse_to_SBS96=bool(sigProfAssignment_parameters["context_type"]),
                    #                    cosmic_version=float(sigProfAssignment_parameters["cosmic_version"]),
                    #                    exome=bool(sigProfAssignment_parameters["exome"]),
                    #                    genome_build=sigProfAssignment_parameters["genome_build"],
                    #                    signature_database=sigProfAssignment_parameters["signature_database"],
                    #                    exclude_signature_subgroups=list(sigProfAssignment_parameters["exclude_signature_subgroups"]),
                    #                    export_probabilities=bool(sigProfAssignment_parameters["export_probabilities"]),
                    #                    export_probabilities_per_mutation=bool(sigProfAssignment_parameters["export_probabilities_per_mutation"]),
                    #                    make_plots=bool(sigProfAssignment_parameters["make_plots"]),
                    #                    sample_reconstruction_plots=bool(sigProfAssignment_parameters["sample_reconstruction_plots"]),
                    #                    verbose=bool(sigProfAssignment_parameters["verbose"]))
                    print("Your run was successful. You can find results in " + sigProfAssignment_parameters["output"] + ".\n")

                elif choice1 == 't':
                    print("""Running *** SigProfilerAssignment *** with these parameters:\n
input_type="vcf",
context_type="96",
collapse_to_SBS96=True,
cosmic_version=3.4,
exome=False,
genome_build="GRCh38",
signature_database=None,
exclude_signature_subgroups=[],
export_probabilities=True,
export_probabilities_per_mutation=True,
make_plots=True,
sample_reconstruction_plots=True,
verbose=False\n""")

                    Analyze.cosmic_fit(samples=output_processed_data_folder,
                                       output=output_directory_sigProfilerAssignment,
                                       input_type="vcf",
                                       context_type="96",
                                       collapse_to_SBS96=True,
                                       cosmic_version=3.4,
                                       exome=False,
                                       genome_build="GRCh38",
                                       nnls_add_penalty=0.05,
                                       nnls_remove_penalty=0.01,
                                       initial_remove_penalty=0.05,
                                       signature_database=None,
                                       exclude_signature_subgroups=[],
                                       export_probabilities=True,
                                       export_probabilities_per_mutation=True,
                                       make_plots=True,
                                       sample_reconstruction_plots=True,
                                       verbose=False)
                    print("Run was successful. You can find results in SigProfilerAssignmentOutput folder.\n")

                elif choice1 == 'b':
                    print("Thank you for using this tool.\n")
                    break
                else:
                    print("Invalid choice. Please choose from available choices.\n")
        if choice == 'e':
            while True:
                print("Options to run SigProfilerExtractor")
                print("h - help and instructions for running SigProfilerExtractor")
                print("r - run SigProfilerExtractor with your own setup of parameters")
                print("t - run SigProfilerExtractor with default setup of parameters")
                print("b - Back to menu\n")
                choice2 = input()
                if choice2 == 'h':
                    print(""""INPUT PARAMETERS AND HYPERPARAMETERS WITH INSTRUCTIONS\n
input_type               String      The type of input:

                                     "vcf": used for vcf format inputs.
                                     "matrix": used for table format inputs using a tab separated file.
                                     "bedpe": used for bedpe files with each SV annotated with its type, size bin, and clustered/non-clustered status. Please check the required format at https://github.com/AlexandrovLab/SigProfilerMatrixGenerator#structural-variant-matrix-generation.
                                     "seg:TYPE": used for a multi-sample segmentation file for copy number analysis. Please check the required format at https://github.com/AlexandrovLab/SigProfilerMatrixGenerator#copy-number-matrix-generation. The accepted callers for TYPE are the following {"ASCAT", "ASCAT_NGS", "SEQUENZA", "ABSOLUTE", "BATTENBERG", "FACETS", "PURPLE", "TCGA"}. For example, when using segmentation file from BATTENBERG then set input_type to "seg:BATTENBERG".

output                   String      The name of the output folder. The output folder will be generated in the current working directory.
input_data               String      Path to input folder for input_type:

                                     vcf
                                     bedpe

                                     Path to file for input_type:

                                     matrix
                                     seg:TYPE

reference_genome         String      The name of the reference genome. The default reference genome is "GRCh37". This parameter is applicable only if the input_type is "vcf".
opportunity_genome       String      The build or version of the reference genome for the reference signatures. The default opportunity genome is GRCh37. If the input_type is "vcf", the opportunity_genome automatically matches the input reference genome value. Only the genomes available in COSMIC are supported (GRCh37, GRCh38, mm9, mm10 and rn6). If a different opportunity genome is selected, the default genome GRCh37 will be used.
context_type             String      A string of mutaion context name/names separated by comma (","). The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is "96,DINUC,ID", where "96" is the SBS96 context, "DINUC" is the DINUCLEOTIDE context and ID is INDEL context.
exome                    Boolean     Defines if the exomes will be extracted. The default value is "False".
NMF Replicates

minimum_signatures       Positive Integer   The minimum number of signatures to be extracted. The default value is 1.
maximum_signatures       Positive Integer   The maximum number of signatures to be extracted. The default value is 25.
nmf_replicates           Positive Integer   The number of iteration to be performed to extract each number signature. The default value is 100.
resample                 Boolean           Default is True. If True, add poisson noise to samples by resampling.
seeds                    String            It can be used to get reproducible resamples for the NMF replicates. A path of a tab separated .txt file containing the replicated id and preset seeds in a two columns dataframe can be passed through this parameter. The Seeds.txt file in the results folder from a previous analysis can be used for the seeds parameter in a new analysis. The Default value for this parameter is "random". When "random", the seeds for resampling will be random for different analysis.
NMF Engines

matrix_normalization     String    Method of normalizing the genome matrix before it is analyzed by NMF. Default is value is "gmm". Other options are, "log2", "custom" or "none".
nmf_init                 String    The initialization algorithm for W and H matrix of NMF. Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar' and 'nndsvd_min'. Default is 'random'.
precision                String    Values should be single or double. Default is single.
min_nmf_iterations       Integer   Value defines the minimum number of iterations to be completed before NMF converges. Default is 10000.
max_nmf_iterations       Integer   Value defines the maximum number of iterations to be completed before NMF converges. Default is 1000000.
nmf_test_conv            Integer   Value defines the number number of iterations to done between checking next convergence. Default is 10000.
nmf_tolerance            Float     Value defines the tolerance to achieve to converge. Default is 1e-15.
Execution

cpu                      Integer   The number of processors to be used to extract the signatures. The default value is -1 which will use all available processors.
gpu                      Boolean   Defines if the GPU resource will used if available. Default is False. If True, the GPU resources will be used in the computation. Note: All available CPU processors are used by default, which may cause a memory error. This error can be resolved by reducing the number of CPU processes through the cpu parameter.
batch_size               Integer   Will be effective only if the GPU is used. Defines the number of NMF replicates to be performed by each CPU during the parallel processing. Default is 1.
Solution Estimation Thresholds

stability                Float     Default is 0.8. The cutoff thresh-hold of the average stability. Solutions with average stabilities below this thresh-hold will not be considered.
min_stability            Float     Default is 0.2. The cutoff thresh-hold of the minimum stability. Solutions with minimum stabilities below this thresh-hold will not be considered.
combined_stability       Float     Default is 1.0. The cutoff thresh-hold of the combined stability (sum of average and minimum stability). Solutions with combined stabilities below this thresh-hold will not be considered.
allow_stability_drop     Boolean   Default is False. Defines if solutions with a drop in stability with respect to the highest stable number of signatures will be considered.
Decomposition

cosmic_version           Float     Takes a positive float among 1, 2, 3, 3.1, 3.2, 3.3, and 3.4. Default is 3.4. Defines the version of the COSMIC reference signatures.
make_decomposition_plots Boolean   Defualt is True. If True, Denovo to Cosmic sigantures decompostion plots will be created as a part the results.
collapse_to_SBS96        Boolean   Defualt is True. If True, SBS288 and SBS1536 Denovo signatures will be mapped to SBS96 reference signatures. If False, those will be mapped to reference signatures of the same context.
Others

get_all_signature_matrices Boolean  If True, the Ws and Hs from all the NMF iterations are generated in the output.
export_probabilities     Boolean  Defualt is True. If False, then doesn't create the probability matrix.\n""""")

                elif choice2 == 'r':
                    print(
                        "*** Enter parameters for *** SigProfilerExtractor *** separated with comma (,) without whitespaces \nExample: "
                        "input_type=vcf,output=outputdir,input_data=input_data_dir,reference_genome=GRCh38,"
                        "opportunity_genome=GRCh38,context_type=SBS96, exome=False,minimum_signatures=1, maximum_signatures=5,nmf_replicates=100,"
                        "resample=True,nmf_init=random, matrix_normalization=gmm,min_nmf_iterations=100, "
                        "max_nmf_iterations=1000, nmf_test_conv=500,nmf_tolerance=1e-15, get_all_signature_matrices=False,nnls_add_penalty=0.01, "
                        "nnls_remove_penalty=0.05,make_decomposition_plots=True,initial_remove_penalty=0.05,stability=0.9, min_stability=0.1\n")
                    print("ENTER PARAMETERS:\n")
                    parameters_string = input()
                    sigProfExtractor_parameters = form_input(parameters_string, "Extractor")
                    print("\nRunning *** SigProfilerExtractor *** with these parameters:\n" + parameters_string)
                    # try catch blok a vypisat zly input
                    # Analyze.cosmic_fit(sigProfExtractor_parameters["samples"], sigProfExtractor_parameters["output"],
                    #                    input_type=sigProfExtractor_parameters["input_type"],
                    #                    context_type=sigProfExtractor_parameters["context_type"],
                    #                    collapse_to_SBS96=bool(sigProfExtractor_parameters["context_type"]),
                    #                    cosmic_version=float(sigProfExtractor_parameters["cosmic_version"]),
                    #                    exome=bool(sigProfExtractor_parameters["context_type"]),
                    #                    genome_build=sigProfExtractor_parameters["genome_build"],
                    #                    signature_database=sigProfExtractor_parameters["signature_database"],
                    #                    exclude_signature_subgroups=list(sigProfExtractor_parameters["exclude_signature_subgroups"]),
                    #                    export_probabilities=bool(sigProfExtractor_parameters["export_probabilities"]),
                    #                    export_probabilities_per_mutation=bool(sigProfExtractor_parameters["export_probabilities_per_mutation"]),
                    #                    make_plots=bool(sigProfExtractor_parameters["make_plots"]),
                    #                    sample_reconstruction_plots=bool(sigProfExtractor_parameters["sample_reconstruction_plots"]),
                    #                    verbose=bool(sigProfExtractor_parameters["verbose"]))
                    print("Your run was successful. You can find results in " + sigProfExtractor_parameters[
                        "output"] + ".\n")

                elif choice2 == 't':
                    print("""Running *** SigProfilerExtractor *** with these parameters:\n
 input_type=vcf,
 output=output_directory_sigProfilerExtractor,
 input_data=output_processed_data_folder,
 reference_genome="GRCh38",
 opportunity_genome="GRCh38",
 context_type="SBS96", exome=False,
 minimum_signatures=1, maximum_signatures=5,
 nmf_replicates=100, resample=True,
 nmf_init="random", matrix_normalization="gmm",
 min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=500,
 nmf_tolerance=1e-15, get_all_signature_matrices=False,
 nnls_add_penalty=0.01, nnls_remove_penalty=0.05,
 make_decomposition_plots=True,
 initial_remove_penalty=0.05,
 stability=0.9, min_stability=0.1\n""")
                    # sig.sigProfilerExtractor("vcf",
                    #                      output_directory_sigProfilerExtractor,
                    #                      output_processed_data_folder,
                    #                      reference_genome="GRCh38",
                    #                      opportunity_genome="GRCh38",
                    #                      context_type="SBS96", exome=False,
                    #                      minimum_signatures=1, maximum_signatures=5,
                    #                      nmf_replicates=100, resample=True,
                    #                      nmf_init="random", matrix_normalization="gmm",
                    #                      min_nmf_iterations=100, max_nmf_iterations=1000, nmf_test_conv=500,
                    #                      nmf_tolerance=1e-15, get_all_signature_matrices=False,
                    #                      nnls_add_penalty=0.01, nnls_remove_penalty=0.05,
                    #                      make_decomposition_plots=True,
                    #                      initial_remove_penalty=0.05,
                    #                      stability=0.9, min_stability=0.1)

                elif choice2 == 'b':
                    print("Thank you for using this tool.\n")
                    break
                else:
                    print("Invalid choice. Please choose from available choices.\n")

        elif choice == 'b':
            print("\nThank you for using this tool.\n")
            break
        else:
            print("\nInvalid choice. Please choose from available choices.\n")


def option_2():
    print("*** SignatureAnalyzer ***")
    SigAnalyzerMafInput(input_raw_data_folder, output_processed_data_folder)
    merged_output_file = "SignatureAnalyzerMerged_output_file.maf"
    merge_files(output_processed_data_folder, merged_output_file)
    print("\nRAW DATA SUCCESSFULLY PROCESSED!\n")

    while True:
        print("Choose your option:")
        print("h - help")
        print("r - run signatureanalyzer with your own options and hyperparameters")
        print("t - run all data files in data folder separately")
        print("x - run merged file formed from every datafile with default hyperparameters")
        print("e - Exit\n")
        choice = input()
        if choice == 'h':
            print("*** HELP & RUN OPTIONS ***\n")
            command = "signatureanalyzer -h"
            try:
                output = subprocess.check_output(command, shell=True, text=True)
                print(output)
            except subprocess.CalledProcessError as e:
                print("Error:", e)

        elif choice == 'r':
            print(
                "*** Enter parameters for *** SignatureAnalyzer *** separated with comma (,) without whitespaces \nExample: \n"
                "maf=input.maf,"
                "outdir=output_directory,"
                "cosmic=cosmic3,"
                "hg_build='hg38.2bit,"
                "nruns=30,"
                "objective=poisson\n")

            if not os.path.exists(output_directory_SignatureAnalyzer):
                os.makedirs(output_directory_SignatureAnalyzer)
            else:
                delete_files_in_folder(output_directory_SignatureAnalyzer)

            print("ENTER PARAMETERS:\n")
            parameters_string = input()
            sigAnalyzer_parameters = form_input(parameters_string, "Analyzer")
            print("\nRunning *** Signature Analyzer *** with these parameters:\n" + parameters_string)
            # sa.run_maf(maf=sigAnalyzer_parameters["maf"],
            #            outdir=sigAnalyzer_parameters["outdir"],
            #            cosmic=sigAnalyzer_parameters["cosmic"],
            #            hg_build=sigAnalyzer_parameters["hg_build"],
            #            nruns=sigAnalyzer_parameters["nruns"],
            #            objective=sigAnalyzer_parameters["objective"],
            #            verbose= sigAnalyzer_parameters["verbose"],
            #            plot_results= sigAnalyzer_parameters["plot_results"],
            #            K0=sigAnalyzer_parameters["K0"],
            #            max_iter=sigAnalyzer_parameters["max_iter"],
            #            del_=sigAnalyzer_parameters["del_"],
            #            tolerance=sigAnalyzer_parameters["tolerance"],
            #            phi=sigAnalyzer_parameters["phi"],
            #            a=sigAnalyzer_parameters["a"],
            #            b= sigAnalyzer_parameters["b"] ,
            #            prior_on_W=sigAnalyzer_parameters["prior_on_W"],
            #            prior_on_H=sigAnalyzer_parameters["prior_on_H"],
            #            report_freq=sigAnalyzer_parameters["report_freq"],
            #            active_thresh=sigAnalyzer_parameters["active_thresh"],
            #            cut_norm=sigAnalyzer_parameters["cut_norm"],
            #            cut_diff=sigAnalyzer_parameters["cut_diff"],
            #            cuda_int=sigAnalyzer_parameters["cuda_int"],
            #            tag=sigAnalyzer_parameters["tag"])
            print("Your run was successful. You can find results in " + sigAnalyzer_parameters["output"] + ".\n")

        elif choice == 't':
            if not os.path.exists(output_directory_SignatureAnalyzer):
                os.makedirs(output_directory_SignatureAnalyzer)
            else:
                delete_files_in_folder(output_directory_SignatureAnalyzer)
            for filename in os.listdir(output_processed_data_folder):
                if os.path.isfile(os.path.join(output_processed_data_folder, filename)):
                    sa.run_maf(os.path.join(output_processed_data_folder, filename),
                               outdir=f"{output_directory_SignatureAnalyzer}/{filename.split('.')[0]}_output",
                               cosmic='cosmic3',
                               hg_build='hg38.2bit',
                               nruns=50,
                               objective="poisson",
                               verbose=True)
                print("Your file: " + filename + " has been successfully analyzed!\n")

        elif choice == 'x':
            if not os.path.exists(output_directory_SignatureAnalyzer):
                os.makedirs(output_directory_SignatureAnalyzer)
            else:
                delete_files_in_folder(output_directory_SignatureAnalyzer)
            # mafi, spectra_sbs = sa.get_spectra_from_maf(pd.read_csv(merged_output_file, sep="\t", header='infer'), cosmic='cosmic3', hgfile='hg38.2bit')

            # print(spectra_sbs)
            sa.run_maf(merged_output_file,
                       outdir=f"{merged_output_file.split('.')[0]}_output",
                       cosmic='cosmic3',
                       hg_build='hg38.2bit',
                       nruns=10,
                       objective="poisson")
            print("Your file: " + merged_output_file + " has been successfully analyzed!\n")

        elif choice == 'e':
            print("Thank you for using this tool.\n")
            break
        else:
            print("Invalid choice. Please choose from available choices.\n")


def option_3():
    r_script_path = "sigminer.R"
    print("*** SigMiner ***")
    print("\nRAW DATA SUCCESSFULLY PROCESSED & READY!\n")

    while True:
        print("Choose your option:")
        print("h - Help and instructions for running Sigminer")
        print("r - Run SigMiner code with your own setup of parameters")
        print("d - Run SigMiner code with default parameters")
        print("b - Back to menu")
        choice = input()
        if choice == 'h':
            print(f"""
            nruns - a numeric giving the number of run to perform for each value in range, nrun set to 30~50 is enough to achieve robust result. Default is 30.
            method - specification of the NMF algorithm. Use 'brunet' as default. Available methods for NMF decompositions are 'brunet', 'lee', 'ls-nmf', 'nsNMF', 'offset'.
            nsigs - number of signatures, if you do not put any value optimum value is automatically determined by default
            data_folder - folder with input data in vcf format. Defaulty set to our data.\n""")

        if choice == 'r':
            print(
                "*** Enter parameters for *** SigMiner *** separated with comma (,) without whitespaces \nExample: "
                "nruns=30,nsigs=5,method=brunet,threshold=0.05\n")
            print("ENTER PARAMETERS:\n")
            parameters_string = input()
            sigMiner_parameters = form_input(parameters_string, "SigMiner")

            print("\nRunning *** SigMiner - De Novo Extraction of Signatures & Signature Exposure Fitting and Optimization *** with these parameters:\n" + parameters_string)
            try:
                nruns = int(sigMiner_parameters["nruns"])
                nsigs = int(sigMiner_parameters["nsigs"])
                method = str(sigMiner_parameters["method"])
                treshold = str(sigMiner_parameters["threshold"])
            except (KeyError, ValueError) as e:
                print(f"Error: Invalid parameter type or missing parameter: {e}")
                return
            try:
                subprocess.run(["Rscript", r_script_path, nruns, nsigs, treshold, method],check=True)
                print("\nR script execution completed successfully.\n")
                print("You can find results in SigMinerOutput.\n")
            except subprocess.CalledProcessError as e:
                print(f"Error: R script execution failed with return code {e.returncode}.\n")

        if choice == 'd':

            sigMiner_parameters = set_default_parameters_SigMiner()
            print("\nRunning *** SigMiner - De Novo Extraction of Signatures & Signature Exposure Fitting and Optimization *** with default parameters:\n")
            try:
                subprocess.run(["Rscript", r_script_path,  str(sigMiner_parameters["nruns"]), str(sigMiner_parameters["nsigs"]), str(sigMiner_parameters["threshold"]), str(sigMiner_parameters["method"])],check=True)
                print("\nR script execution completed successfully.\n")
                print("You can find results in SigMinerOutput.\n")
            except subprocess.CalledProcessError as e:
                print(f"\nError: R script execution failed with return code {e.returncode}.\n")

        elif choice == 'b':
            print("Thank you for using this tool.\n")
            break
        else:
            print("Invalid choice. Please choose from available choices.\n")

def main():
    while True:
        print("\nMutational Signatures Analysis Tools")
        print("1 - SigProfiler")
        print("2 - SignatureAnalyzer")
        print("3 - SigMiner")
        print("4 - Exit")

        choice = input("Choose software for analyzing samples:\n")

        if choice == '1':
            option_1()
        elif choice == '2':
            option_2()
        elif choice == '3':
            option_3()
        elif choice == '4':
            print("Exiting program. Thank you for using this tool.")
            break
        else:
            print("Invalid choice. Please enter a number between 1 and 5.")


if __name__ == "__main__":
    main()
