
"""
Author: Lina Sieverling
Affiliation: DKFZ Heidelberg
Aim: A Snakemake workflow to find telomere insertions
Date: Thu Aug 18 17:46:12 CEST 2016
Run: snakemake -s <Snakefile> --configfile <config.yaml> 
Run Example: 

source activate snakemake_telomere_insertions
snakemake -s /home/sieverli/Code/telomere_insertion_analysis/snakemake_telomere_insertions/Snakefile --configfile /abi/data/sieverling/projects/NB_Telomeres/src/config_snakemake_telomere_insertions.yaml --drmaa " -o /ibios/temp2/lina/cluster_messages/TelomereInsertions/NB -j oe -M l.sieverling@dkfz.de -m a -l walltime={params.walltime},mem={params.mem} -N {params.jobname}" --latency-wait 300 --jobs 100 --printshellcmds --keep-going

"""


#---------------------------------------------------------------------------------------
# get PIDs
#---------------------------------------------------------------------------------------
from os import listdir

if config["pids"] == "all":
    pids = [i for i in listdir(config["results_per_pid_dir"]) if not i.startswith('.')]
else:
    pids = config["pids"].split(' ')


#---------------------------------------------------------------------------------------
# remove PIDs without bam files
#---------------------------------------------------------------------------------------

pids_remove = []

for pid_name in pids:
    for sample_name in config["samples"]:
        if not os.path.exists(config["results_per_pid_dir"] + '/' + pid_name + '/alignment/' + sample_name + '_' + pid_name + config["bam_suffix"]):
            print(pid_name + ": alignment bam file for " + sample_name + " sample is missing, skipping this pid!")
            pids_remove.append(pid_name)
            break


pids = [x for x in pids if x not in pids_remove]


#------------------------------------------------------------------
# rule all 
#------------------------------------------------------------------

if len(config["samples"])==2:
    pids_control = pids
else:
    pids_control = []

rule all:
    input:
        expand(config["telomerehunter_dir"] + '/{pid}/tumor_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam', pid=pids),
        expand(config["telomerehunter_dir"] + '/{pid}/control_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam', pid=pids_control),
        expand(config["telomereinsertion_dir"] + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv', pid=pids),
        expand(config["telomereinsertion_dir"] + '/plots/zoomed_in/{pid}_done.txt', pid=pids)



#------------------------------------------------------------------
# run telomerehunter
#------------------------------------------------------------------

if len(config["samples"])==2:
    input_list = [config["results_per_pid_dir"] + '/{pid}/alignment/' + config["samples"][0] + '_{pid}' + config["bam_suffix"], config["results_per_pid_dir"] + '/{pid}/alignment/' + config["samples"][1] + '_{pid}' + config["bam_suffix"]]
    output_list = [config["telomerehunter_dir"] + '/{pid}/' + config["samples"][0] + '_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam', config["telomerehunter_dir"] + '/{pid}/' + config["samples"][1] + '_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam']
    shell_command_addition = "-ibc {input[1]} -pl "
    node_addition = ",nodes=1:ppn=2"
elif len(config["samples"])==1:
    input_list = [config["results_per_pid_dir"] + '/{pid}/alignment/' + config["samples"][0] + '_{pid}' + config["bam_suffix"], config["results_per_pid_dir"]]
    output_list = [config["telomerehunter_dir"] + '/{pid}/' + config["samples"][0] + '_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam']
    shell_command_addition = ""
    node_addition = ""


rule run_telomerehunter:
    input:
        input_list
    output:
        output_list
    params:
        walltime="24:00:00" + node_addition,
        mem="150m",
        jobname="{pid}_telomerehunter",
        sleep_sec_limit=config["sleep_sec_limit"]
    version: "1.0"
    shell:
        "sleep $((1 + RANDOM % {params.sleep_sec_limit}))s; "
        "set +u; source activate telomereEnv; set -u; "
        "module load R/3.4.2; "
        "time telomerehunter -p {wildcards.pid} -o " + config["telomerehunter_dir"] + " -ibt {input[0]} " + shell_command_addition + "-pff all"



#------------------------------------------------------------------
# find discordant reads
#------------------------------------------------------------------

rule find_discordant_reads:
    input:
        config["telomerehunter_dir"] + '/{pid}/{sample}_TelomerCnt_{pid}/{pid}_filtered_intratelomeric.bam',
    output:
        config["telomereinsertion_dir"] + '/tables/{pid}_{sample}_discordant_reads.tsv'
    params:
        walltime="0:59:00",
        mem="150m",
        jobname="{pid}_find_discordant_reads_{sample}",
        sleep_sec_limit=config["sleep_sec_limit"]
    version: "1.0"
    shell:
        "sleep $((1 + RANDOM % {params.sleep_sec_limit}))s; "
        "set +u; source activate telomereEnv; set -u; "                             # deactivate bash strict mode to set virtual environment: https://bitbucket.org/snakemake/snakemake/wiki/FAQ#markdown-header-my-shell-command-fails-with-with-errors-about-an-unbound-variable-whats-wrong
        "python " + config["src_dir"] + "/find_discordant_reads.py -i {input} -o {output}; "
        "set +u; source deactivate; set -u" 


#------------------------------------------------------------------
# add mate mapping quality
#------------------------------------------------------------------

rule add_mate_mapq:
    input:
        discordant_reads = config["telomereinsertion_dir"] + '/tables/{pid}_{sample}_discordant_reads.tsv',
        bam = config["results_per_pid_dir"] + '/{pid}/alignment/{sample}_{pid}' + config["bam_suffix"]
    output:
        config["telomereinsertion_dir"] + '/tables/{pid}_{sample}_discordant_reads_filtered_with_mapq.tsv'
    params:
        walltime="50:00:00",
        mem="100m",
        jobname="{pid}_add_mate_mapq_{sample}",
        sleep_sec_limit=config["sleep_sec_limit"]
    version: "1.0"
    shell:
        "sleep $((1 + RANDOM % {params.sleep_sec_limit}))s; "
        "set +u; source activate telomereEnv; set -u; "                             # deactivate bash strict mode to set virtual environment: https://bitbucket.org/snakemake/snakemake/wiki/FAQ#markdown-header-my-shell-command-fails-with-with-errors-about-an-unbound-variable-whats-wrong
        "python " + config["src_dir"] + "/add_mate_mapq.py -i {input.discordant_reads} -b {input.bam} -o {output}; "
        "set +u; source deactivate; set -u" 



#------------------------------------------------------------------
# count discordant reads
#------------------------------------------------------------------

paired_t_c_flag = False

if len(config["samples"])==2:
    input_list = [config["telomereinsertion_dir"] + '/tables/{pid}_' + config["samples"][0] + '_discordant_reads_filtered_with_mapq.tsv', config["telomereinsertion_dir"] + '/tables/{pid}_' + config["samples"][1] + '_discordant_reads_filtered_with_mapq.tsv']
    tumor_input = "{input[0]}"
    control_input = "{input[1]}"
    paired_t_c_flag = True
elif len(config["samples"])==1:
    input_list = [config["telomereinsertion_dir"] + '/tables/{pid}_' + config["samples"][0] + '_discordant_reads_filtered_with_mapq.tsv']
    tumor_input = "{input[0]}"
    control_input = "NULL"


if not os.path.exists(config["blacklist"]) and not paired_t_c_flag:
    print("Please provide paired tumor-control samples or a blacklist, otherwise no proper filtering for false positives is possible!")


rule count_discordant_reads:
    input:
        input_list
    output:
        windowTable = config["telomereinsertion_dir"] + '/tables/{pid}_discordant_reads_1_kb_windows.tsv'
    params:
        blacklist=config["blacklist"],
        walltime="0:59:00",
        mem="100m",
        jobname="{pid}_count_discordant_reads"
    version: "1.0"
    shell:
        "R-3.2.2 --no-save --slave --args -t " + tumor_input + " -c " + control_input + " -b {params.blacklist} -o {output.windowTable} " + " -f " + config["R_function_file"] + " < " + config["src_dir"] + "/count_discordant_reads.R"



#------------------------------------------------------------------
# get candidate regions
#------------------------------------------------------------------

rule get_candidate_regions:
    input:
        windowTable = config["telomereinsertion_dir"] + '/tables/{pid}_discordant_reads_1_kb_windows.tsv'
    output:
        candidateRegions = config["telomereinsertion_dir"] + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions.tsv'
    params:
        walltime="0:15:00",
        mem="100m",
        jobname="{pid}_get_candidate_regions"
    version: "1.0"
    shell:
        "R-3.2.2 --no-save --slave --args {input.windowTable} {output.candidateRegions} " + str(config["tumor_discordant_read_lower_limit"]) + " " + str(config["control_discordant_read_upper_limit"]) + " " + str(config["consider_blacklist"]) + " " + config["R_function_file"] + " < " + config["src_dir"] + "/get_candidate_regions.R"



#------------------------------------------------------------------
# predict insertion sites
#------------------------------------------------------------------

if len(config["samples"])==2:
    bam = config["results_per_pid_dir"] + '/{pid}/alignment/{sample}_{pid}' + config["bam_suffix"]
elif len(config["samples"])==1:
    bam = config["results_per_pid_dir"] + '/{pid}/alignment/' + config["samples"][0] + '_{pid}' + config["bam_suffix"]

rule find_fusion_reads:
    input: 
        candidateRegions = config["telomereinsertion_dir"] + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions.tsv',
        bam = bam
    output:
        config["telomereinsertion_dir"] + '/clipped_reads/{pid}_{sample}_clipped_reads.tsv'
    params:
        walltime="100:00:00",
        mem="1g",
        jobname="{pid}_find_fusion_reads"
    version: "1.0"
    shell:
        "R-3.2.2 --no-save --slave --args {input.candidateRegions} {input.bam} {output} " + config["R_function_file"] + " < " + config["src_dir"] + "/find_fusion_reads.R"


rule predict_insertion_sites:
    input: 
        candidateRegions = config["telomereinsertion_dir"] + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions.tsv',
        clippedReads = config["telomereinsertion_dir"] + '/clipped_reads/{pid}_' + config["samples"][0] + '_clipped_reads.tsv',
        discordantReads = config["telomereinsertion_dir"] + '/tables/{pid}_' + config["samples"][0] + '_discordant_reads_filtered_with_mapq.tsv'
    output:
        config["telomereinsertion_dir"] + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended.tsv'
    params:
        walltime="0:10:00",
        mem="100m",
        jobname="{pid}_predict_insertion_sites"
    version: "1.0"
    message: "--- {wildcards.pid}: predict insertion sites ---"
    shell:
        "R-3.2.2 --no-save --slave --args {input.candidateRegions} {input.clippedReads} {input.discordantReads} {output} " + config["R_function_file"] + " < " + config["src_dir"] + "/predict_insertion_sites.R"



#------------------------------------------------------------------
# get consensus sequence of insertion (and bases in sequence microhomology)
#------------------------------------------------------------------

rule get_consensus:
    input: 
        candidateRegions = config["telomereinsertion_dir"] + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended.tsv',
        clippedReads = config["telomereinsertion_dir"] + '/clipped_reads/{pid}_' + config["samples"][0] + '_clipped_reads.tsv'
    output:
        config["telomereinsertion_dir"] + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended_with_consensus.tsv'
    params:
        walltime="0:10:00",
        mem="500m",
        jobname="{pid}_get_consensus"
    version: "1.0"
    message: "--- {wildcards.pid}: get consensus ---"
    shell:
        "R-3.2.2 --no-save --slave --args {input.candidateRegions} {input.clippedReads} {output} " + config["R_function_file"] + " < " + config["src_dir"] + "/get_consensus.R"



#------------------------------------------------------------------
# make plots
#------------------------------------------------------------------

rule make_bed_for_visualization:
    input: 
        candidateRegions = config["telomereinsertion_dir"] + '/candidate_region_tables/{pid}_telomere_insertions_candidate_regions_extended.tsv'
    output:
        outfile1 = config["telomereinsertion_dir"] + '/plots/bedfiles/zoomed_out/{pid}_telomere_insertions.bed',
        outfile2 = config["telomereinsertion_dir"] + '/plots/bedfiles/zoomed_in/{pid}_telomere_insertions.bed'
    params:
        walltime="0:10:00",
        mem="100m",
        jobname="{pid}_make_bed"
    version: "1.0"
    shell:
        "R-3.2.2 --no-save --slave --args {input.candidateRegions} {output.outfile1} {output.outfile2} {wildcards.pid}" + " < " + config["src_dir"] + "/make_bed_for_visualization.R"


if len(config["samples"])==2:

    rule visualize_zoomed_in:
        input: 
            bed = config["telomereinsertion_dir"] + '/plots/bedfiles/zoomed_in/{pid}_telomere_insertions.bed',
            tumor_bam = config["results_per_pid_dir"] + '/{pid}/alignment/' + config["samples"][0] + '_{pid}' + config["bam_suffix"],
            control_bam = config["results_per_pid_dir"] + '/{pid}/alignment/' + config["samples"][1] + '_{pid}' + config["bam_suffix"],
            discordant_reads_tumor = config["telomereinsertion_dir"] + '/tables/{pid}_' + config["samples"][0] + '_discordant_reads_filtered_with_mapq.tsv',
            discordant_reads_control = config["telomereinsertion_dir"] + '/tables/{pid}_' + config["samples"][1] + '_discordant_reads_filtered_with_mapq.tsv',
            clipped_reads_tumor = config["telomereinsertion_dir"] + '/clipped_reads/{pid}_' + config["samples"][0] + '_clipped_reads.tsv',
            clipped_reads_control = config["telomereinsertion_dir"] + '/clipped_reads/{pid}_' + config["samples"][1] + '_clipped_reads.tsv'
        output:
            config["telomereinsertion_dir"] + '/plots/zoomed_in/{pid}_done.txt'
        params:
            walltime="10:00:00",
            mem="3g",
            jobname="{pid}_visualize_zoomed_in",
            sleep_sec_limit=config["sleep_sec_limit"]
        version: "1.0"
        shell:
            "sleep $((1 + RANDOM % {params.sleep_sec_limit}))s; "
            "set +u; source activate telomereEnv; set -u; "
            "python " + config["src_dir"] + "/visualize_telomere_insertions.py \
                   --control {input.control_bam} \
                   --tumor {input.tumor_bam}  \
                   --ref /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef/hs37d5.fa \
                   --bed {input.bed} \
                   --samtoolsbin samtools-1.3.1 \
                   --colored_reads_tumor {input.discordant_reads_tumor} \
                   --colored_reads_control {input.discordant_reads_control} \
                   --clipped_reads_tumor {input.clipped_reads_tumor} \
                   --clipped_reads_control {input.clipped_reads_control} \
                   --prefix " + config["telomereinsertion_dir"] + "/plots/zoomed_in/ \
                   --outfile {output}"


elif len(config["samples"])==1:
    rule visualize_zoomed_in:
        input: 
            bed = config["telomereinsertion_dir"] + '/plots/bedfiles/zoomed_in/{pid}_telomere_insertions.bed',
            tumor_bam = config["results_per_pid_dir"] + '/{pid}/alignment/' + config["samples"][0] + '_{pid}' + config["bam_suffix"],
            discordant_reads_tumor = config["telomereinsertion_dir"] + '/tables/{pid}_' + config["samples"][0] + '_discordant_reads_filtered_with_mapq.tsv',
            clipped_reads_tumor = config["telomereinsertion_dir"] + '/clipped_reads/{pid}_' + config["samples"][0] + '_clipped_reads.tsv' 
        output:
            config["telomereinsertion_dir"] + '/plots/zoomed_in/{pid}_done.txt'
        params:
            walltime="10:00:00",
            mem="3g",
            jobname="{pid}_visualize_zoomed_in",
            sleep_sec_limit=config["sleep_sec_limit"]
        version: "1.0"
        shell:
            "sleep $((1 + RANDOM % {params.sleep_sec_limit}))s; "
            "set +u; source activate telomereEnv; set -u; "
            "python " + config["src_dir"] + "/visualize_telomere_insertions.py \
                   --tumor {input.tumor_bam}  \
                   --ref /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef/hs37d5.fa \
                   --bed {input.bed} \
                   --samtoolsbin samtools-1.3.1 \
                   --colored_reads_tumor {input.discordant_reads_tumor} \
                   --clipped_reads_tumor {input.clipped_reads_tumor} \
                   --prefix " + config["telomereinsertion_dir"] + "/plots/zoomed_in/ \
                   --outfile {output}"

