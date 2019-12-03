#!/usr/bin/env nextflow


if (params.help) {
    log.info """
    ===========================================
      QFO ORTHOLOGY BENCHMARKING PIPELINE
    ===========================================
    Usage:
    Run the pipeline with default parameters:
    nexflow run benchmark.nf

    Run with user parameters:
    nextflow run benchmark.nf --input {orthology.predictions} --participant_id {tool.name} --results_dir {results.dir}

    Mandatory arguments:
        --input                 Predicted orthologs in TSV or orthoxml file

        --participant_id        Name of the tool / method


    Additional options:
        --challenges_ids        List of benchmarks / challenges to be run by the pipeline.
                                Separate challenges with space, and don't forget to quote.
                                Defaults to all available benchmarks:
                                "${params.challenges_ids}"

        --event_year            QfO Reference Proteomes release, defaults to 2018

        --community_id          Name or OEB permanent ID for the benchmarking community

        --go_evidences          Evidence filter of GO annotation used in the GO benchmark
                                Defaults to experimental annotations

        --goldstandard_dir      Dir that contains the benchmarking datasets needed to execute
        --assess_dir            Dir where the result data for the benchmark are stored (e.g. outdir of previous runs)

        --results_dir           Base result for all the following output directories, unless overwritten
        --validation_result     The output directory where the results from validation step will be saved (currently not used)
        --assessment_results    The output directory where the results from the computed metrics step will be saved
        --outdir                The output directory where the consolidation of the benchmark will be saved
        --statsdir              The output directory with nextflow statistics
        --data_model_export_dir The output dir where json file with benchmarking data model contents will be saved
        --otherdir              The output directory where custom results will be saved (no directory inside)

    Flags:
        --help                  Display this message¬
    """.stripIndent()

    exit 1
}

log.info """
         ==============================================
          QFO ORTHOLOGY BENCHMARKING PIPELINE
         ==============================================
         input file: ${params.input}
         method name : ${params.participant_id}
         goldstandard path (refset): ${params.goldstandard_dir}
         benchmarking community = ${params.community_id}
         selected benchmarks: ${params.challenges_ids}
         Evidence filter for GO benchmark: ${params.go_evidences}
         Public Benchmark results: ${params.assess_dir}
         validation results directory: ${params.validation_result}
         assessment results directory: ${params.assessment_results}
         consolidated benchmark results directory: ${params.outdir}
         Statistics results about nextflow run: ${params.statsdir}
         Benchmarking data model file location: ${params.data_model_export_dir}
         Directory with community-specific results: ${params.otherdir}
         """
    .stripIndent()


//input
predictions = file(params.input)
method_name = params.participant_id.replaceAll("\\s","_")
refset_dir = Channel.fromPath(params.goldstandard_dir, type: 'dir')
benchmarks = params.challenges_ids
community_id = params.community_id
benchmark_data = file(params.assess_dir)
go_evidences = params.go_evidences
tree_clades = Channel.from("Luca", "Vertebrata", "Fungi", "Eukaryota")
tree_clades2 = Channel.from("Luca", "Vertebrata", "Fungi", "Eukaryota")
genetree_sets = Channel.from("SwissTrees", "TreeFam-A")
tree_clades0 = Channel.from("Eukaryota", "Fungi", "Bacteria")

//output
//validation_out = file(params.validation_result)
assessment_out = file(params.assessment_results)
result_file_path = file(params.outdir)
data_model_export_dir = file(params.data_model_export_dir)
otherdir = file(params.otherdir)

refset_dir.into { refset_dir_validate; refset_dir_convert; refset_dir_go_test; refset_dir_ec_test; refset_dir_std; refset_dir_g_std; refset_dir_g_std_v2; refset_dir_geneTrees }

/*
 * validate input file
 */
process validate_input_file {
    label "py"

    input:
    file predictions
    file refset_dir from refset_dir_validate
    val benchmarks
    val community_id
    val method_name

    output:
    val task.exitStatus into EXIT_STAT
    file "participant.json" into PARTICIPANT_STUB

    """
    echo "DEBUG: mountpoints"
    mount
    echo "END DEBUG"
    echo
    #ls -l /home/vre/projects/wrapper_test/javi_tmp/input/SonicParanoid_default.rels.raw.gz
    ls -la $refset_dir
	echo
    /benchmark/validate.py --com $community_id --challenges_ids "$benchmarks" --participant "$method_name" --out "participant.json" $refset_dir/mapping.json.gz $predictions
    """
}

db = Channel.create()

/*
 * extract pairwise predictions and store in darwin compatible database
 */
process convertPredictions {

    label "py"

    input:
    val file_validated from EXIT_STAT
    file predictions
    file refset_dir from refset_dir_convert

    output:
    file 'predictions.db' into db_go_test,db_ec_test,db_std,db_g_std,db_g_std_v2,db_geneTrees

    when:
    file_validated == 0

    """
    /benchmark/map_relations.py --out predictions.db $refset_dir/mapping.json.gz $predictions
    """
}



process go_benchmark {

    label "darwin"

    input:
    file db from db_go_test
    val method_name
    file refset_dir from refset_dir_go_test
    val go_evidences
    val community_id
    val result_file_path
    // for mountpoint 
    file predictions

    output:
    file "GO.json" into GO_STUB

    when:
    benchmarks =~ /GO/

    """
    /benchmark/GoTest.sh -o "${result_file_path}/GO" -a GO.json -c "$community_id" -e "$go_evidences" $db "$method_name" $refset_dir
    """
}

process ec_benchmark {

    label "darwin"

    input:
    file db from db_ec_test
    val method_name
    file refset_dir from refset_dir_ec_test
    val community_id
    val result_file_path
    // for mountpoint 
    file predictions

    output:
    file "EC.json" into EC_STUB


    when:
    benchmarks =~ /EC/


    """
    /benchmark/EcTest.sh -o "${result_file_path}/EC" -a EC.json -c "$community_id" $db "$method_name" $refset_dir
    """
}

process speciestree_benchmark {

    label "darwin"
    tag "$clade"

    input:
    file db from db_std
    val method_name
    file refset_dir from refset_dir_std
    val clade from tree_clades0
    val community_id
    val result_file_path
    // for mountpoint 
    file predictions

    output:
    file "STD_${clade}.json" into STD_STUB

    when:
    benchmarks =~ /STD_$clade/


    """
    /benchmark/SpeciesTreeDiscordanceTest.sh -o "${result_file_path}/STD_${clade}" -a "STD_${clade}.json" -c "$community_id" -p $clade -m 0 $db "$method_name" $refset_dir
    """
}


process g_speciestree_benchmark {

    label "darwin"
    tag "$clade"

    input:
    file db from db_g_std
    val method_name
    file refset_dir from refset_dir_g_std
    val clade from tree_clades
    val community_id
    val result_file_path
    // for mountpoint 
    file predictions

    output:
    file "G_STD_${clade}.json" into G_STD_STUB

    when:
    benchmarks =~ /G_STD_$clade/


    """
    /benchmark/SpeciesTreeDiscordanceTest.sh -o "${result_file_path}/G_STD_${clade}" -a "G_STD_${clade}.json" -c "$community_id" -p $clade -m 1 $db "$method_name" $refset_dir
    """
}

process g_speciestree_benchmark_variant2 {
    label "darwin"
    tag "$clade"

    input:
    file db from db_g_std_v2
    val method_name
    file refset_dir from refset_dir_g_std_v2
    val community_id
    val result_file_path
    val clade from tree_clades2
    // for mountpoint 
    file predictions

    output:
    file "G_STD2_${clade}.json" into G_STD2_STUB

    when:
    benchmarks =~ /G_STD2_$clade/


    """
    /benchmark/SpeciesTreeDiscordanceTest.sh -o "${result_file_path}/G_STD2_${clade}" -a "G_STD2_${clade}.json" -c "$community_id" -p $clade -m 2 $db "$method_name" $refset_dir
    """
}


process reference_genetrees_benchmark {
    label "darwin"
    tag "$testset"

    input:
    file db from db_geneTrees
    val method_name
    file refset_dir from refset_dir_geneTrees
    val testset from genetree_sets
    val community_id
    val result_file_path
    // for mountpoint 
    file predictions

    output:
    file "${testset}.json" into REFPHYLO_STUB

    when:
    benchmarks =~ /$testset/


    """
    /benchmark/RefPhyloTest.sh -o "$result_file_path/${testset}" -a "${testset}.json" -t "$testset" $db "$method_name" $refset_dir
    """
}


challenge_assessments = GO_STUB.mix(EC_STUB, STD_STUB, G_STD_STUB, G_STD2_STUB, REFPHYLO_STUB)

process consolidate {
    label "py"

    input:
    file participants from PARTICIPANT_STUB.collect()
    file challenge_stubs from challenge_assessments.collect()
    file benchmark_data
    val assessment_out
    val data_model_export_dir
    val result_file_path
    //for mountpoint
    file predictions

    """
    python /benchmark/manage_assessment_data.py -m $challenge_stubs -b $benchmark_data -o $result_file_path
    python /benchmark/merge_data_model_files.py -p $participants  -m $challenge_stubs -r $result_file_path -a $assessment_out -o $data_model_export_dir
    """
}



workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
