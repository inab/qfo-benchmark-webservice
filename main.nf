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
refset_dir = file(params.goldstandard_dir, type: 'dir')
benchmarks = params.challenges_ids
benchmarks_arr = params.challenges_ids.split(/ +/)
community_id = params.community_id
benchmark_data = file(params.assess_dir, type: 'dir')
go_evidences = params.go_evidences
tree_clades = ["Luca", "Vertebrata", "Fungi", "Eukaryota"]
tree_clades2 = ["Luca", "Vertebrata", "Fungi", "Eukaryota"]
genetree_sets = ["SwissTrees", "TreeFam-A"]
tree_clades0 = ["Eukaryota", "Fungi", "Bacteria"]

//output
//validation_out = file(params.validation_result)
assessment_out = file(params.assessment_results)
result_file_path = file(params.outdir, type: 'dir')
data_model_export_dir = file(params.data_model_export_dir, type: 'dir')
otherdir = file(params.otherdir, type: 'dir')



/*
 * validate input file
 */
process validate_input_file {
    label "py"

    input:
    file predictions
    file refset_dir
    val benchmarks
    val community_id
    val method_name

    output:
    val task.exitStatus into EXIT_STAT
    file "participant.json" into PARTICIPANT_STUB

    """
    /benchmark/validate.py --com $community_id --challenges_ids "$benchmarks" --participant "$method_name" --out "participant.json" $refset_dir/mapping.json.gz $predictions
    """
}

// These channels rule next steps
db_go_test = Channel.create()
db_ec_test = Channel.create()
db_std = Channel.create()
db_g_std = Channel.create()
db_g_std_v2 = Channel.create()
db_geneTrees = Channel.create()

/*
 * extract pairwise predictions and store in darwin compatible database
 */
process convertPredictions {

    label "py"

    input:
    val file_validated from EXIT_STAT
    file predictions
    file refset_dir

    output:
    file 'predictions.db' into predictions_db

    when:
    file_validated == 0

    """
    /benchmark/map_relations.py --out predictions.db $refset_dir/mapping.json.gz $predictions
    """
}

process scheduleMetrics {
    
    input:
    file predictions_db
    val benchmark from benchmarks_arr
    
    // Setting up the cascade of events
    exec:
    
    def m
    switch(benchmark) {
        case "GO":
            predictions_db.into { db_go_test; predictions_db }
            break
        case "EC":
            predictions_db.into { db_ec_test; predictions_db }
            break
        case ~/^STD_/:
            m = benchmark =~ /^STD_(.+)$/
            def clade0 = m[0][1]
            if(tree_clades0.contains(clade0)) {
                predictions_db.into { pred_dup; predictions_db }
                pred_dup.map {pred_db -> tuple(clade0,pred_db)}.set { db_std }
            }
            break
        case ~/^G_STD_/:
            m = benchmark =~ /^G_STD_(.+)$/
            def clade = m[0][1]
            if(tree_clades.contains(clade)) {
                predictions_db.into { pred_dup; predictions_db }
                pred_dup.map {pred_db -> tuple(clade,pred_db)}.set { db_g_std }
            }
            break
        case ~/^G_STD2_/:
            m = benchmark =~ /^G_STD2_(.+)$/
            def clade2 = m[0][1]
            if(tree_clades2.contains(clade2)) {
                predictions_db.into { pred_dup; predictions_db }
                pred_dup.map {pred_db -> tuple(clade2,pred_db)}.set { db_g_std_v2 }
            }
            break
        default:
            if(genetree_sets.contains(benchmark)) {
                predictions_db.into { pred_dup; predictions_db }
                pred_dup.map {pred_db -> tuple(benchmark,pred_db)}.set { db_geneTrees }
            }
            break
    }
}



process go_benchmark {

    label "darwin"

    input:
    file db from db_go_test
    val method_name
    file refset_dir
    val go_evidences
    val benchmarks_arr
    val community_id
    file result_file_path
    // for mountpoint 
    file predictions
    
    output:
    file "GO.json" into GO_STUB
    
    """
    /benchmark/GoTest.sh -o "${result_file_path}/GO" -a GO.json -c "$community_id" -e "$go_evidences" $db "$method_name" $refset_dir
    """
}

process ec_benchmark {

    label "darwin"

    input:
    file db from db_ec_test
    val method_name
    file refset_dir
    val benchmarks_arr
    val community_id
    file result_file_path
    // for mountpoint 
    file predictions

    output:
    file "EC.json" into EC_STUB
    
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
    file refset_dir
    val clade from tree_clades0
    val benchmarks_arr
    val community_id
    file result_file_path
    // for mountpoint 
    file predictions

    output:
    file "STD_${clade}.json" into STD_STUB

    when:
    benchmarks_arr.contains("STD_$clade")


    """
    /benchmark/SpeciesTreeDiscordanceTest.sh -o "${result_file_path}/STD_${clade}" -a "STD_${clade}.json" -c "$community_id" -p $clade -m 0 $db "$method_name" $refset_dir
    """
}


process g_speciestree_benchmark {

    label "darwin"
    tag "$clade"

    input:
    set clade, file(db) from db_g_std
    val method_name
    file refset_dir
    val benchmarks_arr
    val community_id
    file result_file_path
    // for mountpoint 
    file predictions

    output:
    file "G_STD_${clade}.json" into G_STD_STUB

    """
    /benchmark/SpeciesTreeDiscordanceTest.sh -o "${result_file_path}/G_STD_${clade}" -a "G_STD_${clade}.json" -c "$community_id" -p $clade -m 1 $db "$method_name" $refset_dir
    """
}

process g_speciestree_benchmark_variant2 {
    label "darwin"
    tag "$clade"

    input:
    set clade, file(db) from db_g_std_v2
    val method_name
    file refset_dir
    val benchmarks_arr
    val community_id
    file result_file_path
    // for mountpoint 
    file predictions

    output:
    file "G_STD2_${clade}.json" into G_STD2_STUB

    """
    /benchmark/SpeciesTreeDiscordanceTest.sh -o "${result_file_path}/G_STD2_${clade}" -a "G_STD2_${clade}.json" -c "$community_id" -p $clade -m 2 $db "$method_name" $refset_dir
    """
}


process reference_genetrees_benchmark {
    label "darwin"
    tag "$testset"

    input:
    set testset, file(db) from db_geneTrees
    val method_name
    file refset_dir
    val benchmarks_arr
    val community_id
    file result_file_path
    // for mountpoint 
    file predictions

    output:
    file "${testset}.json" into REFPHYLO_STUB

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
    file assessment_out
    file data_model_export_dir
    file result_file_path
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
