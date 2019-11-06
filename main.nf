#!/usr/bin/env nextflow

dataFolder = "PF16657"
leilaDir = "/users/cn/lmansouri/PROJECTS/Phylo3D/NEW_Phylo3D/TMalign/other_distances/d1_unweighted/"
params.phylip = "${leilaDir}/${dataFolder}/*_tmalign.ph"
// params.phylip = "$baseDir/data/${familyFolder}/*_tmalign.ph"

//params for <01>
params.method ="BioNJ"
params.modelProtein = "LG"
params.gammaRate = "1.0"
params.seedValue = "5"
params.replicates = "100"
params.fastme = true

//params for <02>
params.msa = "tmalign"
params.mode = "4"   //1-6
params.exp = "2"    //(2=square; 3=cubic)
params.weight = "1"

params.fasta = "${leilaDir}/${dataFolder}/*_${params.msa}.fa"
params.templates = "${leilaDir}/${dataFolder}/*_ref.template_list2"
params.pdb = "${leilaDir}/${dataFolder}/*.pdb"
params.replicatesNum = "100"
params.evaluate3DVal = "distances"

params.output = "$baseDir/${dataFolder}/results/"

// Channels containing ph files <01>
if ( params.phylip ) {
  Channel
  .fromPath(params.phylip)
  .map { item -> [ item.baseName, item] }
  .view()
  .set { phySeqs }
}
// Channels containing fasta files <02>
if ( params.fasta ) {
  Channel
  .fromPath(params.fasta)
  .map { item -> [ item.baseName, item] }
  .view()
  .set { fastaSeqs }
}
// Channels containing tempaltes files <02>
if ( params.templates ) {
  Channel
  .fromPath(params.templates)
  .map { item -> [ item.baseName, item] }
  .view()
  .set { templates }
}
// Channels containing PDB files <02>
if ( params.pdb ) {
  Channel
  .fromPath(params.pdb)
  .collect()
  .view()
  .set { pdbFiles }
}

process run_fastme_tmalign {
    tag "${id}"
    publishDir "${params.output}", mode: 'copy', overwrite: true        //TODO diff folder for diff outChannel

    input:
      set val(id), file(seqs) from phySeqs

    output:
     file("*.replicates") into replicatesOut
     file("*.nwk") into newickOut

    when:
      params.fastme

    script:
      """
      ${baseDir}/bin/fastme -i ${seqs} -o ${seqs}.nwk -m ${params.method} -p ${params.modelProtein} -g ${params.gammaRate} -s -n -z ${params.seedValue} -b ${params.replicates} -B ${seqs}.replicates
      """
}

process phylo3d_unweighted_d1_ratio {
    tag "${id}"
    publishDir "${params.output}", mode: 'copy', overwrite: true        //TODO diff folder for diff outChannel

    input:
      set val(id), file(fasta) from fastaSeqs
      set val(id2), file(template) from templates
      file(pdb) from pdbFiles

    output:
      set val(id), file("*.trees") into treesOut
      file("*.matrices") into matrixOut

    //TODO <<script>> -output as variable
    script:
    """
    export THREED_TREE_MODE=${params.mode}
    export THREED_TREE_NO_WEIGHTS=${params.weight}
    export THREED_TREE_MODE_EXP=${params.exp}

    t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +tree replicates ${params.replicatesNum} +evaluate3D ${params.evaluate3DVal} +tree2bs first +print_replicates -output newick > ${id}_unw_d1_ratio_${params.msa}_${params.mode}-${params.exp}.trees
    t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +tree replicates ${params.replicatesNum} +evaluate3D ${params.evaluate3DVal} +tree2bs first +print_replicates -output dm > ${id}_unw_d1_ratio_${params.msa}_${params.mode}-${params.exp}.matrices
    """
}
//process extr_fastme_per_family_D1_unweighted{
//    tag "${id}"
//    publishDir "${params.output}", mode: 'copy', overwrite: true

//    input:

//    output:

//    script:
//    """
    //split matrixOut

    //run fastme for each matrix ?? parallel??

    //mode number 1 -> not sure if it makes sense
//    """
//}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}" 
}

