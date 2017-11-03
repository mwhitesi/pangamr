#!/usr/bin/env python

"""Config Settings for Pangenome AMR analysis

Defines file locations of critical datasets


"""

PANSEQ = {
	'pangenome_file': '/media/poolhouse/workspace/l_amr/data/PATRIC/salmonella/panseq/pan_genome.txt'
}

PHENOTYPE = {
	'amr_file': '/home/matt/workspace/l_amr/patric_tools/data/salmonella/amr.csv'
}

ANNOTATION = {
	'blast_file': '/media/poolhouse/workspace/l_amr/data/PATRIC/salmonella/blast/parsed_pangenome_blast_report.txt'
}

SH = {
	'amr': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/amr.jpkl',
	'sample_index': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/sample_index.jpkl',
	'amr_list': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/amr_list.jpkl',
	'locus_list': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/locus_list.jpkl',
	'pg': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/pg.jpkl',
	'tree': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/SVM_SH/SVM_SH_Pasnp/parsnp.tree',
	'test_train_index': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/test_train_index.jpkl',
	'amp_rfc': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/amp_rfc.jpkl',
	'amp_gbc': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/amp_gbc.jpkl',
	'amp_xbc': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/amp_xbc.jpkl',
	'fox_rfc': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/fox_rfc.jpkl',
	'fox_gbc': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/fox_gbc.jpkl',
	'fox_xbc': '/media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/ml/fox_xbc.jpkl',
}