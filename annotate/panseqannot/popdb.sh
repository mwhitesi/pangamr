#!/bin/bash

python annot.py /media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/panseq/percentIdentityCutoff90__fragmentationSize1000__coreGenomeThreshold180/
python loader.py rgi /media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/panseq/percentIdentityCutoff90__fragmentationSize1000__coreGenomeThreshold180/ /media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/annotations/rgi/report.txt
python loader.py resfinder /media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/panseq/percentIdentityCutoff90__fragmentationSize1000__coreGenomeThreshold180/ /media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/annotations/resfinder/resfinder_blast_report.tsv
python loader.py resfams /media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/panseq/percentIdentityCutoff90__fragmentationSize1000__coreGenomeThreshold180/ /media/poolhouse/workspace/l_amr/data/NML/heidelburg_james/annotations/resfams/Resfams.scan
