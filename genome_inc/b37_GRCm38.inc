ifndef B37_GRCM38_INC
include usb-modules-v2/genome_inc/b37.inc
REF = b37_GRCm38
REF_DIR = /scicore/home/terracci/GROUP/ref_nobackup

ANN_DIR := $(REF_DIR)/annotation_data/b37_GRCm38
REF_FASTA := $(REF_DIR)/genomes/human_g1k_v37_decoy_GRCm38.fasta
REF_DICT = $(REF_FASTA:.fasta=.dict)

STAR_GENOME_DIR := $(REF_DIR)/genome_index/star_b37_GRCm38

RSEM_INDEX = $(REF_DIR)/genome_index/rsem_b37_GRCm38/rsem_b37_GRCm38

GENCODE_GTF = $(ANN_DIR)/gencode_b37_GRCm38.gtf
GENCODE_GENE_GTF = $(ANN_DIR)/gencode_b37_GRCm38_genes.gtf

endif
B37_GRCM38_INC = true

