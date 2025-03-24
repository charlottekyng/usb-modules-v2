# run brass
MUT_CALLER = brass

include usb-modules-v2/Makefile.inc

LOGDIR = log/brass.$(NOW)
PHONY += brass

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

HIGH_DEPTH_BED = $(HOME)/share/reference/extremedepth.bed
BRASS_REPEATS = $(HOME)/share/reference/brassRepeats.bed.gz
GENOME_CACHE = $(HOME)/share/reference/Homo_sapiens.GRCh37.74.vagrent.cache.gz
BRASS_PROTOCOL = WGS # WGS|WXS|RNA
BRASS_OPTS = -g $(REF_FASTA) -s HUMAN -as 37 -pr $(BRASS_PROTOCOL) -gc $(GENOME_CACHE) -d $(HIGH_DEPTH_BED) -r $(BRASS_REPEATS) 

brass : $(foreach pair,$(SAMPLE_PAIRS),brass/$(pair).brass_timestamp)

bam/%.bam.bas : bam/%.bam
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(BRASS) bam_stats \
	-i $< \
	-o $@ \
	-r $(REF_FASTA).fai")

brass/genome_window.bed :
	$(call LSCRIPT_MEM,2G,3G,"$(BEDTOOLS) makewindows -g $(REF_FASTA) -w 10000 > $@")
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(BRASS) bedtools makewindows \
	-g $(REF_FASTA) \
	-w 10000 > $@")

define brass-tumor-normal
# Main call
brass/$1_$2.brass_input_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(BRASS) brass.pl \
	-x $$(if $$(findstring hg38,$$(REF)),usb-modules-v2/resources/human.hg38.excl.bed,\
	$$(if $$(findstring b37,$$(REF)),usb-modules-v2/resources/human.hg19.excl.bed,\
	$$(if $$(findstring GRCm38,$$(REF)),usb-modules-v2/resources/mouse.mm10.excl.bed,))) \
	-g $$(REF_FASTA) \
	-o $$@ \
	$$^")

brass/$1_$2.brass_input_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas
	$$(call LSCRIPT_PARALLEL_MEM,2,6G,10G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -c 2 -p input -o brass/$1_$2 && touch $$@")

brass/$1_$2.brass_group_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_input_timestamp
	$$(call LSCRIPT_MEM,7G,9G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -o brass/$1_$2 -p group && touch $$@")

brass/$1_$2.brass_filter_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_group_timestamp
	$$(call LSCRIPT_MEM,7G,9G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -o brass/$1_$2 -p filter && touch $$@")

brass/$1_$2.brass_assemble_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_filter_timestamp
	$$(call LSCRIPT_PARALLEL_MEM,8,2G,3G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -o brass/$1_$2 -c 8 -p assemble && touch $$@")

brass/$1_$2.brass_grass_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_assemble_timestamp
	$$(call LSCRIPT_MEM,7G,9G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -o brass/$1_$2 -p grass && touch $$@")

brass/$1_$2.brass_tabix_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_grass_timestamp
	$$(call LSCRIPT_MEM,2G,3G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -o brass/$1_$2 -p tabix && touch $$@")

brass/$1_$2.brass_timestamp : brass/$1_$2.brass_tabix_timestamp
	$$(INIT) touch $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call brass-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

.PHONY: $(PHONY)