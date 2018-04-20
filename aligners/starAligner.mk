include usb-modules-v2/Makefile.inc
include usb-modules-v2/aligners/align.inc

LOGDIR ?= log/star.$(NOW)

.PHONY: star
.DELETE_ON_ERROR:

ALIGNER := star
override BAM_SORT := false
override BAM_FIX_RG := true

ifeq ($(strip $(PRIMARY_ALIGNER)),star)
STAR_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
else
STAR_BAMS = $(foreach sample,$(SAMPLES),star/$(sample).star.bam)
endif

star : $(STAR_BAMS) $(addsuffix .bai,$(STAR_BAMS)) star/all$(PROJECT_PREFIX).ReadsPerGene.out.tab star/all$(PROJECT_PREFIX).alignment_stats.txt 

star/%.Aligned.sortedByCoord.out.bam : fastq/%.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/%.2.fastq.gz)
	$(call RUN,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(STAR_MODULE),"\
	STAR --runMode alignReads --runThreadN 4 \
	--genomeDir $(STAR_GENOME_DIR) \
	--readFilesIn $< $(if $(findstring true,$(PAIRED_END)),$(word 2,$^)) \
	--readFilesCommand gunzip -c \
	--sjdbScore 2 --sjdbOverhang 100 \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 \
	--outFileNamePrefix $(@D)/$*. --outFilterType BySJout \
	--outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --outFilterMultimapScoreRange 1 \
	--outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 \
	--outSAMstrandField intronMotif --outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMattrIHstart 0 \
	--chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimSegmentReadGapMax 3 \
	--quantMode GeneCounts TranscriptomeSAM --twopassMode Basic && \
	$(RMR) $(@D)/$*._STARgenome $(@D)/$*._STARpass1")

star/%.star.bam : star/%.Aligned.sortedByCoord.out.bam
	$(INIT) ln -f $< $@

star/%.Aligned.toTranscriptome.out.bam : star/%.Aligned.sortedByCoord.out.bam
	
star/%.Chimeric.out.junction : star/%.Aligned.sortedByCoord.out.bam
	
star/%.Chimeric.out.sam : star/%.Aligned.sortedByCoord.out.bam
	
star/%.ReadsPerGene.out.tab : star/%.Aligned.sortedByCoord.out.bam
	
star/%.Log.final.out : star/%.Aligned.sortedByCoord.out.bam
	
star/%.Unmapped.out.mate1 : star/%.Aligned.sortedByCoord.out.bam
	
star/%.Unmapped.out.mate2 : star/%.Aligned.sortedByCoord.out.bam
	
%.gz : %
	$(GZIP) $<

# star/%.Unmapped.out.mate1.gz $(if $(findstring true,$(PAIRED_END)),star/%.Unmapped.out.mate2.gz) star/%.Chimeric.out.bam
bam/%.bam : star/%.star.$(BAM_SUFFIX) 
	$(INIT) ln -f $< $@

star/all$(PROJECT_PREFIX).ReadsPerGene.out.tab : $(foreach sample,$(SAMPLES),star/$(sample).ReadsPerGene.out.tab)
	perl -p -e "s/N_unmapped/GENE\t\t\t\nN_unmapped/;" `ls star/secondpass/*.ReadsPerGene.out.tab|head -1` | cut -f 1 > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; perl -p -e "s/N_unmapped/\t$$sample\t\t\nN_unmapped/;" \
	star/secondpass/$$sample.ReadsPerGene.out.tab | cut -f 2 | paste $@ - > $@.tmp; mv $@.tmp $@; done

star/all$(PROJECT_PREFIX).alignment_stats.txt : $(foreach sample,$(SAMPLES),star/$(sample).Log.final.out)
	$(INIT) \
	{ \
	grep "|" $< | cut -f1 -d '|' | perl -p -e "s/ //g; s/\n/\t/g;"; echo ""; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.Log.final.out}); \
		echo -n $$samplename;\
		grep "|" $$metrics | cut -f2 -d '|' | perl -p -e "s/ //g; s/\n/\t/g;"; echo ""; \
	done; \
	} > $@



include usb-modules-v2/fastq_tools/fastq.mk
include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/aligners/align.mk
