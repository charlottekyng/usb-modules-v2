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
STAR_BAMS = $(foreach sample,$(SAMPLES),star/secondpass/$(sample).star.bam)
endif

star : $(STAR_BAMS) $(addsuffix .bai,$(STAR_BAMS)) star/all$(PROJECT_PREFIX).ReadsPerGene.out.tab star/all$(PROJECT_PREFIX).alignment_stats.txt 

star/firstpass/%.SJ.out.tab : fastq/%.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/%.2.fastq.gz)
	$(call RUN,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(STAR_MODULE),"\
	$(MKDIR) star star/firstpass/; \
	STAR --runMode alignReads \
	--runThreadN 4 --genomeDir $(STAR_GENOME_DIR) --readFilesIn $< $(if $(findstring true,$(PAIRED_END)),$(word 2,$^)) \
	--readFilesCommand gunzip -c \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
	--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
	--outFileNamePrefix $(@D)/$*. \
	--outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMattrIHstart 0")

star/secondpass/%.Aligned.sortedByCoord.out.bam : fastq/%.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/%.2.fastq.gz) star/firstpass/%.SJ.out.tab
	$(call RUN,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(STAR_MODULE),"\
	$(MKDIR) star star/secondpass/; \
	STAR --runMode alignReads \
	--runThreadN 4 --genomeDir $(STAR_GENOME_DIR) --readFilesIn $< $(if $(findstring true,$(PAIRED_END)),$(word 2,$^)) \
	--readFilesCommand gunzip -c \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 10 --alignIntronMax 200000 --alignMatesGapMax 200000 \
	--alignSJstitchMismatchNmax 5 -1 5 5 --limitSjdbInsertNsj 5000000 \
	--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
	--outFileNamePrefix star/secondpass/$*. \
	--sjdbFileChrStartEnd $(filter %.SJ.out.tab,$^) \
	--outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped None --outMultimapperOrder Random --outSAMattrIHstart 0 \
	--chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimSegmentReadGapMax parameter 3 \
	--quantMode GeneCounts TranscriptomeSAM && \
	$(RMR) $(@D)/$*._STARgenome")

#	mv $(@D)/$*.Aligned.sortedByCoord.out.bam $(@D)/$*.star.bam \

star/secondpass/%.star.bam : star/secondpass/%.Aligned.sortedByCoord.out.bam
	$(INIT) ln -f $< $@

star/secondpass/%.Chimeric.out.junction : star/secondpass/%.star.bam
	
star/secondpass/%.ReadsPerGene.out.tab : star/secondpass/%.star.bam
	
star/secondpass/%.Aligned.toTranscriptome.out.bam : star/secondpass/%.star.bam

star/secondpass/%.Log.final.out : star/secondpass/%.star.bam
	

bam/%.bam : star/secondpass/%.star.$(BAM_SUFFIX)
	$(INIT) ln -f $< $@


star/all$(PROJECT_PREFIX).ReadsPerGene.out.tab : $(foreach sample,$(SAMPLES),star/secondpass/$(sample).ReadsPerGene.out.tab)
	perl -p -e "s/N_unmapped/GENE\t\t\t\nN_unmapped/;" `ls star/secondpass/*.ReadsPerGene.out.tab|head -1` | cut -f 1 > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; perl -p -e "s/N_unmapped/\t$$sample\t\t\nN_unmapped/;" \
	star/secondpass/$$sample.ReadsPerGene.out.tab | cut -f 2 | paste $@ - > $@.tmp; mv $@.tmp $@; done

star/all$(PROJECT_PREFIX).alignment_stats.txt : $(foreach sample,$(SAMPLES),star/secondpass/$(sample).Log.final.out)
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
