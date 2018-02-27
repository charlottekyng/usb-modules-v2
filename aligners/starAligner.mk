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

star : $(STAR_BAMS) $(addsuffix .bai,$(STAR_BAMS)) star/all.ReadsPerGene.out.tab

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
	

bam/%.bam : star/secondpass/%.star.$(BAM_SUFFIX)
	$(INIT) ln -f $< $@


#define align-star-secondpass
#star/secondpass/$1.star.bam : fastq/$1.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/$1.2.fastq.gz) $(foreach sample,$(SAMPLES),star/firstpass/$(sample).SJ.out.tab)
#	$$(call LSCRIPT_PARALLEL_MEM,8,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_LONG),"$$(MKDIR) star star/secondpass/; \
#	$$(LOAD_STAR_MODULE); STAR --runMode alignReads \
#	--runThreadN 8 --genomeDir $$(STAR_GENOME_DIR) --readFilesIn $$< $$(if $$(findstring true,$(PAIRED_END)),$$(word 2,$$^)) \
#	--readFilesCommand gunzip -c \
#	--alignSJoverhangMin 8 --alignSJDBoverhangMin 10 --alignIntronMax 200000 --alignMatesGapMax 200000 \
#	--alignSJstitchMismatchNmax 5 -1 5 5 --limitSjdbInsertNsj 5000000 \
#	--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
#	--outFileNamePrefix star/secondpass/$$(basename $$(basename $$(notdir $$@))). \
#	--sjdbFileChrStartEnd $$(filter %.SJ.out.tab,$$^) \
#	--outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate \
#	--outReadsUnmapped None --outMultimapperOrder Random --outSAMattrIHstart 0 \
#	--chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimSegmentReadGapMax parameter 3 \
#	--quantMode GeneCounts TranscriptomeSAM && mv star/secondpass/$1.Aligned.sortedByCoord.out.bam star/secondpass/$1.star.bam")

#star/secondpass/$1.Chimeric.out.junction : star/secondpass/$1.star.bam
#star/secondpass/$1.ReadsPerGene.out.tab : star/secondpass/$1.star.bam
#star/secondpass/$1.Aligned.toTranscriptome.out.bam : star/secondpass/$1.star.bam
#endef
#$(foreach sample,$(SAMPLES),\
#	$(eval $(call align-star-secondpass,$(sample))))

#bam/%.bam : star/secondpass/%.star.bam
#	$(INIT) ln -f $< $@

star/all.ReadsPerGene.out.tab : $(foreach sample,$(SAMPLES),star/secondpass/$(sample).ReadsPerGene.out.tab)
	perl -p -e "s/N_unmapped/GENE\t\t\t\nN_unmapped/;" `ls star/secondpass/*.ReadsPerGene.out.tab|head -1` | cut -f 1 > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; perl -p -e "s/N_unmapped/\t$$sample\t\t\nN_unmapped/;" \
	star/secondpass/$$sample.ReadsPerGene.out.tab | cut -f 2 | paste $@ - > $@.tmp; mv $@.tmp $@; done

include usb-modules-v2/fastq_tools/fastq.mk
include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/aligners/align.mk
