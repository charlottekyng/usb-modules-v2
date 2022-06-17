include usb-modules-v2/Makefile.inc

LOGDIR ?= log/tcell_extrect.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : tcell_extrect


define extrect
tcell_extrect : tcell_extrect/$1_TCRA.txt.gz

tcell_extrect/$1_TCRA.txt.gz : tcell_extrect/$1.resTcellExTRECT.txt
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_SHORT),,"\
	$$(GZIP) tcell_extrect/$1_TCRA.txt")

tcell_extrect/%.resTcellExTRECT.txt : bam/$1.bam facets/cncf/$1_$2.out facets/cncf/$1_$2.cncf.txt
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(R4_MODULE) $$(SAMTOOLS_MODULE),"\
	$$(TCELLEXTRECT) --sample_id $1 \
	--TCRA_exons $$(TCRA_EXONS) \
	--purity `grep Purity $$(filter %.out,$$^) | cut -f2 -d'=' | tr -d ' ' | sed 's/NA/0.1/'` \
	--cncf $$(filter %.cncf.txt,$$^) \
	--outdir tcell_extrect \
	$$(filter %.bam,$$^)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call extrect,$(tumor.$(pair)),$(normal.$(pair)))))
