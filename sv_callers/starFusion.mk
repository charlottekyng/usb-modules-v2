include usb-modules-v2/Makefile.inc

LOGDIR ?= log/star_fusion.$(NOW)

PHONY += star_fusion
.PHONY : $(PHONY)
.DELETE_ON_ERROR:

star_fusion : star_fusion/all.star-fusion.STAR-Fusion.final
define star-fusion
star_fusion/$1/star-fusion.STAR-Fusion.filter.ok : star/$1.Chimeric.out.junction
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(STAR_FUSION_MODULE),"\
	$$(STAR_FUSION) --genome_lib_dir $$(STAR_CTAT_DIR) \
	-J $$< --output_dir star_fusion/$1 --max_promiscuity $$(STAR_FUSION_MAX_PROMISCUITY) \
	--min_novel_junction_support $$(STAR_FUSION_MIN_JUNC_SUPP) --min_alt_pct_junction $$(STAR_FUSION_MIN_ALT_PCT_JUNC)")

star_fusion/$1/star-fusion.fusion_candidates.final : star_fusion/$1/star-fusion.STAR-Fusion.filter.ok
endef
$(foreach sample,$(SAMPLES),$(eval $(call star-fusion,$(sample))))


star_fusion/all.star-fusion.STAR-Fusion.final : $(foreach sample,$(SAMPLES),star_fusion/$(sample)/star-fusion.fusion_candidates.final)
	$(INIT) \
	{ \
	sed "s/^#fusion/SAMPLE\tfusion/" $< | head -1; \
	for fusion in $^; do \
		samplename=`dirname $${fusion}`; \
		samplename=`basename $${samplename}`; \
		sed "/^#/d; s/^/$$samplename\t/" $$fusion; \
	done; \
	} > $@

include usb-modules-v2/aligners/starAligner.mk
