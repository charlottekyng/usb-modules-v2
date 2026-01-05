.DEFAULT_GOAL := conSV

MUT_CALLER = conSV

include usb-modules-v2/Makefile.inc
include usb-modules-v2/vcf_tools/vcftools.mk

LOGDIR = log/conSV.$(NOW)
PHONY += conSV conSV_disruption conSV_fusion

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

SV_CALLERS ?= brass delly gridss manta svaba
SV_CALLERS := $(strip $(SV_CALLERS))
SV_CALLER_STRING := $(subst $(space),.,$(SV_CALLERS))
MIN_CALLERS ?= 2
SLOPE ?= 200

conSV : $(foreach pair,$(SAMPLE_PAIRS),conSV/conSV_$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))/$(tumor.$(pair)).$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))_slope_$(SLOPE).tsv)
conSV_disruption : $(foreach pair,$(SAMPLE_PAIRS),conSV/conSV_$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))/$(tumor.$(pair)).$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))_slope_$(SLOPE)_gene_disruption.tsv)
conSV_fusion : $(foreach pair,$(SAMPLE_PAIRS),conSV/conSV_$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))/$(tumor.$(pair)).$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))_slope_$(SLOPE)_fusion.tsv)

# create folder and organise data
#BRASS no filter implemented so taken as it is
define preprocess_brass
conSV/brass/$1.vcf : brass/$1_$2/$1_vs_$2.annot.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef

#the other callers have PASS filter implemented
define preprocess_delly
conSV/delly/$1.vcf : delly/$1.somatic.SVpass.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef

define preprocess_gridss
conSV/gridss/$1.vcf : gridss/$1.somatic.SVpass.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef

define preprocess_manta
conSV/manta/$1.vcf : manta/$1/results/variants/somaticSV.SVpass.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef

define preprocess_svaba
conSV/svaba/$1.vcf : svaba/$1.svaba.somatic.sv.SVpass.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef

# Preprocess step for SV callers to populate the fodler and rename files appropriately
define preprocess_dispatch
$(if $(filter brass,$3),$(call preprocess_brass,$1,$2))
$(if $(filter delly,$3),$(call preprocess_delly,$1,$2))
$(if $(filter gridss,$3),$(call preprocess_gridss,$1,$2))
$(if $(filter manta,$3),$(call preprocess_manta,$1,$2))
$(if $(filter svaba,$3),$(call preprocess_svaba,$1,$2))
$(if $(filter-out brass delly gridss manta svaba,$3),$(error Unknown SV caller '$3'))
endef
$(foreach pair,$(SAMPLE_PAIRS), \
  $(foreach caller,$(SV_CALLERS), \
    $(eval $(call preprocess_dispatch,$(tumor.$(pair)),$(normal.$(pair)),$(caller))) \
  ) \
)

# Preprocess step to normalize and annotate each SV caller's filtered VCF using VariantExtractor from oncoliner
define variant_preprocess_varExtr
conSV/variantExtr/$1.$2.varExtr.vcf : conSV/$2/$1.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	$$(VARIANT_EXTRACTOR) python3 $$(IMG_DIR)/oncoliner/VariantExtractor.py $$< $$@")
endef
define variant_preprocess_copy
conSV/variantExtr/$1.$2.varExtr.vcf : conSV/$2/$1.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef
define variant_preprocess_dispatch
$(if $(filter delly,$2),\
  $(call variant_preprocess_copy,$1,$2),\
  $(call variant_preprocess_varExtr,$1,$2))
endef
# Loop through each sample pair and caller to preprocess each VCF
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach caller,$(SV_CALLERS), \
		$(eval $(call variant_preprocess_dispatch,$(tumor.$(pair)),$(caller))) \
	) \
)

#making consensus SVs
define conSV
conSV/conSV_$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))/$1.$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))_slope_$(SLOPE).tsv : $(foreach caller,$(SV_CALLERS), conSV/variantExtr/$1.$(caller).varExtr.vcf)
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(CONSV) --min_callers $$(MIN_CALLERS) --slope $$(SLOPE) --sv_callers $$(SV_CALLER_STRING) --output $$@ --input $1")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call conSV,$(tumor.$(pair)))) \
)


#making consensus SVs disruption annotation: 🚧 wip (implemented only for hg38 and for oncoKB anno file)
#also promoter rigidly defined as region 200bp upstream and 50bp downstream genes TSS
define conSV_disruption
conSV/conSV_$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))/$1.$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))_slope_$(SLOPE)_gene_disruption.tsv : conSV/conSV_$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))/$1.$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))_slope_$(SLOPE).tsv
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(CONSV_disruption) --oncokb_file $$(ONCOKB_FILE) --oncokb_minresources $$(ONCOKB_MINRESOURCES) --output $$@ --input $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call conSV_disruption,$(tumor.$(pair)))) \
)

#making consensus SVs fusion gene annotation: 🚧 wip (implemented only for hg38 and naive way to consider fusion)
# genes to consider for fusion are defined in the conSV_fusion.R script in future could be parameterized
define conSV_fusion
conSV/conSV_$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))/$1.$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))_slope_$(SLOPE)_fusion.tsv : conSV/conSV_$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))/$1.$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))_slope_$(SLOPE).tsv
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(CONSV_fusion) --output $$@ --input $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call conSV_fusion,$(tumor.$(pair)))) \
)