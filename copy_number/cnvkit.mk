include usb-modules-v2/Makefile.inc

LOGDIR ?= log/cnvkit.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : cnvkit

cnvkit : cnvkit/cns/all$(PROJECT_PREFIX).heatmap.pdf cnvkit/cns/all$(PROJECT_PREFIX).geneCN.GL_LRR.pdf \
cnvkit/cns/all$(PROJECT_PREFIX).geneCN.log2.pdf cnvkit/cns/all$(PROJECT_PREFIX).cns.pdf.tar.gz

cnvkit/rsem/%.genes.results : rsem/%.genes.results
	$(MKDIR) cnvkit cnvkit/rsem; \
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(CNVKIT_PROCESS_RSEM) --inputRSEMFile $< --gtf $(GENCODE_CODING_GENE_GTF) \
	--geneBiotype protein_coding --outputFile $@")

cnvkit/all$(PROJECT_PREFIX).summary.out : $(foreach sample,$(SAMPLES),cnvkit/rsem/$(sample).genes.results) $(foreach normal,$(PANEL_OF_NORMAL_SAMPLES),cnvkit/rsem/$(normal).genes.results)
	$(MKDIR) cnvkit cnvkit/cnr; \
	ln $(CNVKIT_GENE_RESOURCE) .; ln $(CNVKIT_CORRELATION) .; \
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SINGULARITY_MODULE),"\
	$(SINGULARITY_EXEC) $(CNVKIT_IMG) $(CNVKIT) import-rna -f rsem \
	-n $(foreach normal,$(PANEL_OF_NORMAL_SAMPLES),./cnvkit/rsem/$(normal).genes.results) \
	-g ./$(notdir $(CNVKIT_GENE_RESOURCE)) \
	-c ./$(notdir $(CNVKIT_CORRELATION)) \
	-o ./$@ -d ./cnvkit/cnr \
	$(foreach sample,$(SAMPLES),./cnvkit/rsem/$(sample).genes.results) && \
	$(RM) ./$(notdir $(CNVKIT_GENE_RESOURCE)) ./$(notdir $(CNVKIT_CORRELATION))")

cnvkit/cnr/%.cnr : cnvkit/all$(PROJECT_PREFIX).summary.out
	
	
cnvkit/cns/%.cns : cnvkit/cnr/%.cnr
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SINGULARITY_MODULE),"\
	$(SINGULARITY_EXEC) $(CNVKIT_IMG) $(CNVKIT) segment \
	./$< -o ./$@ --drop-low-coverage")
	
cnvkit/cns/%.cns.pdf : cnvkit/cns/%.cns cnvkit/cnr/%.cnr
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SINGULARITY_MODULE),"\
	$(SINGULARITY_EXEC) $(CNVKIT_IMG) $(CNVKIT) scatter \
	-s ./$< -o ./$@ $(<<)")
	
cnvkit/cns/all$(PROJECT_PREFIX).heatmap.pdf : $(foreach sample,$(SAMPLES),cnvkit/cns/$(sample).cns)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SINGULARITY_MODULE),"\
	$(SINGULARITY_EXEC) $(CNVKIT_IMG) $(CNVKIT) heatmap -d \
	$(foreach sample,$(SAMPLES),./cnvkit/cns/$(sample).cns) -o ./$@")
	
cnvkit/cns/all$(PROJECT_PREFIX).cns.pdf.tar.gz : $(foreach sample,$(SAMPLES),cnvkit/cns/$(sample).cns.pdf)
	$(INIT) tar -czf $@ $^

cnvkit/cns/all$(PROJECT_PREFIX).geneCN.%.pdf : cnvkit/cns/all$(PROJECT_PREFIX).geneCN.%.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(CNVKIT_GENE_CN_PLOT) $(CNVKIT_GENE_CN_PLOT_OPTS) $< $@")
	
cnvkit/cns/all$(PROJECT_PREFIX).geneCN.GL_LRR.txt : $(foreach sample,$(SAMPLES),cnvkit/cns/$(sample).cns)
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_SHORT),$(R_MODULE),"\
	$(CNVKIT_GENE_CN) $(CNVKIT_GENE_CN_OPTS) \
	--outFile $(@D)/all$(PROJECT_PREFIX).geneCN $^")

cnvkit/cns/all$(PROJECT_PREFIX).geneCN.log2.txt : cnvkit/cns/all$(PROJECT_PREFIX).geneCN.GL_LRR.txt
	
	



#include usb-modules-v2/rnaseq/rsem.mk
