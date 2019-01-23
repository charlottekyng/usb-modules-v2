include usb-modules-v2/Makefile.inc

GISTIC_OPTS = -genegistic 0 -smallmem 1 -maxseg 5000 -savegene 1 -saveseg 1 -savedata 0 -v 30
GISTIC_OPTS += -qvt 0.25 -conf 0.99 -broad 1 -brlen 0.5 -rx 0
GISTIC_OPTS += -ta $(GISTIC_THRESHOLD) -td $(GISTIC_THRESHOLD) -js $(GISTIC_JS)

LOGDIR = log/gisticFacets.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

MEM := 2G
PE := 1

CNV_SIZE =  300000

#all : gistic/segmentationfile.txt
all : $(foreach size,$(CNV_SIZE),gistic/$(PROJECT_PREFIX)gistic_cnv$(size)/timestamp)

gistic/$(PROJECT_PREFIX)markersfile.txt : gistic/$(PROJECT_PREFIX)segmentationfile.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(GISTIC_MAKE_MARKERS_FILE) --outFile $@ --targetsFile $(TARGETS_FILE_INTERVALS) $<")
	
gistic/$(PROJECT_PREFIX)segmentationfile.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	$(call RUN,8,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(GISTIC_MAKE_SEG_FILE) --outFile $@ --targetsFile $(TARGETS_FILE_INTERVALS) $^")

gistic/$(PROJECT_PREFIX)cnv.$(CNV_SIZE).txt : gistic/$(PROJECT_PREFIX)markersfile.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(GISTIC_MAKE_CNV_FILE) --outFile $@ --dgvFile $(DGV_FILE) --cnvSize $(CNV_SIZE) $<")

gistic/$(PROJECT_PREFIX)gistic_cnv%/timestamp : gistic/$(PROJECT_PREFIX)segmentationfile.txt gistic/$(PROJECT_PREFIX)markersfile.txt gistic/$(PROJECT_PREFIX)cnv.%.txt
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_VSHORT),,"\
	export MCR_DIR=/scicore/home/pissal00/GROUP/usr_nobackup/local/MCR_R2014a/; \
	export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH); \
	umask 002; $(MKDIR) $(@D); $(GISTIC) -b $(@D) -seg $< -mk $(<<) -refgene $(GISTIC_REF) \
	-cnv $(<<<) $(GISTIC_OPTS) 2>&1 && touch $@")
