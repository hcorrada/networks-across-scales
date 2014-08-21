SOURCES=$(wildcard *.text)

HTMLOUT=$(patsubst %.text,%.html,$(SOURCES))
PDFOUT=$(patsubst %.text,%.pdf,$(SOURCES))

SUBDIRS=hidden_text homeworks lectures/classification

CLEAN_SUBDIRS=$(patsubst %,clean-%,$(SUBDIRS))

.PHONY: all clean subdirs clean-subdirs $(SUBDIRS) $(CLEAN_SUBDIRS)

all: subdirs $(HTMLOUT) $(PDFOUT)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean: clean-subdirs
	rm -f *.html *.pdf

clean-subdirs: $(CLEAN_SUBDIRS)

$(CLEAN_SUBDIRS):
	make -C $(patsubst clean-%,%,$@) clean

%.html: %.text
	pandoc -s --toc -c css/style.css -o $@ $<

%.pdf: %.text
	markdown2pdf -H latex_include.tex -o $@ $<
