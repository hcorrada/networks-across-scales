SUBDIRS=text 
CLEAN_SUBDIRS=$(patsubst %,clean-%,$(SUBDIRS))
COURSE=CMSC423

LOCALHOST=$(HOME)/Sites/$(COURSE)
REMOTEHOST=willow.umiacs.umd.edu:/fs/www-umiacs-users/hcorrada/$(COURSE)

.PHONY: all clean subdirs $(SUBDIRS) clean-subdirs $(CLEAN_SUBDIRS)

all: pushlocal

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

#text: src

clean: clean-subdirs

clean-subdirs: $(CLEAN_SUBDIRS)

$(CLEAN_SUBDIRS):
	$(MAKE) -C $(patsubst clean-%,%,$@) clean

pushlocal: subdirs
	rsync -avz text/*.html $(LOCALHOST)
	rsync -avz ../css $(LOCALHOST)
	rsync -avz images $(LOCALHOST)
	rsync -avz text/*.pdf $(LOCALHOST)/pdf
	rsync -avz readings/ $(LOCALHOST)/readings
	rsync -avz text/lectures/ --include='*.pdf' --include='*.html' --include='*.R' --include='*.Rmd' --exclude='*.*' $(LOCALHOST)/lectures
#	rsync -avz src/*.R $(LOCALHOST)/src
#	rsync -avz ../Data $(LOCALHOST)
	rsync -avz text/homeworks/ --include='*.pdf' --include='*.html' --exclude='*.*' $(LOCALHOST)/homeworks

pushremote: subdirs
	rsync -avz text/*.html $(REMOTEHOST)
	rsync -avz ../css $(REMOTEHOST)
	rsync -avz images/ $(REMOTEHOST)/images
	rsync -avz text/*.pdf $(REMOTEHOST)/pdf/
	rsync -avz readings/ $(REMOTEHOST)/readings
	rsync -avz text/lectures/ --include='*.pdf' --include='*.html' --include='*.R' --include='*.Rmd' --include='*.pptx' --exclude='*.*' $(REMOTEHOST)/lectures
	rsync -avz ../Data $(REMOTEHOST)
#	rsync -avz src/*.R $(REMOTEHOST)/src/
	rsync -avz text/homeworks/ --include='*.pdf' --include='*.R' --include='*.html' --include='*.Rmd' --exclude='*.*' $(REMOTEHOST)/homeworks



