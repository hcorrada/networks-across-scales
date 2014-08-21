SOURCES=projectI.text projectII.text projectIII.text projectIV.text classificationNotes.text
HTMLOUT=$(patsubst %.text,%.html,$(SOURCES))
PDFOUT=$(patsubst %.text,%.pdf,$(SOURCES))

.PHONY: all clean

#all: $(HTMLOUT) $(PDFOUT)
all: $(HTMLOUT)

clean:
	rm -f *.html *.pdf

%.html: %.text
	pandoc -s --mathjax --toc -c ../css/lecture_style.css -o $@ $<

%.pdf: %.text 
	cd .tmp && \
	pandoc -s -w latex -H ../../latex_include.tex -o tmp.tex ../$< && \
	sed -e 's/\(.*includegraphics\)\({.*}\)/\\begin{figure}[htp]\\centering \1[width=4in]\2\\end{figure}/g' tmp.tex > $(basename $@).tex && \
	pdflatex $(basename $@).tex && \
	mv $@ ../$@ 





