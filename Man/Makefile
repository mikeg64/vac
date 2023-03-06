# This Makefile can produce a single Postscript file MANUAL.ps.gz 
# for PRINTING the manual. Users do not need to remake this file
# unless there is a need to change something relative to the distributed
# MANUAL.ps.gz file. To print out the manual type
#
#    make print
#
# or
#
#    gunzip MANUAL.ps.gz ; lpr -s MANUAL.ps 
#
# To reproduce MANUAL.ps.gz, simply type 
#
#    make
#
# in the Man directory. First you will get the "ALL.html" file
# which is then translated to ALL.tex by the Perl script "html2latex.pl" 
# (requires PERL5), then Latex2e produces ALL.dvi (it runs three times 
# in order to get the table of contents right), dvips produces MANUAL.ps, 
# and finally gzip compresses it to MANUAL.ps.gz
#
# These steps may take a while. There are many warnings from Latex2e 
# about underfull hboxes, but these can be ignored.
#
# You can edit this Makefile (e.g. to omit some Manual pages, 
# to change the name of "dvips" or other non-standard codes)
# and the "MANhead.tex" file to change the style (e.g. onesided).
#
# You are NOT allowed to modify the content of the manuals or to
# omit the title page showing the title, author, and date.
#
# After producing a satisfactory MANUAL.ps.gz file type
#
#    make clean
#
# to delete the intermediate files.

HTML= WELCOME.html Features.html USAGE.html \
      mpi.html Weblink.html Examples.html equations.html vacusr.html \
      discretization.html axial.html inipar.html par.html methods.html \
      idl.html matlab.html sm.html gnuplot.html java.html \
      convert.html fileformat.html \
      source.html srclimits.html vacpp.html varnames.html \
      bugs.html Diary.html Credits.html REGISTER.html 

manual: MANUAL.ps.gz

ALL.html: $(HTML) Makefile
	/bin/rm -f ALL.html
	cat $(HTML) > ALL.html

ALL.tex: ALL.html ./html2latex.pl
	./html2latex.pl ALL.html > ALL.tex

ALL.dvi: ALL.tex MANhead.tex
	latex ALL.tex
	latex ALL.tex
	latex ALL.tex

MANUAL.ps.gz:	ALL.dvi
	dvips -o MANUAL.ps ALL.dvi
	/bin/rm -f MANUAL.ps.gz
	gzip MANUAL.ps

print:  MANUAL.ps.gz
	gunzip -c MANUAL.ps.gz | lpr

clean:
	/bin/rm -f ALL.* *~ *.aux *.log *.dvi
