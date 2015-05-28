#!/bin/sh


getModifiedDates() {
  sourceFile=$1
  targetFile=$2

  sourceDate=`stat -s $sourceFile | sed -n -E "s/.*st_mtime=([0-9]+).*/\1/p"`
  targetDate=`stat -s $targetFile | sed -n -E "s/.*st_mtime=([0-9]+).*/\1/p"`
}


getModifiedDates article.Rnw article.tex

if [ $sourceDate -ge $targetDate ]; then
  R --vanilla CMD Sweave article.Rnw
  [ $? -gt 0 ] && exit 1
fi

buildBib=0

getModifiedDates blme.bib article.pdf

[ $sourceDate -ge $targetDate ] && buildBib=1

getModifiedDates article.tex article.pdf

if [ $buildBib -gt 0 ]; then
  R CMD pdflatex article.tex && \
  R CMD bibtex article.aux && \
  R CMD bibtex article.aux
  [ $? -gt 0 ] && exit 1
fi

if [ $sourceDate -ge $targetDate ]; then
  R CMD pdflatex article.tex && \
  R CMD pdflatex article.tex
fi