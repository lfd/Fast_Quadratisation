
.PHONY: all clean

all: ext data plots

clean:
	rm -rf RPlots/
	rm -rf build/
	rm -rf img-gen/

pre_compile:
	cp -r pics/* RPlots/ 
	./pre_compile_pics.sh

plots: 
	Rscript plots.r

data:
	python Test.py

extract:
	7z x -y data.zip.001
