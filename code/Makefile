
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
	mkdir -p Experiments Experiments_Industrial Experiments_Industrial_paperWithoutSA Experiments_Industrial_Test
	7z x -y data.zip.001
