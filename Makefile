# Makefile for the lazy

build_directory := build
install_directory := $(PWD)/install

all:
	ln -sf $(PWD)/cformat.sh $(shell git rev-parse --git-common-dir)/hooks/pre-commit
	cmake . -B$(build_directory) -DCMAKE_INSTALL_PREFIX=$(install_directory)
	cmake --build $(build_directory) -- $(MAKEFLAGS)
	cmake --install $(build_directory) | grep -v "Up-to-date"

clean:
	rm -rf $(build_directory) $(install_directory) ratpac.sh ratpac.csh RatpacConfig.cmake
