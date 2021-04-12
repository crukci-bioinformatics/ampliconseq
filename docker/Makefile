
version ?= latest

build:
	docker build --tag crukcibioinformatics/ampliconseq:${version} .

release: build
	docker push crukcibioinformatics/ampliconseq:${version}

