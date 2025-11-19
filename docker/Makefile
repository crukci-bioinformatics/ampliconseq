
version ?= latest

build:
	docker build --tag crukcibioinformatics/ampliconseq:${version} .

release: build
	docker push crukcibioinformatics/ampliconseq:${version}

singularity: build
	singularity build ampliconseq-${version}.sif docker-daemon://crukcibioinformatics/ampliconseq:${version}

