# Definitions
repository = gcr.io/ucsd-medicine-cast
identifier = plink2
version = 0.0.0
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = ${repository}/${identifier}
tag = ${version}--${git_commit}

# Steps
build:
	# do the docker build
	docker build --no-cache -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest

push:
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest