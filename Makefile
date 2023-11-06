include set_env.mk

## build-image                                 : Build the docker image for notebooks
build-image: Dockerfile
	docker build -t agbiome/rna_seq_test . \
	&& touch $@

## run-notebook                                : Run the notebook server
.PHONY: run-notebook
run-notebook: build-image stop-notebook
	DATA_DIR=$(DATA_DIR) ./jupyter-lab.sh

## stop-notebook                               : Stop the notebook server
.PHONY: stop-notebook
stop-notebook:
	docker stop rna_seq_test_$(HOSTNAME)_$(UID) || echo "already stopped"
	sleep 1

## get-url                                     : Get the URL of the notebook server
.PHONY: get-url
get-url:
	@docker exec --user jovyan -it rna_seq_test_$(HOSTNAME)_$(UID) jupyter lab list |     sed -e 's/http:\/\/[a-zA-Z0-9]*/http:\/\/localhost/'
