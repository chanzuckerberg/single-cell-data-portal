IMAGE=corpora-upload-${DEPLOYMENT_STAGE}

# make build
#	build a local image for testing
# make push
#	push the latest local image to the ECR repository

bump_version:
	expr `cat VERSION` + 1 > VERSION

build:
	mkdir -p code/chalice/api_server
	cp -R ../corpora/dataset_processing code/
	cp -R ../corpora/common code/
	cp -R ../../requirements.txt code/
	cp -R ../chalice/api_server/requirements.txt code/chalice/api_server/requirements.txt
	docker build -t $(IMAGE) .
	rm -rf code/

login:
	$(eval REPOSITORY_URI = $(shell aws ecr describe-repositories --repository-name $(IMAGE) | jq -r '.repositories[].repositoryUri'))
	$(eval ECR_URI = $(shell dirname $(REPOSITORY_URI)))
	aws ecr get-login-password --region us-west-2 | docker login --username AWS --password-stdin $(ECR_URI)

push: login 
	$(eval REPOSITORY_URI = $(shell aws ecr describe-repositories --repository-name $(IMAGE) | jq -r '.repositories[].repositoryUri'))
	$(eval VERSION=$(shell cat VERSION))
	docker tag $(IMAGE):latest $(REPOSITORY_URI):$(VERSION)
	docker push $(REPOSITORY_URI):$(VERSION)

promote: login
	$(eval REPOSITORY_URI = $(shell aws ecr describe-repositories --repository-name $(IMAGE) | jq -r '.repositories[].repositoryUri'))
	$(eval VERSION=$(shell cat VERSION))
	docker tag $(REPOSITORY_URI):$(VERSION) $(REPOSITORY_URI):active
	docker push $(REPOSITORY_URI):active

.phony: build bump_version promote push login
