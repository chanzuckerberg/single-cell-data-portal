default: plan

.PHONY: init
init:
	cd envs/$(DEPLOYMENT_STAGE) && make init

.PHONY: plan
plan:
	cd envs/$(DEPLOYMENT_STAGE) && make plan

.PHONY: apply
apply:
	cd envs/$(DEPLOYMENT_STAGE) && make apply
