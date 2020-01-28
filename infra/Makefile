DIRS=${shell find . -name "*.tf" -exec dirname {} \; | sort --unique}
COMPONENTS=${shell for d in $(DIRS); do basename $$d; done}

all: init-all

init-all:
	@for c in $(COMPONENTS); do \
		$(MAKE)	init COMPONENT=$$c || exit 1; \
	done

plan-all:
	@for c in $(COMPONENTS); do \
		$(MAKE)	plan COMPONENT=$$c || exit 1; \
	done

apply-all:
	@for c in $(COMPONENTS); do \
		$(MAKE)	apply COMPONENT=$$c || exit 1; \
	done

destroy-all:
	@for c in $(COMPONENTS); do \
		$(MAKE)	destroy COMPONENT=$$c || exit 1; \
	done

clean-all:
	@for c in $(COMPONENTS); do \
		$(MAKE)	clean COMPONENT=$$c || exit 1; \
	done

plan: init
	@echo $(COMPONENTS)
	cd $(COMPONENT); terraform plan -detailed-exitcode

apply: init
	@echo $(COMPONENTS)
	cd $(COMPONENT); terraform apply

destroy: init
	cd $(COMPONENT); terraform destroy

init:
	rm -rf $(COMPONENT)/.terraform/*.tfstate
	./build_deploy_config.py $(COMPONENT)
	cd $(COMPONENT); terraform init;

clean:
	cd $(COMPONENT); rm -rf .terraform

.PHONY: init-all plan-all apply-all clean-all plan apply destroy init clean
