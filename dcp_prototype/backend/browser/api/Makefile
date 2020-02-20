build:
	mkdir -p chalicelib/browser/code
	mkdir -p chalicelib/config
	cp -R ../code chalicelib/browser
	cp browser-api.yml chalicelib/config/browser-api.yml

deploy: build
	chalice deploy --no-autogen-policy --stage $(DEPLOYMENT_STAGE) --api-gateway-stage $(DEPLOYMENT_STAGE) --connection-timeout 240
