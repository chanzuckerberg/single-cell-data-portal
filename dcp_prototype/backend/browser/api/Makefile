build:
	mkdir -p chalicelib/config
	mkdir -p chalicelib/browser/rds
	mkdir -p chalicelib/browser/config
	cp -R ../rds chalicelib/browser
	cp -R ../config chalicelib/browser
	cp browser-api.yml chalicelib/config/browser-api.yml

deploy: build
	chalice deploy --no-autogen-policy --stage $(DEPLOYMENT_STAGE) --api-gateway-stage $(DEPLOYMENT_STAGE) --connection-timeout 240
