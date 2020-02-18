build:
	mkdir -p chalicelib/browser/rds
	cp ../rds/__init__.py chalicelib/browser/rds
	cp ../rds/browser_orm.py chalicelib/browser/rds
	cp ../rds/db_config.py chalicelib/browser/rds
	mkdir -p chalicelib/config
	cp browser-api.yml chalicelib/config/browser-api.yml
	shopt -s nullglob; for wheel in vendor.in/*/*.whl; do unzip -q -o -d vendor $$wheel; done

deploy: build
	chalice deploy --no-autogen-policy --stage $(DEPLOYMENT_STAGE) --api-gateway-stage $(DEPLOYMENT_STAGE) --connection-timeout 240
