include ../common.mk

.PHONY: clean
clean:
	rm -f common/web/templates/index.html
	rm -rf common/web/static
	rm -f common/web/csp-hashes.json

.PHONY: unit-test
unit-test:
	PYTHONWARNINGS=ignore:ResourceWarning coverage run \
		--source=app,common,compute,dataset \
		--omit=.coverage,venv \
		-m unittest discover \
		--start-directory tests/unit \
		--top-level-directory .. \
		--verbose; test_result=$$?; \
	exit $$test_result \


.PHONY: test-annotations-performance
test-annotations-performance:
	python tests/performance/performance_test_annotations_backend.py

.PHONY: test-annotations-scale
test-annotations-scale:
	locust -f tests/performance/scale_test_annotations.py --headless -u 30 -r 10 --host https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/e/ --run-time 5m 2>&1 | tee locust_dev_stats.txt

.PHONY: e2e
e2e:
	python -m unittest discover \
	--start-directory tests/e2e \
	--top-level-directory .. \
	--verbose; test_result=$$?; \
