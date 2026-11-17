.DEFAULT_GOAL := library

.PHONY: format format-check

format format-check:
	$(MAKE) -C metal-bsr-prototype $@

.PHONY: library install
library install:
	$(MAKE) -C metal-bsr-prototype $@

.PHONY: api-example api-example-run
api-example api-example-run:
	$(MAKE) -C metal-bsr-prototype $@
