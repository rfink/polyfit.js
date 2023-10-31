# Copyright © 2023, P. Lutus and Ryan Fink. Released under the GPL https://www.gnu.org/licenses/gpl-3.0.en.html
REPORTER = spec
TESTS = $(wildcard test/*.test.js)
DROPAFTER=false

components:
	@component install --dev

build: components
	@component build

test:
	@NODE_ENV=test ./node_modules/.bin/mocha $(TESTS) \
		--require "should" \
		--growl \
		--reporter $(REPORTER)

.PHONY: test
