## Add your own custom Makefile targets here

data/gocam.yaml:
	poetry run gocam fetch > $@.tmp && mv $@.tmp $@
