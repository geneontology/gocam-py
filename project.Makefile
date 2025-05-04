## Add your own custom Makefile targets here

RUN = poetry run
Q = linkml-store  -d gocams::main

data/gocam.yaml:
	$(RUN) gocam fetch -f yaml > $@.tmp && mv $@.tmp $@


MODEL_COUNTS_BY = taxon term enabled-by mf occurs-in part-of provided-by causal-edge-predicate

all_reports: $(patsubst %, reports/model-counts-by-%.csv, $(MODEL_COUNTS_BY))

reports/model-counts-by-taxon.csv:
	$(Q) fq -S taxon -O csv -o $@

reports/model-counts-by-term.csv:
	$(Q) fq -S objects -O csv -o $@

reports/model-counts-by-enabled-by.csv:
	$(Q) fq -S activities.enabled_by.term -O csv -o $@

reports/model-counts-by-mf.csv:
	$(Q) fq -S activities.molecular_function.term -O csv -o $@

reports/model-counts-by-occurs-in.csv:
	$(Q) fq -S activities.occurs_in.term -O csv -o $@

reports/model-counts-by-part-of.csv:
	$(Q) fq -S activities.part_of.term -O csv -o $@

reports/model-counts-by-provided-by.csv:
	$(Q) fq -S provenances.provided_by -O csv -o $@

reports/model-counts-by-causal-edge-predicate.csv:
	$(Q) fq -S activities.causal_associations.predicate -O csv -o $@

reports/model-pmid.csv:
	$(Q) query |  jq -r '.[] | .id as $$id | [.. | objects | .reference? | select(. != null and startswith("PMID:"))] | .[] | [$$id, .] | @tsv' > $@

reports/describe.txt:
	$(Q) describe -o $@
