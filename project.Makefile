## Add your own custom Makefile targets here

RUN = uv run
Q = linkml-store  -d gocams::main
QI = linkml-store  -d gocams::indexed

HM = linkml-store plot heatmap
HMCY = $(HM) --cluster y --cluster-method ward --cluster-metric cosine

data/gocam.yaml:
	$(RUN) gocam fetch -f yaml > $@.tmp && mv $@.tmp $@
.PRECIOUS: data/gocam.yaml

data/gocam-indexed.yaml: data/gocam.yaml
	$(RUN) gocam index-models $< -o $@.tmp && mv $@.tmp $@
.PRECIOUS: data/gocam.yaml

data/gocam-indexed.json: data/gocam.yaml
	$(RUN) gocam index-models -O json $< -o $@.tmp && mv $@.tmp $@
.PRECIOUS: data/gocam.yaml

data/gocam-flattened.jsonl: data/gocam-indexed.json
	$(RUN) gocam flatten-models -O jsonl $< -o $@.tmp && mv $@.tmp $@
.PRECIOUS: data/gocam-flattened.yaml

data/gocam.owl: data/gocam.yaml
	$(RUN) gocam convert -O owl $< -o $@.tmp && mv $@.tmp $@
.PRECIOUS: data/gocam.owl

data/gocam.cx2: data/gocam.yaml
	$(RUN) gocam convert -O cx2 $< -o $@.tmp && mv $@.tmp $@

mongodb-load:
	linkml-store -d gocams insert -f yamll --replace data/gocam.yaml
	linkml-store -d gocams::flattened insert --replace data/gocam-flattened.jsonl

mongodb-load-flattened: data/gocam-flattened.jsonl
	linkml-store -d gocams -c flattened insert --replace $<

data/gocam-flattened.slim.json: data/gocam-flattened.jsonl
	linkml-store -d gocams::flattened query -s "[id, title,taxon,status,model_activity_part_of_rollup_label,model_activity_occurs_in_rollup_label,model_activity_enabled_by_terms_id,number_of_activities,length_of_longest_causal_association_path,number_of_strongly_connected_components]" -O json -o $@.tmp && mv $@.tmp $@
.PRECIOUS: data/gocam-flattened.slim.json

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

reports/activity-counts-by-provided-by.csv:
	$(QI) fq -S provenances.provided_by+query_index.number_of_activities -O csv -o $@

reports/rows.csv:
	$(QI) query -s "[taxon,query_index.number_of_activities,query_index.number_of_enabled_by_terms,query_index.number_of_causal_associations,query_index.number_of_strongly_connected_components,query_index.length_of_longest_causal_association_path]" -O csv -o $@
.PRECIOUS: reports/rows.csv

reports/heatmap-taxon-by-activities.png: reports/rows.csv
	$(HM) -f csv -x query_index.number_of_activities -y taxon $< -o $@

reports/clustered-heatmap-taxon-by-activities.png: reports/rows.csv
	$(HMCY) -f csv -x query_index.number_of_activities -y taxon $< -o $@

reports/clustered-heatmap-taxon-by-enablers.png: reports/rows.csv
	$(HMCY) -f csv -x query_index.number_of_enabled_by_terms -y taxon $< -o $@

reports/clustered-heatmap-taxon-by-sccs.png: reports/rows.csv
	$(HMCY) -f csv -x query_index.number_of_strongly_connected_components -y taxon $< -o $@

reports/clustered-heatmap-taxon-by-longest-path.png: reports/rows.csv
	$(HMCY) -f csv -x query_index.length_of_longest_causal_association_path -y taxon $< -o $@

reports/model-pmid.csv:
	$(Q) query |  jq -r '.[] | .id as $$id | [.. | objects | .reference? | select(. != null and startswith("PMID:"))] | .[] | [$$id, .] | @tsv' > $@

reports/describe.txt:
	$(Q) describe -o $@

reports/model-counts-by-annoton.csv:
	$(QI) fq -S query_index.annoton_terms.label -O csv -o $@
