## Add your own custom Makefile targets here

RUN = uv run

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
