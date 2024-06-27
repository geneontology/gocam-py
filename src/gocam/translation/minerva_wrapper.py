import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import ClassVar, Dict, Iterable, Iterator, Optional, Tuple, List

import requests
import requests_cache
import yaml
from oaklib import BasicOntologyInterface, get_adapter

from gocam.datamodel import Model, Activity, MolecularFunctionAssociation, BiologicalProcessAssociation, \
    CellularAnatomicalEntityAssociation, CausalAssociation, Object, EvidenceItem, ProvenanceInfo

INDEX_URL = "https://go-public.s3.amazonaws.com/files/gocam-models.json"
GOCAM_ENDPOINT = "https://api.geneontology.org/api/go-cam/"

ENABLED_BY = "RO:0002333"
PART_OF = "BFO:0000050"
OCCURS_IN = "BFO:0000066"

logger = logging.getLogger(__name__)


def _normalize_property(prop: str) -> str:
    """
    Normalize a property.

    Sometimes the JSON will use full URIs, sometimes just the local part
    """
    if "/" in prop:
        return prop.split("/")[-1]
    return prop

def _annotations(obj: Dict) -> Dict:
    """
    Extract annotations from an object (assumes single-valued).

    Annotations are lists of objects with keys "key" and "value".
    """
    return {_normalize_property(a["key"]): a["value"] for a in obj["annotations"]}

def _annotations_multivalued(obj: Dict) -> Dict:
    """
    Extract annotations from an object (assumes multi-valued).

    Annotations are lists of objects with keys "key" and "value".
    """
    anns = defaultdict(list)
    for a in obj["annotations"]:
        anns[a["key"]].append(a["value"])
    return anns

MAIN_TYPES = [
    "molecular_function",
    "biological_process",
    "cellular_component",
    "information biomacromolecule",
    "evidence",
    "chemical entity",
    "anatomical entity",
]

TYPE_MAPPING = {
    "molecular_function": "MolecularFunctionTermObject",
    "biological_process": "BiologicalProcessTermObject",
    "cellular_component": "CellularComponentTermObject",
    "information biomacromolecule": "InformationBiomacromoleculeTermObject",
    "evidence": "EvidenceTermObject",
    "chemical entity": "ChemicalEntityTermObject",
    "anatomical entity": "AnatomicalEntityTermObject",
}


@dataclass
class MinervaWrapper:
    """
    An Wrapper over the current GO API which returns "Minerva" JSON objects.

    TODO: Implement a fact counter to ensure all facts are encountered for and nothing dropped on floor
    """

    name: ClassVar[str] = "gocam"

    default_object_type = "Activity"

    _label_adapter: BasicOntologyInterface = None

    session: requests.Session = field(default_factory=lambda: requests.Session())


    def models(
        self, **kwargs
    ) -> Iterator[Model]:
        """
        Fetch all models from the API.

        :param kwargs:
        :return:
        """

        models = requests.get(INDEX_URL).json()
        for model in models:
            if "gocam" not in model:
                raise ValueError(f"Missing gocam in {model}")
            gocam = model["gocam"]
            yield self.object_by_id(gocam.replace("http://model.geneontology.org/", ""))

    def models_ids(
        self, **kwargs
    ) -> Iterator[Model]:
        models = requests.get(INDEX_URL).json()
        for model in models:
            if "gocam" not in model:
                raise ValueError(f"Missing gocam in {model}")
            gocam = model["gocam"]
            yield gocam.replace("http://model.geneontology.org/", "")


    def internal_object_by_id(self, object_id: str) -> Optional[Dict]:
        session = self.session
        if not object_id:
            raise ValueError(f"Missing object ID: {object_id}")
        logger.info(f"Fetch object: {object_id}")
        # response = session.get(MODEL_URL_TEMPLATE.format(model_id=object_id))
        local_id = object_id.replace("gocam:", "")
        response = session.get(f"{GOCAM_ENDPOINT}/{local_id}")
        response.raise_for_status()
        obj = response.json()
        return obj

    def object_by_id(self, object_id: str) -> Optional[Model]:
        session = self.session
        if not object_id:
            raise ValueError(f"Missing object ID: {object_id}")
        logger.info(f"Fetch object: {object_id}")
        local_id = object_id.replace("gocam:", "")
        response = session.get(f"{GOCAM_ENDPOINT}/{local_id}")
        response.raise_for_status()
        obj = response.json()
        return self.object_from_dict(obj)

    def object_from_dict(self, obj: Dict) -> Optional[Model]:
        id = obj["id"]
        logger.info(f"Processing model id: {id}")
        logger.debug(f"Model: {yaml.dump(obj)}")
        individuals = obj["individuals"]
        individual_to_type = {}  # individual ID to "root" type / category, e.g Evidence, BP, ..
        individual_to_term = {}
        individual_to_annotations = {}
        sources = set()

        id2obj = {}

        def _cls(obj: Dict) -> Optional[str]:
            if obj.get("type", None) == "complement":
                logger.warning(f"Ignoring Complement: {obj}")
                # class expression representing NOT
                return None
            if "id" not in obj:
                raise ValueError(f"No ID for {obj}")
            id = obj["id"]
            id2obj[id] = obj
            return id

        def _evidence_from_fact(fact: Dict) -> Tuple[EvidenceItem, List[ProvenanceInfo]]:
            anns_mv = _annotations_multivalued(fact)
            evidence_inst_ids = anns_mv.get("evidence", [])
            evs = []
            provs = []
            for evidence_inst_id in evidence_inst_ids:
                with_obj = individual_to_annotations.get(evidence_inst_id, {}).get("with", None)
                if with_obj:
                    with_objs = with_obj.split(" | ")
                else:
                    with_objs = None
                prov = ProvenanceInfo(
                    contributor=individual_to_annotations.get(evidence_inst_id, {}).get("contributor", None),
                    date=individual_to_annotations.get(evidence_inst_id, {}).get("date", None),
                )
                ev = EvidenceItem(
                    term=individual_to_term.get(evidence_inst_id, None),
                    reference=individual_to_annotations.get(evidence_inst_id, {}).get("source", None),
                    with_objects=with_objs,
                    provenances=[prov]
                )
                evs.append(ev)
            return evs, provs

        for i in individuals:
            typs = [x["label"] for x in i.get("root-type", []) if x]
            typ = None
            for t in typs:
                if t in MAIN_TYPES:
                    typ = t
                    break
            if not typ:
                logger.warning(f"Could not find type for {i}")
                continue
            individual_to_type[i["id"]] = typ
            terms = [_cls(x) for x in i.get("type", [])]
            terms = [x for x in terms if x]
            if len(terms) > 1:
                logger.warning(f"Multiple terms for {i}: {terms}")
            if not terms:
                logger.warning(f"No terms for {i}")
                continue
            individual_to_term[i["id"]] = terms[0]
            anns = _annotations(i)
            individual_to_annotations[i["id"]] = anns
            if "source" in anns:
                sources.add(anns["source"])
        activities = []
        relationships = []
        gene_ids = set()
        process_ids = set()
        activities_by_mf_id: Dict[str, Activity] = defaultdict(list)
        facts_by_property = defaultdict(list)
        for fact in obj["facts"]:
            facts_by_property[fact["property"]].append(fact)
        enabled_by_facts = facts_by_property.get(ENABLED_BY, [])
        if not enabled_by_facts:
            raise ValueError(f"Missing {ENABLED_BY} in {facts_by_property}")
        for fact in enabled_by_facts:
            s, o = fact["subject"], fact["object"]
            if s not in individual_to_term:
                logger.warning(f"Missing {s} in {individual_to_term}")
                continue
            if o not in individual_to_term:
                logger.warning(f"Missing {o} in {individual_to_term}")
                continue
            gene_id = individual_to_term[o]

            evs, provs = _evidence_from_fact(fact)
            activity = Activity(
                id=s,
                enabled_by=gene_id,
                molecular_function=MolecularFunctionAssociation(term=individual_to_term[s],
                                                                provenances=provs,
                                                                evidence=evs),
            )
            gene_ids.add(gene_id)
            activities.append(activity)
            activities_by_mf_id[s].append(activity)

        for fact in facts_by_property.get(PART_OF, []):
            s, o = fact["subject"], fact["object"]
            if o not in individual_to_term:
                logger.warning(f"Missing {o} in {individual_to_term}")
                continue
            for a in activities_by_mf_id.get(s, []):
                evs, _provs = _evidence_from_fact(fact)
                a.part_of = BiologicalProcessAssociation(term=individual_to_term[o],
                                                         evidence=evs)
                process_ids.add(individual_to_term[o])

        for fact in facts_by_property.get(OCCURS_IN, []):
            s, o = fact["subject"], fact["object"]
            if o not in individual_to_term:
                logger.warning(f"Missing {o} in {individual_to_term}")
                continue
            for a in activities_by_mf_id.get(s, []):
                evs, _provs = _evidence_from_fact(fact)
                a.occurs_in = CellularAnatomicalEntityAssociation(term=individual_to_term[o],
                                                                  evidence=evs)

        for p, facts in facts_by_property.items():
            for fact in facts:
                s, o = fact["subject"], fact["object"]
                sas = activities_by_mf_id.get(s, [])
                oas = activities_by_mf_id.get(o, [])
                if not sas or not oas:
                    continue
                if individual_to_type.get(s, None) != "molecular_function":
                    continue
                if individual_to_type.get(o, None) != "molecular_function":
                    continue
                sa = sas[0]
                oa = oas[0]
                evs, _provs = _evidence_from_fact(fact)
                rel = CausalAssociation(
                    predicate=p,
                    downstream_activity=oa.id,
                    evidence=evs,
                )
                if not sa.causal_associations:
                    sa.causal_associations = []
                sa.causal_associations.append(rel)
                # raise ValueError(f"{rel} in {sa}")
                relationships.append(rel)
        #pmids = {a.reference for a in activities if "reference" in a}
        #pmids = [p for p in pmids if p and p.startswith("PMID")]
        #pubs = self.pubmed_wrapper.objects_by_ids(pmids)
        #pubs_by_id = {p["id"]: p for p in pubs}
        #for a in activities:
        #    if "reference" in a:
        #        ref = a["reference"]
        #        if ref in pubs_by_id:
        #            a["reference_title"] = pubs_by_id[ref]["title"]
        annotations = _annotations(obj)
        annotations_mv = _annotations_multivalued(obj)
        objs = [
            Object(id=obj["id"], label=obj["label"]) for obj in id2obj.values()
        ]
        cam = Model(
            id=id,
            title=annotations["title"],
            status=annotations.get("state", None),
            comments=annotations_mv.get("comment", None),
            taxon=annotations.get("in_taxon", None),
            activities=activities,
            objects=objs,
        )
        return cam
