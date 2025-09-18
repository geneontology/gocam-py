"""
Detection of diamond patterns in GO-CAM causal flows.

Diamond patterns typically represent cases where reactions can be catalyzed by multiple paralogs.
A diamond pattern consists of:
- A common upstream activity (the "source")
- Multiple parallel activities in the middle (the "paralogs")
- A common downstream activity (the "sink")
"""

from collections import defaultdict
from typing import List, Dict, Set, Tuple, Optional
from dataclasses import dataclass, field
import logging
from gocam.datamodel.gocam import Model, Activity, CausalAssociation

logger = logging.getLogger(__name__)


@dataclass
class DiamondPattern:
    """Represents a diamond pattern in a causal flow."""
    source_activity_id: str
    parallel_activity_ids: List[str]
    sink_activity_id: str
    model_id: str
    model_title: str
    # Store activity details for enriched reporting
    activity_details: Optional[Dict[str, Dict[str, str]]] = None
    # Store path lengths and classification
    upstream_hops: int = 0  # Number of hops from source to parallel activities
    downstream_hops: int = 0  # Number of hops from parallel activities to sink
    intermediate_nodes: Optional[Dict[str, List[str]]] = None  # Intermediate nodes on each path
    classification: str = ""  # Classification of diamond type
    is_pure: bool = True  # Whether this is a pure diamond (no extra edges)

    @property
    def num_parallel_paths(self) -> int:
        """Number of parallel paths in the diamond."""
        return len(self.parallel_activity_ids)

    @property
    def diamond_type(self) -> str:
        """Classify the diamond based on its structure."""
        if self.upstream_hops == 1 and self.downstream_hops == 1:
            return "simple"  # Direct diamond: source -> parallel -> sink
        elif self.upstream_hops > 1 and self.downstream_hops == 1:
            return "upstream_complex"  # Multiple hops before parallelization
        elif self.upstream_hops == 1 and self.downstream_hops > 1:
            return "downstream_complex"  # Multiple hops after parallelization
        else:
            return "complex"  # Multiple hops on both sides

    @property
    def total_path_length(self) -> int:
        """Total path length from source to sink."""
        return self.upstream_hops + self.downstream_hops

    def __str__(self) -> str:
        return (f"Diamond in {self.model_id}: "
                f"{self.source_activity_id} -> "
                f"[{', '.join(self.parallel_activity_ids)}] -> "
                f"{self.sink_activity_id} "
                f"(type: {self.diamond_type}, hops: {self.upstream_hops}+{self.downstream_hops})")


class DiamondDetector:
    """Detects diamond patterns in GO-CAM models."""

    def __init__(self):
        self.diamonds: List[DiamondPattern] = []

    def _extract_activity_details(self, activity: Activity) -> Dict[str, str]:
        """Extract gene name and molecular function from an activity."""
        details = {}

        # Extract gene/complex information
        if activity.enabled_by:
            if hasattr(activity.enabled_by, 'term'):
                details['gene'] = activity.enabled_by.term
            else:
                details['gene'] = 'Unknown'
        else:
            details['gene'] = 'Unknown'

        # Extract molecular function
        if activity.molecular_function:
            if hasattr(activity.molecular_function, 'term'):
                details['molecular_function'] = activity.molecular_function.term
            else:
                details['molecular_function'] = 'Unknown'
        else:
            details['molecular_function'] = 'Unknown'

        return details

    def _check_diamond_purity(self, source: str, parallel_activities: List[str],
                               sink: str, downstream_map: Dict[str, Set[str]],
                               upstream_map: Dict[str, Set[str]]) -> bool:
        """
        Check if a diamond is pure (no extra edges beyond the diamond pattern).

        A pure diamond has:
        - Source only connects to parallel activities (no other outgoing edges)
        - Parallel activities only connect to sink (no other outgoing edges)
        - Parallel activities only receive from source (no other incoming edges)
        - Sink only receives from parallel activities (no other incoming edges)
        """
        # Check source: should only connect to parallel activities
        source_targets = downstream_map.get(source, set())
        if source_targets != set(parallel_activities):
            # Source has edges to nodes outside the diamond
            return False

        # Check each parallel activity
        for parallel in parallel_activities:
            # Parallel should only connect to sink
            parallel_targets = downstream_map.get(parallel, set())
            if parallel_targets != {sink}:
                return False

            # Parallel should only receive from source
            parallel_sources = upstream_map.get(parallel, set())
            if parallel_sources != {source}:
                return False

        # Check sink: should only receive from parallel activities
        sink_sources = upstream_map.get(sink, set())
        if sink_sources != set(parallel_activities):
            return False

        return True

    def _find_path_length(self, start: str, end: str, graph: Dict[str, Set[str]]) -> Tuple[int, List[str]]:
        """
        Find the shortest path length and intermediate nodes between two activities.
        Uses BFS to find shortest path.

        Returns:
            Tuple of (path_length, list_of_intermediate_nodes)
        """
        from collections import deque

        if start == end:
            return 0, []

        visited = {start}
        queue = deque([(start, 0, [])])

        while queue:
            node, dist, path = queue.popleft()

            if node in graph:
                for neighbor in graph[node]:
                    if neighbor == end:
                        return dist + 1, path

                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append((neighbor, dist + 1, path + [neighbor]))

        return -1, []  # No path found

    def find_diamonds_in_model(self, model: Model) -> List[DiamondPattern]:
        """
        Find all diamond patterns in a single GO-CAM model.

        Args:
            model: The GO-CAM model to analyze

        Returns:
            List of diamond patterns found in the model
        """
        if not model.activities:
            return []

        # Build activity ID to Activity object mapping
        activity_map = {activity.id: activity for activity in model.activities}

        # Build adjacency lists for the causal graph
        downstream_map = defaultdict(set)  # activity_id -> set of downstream activity_ids
        upstream_map = defaultdict(set)    # activity_id -> set of upstream activity_ids

        for activity in model.activities:
            if activity.causal_associations:
                for causal_assoc in activity.causal_associations:
                    if causal_assoc.downstream_activity:
                        downstream_map[activity.id].add(causal_assoc.downstream_activity)
                        upstream_map[causal_assoc.downstream_activity].add(activity.id)

        diamonds_found = []

        # For each activity, check if it's a source of a diamond
        for source_id in list(downstream_map.keys()):
            downstream_activities = downstream_map[source_id]

            # Need at least 2 downstream activities for potential diamond
            if len(downstream_activities) < 2:
                continue

            # Check if any pair of downstream activities reconverge
            downstream_list = list(downstream_activities)
            for i in range(len(downstream_list)):
                for j in range(i + 1, len(downstream_list)):
                    activity1 = downstream_list[i]
                    activity2 = downstream_list[j]

                    # Find common downstream activities (sinks)
                    common_downstream = downstream_map[activity1] & downstream_map[activity2]

                    for sink_id in common_downstream:
                        # Check if this is a true diamond (no direct connection from source to sink)
                        if sink_id not in downstream_map[source_id]:
                            # Found a diamond pattern with 2 parallel paths
                            parallel_activities = [activity1, activity2]

                            # Check if there are more parallel paths
                            for k in range(len(downstream_list)):
                                if k != i and k != j:
                                    activity_k = downstream_list[k]
                                    if sink_id in downstream_map[activity_k]:
                                        parallel_activities.append(activity_k)

                            # Calculate path lengths
                            # Check upstream path (source to first parallel activity)
                            upstream_hop, _ = self._find_path_length(source_id, parallel_activities[0], downstream_map)

                            # Check downstream path (first parallel activity to sink)
                            downstream_hop, _ = self._find_path_length(parallel_activities[0], sink_id, downstream_map)

                            # Check if diamond is pure (no extra edges)
                            is_pure = self._check_diamond_purity(
                                source_id, parallel_activities, sink_id,
                                downstream_map, upstream_map
                            )

                            # Extract activity details for enriched reporting
                            activity_details = {}
                            for act_id in [source_id] + parallel_activities + [sink_id]:
                                if act_id in activity_map:
                                    activity_details[act_id] = self._extract_activity_details(activity_map[act_id])
                                else:
                                    activity_details[act_id] = {'gene': 'Unknown', 'molecular_function': 'Unknown'}

                            diamond = DiamondPattern(
                                source_activity_id=source_id,
                                parallel_activity_ids=sorted(parallel_activities),
                                sink_activity_id=sink_id,
                                model_id=model.id,
                                model_title=model.title,
                                activity_details=activity_details,
                                upstream_hops=upstream_hop,
                                downstream_hops=downstream_hop,
                                is_pure=is_pure
                            )

                            # Avoid duplicates
                            if not self._is_duplicate_diamond(diamond, diamonds_found):
                                diamonds_found.append(diamond)

        return diamonds_found

    def _is_duplicate_diamond(self, diamond: DiamondPattern, existing_diamonds: List[DiamondPattern]) -> bool:
        """Check if a diamond pattern already exists in the list."""
        for existing in existing_diamonds:
            if (existing.source_activity_id == diamond.source_activity_id and
                existing.sink_activity_id == diamond.sink_activity_id and
                set(existing.parallel_activity_ids) == set(diamond.parallel_activity_ids)):
                return True
        return False

    def find_diamonds_in_models(self, models: List[Model]) -> List[DiamondPattern]:
        """
        Find all diamond patterns across multiple GO-CAM models.

        Args:
            models: List of GO-CAM models to analyze

        Returns:
            List of all diamond patterns found
        """
        all_diamonds = []
        for model in models:
            try:
                diamonds = self.find_diamonds_in_model(model)
                all_diamonds.extend(diamonds)
                if diamonds:
                    logger.info(f"Found {len(diamonds)} diamond(s) in model {model.id}")
            except Exception as e:
                logger.error(f"Error processing model {model.id}: {e}")

        self.diamonds = all_diamonds
        return all_diamonds

    def get_paralog_candidates(self, models: List[Model]) -> Dict[str, List[DiamondPattern]]:
        """
        Group diamond patterns by their parallel activities to identify potential paralogs.

        Returns:
            Dictionary mapping activity descriptions to diamond patterns they participate in
        """
        if not self.diamonds:
            self.find_diamonds_in_models(models)

        paralog_groups = defaultdict(list)

        # Build activity ID to enabled_by mapping
        activity_to_gene = {}
        for model in models:
            if model.activities:
                for activity in model.activities:
                    if activity.enabled_by and hasattr(activity.enabled_by, 'term'):
                        activity_to_gene[activity.id] = activity.enabled_by.term

        # Group diamonds by their parallel gene products
        for diamond in self.diamonds:
            parallel_genes = []
            for activity_id in diamond.parallel_activity_ids:
                if activity_id in activity_to_gene:
                    parallel_genes.append(activity_to_gene[activity_id])

            if len(parallel_genes) >= 2:
                key = tuple(sorted(parallel_genes))
                paralog_groups[str(key)].append(diamond)

        return dict(paralog_groups)

    def generate_report(self, models: List[Model]) -> str:
        """
        Generate a human-readable report of diamond patterns found.

        Args:
            models: List of GO-CAM models that were analyzed

        Returns:
            Formatted report as a string
        """
        if not self.diamonds:
            self.find_diamonds_in_models(models)

        report = []
        report.append("GO-CAM Diamond Pattern Analysis Report")
        report.append("=" * 50)
        report.append(f"\nTotal models analyzed: {len(models)}")
        report.append(f"Total diamond patterns found: {len(self.diamonds)}")

        if not self.diamonds:
            report.append("\nNo diamond patterns detected.")
            return "\n".join(report)

        # Analyze diamond types and purity
        type_counts = defaultdict(int)
        width_counts = defaultdict(int)
        pure_count = 0
        impure_count = 0

        for diamond in self.diamonds:
            type_counts[diamond.diamond_type] += 1
            width_counts[diamond.num_parallel_paths] += 1
            if diamond.is_pure:
                pure_count += 1
            else:
                impure_count += 1

        report.append("\nDiamond Purity:")
        report.append(f"  Pure diamonds: {pure_count} ({round(pure_count * 100 / len(self.diamonds)) if self.diamonds else 0}%)")
        report.append(f"  Impure diamonds: {impure_count} ({round(impure_count * 100 / len(self.diamonds)) if self.diamonds else 0}%)")

        report.append("\nDiamond Type Distribution:")
        for dtype, count in sorted(type_counts.items()):
            report.append(f"  {dtype}: {count}")

        report.append("\nParallel Path Width Distribution:")
        for width, count in sorted(width_counts.items()):
            report.append(f"  {width} paths: {count} diamonds")

        # Group diamonds by model
        diamonds_by_model = defaultdict(list)
        for diamond in self.diamonds:
            diamonds_by_model[diamond.model_id].append(diamond)

        report.append(f"\nModels with diamond patterns: {len(diamonds_by_model)}")
        report.append("\n" + "-" * 50)

        for model_id, model_diamonds in diamonds_by_model.items():
            model_title = model_diamonds[0].model_title if model_diamonds else "Unknown"
            report.append(f"\nModel: {model_id}")
            report.append(f"Title: {model_title}")
            report.append(f"Diamond patterns: {len(model_diamonds)}")

            for i, diamond in enumerate(model_diamonds, 1):
                report.append(f"  Diamond {i}:")
                report.append(f"    Type: {diamond.diamond_type}")
                report.append(f"    Purity: {'Pure' if diamond.is_pure else 'Impure (has extra edges)'}")
                report.append(f"    Structure: {diamond.upstream_hops} hop(s) upstream, {diamond.downstream_hops} hop(s) downstream")
                report.append(f"    Width: {diamond.num_parallel_paths} parallel paths")

                # Source activity details
                source_details = diamond.activity_details.get(diamond.source_activity_id, {}) if diamond.activity_details else {}
                source_gene = source_details.get('gene', 'Unknown')
                source_mf = source_details.get('molecular_function', 'Unknown')
                report.append(f"    Source: {diamond.source_activity_id}")
                report.append(f"      Gene: {source_gene}")
                report.append(f"      MF: {source_mf}")

                # Parallel activities details
                report.append(f"    Parallel paths ({diamond.num_parallel_paths}):")
                for parallel_id in diamond.parallel_activity_ids:
                    parallel_details = diamond.activity_details.get(parallel_id, {}) if diamond.activity_details else {}
                    parallel_gene = parallel_details.get('gene', 'Unknown')
                    parallel_mf = parallel_details.get('molecular_function', 'Unknown')
                    report.append(f"      - {parallel_id}")
                    report.append(f"        Gene: {parallel_gene}")
                    report.append(f"        MF: {parallel_mf}")

                # Sink activity details
                sink_details = diamond.activity_details.get(diamond.sink_activity_id, {}) if diamond.activity_details else {}
                sink_gene = sink_details.get('gene', 'Unknown')
                sink_mf = sink_details.get('molecular_function', 'Unknown')
                report.append(f"    Sink: {diamond.sink_activity_id}")
                report.append(f"      Gene: {sink_gene}")
                report.append(f"      MF: {sink_mf}")

        # Add paralog analysis
        paralog_groups = self.get_paralog_candidates(models)
        if paralog_groups:
            report.append("\n" + "=" * 50)
            report.append("Potential Paralog Groups")
            report.append("-" * 50)

            for genes, diamonds in paralog_groups.items():
                report.append(f"\nGene products: {genes}")
                report.append(f"Appears in {len(diamonds)} diamond pattern(s)")

        return "\n".join(report)