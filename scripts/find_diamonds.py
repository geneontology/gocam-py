#!/usr/bin/env python
"""
Script to find diamond patterns in GO-CAM models.

Diamond patterns typically indicate reactions that can be catalyzed by multiple paralogs.
"""

import sys
import yaml
import json
import logging
from pathlib import Path
from typing import List
import click

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from gocam.datamodel.gocam import Model
from gocam.analysis.diamond_detector import DiamondDetector

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def load_models(file_path: Path) -> List[Model]:
    """Load GO-CAM models from a file."""
    logger.info(f"Loading models from {file_path}")

    models_data = []
    with open(file_path, 'r') as f:
        if file_path.suffix == '.yaml':
            # Handle multi-document YAML
            for doc in yaml.safe_load_all(f):
                if doc:
                    models_data.append(doc)
        elif file_path.suffix == '.json':
            data = json.load(f)
            # Handle both single model and list of models
            if isinstance(data, dict):
                if 'models' in data:
                    models_data = data['models']
                else:
                    models_data = [data]
            elif isinstance(data, list):
                models_data = data
            else:
                raise ValueError(f"Unexpected data format in {file_path}")
        else:
            raise ValueError(f"Unsupported file format: {file_path.suffix}")

    models = []
    for model_data in models_data:
        try:
            model = Model(**model_data)
            models.append(model)
        except Exception as e:
            logger.warning(f"Failed to parse model: {e}")

    logger.info(f"Loaded {len(models)} models")
    return models


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file for the report')
@click.option('--format', '-f', type=click.Choice(['text', 'json', 'yaml']), default='text', help='Output format')
@click.option('--min-parallel', type=int, default=2, help='Minimum number of parallel paths to report')
def main(input_file, output, format, min_parallel):
    """Find diamond patterns in GO-CAM models."""

    input_path = Path(input_file)

    # Load models
    models = load_models(input_path)

    if not models:
        logger.error("No models loaded")
        return 1

    # Detect diamonds
    detector = DiamondDetector()
    diamonds = detector.find_diamonds_in_models(models)

    # Filter by minimum parallel paths
    if min_parallel > 2:
        diamonds = [d for d in diamonds if d.num_parallel_paths >= min_parallel]
        logger.info(f"Filtered to {len(diamonds)} diamonds with >= {min_parallel} parallel paths")

    # Generate output
    if format == 'text':
        output_content = detector.generate_report(models)
    elif format == 'json':
        output_data = []
        for diamond in diamonds:
            diamond_data = {
                'model_id': diamond.model_id,
                'model_title': diamond.model_title,
                'diamond_type': diamond.diamond_type,
                'is_pure': diamond.is_pure,
                'structure': {
                    'upstream_hops': diamond.upstream_hops,
                    'downstream_hops': diamond.downstream_hops,
                    'total_path_length': diamond.total_path_length,
                    'width': diamond.num_parallel_paths
                },
                'source': {
                    'id': diamond.source_activity_id,
                    'gene': diamond.activity_details.get(diamond.source_activity_id, {}).get('gene', 'Unknown') if diamond.activity_details else 'Unknown',
                    'molecular_function': diamond.activity_details.get(diamond.source_activity_id, {}).get('molecular_function', 'Unknown') if diamond.activity_details else 'Unknown'
                },
                'parallel_activities': [],
                'sink': {
                    'id': diamond.sink_activity_id,
                    'gene': diamond.activity_details.get(diamond.sink_activity_id, {}).get('gene', 'Unknown') if diamond.activity_details else 'Unknown',
                    'molecular_function': diamond.activity_details.get(diamond.sink_activity_id, {}).get('molecular_function', 'Unknown') if diamond.activity_details else 'Unknown'
                }
            }

            # Add parallel activity details
            for parallel_id in diamond.parallel_activity_ids:
                parallel_details = {
                    'id': parallel_id,
                    'gene': diamond.activity_details.get(parallel_id, {}).get('gene', 'Unknown') if diamond.activity_details else 'Unknown',
                    'molecular_function': diamond.activity_details.get(parallel_id, {}).get('molecular_function', 'Unknown') if diamond.activity_details else 'Unknown'
                }
                diamond_data['parallel_activities'].append(parallel_details)

            output_data.append(diamond_data)
        output_content = json.dumps(output_data, indent=2)
    elif format == 'yaml':
        output_data = []
        for diamond in diamonds:
            output_data.append({
                'model_id': diamond.model_id,
                'model_title': diamond.model_title,
                'source': diamond.source_activity_id,
                'parallel_activities': diamond.parallel_activity_ids,
                'sink': diamond.sink_activity_id,
                'num_parallel_paths': diamond.num_parallel_paths
            })
        output_content = yaml.dump(output_data, default_flow_style=False)

    # Write output
    if output:
        with open(output, 'w') as f:
            f.write(output_content)
        logger.info(f"Report written to {output}")
    else:
        print(output_content)

    # Summary
    logger.info(f"Found {len(diamonds)} diamond patterns across {len(models)} models")

    # Get paralog insights
    paralog_groups = detector.get_paralog_candidates(models)
    if paralog_groups:
        logger.info(f"Identified {len(paralog_groups)} potential paralog groups")

    return 0


if __name__ == '__main__':
    sys.exit(main())