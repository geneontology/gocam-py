from dataclasses import dataclass
from enum import Enum


class Color(str, Enum):
    BLUE = "#6495ED"
    PURPLE = "#800080"
    DARK_SLATE = "#2F4F4F"
    LIGHT_RED = "#fF9999"
    RED = "#FF0000"
    GREEN = "#008800"
    LIGHT_GREEN = "#95e095"
    PINK = "#ED6495"
    DARK_BLUE = "#483D8B"
    LIGHT_BLUE = "#add8e6"
    WHITE = "#FFFAFA"
    BLACK = "#000000"
    GREY = "#CCCCCC"
    PALE_AQUA = "#B2DFDB"
    LAVENDER_PINK = "#E2BDE7"
    MINT_GREEN = "#C8E6C9"


class Width(int, Enum):
    DEFAULT = 5


class ArrowShape(str, Enum):
    CIRCLE = "circle"
    DIAMOND = "diamond"
    SQUARE = "square"
    TEE = "tee"
    TRIANGLE = "triangle"


class LineStyle(str, Enum):
    DASHED = "dashed"
    SOLID = "solid"


class NodeType(str, Enum):
    GENE = "gene"
    COMPLEX = "complex"
    MOLECULE = "molecule"


class RelationType(str, Enum):
    CAUSALLY_UPSTREAM_OF_POSITIVE_EFFECT = "RO:0002304"
    CAUSALLY_UPSTREAM_OF_NEGATIVE_EFFECT = "RO:0002305"
    CONSTITUTIVELY_UPSTREAM_OF = "RO:0012009"
    DIRECTLY_NEGATIVELY_REGULATES = "RO:0002630"
    DIRECTLY_POSITIVELY_REGULATES = "RO:0002629"
    HAS_INPUT = "RO:0002233"
    HAS_OUTPUT = "RO:0002234"
    INDIRECTLY_NEGATIVELY_REGULATES = "RO:0002407"
    INDIRECTLY_POSITIVELY_REGULATES = "RO:0002409"
    IS_SMALL_MOLECULE_INHIBITOR_OF = "RO:0012006"
    IS_SMALL_MOLECULE_ACTIVATOR_OF = "RO:0012005"
    NEGATIVELY_REGULATES = "RO:0002212"
    POSITIVELY_REGULATES = "RO:0002213"
    PROVIDES_INPUT_FOR = "RO:0002413"
    REMOVES_INPUT_FOR = "RO:0012010"


class NodeShape(str, Enum):
    ELLIPSE = "ellipse"
    OCTAGON = "octagon"
    RECTANGLE = "rectangle"
    ROUND_RECTANGLE = "round-rectangle"


@dataclass
class RelationStyle:
    line_style: LineStyle
    arrow_shape: ArrowShape
    label: str
    color: Color
    width: Width


@dataclass
class NodeStyle:
    color: Color
    height: int
    label_font_size: int
    label_max_width: int
    shape: NodeShape
    width: int


RELATIONS = {
    RelationType.CAUSALLY_UPSTREAM_OF_POSITIVE_EFFECT: RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.TRIANGLE,
        label="causally upstream of, positive effect",
        color=Color.LIGHT_GREEN,
        width=Width.DEFAULT,
    ),
    RelationType.CAUSALLY_UPSTREAM_OF_NEGATIVE_EFFECT: RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.TEE,
        label="causally upstream of, negative effect",
        color=Color.LIGHT_RED,
        width=Width.DEFAULT,
    ),
    RelationType.CONSTITUTIVELY_UPSTREAM_OF: RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.CIRCLE,
        label="constitutively upstream of",
        color=Color.LIGHT_GREEN,
        width=Width.DEFAULT,
    ),
    RelationType.DIRECTLY_NEGATIVELY_REGULATES: RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.TEE,
        label="directly negatively regulates",
        color=Color.RED,
        width=Width.DEFAULT,
    ),
    RelationType.DIRECTLY_POSITIVELY_REGULATES: RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.TRIANGLE,
        label="directly positively regulates",
        color=Color.GREEN,
        width=Width.DEFAULT,
    ),
    RelationType.HAS_INPUT: RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="has input",
        color=Color.BLUE,
        width=Width.DEFAULT,
    ),
    RelationType.HAS_OUTPUT: RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="has output",
        color=Color.PINK,
        width=Width.DEFAULT,
    ),
    RelationType.INDIRECTLY_NEGATIVELY_REGULATES: RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.TEE,
        label="indirectly negatively regulates",
        color=Color.RED,
        width=Width.DEFAULT,
    ),
    RelationType.INDIRECTLY_POSITIVELY_REGULATES: RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.TRIANGLE,
        label="indirectly positively regulates",
        color=Color.GREEN,
        width=Width.DEFAULT,
    ),
    RelationType.IS_SMALL_MOLECULE_INHIBITOR_OF: RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.TEE,
        label="is small molecule inhibitor of",
        color=Color.RED,
        width=Width.DEFAULT,
    ),
    RelationType.IS_SMALL_MOLECULE_ACTIVATOR_OF: RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.TRIANGLE,
        label="is small molecule activator of",
        color=Color.GREEN,
        width=Width.DEFAULT,
    ),
    RelationType.NEGATIVELY_REGULATES: RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.TEE,
        label="negatively regulates",
        color=Color.RED,
        width=Width.DEFAULT,
    ),
    RelationType.POSITIVELY_REGULATES: RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.TRIANGLE,
        label="positively regulates",
        color=Color.GREEN,
        width=Width.DEFAULT,
    ),
    RelationType.PROVIDES_INPUT_FOR: RelationStyle(
        line_style=LineStyle.SOLID,
        # The widget uses `circle-triangle` in this case, but that shape is not part of the CX2
        # spec. So we need to intentionally deviate and use a supported shape here.
        # See: https://cytoscape.org/cx/cx2/cx2-visual-styles/#2361-edge-source-arrow-shape
        arrow_shape=ArrowShape.CIRCLE,
        label="provides input for",
        color=Color.PURPLE,
        width=Width.DEFAULT,
    ),
    RelationType.REMOVES_INPUT_FOR: RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.SQUARE,
        label="removes input for",
        color=Color.LIGHT_RED,
        width=Width.DEFAULT,
    ),
}

NODE_STYLES = {
    NodeType.GENE: NodeStyle(
        color=Color.MINT_GREEN,
        height=40,
        label_font_size=12,
        label_max_width=80,
        shape=NodeShape.RECTANGLE,
        width=85,
    ),
    NodeType.COMPLEX: NodeStyle(
        color=Color.LAVENDER_PINK,
        height=70,
        label_font_size=15,
        label_max_width=100,
        shape=NodeShape.RECTANGLE,
        width=110,
    ),
    NodeType.MOLECULE: NodeStyle(
        color=Color.PALE_AQUA,
        height=40,
        label_font_size=10,
        label_max_width=70,
        shape=NodeShape.ELLIPSE,
        width=70,
    ),
}

VISUAL_PROPERTIES = {
    "default": {
        "edge": {
            "EDGE_CURVED": True,
            "EDGE_LABEL_AUTOROTATE": False,
            "EDGE_LABEL_BACKGROUND_COLOR": "#FFFFFF",
            "EDGE_LABEL_BACKGROUND_OPACITY": 1,
            "EDGE_LABEL_BACKGROUND_SHAPE": "rectangle",
            "EDGE_LABEL_COLOR": "#000000",
            "EDGE_LABEL_FONT_FACE": {
                "FONT_FAMILY": "sans-serif",
                "FONT_NAME": "Dialog",
                "FONT_STYLE": "normal",
                "FONT_WEIGHT": "normal",
            },
            "EDGE_LABEL_FONT_SIZE": 10,
            "EDGE_LABEL_MAX_WIDTH": 200,
            "EDGE_LABEL_OPACITY": 1,
            "EDGE_LABEL_POSITION": {
                "EDGE_ANCHOR": "C",
                "JUSTIFICATION": "center",
                "LABEL_ANCHOR": "C",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
            },
            "EDGE_LABEL_ROTATION": 0,
            "EDGE_LINE_COLOR": "#848484",
            "EDGE_LINE_STYLE": "solid",
            "EDGE_OPACITY": 1,
            "EDGE_SELECTED": "false",
            "EDGE_SELECTED_PAINT": "#FFFF00",
            "EDGE_SOURCE_ARROW_COLOR": "#000000",
            "EDGE_SOURCE_ARROW_SELECTED_PAINT": "#FFFF00",
            "EDGE_SOURCE_ARROW_SHAPE": "none",
            "EDGE_SOURCE_ARROW_SIZE": 6,
            "EDGE_STACKING": "AUTO_BEND",
            "EDGE_STACKING_DENSITY": 0.5,
            "EDGE_STROKE_SELECTED_PAINT": "#FFFF00",
            "EDGE_TARGET_ARROW_COLOR": "#000000",
            "EDGE_TARGET_ARROW_SELECTED_PAINT": "#FFFF00",
            "EDGE_TARGET_ARROW_SHAPE": "none",
            "EDGE_TARGET_ARROW_SIZE": 6,
            "EDGE_VISIBILITY": "element",
            "EDGE_WIDTH": 2,
            "EDGE_Z_ORDER": 0,
        },
        "network": {"NETWORK_BACKGROUND_COLOR": "#FFFFFF"},
        "node": {
            "COMPOUND_NODE_PADDING": "10.0",
            "COMPOUND_NODE_SHAPE": "ROUND_RECTANGLE",
            "NODE_BACKGROUND_COLOR": "#FFFFFF",
            "NODE_BACKGROUND_OPACITY": 1,
            "NODE_BORDER_COLOR": "#CCCCCC",
            "NODE_BORDER_OPACITY": 1,
            "NODE_BORDER_STYLE": "solid",
            "NODE_BORDER_WIDTH": 1,
            "NODE_CUSTOMGRAPHICS_POSITION_1": {
                "ENTITY_ANCHOR": "C",
                "GRAPHICS_ANCHOR": "C",
                "JUSTIFICATION": "center",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
            },
            "NODE_CUSTOMGRAPHICS_POSITION_2": {
                "ENTITY_ANCHOR": "C",
                "GRAPHICS_ANCHOR": "C",
                "JUSTIFICATION": "center",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
            },
            "NODE_CUSTOMGRAPHICS_POSITION_3": {
                "ENTITY_ANCHOR": "C",
                "GRAPHICS_ANCHOR": "C",
                "JUSTIFICATION": "center",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
            },
            "NODE_CUSTOMGRAPHICS_POSITION_4": {
                "ENTITY_ANCHOR": "C",
                "GRAPHICS_ANCHOR": "C",
                "JUSTIFICATION": "center",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
            },
            "NODE_CUSTOMGRAPHICS_POSITION_5": {
                "ENTITY_ANCHOR": "C",
                "GRAPHICS_ANCHOR": "C",
                "JUSTIFICATION": "center",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
            },
            "NODE_CUSTOMGRAPHICS_POSITION_6": {
                "ENTITY_ANCHOR": "C",
                "GRAPHICS_ANCHOR": "C",
                "JUSTIFICATION": "center",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
            },
            "NODE_CUSTOMGRAPHICS_POSITION_7": {
                "ENTITY_ANCHOR": "C",
                "GRAPHICS_ANCHOR": "C",
                "JUSTIFICATION": "center",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
            },
            "NODE_CUSTOMGRAPHICS_POSITION_8": {
                "ENTITY_ANCHOR": "C",
                "GRAPHICS_ANCHOR": "C",
                "JUSTIFICATION": "center",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
            },
            "NODE_CUSTOMGRAPHICS_POSITION_9": {
                "ENTITY_ANCHOR": "C",
                "GRAPHICS_ANCHOR": "C",
                "JUSTIFICATION": "center",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
            },
            "NODE_CUSTOMGRAPHICS_SIZE_1": 50,
            "NODE_CUSTOMGRAPHICS_SIZE_2": 50,
            "NODE_CUSTOMGRAPHICS_SIZE_3": 50,
            "NODE_CUSTOMGRAPHICS_SIZE_4": 50,
            "NODE_CUSTOMGRAPHICS_SIZE_5": 50,
            "NODE_CUSTOMGRAPHICS_SIZE_6": 50,
            "NODE_CUSTOMGRAPHICS_SIZE_7": 50,
            "NODE_CUSTOMGRAPHICS_SIZE_8": 50,
            "NODE_CUSTOMGRAPHICS_SIZE_9": 50,
            "NODE_HEIGHT": 35,
            "NODE_LABEL_BACKGROUND_COLOR": "#B6B6B6",
            "NODE_LABEL_BACKGROUND_OPACITY": 1,
            "NODE_LABEL_BACKGROUND_SHAPE": "none",
            "NODE_LABEL_COLOR": "#000000",
            "NODE_LABEL_FONT_FACE": {
                "FONT_FAMILY": "sans-serif",
                "FONT_NAME": "SansSerif",
                "FONT_STYLE": "normal",
                "FONT_WEIGHT": "normal",
            },
            "NODE_LABEL_FONT_SIZE": 12,
            "NODE_LABEL_MAX_WIDTH": 200,
            "NODE_LABEL_OPACITY": 1,
            "NODE_LABEL_POSITION": {
                "HORIZONTAL_ALIGN": "center",
                "HORIZONTAL_ANCHOR": "center",
                "JUSTIFICATION": "center",
                "MARGIN_X": 0,
                "MARGIN_Y": 0,
                "VERTICAL_ALIGN": "center",
                "VERTICAL_ANCHOR": "center",
            },
            "NODE_LABEL_ROTATION": 0,
            "NODE_SELECTED": False,
            "NODE_SELECTED_PAINT": "#FFFF00",
            "NODE_SHAPE": "round-rectangle",
            "NODE_VISIBILITY": "element",
            "NODE_WIDTH": 75,
            "NODE_X_LOCATION": 0,
            "NODE_Y_LOCATION": 0,
            "NODE_Z_LOCATION": 0,
        },
    },
    "edgeMapping": {
        "EDGE_LINE_COLOR": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "represents",
                "map": [
                    {"v": key, "vp": value.color} for key, value in RELATIONS.items()
                ],
                "type": "string",
            },
        },
        "EDGE_LINE_STYLE": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "represents",
                "map": [
                    {"v": key, "vp": value.line_style}
                    for key, value in RELATIONS.items()
                ],
                "type": "string",
            },
        },
        "EDGE_TARGET_ARROW_COLOR": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "represents",
                "map": [
                    {"v": key, "vp": value.color} for key, value in RELATIONS.items()
                ],
                "type": "string",
            },
        },
        "EDGE_TARGET_ARROW_SHAPE": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "represents",
                "map": [
                    {"v": key, "vp": value.arrow_shape}
                    for key, value in RELATIONS.items()
                ],
                "type": "string",
            },
        },
        "EDGE_WIDTH": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "represents",
                "map": [
                    {"v": key, "vp": value.width} for key, value in RELATIONS.items()
                ],
                "type": "integer",
            },
        },
    },
    "nodeMapping": {
        "NODE_BACKGROUND_COLOR": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "type",
                "map": [
                    {"v": key, "vp": value.color} for key, value in NODE_STYLES.items()
                ],
                "type": "string",
            },
        },
        "NODE_HEIGHT": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "type",
                "map": [
                    {"v": key, "vp": value.height} for key, value in NODE_STYLES.items()
                ],
                "type": "string",
            },
        },
        "NODE_LABEL": {
            "type": "PASSTHROUGH",
            "definition": {"attribute": "name", "type": "string"},
        },
        "NODE_LABEL_FONT_SIZE": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "type",
                "map": [
                    {"v": key, "vp": value.label_font_size}
                    for key, value in NODE_STYLES.items()
                ],
                "type": "integer",
            },
        },
        "NODE_LABEL_MAX_WIDTH": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "type",
                "map": [
                    {"v": key, "vp": value.label_max_width}
                    for key, value in NODE_STYLES.items()
                ],
                "type": "integer",
            },
        },
        "NODE_SHAPE": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "type",
                "map": [
                    {"v": key, "vp": value.shape} for key, value in NODE_STYLES.items()
                ],
                "type": "string",
            },
        },
        "NODE_WIDTH": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "type",
                "map": [
                    {"v": key, "vp": value.width} for key, value in NODE_STYLES.items()
                ],
                "type": "integer",
            },
        },
    },
}

VISUAL_EDITOR_PROPERTIES = {
    "properties": {
        "nodeSizeLocked": False,
        "arrowColorMatchesEdge": False,
        "nodeCustomGraphicsSizeSync": True,
    }
}
