from dataclasses import dataclass
from enum import Enum


class Color(str, Enum):
    LIGHT_BLUE = "#ADD8E6"
    CORNFLOWER_BLUE = "#6495ED"
    MEDIUM_AQUAMARINE = "#66CDAA"
    DARK_SLATE_GRAY = "#2F4F4F"
    RED = "#FF0000"
    GREEN = "#008000"
    HOT_PINK = "#ED6495"
    DARK_SALMON = "#E9967A"
    DARK_GOLDENROD = "#B8860B"
    DARK_SLATE_BLUE = "#483D8B"
    LIGHT_GOLD = "#FFCC66"
    PALE_LAVENDER = "#EBE3F9"
    LIGHT_GRAY = "#DDDDDD"
    PALE_AQUA = "#E0F2F1"


class Width(int, Enum):
    DEFAULT = 4
    SMALL = 2
    INDIRECT = 3
    DIRECT = 4


class ArrowShape(str, Enum):
    CIRCLE = "circle"
    DIAMOND = "diamond"
    TEE = "tee"
    TRIANGLE = "triangle"


class LineStyle(str, Enum):
    DASHED = "dashed"
    SOLID = "solid"


class NodeType(str, Enum):
    GENE = "gene"
    COMPLEX = "complex"
    MOLECULE = "molecule"


@dataclass
class RelationStyle:
    line_style: LineStyle
    arrow_shape: ArrowShape
    label: str
    color: Color
    width: Width


RELATIONS = {
    "BFO:0000050": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="part of",
        color=Color.LIGHT_BLUE,
        width=Width.DEFAULT,
    ),
    "BFO:0000051": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="has part",
        color=Color.CORNFLOWER_BLUE,
        width=Width.DEFAULT,
    ),
    "BFO:0000066": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="occurs in",
        color=Color.MEDIUM_AQUAMARINE,
        width=Width.DEFAULT,
    ),
    "RO:0002211": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.CIRCLE,
        label="regulates",
        color=Color.DARK_SLATE_GRAY,
        width=Width.INDIRECT,
    ),
    "RO:0002212": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.TEE,
        label="negatively regulates",
        color=Color.RED,
        width=Width.INDIRECT,
    ),
    "RO:0002630": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.TEE,
        label="directly negatively regulates",
        color=Color.RED,
        width=Width.DIRECT,
    ),
    "RO:0002213": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.TRIANGLE,
        label="positively regulates",
        color=Color.GREEN,
        width=Width.INDIRECT,
    ),
    "RO:0002629": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.TRIANGLE,
        label="directly positively regulates",
        color=Color.GREEN,
        width=Width.DIRECT,
    ),
    "RO:0002233": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="has input",
        color=Color.CORNFLOWER_BLUE,
        width=Width.DEFAULT,
    ),
    "RO:0002234": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="has output",
        color=Color.HOT_PINK,
        width=Width.DEFAULT,
    ),
    "RO:0002331": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="involved in",
        color=Color.DARK_SALMON,
        width=Width.DEFAULT,
    ),
    "RO:0002333": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="enabled by",
        color=Color.DARK_GOLDENROD,
        width=Width.DEFAULT,
    ),
    "RO:0002411": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.CIRCLE,
        label="causally upstream of",
        color=Color.DARK_SLATE_BLUE,
        width=Width.INDIRECT,
    ),
    "RO:0002418": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.CIRCLE,
        label="causally upstream of or within",
        color=Color.DARK_SLATE_BLUE,
        width=Width.INDIRECT,
    ),
    "RO:0002408": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.TEE,
        label="directly inhibits",
        color=Color.RED,
        width=Width.DIRECT,
    ),
    "RO:0002406": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.TRIANGLE,
        label="directly activates",
        color=Color.GREEN,
        width=Width.DIRECT,
    ),
    "RO:0002305": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.CIRCLE,
        label="causally upstream of, negative effect",
        color=Color.RED,
        width=Width.INDIRECT,
    ),
    "RO:0004046": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.CIRCLE,
        label="causally upstream of or within, negative effect",
        color=Color.RED,
        width=Width.INDIRECT,
    ),
    "RO:0002304": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.CIRCLE,
        label="causally upstream of, positive effect",
        color=Color.GREEN,
        width=Width.INDIRECT,
    ),
    "RO:0004047": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.CIRCLE,
        label="causally upstream of or within, positive effect",
        color=Color.GREEN,
        width=Width.INDIRECT,
    ),
    "annotation": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.DIAMOND,
        label="annotation",
        color=Color.DARK_SLATE_BLUE,
        width=Width.DEFAULT,
    ),
    "instance_of": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="activity",
        color=Color.LIGHT_GOLD,
        width=Width.DEFAULT,
    ),
    "RO:0002413": RelationStyle(
        line_style=LineStyle.SOLID,
        arrow_shape=ArrowShape.CIRCLE,
        label="directly provides input for",
        color=Color.LIGHT_BLUE,
        width=Width.SMALL,
    ),
    "RO:0002407": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.TRIANGLE,
        label="indirectly positively regulates",
        color=Color.GREEN,
        width=Width.INDIRECT,
    ),
    "RO:0002409": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.TEE,
        label="indirectly negatively regulates",
        color=Color.RED,
        width=Width.INDIRECT,
    ),
    "RO:0012009": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.CIRCLE,
        label="constitutively upstream of",
        color=Color.DARK_SLATE_BLUE,
        width=Width.INDIRECT,
    ),
    "RO:0012010": RelationStyle(
        line_style=LineStyle.DASHED,
        arrow_shape=ArrowShape.CIRCLE,
        label="removes input for",
        color=Color.DARK_SLATE_BLUE,
        width=Width.INDIRECT,
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
                    {"v": NodeType.GENE, "vp": Color.PALE_LAVENDER},
                    {"v": NodeType.COMPLEX, "vp": Color.LIGHT_GRAY},
                    {"v": NodeType.MOLECULE, "vp": Color.PALE_AQUA},
                ],
                "type": "string",
            },
        },
        "NODE_HEIGHT": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "type",
                "map": [
                    {"v": NodeType.GENE, "vp": 40},
                    {"v": NodeType.COMPLEX, "vp": 70},
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
                    {"v": NodeType.COMPLEX, "vp": 15},
                    {"v": NodeType.MOLECULE, "vp": 10},
                ],
                "type": "integer",
            },
        },
        "NODE_LABEL_MAX_WIDTH": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "type",
                "map": [
                    {"v": NodeType.GENE, "vp": 80},
                    {"v": NodeType.COMPLEX, "vp": 100},
                    {"v": NodeType.MOLECULE, "vp": 70},
                ],
                "type": "integer",
            },
        },
        "NODE_SHAPE": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "type",
                "map": [
                    {"v": NodeType.GENE, "vp": "ellipse"},
                    {"v": NodeType.COMPLEX, "vp": "octagon"},
                ],
                "type": "string",
            },
        },
        "NODE_WIDTH": {
            "type": "DISCRETE",
            "definition": {
                "attribute": "type",
                "map": [
                    {"v": NodeType.GENE, "vp": 85},
                    {"v": NodeType.COMPLEX, "vp": 110},
                    {"v": NodeType.MOLECULE, "vp": 70},
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
