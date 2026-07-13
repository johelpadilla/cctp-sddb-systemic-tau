"""
Native ordinal alarm detectors (methodological alternatives to abs-z).

Option 1 (OPC) and Option 2 (SDD) are independent pure detectors.
Exploratory cascade fusion (SDD → OPC confirm) is a separate post-process
module; it does not alter opc_detect / sdd_detect alarm math.
"""

from .opc_detector import opc_detect, opc_alarm_at
from .sdd_detector import sdd_detect, total_variation, sdd_alarm_at
from .cascade_fusion import (
    CONFIRM_WINDOW_H,
    CONFIRM_WINDOW_MIN,
    cascade_sdd_confirm_opc,
    cascade_first_alarm_index,
    cascade_first_causal_detection,
)

__all__ = [
    "opc_detect",
    "opc_alarm_at",
    "sdd_detect",
    "sdd_alarm_at",
    "total_variation",
    "CONFIRM_WINDOW_H",
    "CONFIRM_WINDOW_MIN",
    "cascade_sdd_confirm_opc",
    "cascade_first_alarm_index",
    "cascade_first_causal_detection",
]
