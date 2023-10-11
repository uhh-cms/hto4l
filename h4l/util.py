# coding: utf-8

"""
Helpful utils.
"""

from __future__ import annotations

__all__ = ["IF_NANO_V9", "IF_NANO_V10"]

from typing import Any

from columnflow.columnar_util import ArrayFunction, deferred_column


@deferred_column
def IF_NANO_V9(self, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 9 else None


@deferred_column
def IF_NANO_V10(self, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version >= 10 else None
