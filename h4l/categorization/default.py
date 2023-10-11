# coding: utf-8

"""
Example selection methods for categories.
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import


ak = maybe_import("awkward")


#
# categorizer functions used by categories definitions
#

@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # fully inclusive selection
    return events, ak.ones_like(events.event) == 1
