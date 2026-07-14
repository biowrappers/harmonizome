"""Python wrapper for the Ma'ayan Lab Harmonizome API and dataset downloads."""

from .harmonizome import VERSION, Entity, Harmonizome

__version__ = VERSION
__all__ = ["Harmonizome", "Entity", "VERSION"]
