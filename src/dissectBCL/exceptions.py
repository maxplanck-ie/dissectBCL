"""Exceptions specific to dissectBCL.

FatalError: raised by pure helpers when they detect a condition that
requires the caller to mail the relevant team and stop dissect. The
caller catches it, sends a fatal mail and calls sys.exit(1).
"""


class FatalError(Exception):
    """A condition that stops dissect and requires human intervention."""
