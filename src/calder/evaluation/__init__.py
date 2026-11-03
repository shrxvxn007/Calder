"""Offline evaluation. Hidden labels are permitted here, never in acquisition."""

from calder.evaluation.reports import build_context_comparison, build_report

__all__ = ["build_report", "build_context_comparison"]
