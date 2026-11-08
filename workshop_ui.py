"""Redirect old workshop shortcuts to the practical docking workspace.

Saved workshop projects and notes remain on disk. The interface no longer
requires quiz answers, predictions or discussion checkboxes.
"""
from nicegui import ui


def register_workshop() -> None:
    @ui.page("/workshop")
    def workshop_redirect():
        ui.navigate.to("/study")
