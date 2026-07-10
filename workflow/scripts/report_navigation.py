#!/usr/bin/env python3
"""Shared static report navigation helpers for ASPIS HTML reports."""

from __future__ import annotations

import html
import os
from pathlib import Path
from typing import Any


ReportMapItem = dict[str, Any]


def report_map_item(
    label: str,
    target: str | Path = "",
    *,
    children: list[ReportMapItem] | None = None,
    planned: bool = False,
) -> ReportMapItem:
    return {
        "label": label,
        "target": target,
        "children": children or [],
        "planned": planned,
    }


def report_map_css() -> str:
    return """
    .report-shell { align-items: start; display: grid; gap: 24px; grid-template-columns: 260px minmax(0, 1fr); }
    .report-content { min-width: 0; }
    .report-map { background: #fff; border: 1px solid #d0d7de; border-radius: 6px; max-height: calc(100vh - 32px); overflow: auto; padding: 0.85rem; position: sticky; top: 16px; }
    .report-map-title { border-bottom: 1px solid #d0d7de; color: #24292f; font-size: 1rem; font-weight: 700; margin: 0 0 0.65rem; padding-bottom: 0.45rem; }
    .report-map ul { list-style: none; margin: 0; padding-left: 0; }
    .report-map ul ul { border-left: 1px solid #d0d7de; margin: 0.35rem 0 0.45rem 0.45rem; padding-left: 0.75rem; }
    .report-map li { margin: 0.35rem 0; }
    .report-map a { border-radius: 4px; color: #0969da; display: block; padding: 0.15rem 0.35rem; text-decoration: none; }
    .report-map a:hover { background: #f6f8fa; text-decoration: none; }
    .report-map a.nav-active { background: #0969da; color: #fff; font-weight: 700; }
    .report-map a.nav-active:hover { background: #0969da; color: #fff; }
    .report-map .nav-missing { color: #57606a; }
    .report-map .nav-current { color: #24292f; font-weight: 700; }
    @media (max-width: 1050px) {
      .report-shell { display: block; }
      .report-map { max-height: none; margin-bottom: 1rem; position: static; }
    }
    """


def report_map_script() -> str:
    return """
  <script>
    (function () {
      const links = Array.from(document.querySelectorAll('.report-map a[data-report-nav-target]'));
      if (!links.length) {
        return;
      }
      const sections = links
        .map(link => document.getElementById(link.dataset.reportNavTarget))
        .filter(Boolean);
      function setActive(sectionId) {
        links.forEach(link => {
          link.classList.toggle('nav-active', link.dataset.reportNavTarget === sectionId);
        });
      }
      function updateActive() {
        let current = sections.length ? sections[0].id : '';
        sections.forEach(section => {
          const rect = section.getBoundingClientRect();
          if (rect.top <= 140) {
            current = section.id;
          }
        });
        if (window.location.hash) {
          const hashed = window.location.hash.slice(1);
          if (document.getElementById(hashed)) {
            const rect = document.getElementById(hashed).getBoundingClientRect();
            if (rect.top >= -220 && rect.top <= window.innerHeight * 0.6) {
              current = hashed;
            }
          }
        }
        if (current) {
          setActive(current);
        }
      }
      links.forEach(link => {
        link.addEventListener('click', () => setActive(link.dataset.reportNavTarget));
      });
      window.addEventListener('scroll', updateActive, { passive: true });
      window.addEventListener('resize', updateActive);
      window.addEventListener('hashchange', updateActive);
      updateActive();
    }());
  </script>
    """


def _rel_href(path: Path, base_dir: Path) -> str:
    if path.is_absolute():
        return os.path.relpath(path, start=base_dir).replace(os.sep, "/")
    return path.as_posix()


def _relative_target_href(path: Path, base_dir: Path) -> tuple[str, bool]:
    base_relative = base_dir / path
    if base_relative.exists():
        return path.as_posix(), True
    if path.exists():
        return os.path.relpath(path, start=base_dir).replace(os.sep, "/"), True
    return path.as_posix(), False


def _target_href(target: str | Path, base_dir: Path) -> tuple[str, bool, bool]:
    if not target:
        return "", False, False
    if isinstance(target, Path):
        if target.is_absolute():
            return _rel_href(target, base_dir), target.exists(), True
        href, exists = _relative_target_href(target, base_dir)
        return href, exists, True
    target_text = str(target)
    if target_text.startswith("#") or "://" in target_text:
        return target_text, True, False
    path = Path(target_text)
    if path.is_absolute():
        return _rel_href(path, base_dir), path.exists(), True
    href, exists = _relative_target_href(path, base_dir)
    return href, exists, True


def _render_item(item: ReportMapItem, base_dir: Path) -> str:
    label = html.escape(str(item.get("label", "")))
    href, exists, filesystem_target = _target_href(item.get("target", ""), base_dir)
    planned = bool(item.get("planned", False))
    children = item.get("children") or []
    if href and (exists or planned or not filesystem_target):
        data_attr = ""
        if href.startswith("#") and len(href) > 1:
            data_attr = f' data-report-nav-target="{html.escape(href[1:])}"'
        label_html = f'<a href="{html.escape(href)}"{data_attr}>{label}</a>'
    elif href:
        label_html = f'<span class="nav-missing">{label}</span>'
    else:
        label_html = f'<span class="nav-current">{label}</span>'
    child_html = ""
    if children:
        child_html = "<ul>" + "".join(_render_item(child, base_dir) for child in children) + "</ul>"
    return f"<li>{label_html}{child_html}</li>"


def report_map_sidebar(title: str, items: list[ReportMapItem], base_dir: Path) -> str:
    return (
        '<aside class="report-map" aria-label="Report map">'
        f'<div class="report-map-title">{html.escape(title)}</div>'
        "<ul>"
        + "".join(_render_item(item, base_dir) for item in items)
        + "</ul></aside>"
    )


def report_shell_open(title: str, items: list[ReportMapItem], base_dir: Path) -> str:
    return f'<div class="report-shell">{report_map_sidebar(title, items, base_dir)}<main class="report-content">'


def report_shell_close() -> str:
    return "</main></div>"
