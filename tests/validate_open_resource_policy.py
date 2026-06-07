#!/usr/bin/env python3
"""Validate ASPIS open-resource source policy."""

from __future__ import annotations

from pathlib import Path

import yaml


POLICY = Path("config/aspis_open_resource_sources.example.yaml")
CONTROLLED_LICENSE_STATUS = {"open", "user_provided", "restricted", "unknown"}
RESTRICTED_SOURCES = {"KEGG", "MSigDB", "SignalP", "TMHMM", "DeepTMHMM"}
REQUIRED_KEYS = {
    "resource_id",
    "source",
    "resource_kind",
    "license_status",
    "recommended_for_default_bundle",
    "auto_download",
    "notes",
}


def flatten_open_sources(data: dict) -> list[dict]:
    rows: list[dict] = []
    for section, entries in data.get("open_resource_sources", {}).items():
        if not isinstance(entries, list) or not entries:
            raise AssertionError(f"open_resource_sources.{section} must be a non-empty list")
        for entry in entries:
            if not isinstance(entry, dict):
                raise AssertionError(f"open_resource_sources.{section} contains a non-object entry")
            row = dict(entry)
            row["section"] = section
            rows.append(row)
    return rows


def validate_entry(entry: dict, *, restricted_section: bool) -> None:
    missing = REQUIRED_KEYS - set(entry)
    if missing:
        raise AssertionError(f"{entry.get('resource_id', '<unknown>')} is missing keys: {sorted(missing)}")
    resource_id = str(entry["resource_id"])
    source = str(entry["source"])
    license_status = str(entry["license_status"])
    if license_status not in CONTROLLED_LICENSE_STATUS:
        raise AssertionError(f"{resource_id} has uncontrolled license_status: {license_status!r}")
    if entry["auto_download"] is not False:
        raise AssertionError(f"{resource_id} must not enable automatic downloads in the policy example")
    if not str(entry.get("notes", "")).strip():
        raise AssertionError(f"{resource_id} needs a note explaining use constraints")
    if source in RESTRICTED_SOURCES and license_status != "restricted":
        raise AssertionError(f"{resource_id} is a restricted source but is labeled {license_status!r}")
    if restricted_section:
        if license_status != "restricted":
            raise AssertionError(f"manual-only resource {resource_id} must have license_status restricted")
        if entry["recommended_for_default_bundle"] is not False:
            raise AssertionError(f"manual-only resource {resource_id} cannot be recommended for the default bundle")
    if license_status == "restricted" and entry["recommended_for_default_bundle"] is not False:
        raise AssertionError(f"restricted resource {resource_id} cannot be recommended for the default bundle")


def main() -> int:
    data = yaml.safe_load(POLICY.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise AssertionError(f"{POLICY} must contain a YAML mapping")
    open_entries = flatten_open_sources(data)
    if not any(entry["source"] == "GO" and entry["recommended_for_default_bundle"] for entry in open_entries):
        raise AssertionError("policy must recommend GO for the default open RNA-seq bundle")
    if not any(entry["source"] == "Reactome" and entry["recommended_for_default_bundle"] for entry in open_entries):
        raise AssertionError("policy must recommend Reactome for the default open RNA-seq bundle")
    for entry in open_entries:
        validate_entry(entry, restricted_section=False)
        if entry["source"] in RESTRICTED_SOURCES:
            raise AssertionError(f"restricted source {entry['source']} must not appear under open_resource_sources")
    restricted_entries = data.get("restricted_or_manual_only", [])
    if not isinstance(restricted_entries, list) or not restricted_entries:
        raise AssertionError("restricted_or_manual_only must be a non-empty list")
    observed_restricted = {str(entry.get("source", "")) for entry in restricted_entries if isinstance(entry, dict)}
    missing_restricted = RESTRICTED_SOURCES - observed_restricted
    if missing_restricted:
        raise AssertionError(f"manual-only policy misses restricted sources: {sorted(missing_restricted)}")
    for entry in restricted_entries:
        if not isinstance(entry, dict):
            raise AssertionError("restricted_or_manual_only contains a non-object entry")
        validate_entry(entry, restricted_section=True)
    print("open resource source policy ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())