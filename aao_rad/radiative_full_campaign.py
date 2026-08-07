#!/usr/bin/env python3
"""Plan and execute full-catalog radiative mode-4 campaigns.

The per-stratum mode-4 implementation remains the authority for preparing,
running, validating, and finalizing AAO jobs.  This module adds the missing
campaign layer: it turns an audited cumulative stratum queue into immutable
mode-4 manifests, scheduler-sized tasks, resumable receipts, pooled envelope
reports, automatically sized follow-up campaigns, and production manifests.

No scheduler command is executed by this program.  ``emit-swif`` writes an
auditable submission script which the user may inspect before submitting it.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import shlex
import socket
import stat
import traceback
from pathlib import Path
from typing import Iterable, Optional

import radiative_mode4


CAMPAIGN_SCHEMA = "aao-rad-full-campaign-v1"
QUEUE_SCHEMA = "aao-rad-cumulative-stratum-queue-v1"
CALIBRATION_SCHEMA = radiative_mode4.CALIBRATION_SCHEMA
RUN_SCHEMA = radiative_mode4.RUN_SCHEMA
READY_STATES = {"ready", "ready_provisional_zero_complement"}
SELECTABLE_WORK_CATEGORIES = {
    "supported_calibration",
    "guard_refinement",
    "targeted_discovery",
}


class CampaignError(RuntimeError):
    """Raised when a full-campaign artifact violates an invariant."""


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _load_json(path: Path) -> dict:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise CampaignError(f"{path}: invalid JSON: {error}") from error
    if not isinstance(payload, dict):
        raise CampaignError(f"{path}: expected a JSON object")
    return payload


def _source_revision() -> str:
    return radiative_mode4._source_revision()


def _require_file(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    return resolved


def _queue_records(
    queue_path: Path,
    *,
    selection: str,
    work_categories: Iterable[str],
) -> tuple[dict, list[dict[str, object]]]:
    queue = _load_json(queue_path)
    if queue.get("schema") != QUEUE_SCHEMA:
        raise CampaignError(f"{queue_path}: expected schema {QUEUE_SCHEMA}")
    raw = queue.get("strata")
    if not isinstance(raw, list) or not raw:
        raise CampaignError(f"{queue_path}: queue strata are empty")
    categories = set(work_categories)
    unknown = categories - SELECTABLE_WORK_CATEGORIES
    if unknown:
        raise ValueError(
            "unknown work categories: " + ", ".join(sorted(unknown))
        )
    if selection == "selected":
        records = [
            item
            for item in raw
            if item.get("selected") is True
            and str(item.get("work_category")) in categories
        ]
    elif selection == "data-occupied":
        records = [
            item
            for item in raw
            if int(item.get("data_events") or 0) > 0
            and str(item.get("work_category")) in categories
        ]
    elif selection == "all-catalog":
        records = list(raw)
    else:
        raise ValueError(f"unsupported queue selection {selection!r}")
    if not records:
        raise CampaignError("queue selection produced no strata")
    flat_indices = [int(item["flat_index"]) for item in records]
    stratum_ids = [str(item["stratum_id"]) for item in records]
    if len(set(flat_indices)) != len(flat_indices):
        raise CampaignError("queue selection contains duplicate flat indices")
    if len(set(stratum_ids)) != len(stratum_ids):
        raise CampaignError("queue selection contains duplicate strata")
    return queue, sorted(records, key=lambda item: int(item["flat_index"]))


def _verify_queue_config(queue: dict, config_path: Path) -> None:
    expected = queue.get("analysis_config_sha256")
    if expected is not None and expected != _sha256(config_path):
        raise CampaignError(
            "queue analysis configuration hash differs from --config"
        )


def _write_flat_indices(path: Path, records: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f"{int(item['flat_index'])}\n" for item in records),
        encoding="utf-8",
    )


def _prepare_namespace(
    *,
    command: str,
    config: Path,
    recipes: Path,
    legacy_input: Path,
    output: Path,
    candidate: str,
    core_fraction: float,
    replicas: int,
    seed_base: int,
    flat_index_file: Path,
    apply_y_max: bool,
    heartbeat_interval: int,
    generator_revision: str,
    trials: Optional[int] = None,
    inside_fraction: Optional[float] = None,
    envelope_report: Optional[Path] = None,
    events_per_stratum: Optional[int] = None,
    refinements: Optional[Path] = None,
) -> argparse.Namespace:
    return argparse.Namespace(
        command=command,
        config=config,
        recipes=recipes,
        refinements=refinements,
        input=legacy_input,
        output=output,
        candidate=candidate,
        core_fraction=core_fraction,
        replicas=replicas,
        seed_base=seed_base,
        bin_start=0,
        bin_stop=None,
        flat_indices=None,
        flat_index_file=flat_index_file,
        apply_y_max=apply_y_max,
        heartbeat_interval=heartbeat_interval,
        generator_revision=generator_revision,
        overwrite=False,
        trials=trials,
        calibration_inside_guard_fraction=inside_fraction,
        sigr_max=None,
        envelope_report=envelope_report,
        events_per_stratum=events_per_stratum,
        allow_envelope_revision_mismatch=False,
        envelope_revision_compatibility_rationale=None,
    )


def _manifest_tasks(
    manifest_path: Path,
    *,
    stage: str,
    first_task_id: int = 1,
) -> list[dict[str, object]]:
    manifest = _load_json(manifest_path)
    tasks: list[dict[str, object]] = []
    for offset, record in enumerate(manifest.get("runs") or []):
        output_stem = manifest_path.parent / str(record["output_stem"])
        tasks.append(
            {
                "task_id": first_task_id + offset,
                "stage": stage,
                "operation": manifest["operation"],
                "manifest": str(manifest_path),
                "manifest_sha256": _sha256(manifest_path),
                "flat_index": int(record["flat_index"]),
                "stratum_id": str(record["stratum_id"]),
                "replica_index": int(record["replica_index"]),
                "seed": int(record["seed"]),
                "requested_trials": int(record.get("trials_requested") or 0),
                "requested_events": int(record.get("events_requested") or 0),
                "expected_run": str(Path(str(output_stem) + ".json")),
            }
        )
    if not tasks:
        raise CampaignError(f"{manifest_path}: manifest has no runs")
    return tasks


TASK_FIELDS = (
    "task_id",
    "stage",
    "operation",
    "manifest",
    "manifest_sha256",
    "flat_index",
    "stratum_id",
    "replica_index",
    "seed",
    "requested_trials",
    "requested_events",
    "expected_run",
)


def _write_tasks(root: Path, tasks: list[dict[str, object]]) -> None:
    with (root / "tasks.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.DictWriter(
            destination, fieldnames=TASK_FIELDS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(tasks)
    (root / "task_ids.txt").write_text(
        "".join(f"{int(task['task_id'])}\n" for task in tasks),
        encoding="utf-8",
    )


def _campaign_payload(
    *,
    kind: str,
    root: Path,
    tasks: list[dict[str, object]],
    manifests: list[Path],
    pool_manifests: list[Path],
    frozen_inputs: dict[str, object],
    selection: dict[str, object],
    parent: Optional[Path] = None,
    calibration_report: Optional[Path] = None,
    followup: Optional[dict[str, object]] = None,
) -> dict[str, object]:
    return {
        "schema": CAMPAIGN_SCHEMA,
        "created_utc": _now(),
        "driver_revision": _source_revision(),
        "driver_source_sha256": _sha256(Path(__file__).resolve()),
        "kind": kind,
        "root": str(root),
        "parent_campaign": str(parent) if parent is not None else None,
        "parent_campaign_sha256": (
            _sha256(parent) if parent is not None else None
        ),
        "calibration_report": (
            str(calibration_report)
            if calibration_report is not None
            else None
        ),
        "calibration_report_sha256": (
            _sha256(calibration_report)
            if calibration_report is not None
            else None
        ),
        "frozen_inputs": frozen_inputs,
        "selection": selection,
        "manifests": [
            {"path": str(path), "sha256": _sha256(path)}
            for path in manifests
        ],
        "pool_manifests": [
            {"path": str(path), "sha256": _sha256(path)}
            for path in pool_manifests
        ],
        "task_count": len(tasks),
        "tasks_tsv": str(root / "tasks.tsv"),
        "tasks_tsv_sha256": None,
        "task_ids": str(root / "task_ids.txt"),
        "task_ids_sha256": None,
        "followup": followup,
    }


def _finish_campaign(root: Path, payload: dict, tasks: list[dict]) -> Path:
    _write_tasks(root, tasks)
    payload["tasks_tsv_sha256"] = _sha256(root / "tasks.tsv")
    payload["task_ids_sha256"] = _sha256(root / "task_ids.txt")
    path = root / "campaign.json"
    _write_json(path, payload)
    return path


def plan_calibration(args: argparse.Namespace) -> Path:
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(output)
    queue_path = _require_file(args.queue)
    config_path = _require_file(args.config)
    recipes_path = _require_file(args.recipes)
    input_path = _require_file(args.input)
    refinements_path = (
        _require_file(args.refinements) if args.refinements else None
    )
    categories = args.work_category or sorted(SELECTABLE_WORK_CATEGORIES)
    queue, records = _queue_records(
        queue_path,
        selection=args.selection,
        work_categories=categories,
    )
    _verify_queue_config(queue, config_path)
    output.mkdir(parents=True)
    selection_path = output / "selection" / "flat_indices.txt"
    _write_flat_indices(selection_path, records)
    stage = output / "stages" / "calibration_initial"
    manifest_path = radiative_mode4.prepare(
        _prepare_namespace(
            command="prepare-calibration",
            config=config_path,
            recipes=recipes_path,
            legacy_input=input_path,
            refinements=refinements_path,
            output=stage,
            candidate=args.candidate,
            core_fraction=args.core_fraction,
            replicas=args.replicas,
            seed_base=args.seed_base,
            flat_index_file=selection_path,
            apply_y_max=args.apply_y_max,
            heartbeat_interval=args.heartbeat_interval,
            generator_revision=args.generator_revision,
            trials=args.trials,
            inside_fraction=args.inside_guard_trial_fraction,
        )
    )
    manifest_path = manifest_path.resolve()
    tasks = _manifest_tasks(manifest_path, stage="calibration_initial")
    category_counts = {
        category: sum(
            str(record.get("work_category")) == category
            for record in records
        )
        for category in sorted(SELECTABLE_WORK_CATEGORIES)
    }
    payload = _campaign_payload(
        kind="calibration",
        root=output,
        tasks=tasks,
        manifests=[manifest_path],
        pool_manifests=[manifest_path],
        frozen_inputs={
            "queue": str(queue_path),
            "queue_sha256": _sha256(queue_path),
            "config": str(stage / "analysis_config.json"),
            "config_sha256": _sha256(stage / "analysis_config.json"),
            "recipes": str(stage / "continuous_guard_recipes.json"),
            "recipes_sha256": _sha256(
                stage / "continuous_guard_recipes.json"
            ),
            "legacy_input": str(stage / "legacy_input.inp"),
            "legacy_input_sha256": _sha256(stage / "legacy_input.inp"),
            "refinements": (
                str(stage / "guard_refinements.json")
                if refinements_path is not None
                else None
            ),
            "refinements_sha256": (
                _sha256(stage / "guard_refinements.json")
                if refinements_path is not None
                else None
            ),
        },
        selection={
            "policy": args.selection,
            "work_categories": categories,
            "category_counts": category_counts,
            "selected_strata": len(records),
            "flat_indices": str(selection_path),
            "flat_indices_sha256": _sha256(selection_path),
            "queue_summary": queue.get("summary"),
        },
    )
    return _finish_campaign(output, payload, tasks)


def _load_campaign(path: Path) -> tuple[Path, dict]:
    campaign_path = _require_file(path)
    payload = _load_json(campaign_path)
    if payload.get("schema") != CAMPAIGN_SCHEMA:
        raise CampaignError(
            f"{campaign_path}: expected schema {CAMPAIGN_SCHEMA}"
        )
    root = Path(str(payload["root"])).resolve()
    if campaign_path != root / "campaign.json":
        raise CampaignError("campaign path differs from its frozen root")
    for key in ("tasks_tsv", "task_ids"):
        artifact = _require_file(Path(str(payload[key])))
        if _sha256(artifact) != payload[f"{key}_sha256"]:
            raise CampaignError(f"campaign {key} artifact changed")
    for item in payload.get("manifests") or []:
        manifest = _require_file(Path(str(item["path"])))
        if _sha256(manifest) != item["sha256"]:
            raise CampaignError(f"{manifest}: manifest changed")
    for item in payload.get("pool_manifests") or []:
        manifest = _require_file(Path(str(item["path"])))
        if _sha256(manifest) != item["sha256"]:
            raise CampaignError(f"{manifest}: pooled manifest changed")
    parent = payload.get("parent_campaign")
    if parent is not None:
        parent_path = _require_file(Path(str(parent)))
        if _sha256(parent_path) != payload.get("parent_campaign_sha256"):
            raise CampaignError("parent campaign changed")
    calibration = payload.get("calibration_report")
    if calibration is not None:
        calibration_path = _require_file(Path(str(calibration)))
        if _sha256(calibration_path) != payload.get(
            "calibration_report_sha256"
        ):
            raise CampaignError("campaign calibration report changed")
    frozen = payload.get("frozen_inputs") or {}
    for name in ("queue", "config", "recipes", "legacy_input", "refinements"):
        raw_path = frozen.get(name)
        if raw_path is None:
            continue
        artifact = _require_file(Path(str(raw_path)))
        if _sha256(artifact) != frozen.get(f"{name}_sha256"):
            raise CampaignError(f"frozen {name} artifact changed")
    return campaign_path, payload


def _read_tasks(campaign: dict) -> list[dict[str, str]]:
    path = Path(str(campaign["tasks_tsv"]))
    with path.open(encoding="utf-8", newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    if len(rows) != int(campaign["task_count"]):
        raise CampaignError("task table length differs from campaign")
    return rows


def _task_run_valid(task: dict[str, str]) -> bool:
    path = Path(task["expected_run"])
    if not path.is_file():
        return False
    run = _load_json(path)
    return (
        run.get("schema") == RUN_SCHEMA
        and int(run.get("flat_index", -1)) == int(task["flat_index"])
        and int(run.get("replica_index", -1))
        == int(task["replica_index"])
        and str(run.get("source_manifest_sha256"))
        == task["manifest_sha256"]
    )


def run_task(args: argparse.Namespace) -> Path:
    _campaign_path, campaign = _load_campaign(args.campaign)
    matches = [
        task
        for task in _read_tasks(campaign)
        if int(task["task_id"]) == args.task_id
    ]
    if len(matches) != 1:
        raise CampaignError(
            f"campaign has {len(matches)} tasks numbered {args.task_id}"
        )
    task = matches[0]
    run_path = Path(task["expected_run"])
    receipt = Path(str(campaign["root"])) / "receipts" / (
        f"task_{args.task_id:06d}.json"
    )
    if _task_run_valid(task) and not args.overwrite:
        _write_json(
            receipt,
            {
                "schema": "aao-rad-full-campaign-task-receipt-v1",
                "task_id": args.task_id,
                "status": "already_complete",
                "checked_utc": _now(),
                "run": str(run_path),
                "run_sha256": _sha256(run_path),
                "host": socket.gethostname(),
            },
        )
        return run_path
    started = _now()
    try:
        result = radiative_mode4.run(
            argparse.Namespace(
                manifest=Path(task["manifest"]),
                flat_index=int(task["flat_index"]),
                replica_index=int(task["replica_index"]),
                executable=args.executable,
                overwrite=args.overwrite,
            )
        )
        _write_json(
            receipt,
            {
                "schema": "aao-rad-full-campaign-task-receipt-v1",
                "task_id": args.task_id,
                "status": "complete",
                "started_utc": started,
                "completed_utc": _now(),
                "run": str(result),
                "run_sha256": _sha256(result),
                "host": socket.gethostname(),
            },
        )
        return result
    except Exception as error:
        _write_json(
            receipt,
            {
                "schema": "aao-rad-full-campaign-task-receipt-v1",
                "task_id": args.task_id,
                "status": "failed",
                "started_utc": started,
                "failed_utc": _now(),
                "error_type": type(error).__name__,
                "error": str(error),
                "traceback": traceback.format_exc(),
                "host": socket.gethostname(),
            },
        )
        raise


def status(args: argparse.Namespace) -> Path:
    _campaign_path, campaign = _load_campaign(args.campaign)
    root = Path(str(campaign["root"]))
    records: list[dict[str, object]] = []
    counts = {"complete": 0, "failed": 0, "pending": 0}
    for task in _read_tasks(campaign):
        task_id = int(task["task_id"])
        receipt = root / "receipts" / f"task_{task_id:06d}.json"
        if _task_run_valid(task):
            state = "complete"
        elif receipt.is_file() and _load_json(receipt).get("status") == "failed":
            state = "failed"
        else:
            state = "pending"
        counts[state] += 1
        records.append({**task, "status": state, "receipt": str(receipt)})
    payload = {
        "schema": "aao-rad-full-campaign-status-v1",
        "created_utc": _now(),
        "campaign": str(args.campaign.expanduser().resolve()),
        "campaign_kind": campaign["kind"],
        "task_count": len(records),
        "status_counts": counts,
        "complete": counts["complete"] == len(records),
        "tasks": records,
    }
    output = (
        args.output.expanduser().resolve()
        if args.output is not None
        else root / "status.json"
    )
    _write_json(output, payload)
    fields = (*TASK_FIELDS, "status", "receipt")
    with output.with_suffix(".tsv").open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.DictWriter(
            destination, fieldnames=fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(records)
    return output


def _assert_complete(campaign_path: Path, campaign: dict) -> None:
    missing = [
        int(task["task_id"])
        for task in _read_tasks(campaign)
        if not _task_run_valid(task)
    ]
    if missing:
        preview = ", ".join(str(value) for value in missing[:10])
        suffix = "..." if len(missing) > 10 else ""
        raise CampaignError(
            f"{campaign_path}: {len(missing)} tasks are incomplete: "
            f"{preview}{suffix}"
        )


def finalize(args: argparse.Namespace) -> Path:
    campaign_path, campaign = _load_campaign(args.campaign)
    _assert_complete(campaign_path, campaign)
    manifests = [
        Path(str(item["path"])) for item in campaign["pool_manifests"]
    ]
    output = (
        args.output.expanduser().resolve()
        if args.output is not None
        else Path(str(campaign["root"]))
        / (
            "envelope_calibration.json"
            if campaign["kind"] in ("calibration", "calibration_followup")
            else "campaign_weights.json"
        )
    )
    result = radiative_mode4.finalize(
        argparse.Namespace(
            manifests=manifests,
            output=output,
            envelope_safety_factor=args.envelope_safety_factor,
            maximum_duplicate_fraction=args.maximum_duplicate_fraction,
            minimum_component_targets=args.minimum_component_targets,
            minimum_provisional_inside_targets=(
                args.minimum_provisional_inside_targets
            ),
            additional_inside_pilot_run=None,
            allow_zero_complement=args.allow_zero_complement,
            zero_complement_confidence=args.zero_complement_confidence,
            maximum_zero_complement_target_rate=(
                args.maximum_zero_complement_target_rate
            ),
            allow_calibration_revision_mismatch=False,
            revision_compatibility_rationale=None,
        )
    )
    return result


def _zero_required_trials(confidence: float, maximum_rate: float) -> int:
    return int(
        math.ceil(math.log(1.0 - confidence) / math.log(1.0 - maximum_rate))
    )


def _round_trial_bucket(required: float, quantum: int, maximum: int) -> int:
    if required <= quantum:
        return quantum
    ratio = required / quantum
    bucket = quantum * (2 ** int(math.ceil(math.log(ratio, 2.0))))
    return min(maximum, int(bucket))


def _support_followup(
    stratum: dict,
    *,
    minimum_targets: int,
    minimum_inside_targets: int,
    confidence: float,
    maximum_zero_rate: float,
    safety: float,
    quantum: int,
    maximum_trials: int,
    discovery_trials: int,
) -> Optional[dict[str, object]]:
    if str(stratum.get("pilot_readiness")) in READY_STATES:
        return None
    inside = stratum["inside_guard"]
    complement = stratum["guard_complement"]
    inside_targets = int(inside["target_candidates"])
    inside_trials = int(inside["trials"])
    complement_targets = int(complement["target_candidates"])
    complement_trials = int(complement["trials"])
    status_name = str(stratum["recommendation_status"])

    if status_name == "no_candidate_meets_duplicate_limit":
        return {
            "kind": "manual_envelope_review",
            "schedulable": False,
            "reason": status_name,
        }

    inside_missing = max(0, minimum_inside_targets - inside_targets)
    complement_missing = max(0, minimum_targets - complement_targets)
    if inside_targets == 0 and complement_targets == 0:
        return {
            "kind": "discovery",
            "schedulable": True,
            "inside_fraction": 0.5,
            "trials": discovery_trials,
            "reason": status_name,
        }

    required_options: list[tuple[str, float, float]] = []
    if inside_missing > 0:
        if inside_targets > 0 and inside_trials > 0:
            required_inside = safety * inside_missing / (
                inside_targets / inside_trials
            )
        else:
            required_inside = discovery_trials * 0.5
        required_options.append(("inside_support", required_inside, 0.9))

    if complement_missing > 0 and complement_targets > 0:
        required_complement = safety * complement_missing / (
            complement_targets / complement_trials
        )
        required_options.append(
            ("complement_support", required_complement, 0.1)
        )
    elif complement_targets == 0:
        total_required = _zero_required_trials(confidence, maximum_zero_rate)
        additional = max(0, total_required - complement_trials)
        if additional > 0:
            required_options.append(
                ("zero_complement_exposure", float(additional), 0.1)
            )

    if not required_options:
        return {
            "kind": "manual_support_review",
            "schedulable": False,
            "reason": status_name,
        }
    # A single mixed follow-up handles simultaneous deficiencies.  Otherwise
    # use the component-focused allocation and translate component proposals
    # into total proposals.
    kinds = {item[0] for item in required_options}
    if len(required_options) > 1:
        beta = 0.5
        total_required = max(
            required / (beta if "inside" in kind else (1.0 - beta))
            for kind, required, _fraction in required_options
        )
        kind = "both_component_support"
    else:
        kind, component_required, beta = required_options[0]
        component_fraction = beta if "inside" in kind else (1.0 - beta)
        total_required = component_required / component_fraction
    capped = total_required > maximum_trials
    trials = _round_trial_bucket(total_required, quantum, maximum_trials)
    return {
        "kind": kind,
        "schedulable": True,
        "inside_fraction": beta,
        "trials": trials,
        "required_trials_estimate": total_required,
        "capped_at_maximum": capped,
        "reason": status_name,
    }


def plan_followup(args: argparse.Namespace) -> Path:
    parent_path, parent = _load_campaign(args.campaign)
    if parent["kind"] not in ("calibration", "calibration_followup"):
        raise CampaignError("follow-ups require a calibration campaign")
    calibration_path = _require_file(args.calibration)
    calibration = _load_json(calibration_path)
    if calibration.get("schema") != CALIBRATION_SCHEMA:
        raise CampaignError(
            f"{calibration_path}: expected schema {CALIBRATION_SCHEMA}"
        )
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    frozen = parent["frozen_inputs"]
    config = _require_file(Path(str(frozen["config"])))
    recipes = _require_file(Path(str(frozen["recipes"])))
    legacy_input = _require_file(Path(str(frozen["legacy_input"])))
    refinements = (
        _require_file(Path(str(frozen["refinements"])))
        if frozen.get("refinements")
        else None
    )
    source_manifest = _load_json(
        Path(str(parent["pool_manifests"][0]["path"]))
    )
    groups: dict[tuple[str, int, float], list[dict[str, object]]] = {}
    decisions: list[dict[str, object]] = []
    for stratum in calibration.get("strata") or []:
        decision = _support_followup(
            stratum,
            minimum_targets=args.minimum_component_targets,
            minimum_inside_targets=args.minimum_provisional_inside_targets,
            confidence=args.zero_complement_confidence,
            maximum_zero_rate=args.maximum_zero_complement_target_rate,
            safety=args.followup_target_safety_factor,
            quantum=args.trial_quantum,
            maximum_trials=args.maximum_followup_trials,
            discovery_trials=args.discovery_trials,
        )
        if decision is None:
            decisions.append(
                {
                    "flat_index": int(stratum["flat_index"]),
                    "stratum_id": str(stratum["stratum_id"]),
                    "action": "none_ready",
                    "source_status": stratum["recommendation_status"],
                }
            )
            continue
        record = {
            "flat_index": int(stratum["flat_index"]),
            "stratum_id": str(stratum["stratum_id"]),
            "action": decision["kind"],
            "source_status": stratum["recommendation_status"],
            **decision,
        }
        decisions.append(record)
        if decision["schedulable"]:
            key = (
                str(decision["kind"]),
                int(decision["trials"]),
                float(decision["inside_fraction"]),
            )
            groups.setdefault(key, []).append(record)

    manifests: list[Path] = []
    tasks: list[dict[str, object]] = []
    seed_stride = 20_000_000
    for group_index, (key, records) in enumerate(
        sorted(groups.items()), start=1
    ):
        kind, trials, inside_fraction = key
        tag = f"{group_index:03d}_{kind}_{trials}"
        selection = output / "selection" / f"{tag}.txt"
        _write_flat_indices(selection, records)
        stage = output / "stages" / tag
        seed_base = args.seed_base + seed_stride * group_index
        if seed_base + 1000 * 12960 >= 2_147_483_647:
            raise CampaignError("follow-up seed schedule exceeds 32-bit range")
        manifest = radiative_mode4.prepare(
            _prepare_namespace(
                command="prepare-calibration",
                config=config,
                recipes=recipes,
                legacy_input=legacy_input,
                refinements=refinements,
                output=stage,
                candidate=str(source_manifest["guard_candidate"]),
                core_fraction=float(source_manifest["core_fraction"]),
                replicas=1,
                seed_base=seed_base,
                flat_index_file=selection,
                apply_y_max=bool(
                    source_manifest["analysis_selection"].get(
                        "apply_y_max", False
                    )
                ),
                heartbeat_interval=args.heartbeat_interval,
                generator_revision=str(source_manifest["generator_revision"]),
                trials=trials,
                inside_fraction=inside_fraction,
            )
        ).resolve()
        manifests.append(manifest)
        tasks.extend(
            _manifest_tasks(
                manifest,
                stage=tag,
                first_task_id=len(tasks) + 1,
            )
        )

    if not tasks:
        raise CampaignError(
            "calibration has no automatically schedulable follow-ups"
        )
    pool = [
        Path(str(item["path"])) for item in parent["pool_manifests"]
    ] + manifests
    fields = (
        "flat_index",
        "stratum_id",
        "source_status",
        "action",
        "schedulable",
        "inside_fraction",
        "trials",
        "required_trials_estimate",
        "capped_at_maximum",
    )
    with (output / "followup_plan.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        writer = csv.DictWriter(
            destination, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(decisions)
    manual = [
        item
        for item in decisions
        if item.get("schedulable") is False
    ]
    payload = _campaign_payload(
        kind="calibration_followup",
        root=output,
        tasks=tasks,
        manifests=manifests,
        pool_manifests=pool,
        frozen_inputs=frozen,
        selection=parent["selection"],
        parent=parent_path,
        calibration_report=calibration_path,
        followup={
            "groups": len(groups),
            "scheduled_strata": sum(len(items) for items in groups.values()),
            "ready_strata": sum(
                item["action"] == "none_ready" for item in decisions
            ),
            "manual_review_strata": len(manual),
            "manual_review": manual,
            "plan_tsv": str(output / "followup_plan.tsv"),
            "plan_tsv_sha256": _sha256(output / "followup_plan.tsv"),
        },
    )
    return _finish_campaign(output, payload, tasks)


def plan_production(args: argparse.Namespace) -> Path:
    parent_path, parent = _load_campaign(args.campaign)
    calibration_path = _require_file(args.calibration)
    calibration = _load_json(calibration_path)
    if calibration.get("schema") != CALIBRATION_SCHEMA:
        raise CampaignError(
            f"{calibration_path}: expected schema {CALIBRATION_SCHEMA}"
        )
    ready = {
        int(item["flat_index"]): item
        for item in calibration.get("strata") or []
        if str(item.get("pilot_readiness")) in READY_STATES
        and item.get("recommended_envelope") is not None
    }
    expected_file = Path(str(parent["selection"]["flat_indices"]))
    expected = {
        int(line.strip())
        for line in expected_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    missing = sorted(expected - set(ready))
    if missing and not args.allow_incomplete:
        preview = ", ".join(str(value) for value in missing[:20])
        suffix = "..." if len(missing) > 20 else ""
        raise CampaignError(
            f"{len(missing)} selected strata are not pilot-ready: "
            f"{preview}{suffix}"
        )
    selected = sorted(expected & set(ready))
    if not selected:
        raise CampaignError("production selection contains no ready strata")
    output = args.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    selection = output / "selection" / "ready_flat_indices.txt"
    _write_flat_indices(
        selection,
        [{"flat_index": value} for value in selected],
    )
    frozen = parent["frozen_inputs"]
    config = _require_file(Path(str(frozen["config"])))
    recipes = _require_file(Path(str(frozen["recipes"])))
    legacy_input = _require_file(Path(str(frozen["legacy_input"])))
    refinements = (
        _require_file(Path(str(frozen["refinements"])))
        if frozen.get("refinements")
        else None
    )
    source_manifest = _load_json(
        Path(str(parent["pool_manifests"][0]["path"]))
    )
    stage = output / "stages" / "production"
    manifest = radiative_mode4.prepare(
        _prepare_namespace(
            command="prepare",
            config=config,
            recipes=recipes,
            legacy_input=legacy_input,
            refinements=refinements,
            output=stage,
            candidate=str(source_manifest["guard_candidate"]),
            core_fraction=float(source_manifest["core_fraction"]),
            replicas=args.replicas,
            seed_base=args.seed_base,
            flat_index_file=selection,
            apply_y_max=bool(
                source_manifest["analysis_selection"].get(
                    "apply_y_max", False
                )
            ),
            heartbeat_interval=args.heartbeat_interval,
            generator_revision=str(source_manifest["generator_revision"]),
            envelope_report=calibration_path,
            events_per_stratum=args.events_per_stratum,
        )
    ).resolve()
    tasks = _manifest_tasks(manifest, stage="production")
    payload = _campaign_payload(
        kind="production",
        root=output,
        tasks=tasks,
        manifests=[manifest],
        pool_manifests=[manifest],
        frozen_inputs=frozen,
        selection={
            **parent["selection"],
            "ready_strata": len(selected),
            "not_ready_strata": len(missing),
            "allow_incomplete": args.allow_incomplete,
            "flat_indices": str(selection),
            "flat_indices_sha256": _sha256(selection),
            "events_per_stratum": args.events_per_stratum,
            "replicas": args.replicas,
        },
        parent=parent_path,
        calibration_report=calibration_path,
    )
    return _finish_campaign(output, payload, tasks)


def validate_pilots(args: argparse.Namespace) -> Path:
    campaign_path, campaign = _load_campaign(args.campaign)
    if campaign["kind"] != "production":
        raise CampaignError("pilot validation requires a production campaign")
    _assert_complete(campaign_path, campaign)
    runs = [Path(task["expected_run"]) for task in _read_tasks(campaign)]
    output = (
        args.output.expanduser().resolve()
        if args.output is not None
        else Path(str(campaign["root"])) / "pilot_validation.json"
    )
    return radiative_mode4.validate_pilots(
        argparse.Namespace(
            calibration=args.calibration,
            runs=runs,
            output=output,
            minimum_runs=args.minimum_runs,
            minimum_events=args.minimum_events,
            maximum_duplicate_fraction=args.maximum_duplicate_fraction,
            maximum_guard_complement_fraction=(
                args.maximum_guard_complement_fraction
            ),
            maximum_relative_cross_section_difference=(
                args.maximum_relative_cross_section_difference
            ),
            maximum_cross_section_z_score=args.maximum_cross_section_z_score,
            confidence=args.confidence,
            allow_pilot_revision_mismatch=False,
            revision_compatibility_rationale=None,
        )
    )


def emit_swif(args: argparse.Namespace) -> Path:
    campaign_path, campaign = _load_campaign(args.campaign)
    executable = args.executable.expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise FileNotFoundError(executable)
    workflow = args.workflow.strip()
    if not workflow or any(
        character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-"
        for character in workflow
    ):
        raise ValueError("--workflow contains unsupported characters")
    root = Path(str(campaign["root"]))
    wrapper = root / "run_swif_task.sh"
    driver = Path(__file__).resolve()
    executable_sha256 = _sha256(executable)
    wrapper.write_text(
        "#!/bin/bash\n"
        "set -euo pipefail\n"
        "if [ -f /etc/profile.d/modules.sh ]; then\n"
        "  source /etc/profile.d/modules.sh\n"
        "fi\n"
        "module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles 2>/dev/null || true\n"
        "module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles 2>/dev/null || true\n"
        "module load clas12/5.4 2>/dev/null || true\n"
        f"expected_sha256={shlex.quote(executable_sha256)}\n"
        'actual_sha256=$(sha256sum "$3" | awk \'{print $1}\')\n'
        'if [ "$actual_sha256" != "$expected_sha256" ]; then\n'
        '  echo "ERROR: AAO executable hash changed" >&2\n'
        "  exit 3\n"
        "fi\n"
        f"python3 {shlex.quote(str(driver))} run-task "
        '  --campaign "$1" --task-id "$2" --executable "$3"\n',
        encoding="utf-8",
    )
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    output = (
        args.output.expanduser().resolve()
        if args.output is not None
        else root / "submit_swif.sh"
    )
    task_table = Path(str(campaign["tasks_tsv"]))
    script = f"""#!/bin/bash
set -euo pipefail

WORKFLOW={shlex.quote(workflow)}
CAMPAIGN={shlex.quote(str(campaign_path))}
EXECUTABLE={shlex.quote(str(executable))}
TASKS={shlex.quote(str(task_table))}
WRAPPER={shlex.quote(str(wrapper))}

swif2 create -workflow "$WORKFLOW"

tail -n +2 "$TASKS" | while IFS=$'\\t' read -r task_id stage operation manifest manifest_sha256 flat_index stratum_id replica_index seed requested_trials requested_events expected_run; do
  job_name="${{WORKFLOW}}_t$(printf '%06d' "$task_id")_${{stratum_id}}_r$(printf '%03d' "$replica_index")"
  swif2 add-job \\
    -workflow "$WORKFLOW" \\
    -name "$job_name" \\
    -cores {args.cores} \\
    -disk {shlex.quote(args.disk)} \\
    -ram {shlex.quote(args.ram)} \\
    -time {shlex.quote(args.walltime)} \\
    -os el9 \\
    -- /bin/bash "$WRAPPER" "$CAMPAIGN" "$task_id" "$EXECUTABLE"
done

swif2 run "$WORKFLOW"

echo "Submitted $WORKFLOW"
echo "Status: swif2 status $WORKFLOW -summary -problems"
"""
    output.write_text(script, encoding="utf-8")
    output.chmod(output.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    _write_json(
        root / "swif_execution.json",
        {
            "schema": "aao-rad-full-campaign-swif-execution-v1",
            "created_utc": _now(),
            "campaign": str(campaign_path),
            "campaign_sha256": _sha256(campaign_path),
            "workflow": workflow,
            "executable": str(executable),
            "executable_sha256": executable_sha256,
            "wrapper": str(wrapper),
            "wrapper_sha256": _sha256(wrapper),
            "submission_script": str(output),
            "submission_script_sha256": _sha256(output),
            "resources": {
                "cores": args.cores,
                "disk": args.disk,
                "ram": args.ram,
                "walltime": args.walltime,
            },
            "submitted": False,
            "submission_is_external": True,
        },
    )
    return output


def _add_finalize_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--output", type=Path)
    parser.add_argument("--envelope-safety-factor", type=float, default=1.2)
    parser.add_argument("--maximum-duplicate-fraction", type=float, default=0.05)
    parser.add_argument("--minimum-component-targets", type=int, default=20)
    parser.add_argument(
        "--minimum-provisional-inside-targets", type=int, default=1000
    )
    complement = parser.add_mutually_exclusive_group()
    complement.add_argument(
        "--allow-zero-complement",
        dest="allow_zero_complement",
        action="store_true",
        help="allow the audited provisional zero-complement stopping rule",
    )
    complement.add_argument(
        "--strict-complement",
        dest="allow_zero_complement",
        action="store_false",
        help="require observed support in both calibration components",
    )
    parser.set_defaults(allow_zero_complement=True)
    parser.add_argument(
        "--zero-complement-confidence", type=float, default=0.95
    )
    parser.add_argument(
        "--maximum-zero-complement-target-rate", type=float, default=1.0e-6
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    plan = subparsers.add_parser(
        "plan-calibration",
        help="prepare a complete queue calibration and farm task table",
    )
    plan.add_argument("--queue", type=Path, required=True)
    plan.add_argument("--config", type=Path, required=True)
    plan.add_argument("--recipes", type=Path, required=True)
    plan.add_argument("--refinements", type=Path)
    plan.add_argument("--input", type=Path, required=True)
    plan.add_argument("--output", type=Path, required=True)
    plan.add_argument(
        "--selection",
        choices=("selected", "data-occupied", "all-catalog"),
        default="selected",
    )
    plan.add_argument(
        "--work-category",
        action="append",
        choices=sorted(SELECTABLE_WORK_CATEGORIES),
    )
    plan.add_argument("--candidate", default="padding_0p035")
    plan.add_argument("--core-fraction", type=float, default=0.90)
    plan.add_argument("--inside-guard-trial-fraction", type=float, default=0.50)
    plan.add_argument("--trials", type=int, default=7_000_000)
    plan.add_argument("--replicas", type=int, default=1)
    plan.add_argument("--seed-base", type=int, default=907_001)
    plan.add_argument("--heartbeat-interval", type=int, default=100_000)
    plan.add_argument("--apply-y-max", action="store_true")
    plan.add_argument("--generator-revision", required=True)

    run_parser = subparsers.add_parser(
        "run-task", help="run one immutable task from a campaign"
    )
    run_parser.add_argument("--campaign", type=Path, required=True)
    run_parser.add_argument("--task-id", type=int, required=True)
    run_parser.add_argument("--executable", type=Path, required=True)
    run_parser.add_argument("--overwrite", action="store_true")

    status_parser = subparsers.add_parser(
        "status", help="scan task receipts and validated run artifacts"
    )
    status_parser.add_argument("--campaign", type=Path, required=True)
    status_parser.add_argument("--output", type=Path)

    finalize_parser = subparsers.add_parser(
        "finalize", help="pool every completed campaign manifest"
    )
    finalize_parser.add_argument("--campaign", type=Path, required=True)
    _add_finalize_arguments(finalize_parser)

    followup = subparsers.add_parser(
        "plan-followup",
        help="automatically size and prepare only non-ready strata",
    )
    followup.add_argument("--campaign", type=Path, required=True)
    followup.add_argument("--calibration", type=Path, required=True)
    followup.add_argument("--output", type=Path, required=True)
    followup.add_argument("--minimum-component-targets", type=int, default=20)
    followup.add_argument(
        "--minimum-provisional-inside-targets", type=int, default=1000
    )
    followup.add_argument(
        "--zero-complement-confidence", type=float, default=0.95
    )
    followup.add_argument(
        "--maximum-zero-complement-target-rate", type=float, default=1.0e-6
    )
    followup.add_argument(
        "--followup-target-safety-factor", type=float, default=1.5
    )
    followup.add_argument("--trial-quantum", type=int, default=1_000_000)
    followup.add_argument(
        "--maximum-followup-trials", type=int, default=100_000_000
    )
    followup.add_argument("--discovery-trials", type=int, default=20_000_000)
    followup.add_argument("--heartbeat-interval", type=int, default=100_000)
    followup.add_argument("--seed-base", type=int, default=107_000_001)

    production = subparsers.add_parser(
        "plan-production",
        help="prepare per-stratum generation from a ready envelope catalog",
    )
    production.add_argument("--campaign", type=Path, required=True)
    production.add_argument("--calibration", type=Path, required=True)
    production.add_argument("--output", type=Path, required=True)
    production.add_argument("--events-per-stratum", type=int, required=True)
    production.add_argument("--replicas", type=int, default=1)
    production.add_argument("--seed-base", type=int, default=307_000_001)
    production.add_argument("--heartbeat-interval", type=int, default=100_000)
    production.add_argument("--allow-incomplete", action="store_true")

    validation = subparsers.add_parser(
        "validate-pilots", help="validate all completed production pilot runs"
    )
    validation.add_argument("--campaign", type=Path, required=True)
    validation.add_argument("--calibration", type=Path, required=True)
    validation.add_argument("--output", type=Path)
    validation.add_argument("--minimum-runs", type=int, default=2)
    validation.add_argument("--minimum-events", type=int, default=400)
    validation.add_argument("--maximum-duplicate-fraction", type=float, default=0.05)
    validation.add_argument(
        "--maximum-guard-complement-fraction", type=float, default=0.02
    )
    validation.add_argument(
        "--maximum-relative-cross-section-difference", type=float, default=0.10
    )
    validation.add_argument(
        "--maximum-cross-section-z-score", type=float, default=3.0
    )
    validation.add_argument("--confidence", type=float, default=0.95)

    swif = subparsers.add_parser(
        "emit-swif", help="write but do not execute a SWIF2 submission script"
    )
    swif.add_argument("--campaign", type=Path, required=True)
    swif.add_argument("--workflow", required=True)
    swif.add_argument("--executable", type=Path, required=True)
    swif.add_argument("--output", type=Path)
    swif.add_argument("--cores", type=int, default=1)
    swif.add_argument("--disk", default="2gb")
    swif.add_argument("--ram", default="1gb")
    swif.add_argument("--walltime", default="8hr")
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.command == "plan-calibration":
        result = plan_calibration(args)
    elif args.command == "run-task":
        result = run_task(args)
    elif args.command == "status":
        result = status(args)
    elif args.command == "finalize":
        result = finalize(args)
    elif args.command == "plan-followup":
        result = plan_followup(args)
    elif args.command == "plan-production":
        result = plan_production(args)
    elif args.command == "validate-pilots":
        result = validate_pilots(args)
    else:
        result = emit_swif(args)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
