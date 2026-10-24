"""Exercise attempt integrity without requiring the scientific runtime."""

import ast
import hashlib
from pathlib import Path
from types import CodeType, FunctionType, SimpleNamespace

import pytest


@pytest.fixture
def integrity():
    # These helpers use only the standard library; avoid importing Gammapy.
    source = Path(__file__).parents[2] / "v2dl5/scripts/reflected_region_analysis.py"
    tree = ast.parse(source.read_text())
    names = {"_observation_inputs", "_check_resume"}
    tree.body = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    namespace = {"hashlib": hashlib, "Path": Path}
    code = compile(tree, str(source), "exec")
    return SimpleNamespace(
        **{
            item.co_name: FunctionType(item, namespace)
            for item in code.co_consts
            if isinstance(item, CodeType)
        }
    )


def test_missing_selected_input_is_recorded_and_unselected_input_ignored(integrity, tmp_path):
    available = tmp_path / "available.fits"
    available.write_bytes(b"events")
    rows = [
        {"OBS_ID": run, "FILE_DIR": "", "FILE_NAME": name}
        for run, name in [(1, available.name), (2, "missing.fits"), (3, "outside.fits")]
    ]
    data = SimpleNamespace(runs=[1, 2], get_data_store=lambda: SimpleNamespace(hdu_table=rows))
    inputs = integrity._observation_inputs(data, tmp_path)
    assert inputs == {
        str(available.resolve()): hashlib.sha256(b"events").hexdigest(),
        str((tmp_path / "missing.fits").resolve()): None,
    }
    (tmp_path / "missing.fits").write_bytes(b"restored")
    assert integrity._observation_inputs(data, tmp_path) != inputs

def test_exact_resume_verifies_products(integrity, tmp_path):
    product = tmp_path / "flux.ecsv"
    product.write_bytes(b"flux")
    previous = {
        "status": "partial",
        "fingerprint": "same",
        "products": {product.name: hashlib.sha256(b"flux").hexdigest()},
    }
    for state in ("complete", "partial"):
        previous["status"] = state
        integrity._check_resume(previous, "same", tmp_path)
    product.write_bytes(b"corrupt")
    with pytest.raises(FileExistsError, match="product"):
        integrity._check_resume(previous, "same", tmp_path)
    product.unlink()
    with pytest.raises(FileExistsError, match="product"):
        integrity._check_resume(previous, "same", tmp_path)

def test_changed_attempt_is_preserved(integrity, tmp_path):
    stale = tmp_path / "flux_points.ecsv"
    stale.write_text("old result")
    for state in ("complete", "partial", "analysis_failed"):
        previous = {"status": state, "fingerprint": "old", "products": {}}
        with pytest.raises(FileExistsError, match="new output directory"):
            integrity._check_resume(previous, "new", tmp_path)
    assert stale.read_text() == "old result"

def test_resume_requires_completed_manifest(integrity, tmp_path):
    for previous in [
        {"status": "analysis_failed", "fingerprint": "same"},
        {"status": "complete", "fingerprint": "same"},
        {"status": "complete", "fingerprint": "same", "products": {}},
    ]:
        with pytest.raises(FileExistsError, match="new output directory"):
            integrity._check_resume(previous, "same", tmp_path)
