"""Figure loft dumps: OBJ n-gons, VTK POLYDATA, NPZ polygons."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
from _shared.debug_mesh import _te_xyz, dump_debug_mesh, write_loft_obj, write_loft_vtk


def _faces():
    quad = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0], [0.0, 0.0, 1.0]])
    rib = np.array([
        [0.0, 0.0, 0.5], [0.4, 0.0, 0.5], [0.7, 0.1, 0.5],
        [0.4, 0.2, 0.5], [0.0, 0.2, 0.5],
    ])
    return [quad, rib], ["le_su", "rib_0-cell_0"]


def test_te_xyz_picks_max_chord_frac_not_blade_x():
    samples = {
        "le_su": np.array([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]]),
        "te": np.array([[0.6, 0.1, 0.0], [1.0, 0.0, 0.0], [0.6, -0.1, 0.0]]),
        "te_pr": np.array([[0.95, 0.0, 0.0], [0.4, 0.0, 0.0]]),
    }
    pt = _te_xyz(samples)
    assert np.allclose(pt[:2], [1.0, 0.0])
    pt = _te_xyz(samples, oml_te=[1.0, 0.0])
    assert np.allclose(pt[:2], [1.0, 0.0])


def test_write_loft_obj_keeps_rib_ngon(tmp_path: Path):
    faces, kinds = _faces()
    text = Path(write_loft_obj(tmp_path / "mids.obj", faces, kinds)).read_text(encoding="utf-8")
    assert "g le_su" in text and "g rib_0-cell_0" in text
    ngons = [ln for ln in text.splitlines() if ln.startswith("f ") and len(ln.split()) == 6]
    assert ngons == ["f 5 6 7 8 9"]


def test_write_loft_vtk_polydata_roles(tmp_path: Path):
    faces, kinds = _faces()
    text = Path(write_loft_vtk(tmp_path / "mids.vtk", faces, kinds)).read_text(encoding="utf-8")
    assert text.startswith("# vtk DataFile Version 3.0")
    assert "DATASET POLYDATA" in text
    assert "POLYGONS 2 " in text
    assert "COLOR_SCALARS rgb 3" in text
    assert "role 1 2 string" in text
    assert "rib_0-cell_0" in text
    assert "5 4 5 6 7 8" in text


def test_dump_debug_mesh_obj_vtk_npz(tmp_path: Path):
    faces, kinds = _faces()
    tmpl = SimpleNamespace(n_chord=4, loop=(), webs=())
    a, b = faces[0][:2], faces[0][2:]
    walls = [
        (0.0, tmpl, {"le_su": a, "te_pr": b}),
        (1.0, tmpl, {"le_su": a + np.array([0.0, 0.0, 1.0]), "te_pr": b + np.array([0.0, 0.0, 1.0])}),
    ]
    dump_debug_mesh(tmp_path, walls, faces, kinds, label="unit", repo=tmp_path, stem="mids")
    dump_debug_mesh(
        tmp_path, walls, [faces[1]], [kinds[1]],
        label="unit cutaway", repo=tmp_path, stem="cutaway", write_samples=False,
    )
    assert (tmp_path / "mids.obj").is_file()
    assert (tmp_path / "mids.vtk").is_file()
    assert (tmp_path / "cutaway.obj").is_file()
    assert (tmp_path / "cutaway.vtk").is_file()
    data = np.load(tmp_path / "samples.npz", allow_pickle=True)
    assert len(data["polygons"]) == 2
    assert data["polygons"][1].shape == (5, 3)
    assert "rib_0-cell_0" in data["polygon_roles"].tolist()
    assert data["faces"].shape[1:] == (4, 3)
