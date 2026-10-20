import pytest


@pytest.mark.parametrize(
    "script", ["02_scanners", "03_diffusion", "04_structural", "05_functional"]
)
def test_prepare_import_does_not_execute_pipeline(modules, tmp_path, script):
    module = modules(f"prepare/{script}.py")
    assert callable(module.main)
    assert not (tmp_path / "data").exists()
