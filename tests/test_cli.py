import subprocess
import filecmp
from pathlib import Path

def test_2fast2q_cli_creates_expected_output(tmp_path):

    result = subprocess.run(
        ["2fast2q", "-c", "-t"],
        cwd=tmp_path,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    assert result.returncode == 0, f"CLI failed: {result.stderr}"

    subdirs = [d for d in tmp_path.iterdir() if d.is_dir()]
    assert len(subdirs) == 1, f"Expected exactly one output folder, found: {len(subdirs)}"

    output_dir = subdirs[0]
    files_in_output = list(output_dir.iterdir())
    assert len(files_in_output) == 6, f"Expected 6 files in output folder, found: {len(files_in_output)}"

    output_file = output_dir / "compiled.csv"
    assert output_file.exists(), f"compiled.csv not found in output folder {output_dir}"

    import difflib

    if not filecmp.cmp(output_file, expected_file, shallow=False):
        with open(output_file) as f1, open(expected_file) as f2:
            diff = difflib.unified_diff(
                f1.readlines(), f2.readlines(),
                fromfile='actual compiled.csv',
                tofile='expected compiled.csv',
            )
            print(''.join(diff))
        assert False, "compiled.csv does not match expected output"

    expected_file = Path("tests/compiled.csv")
    assert filecmp.cmp(output_file, expected_file), "Output file does not match expected result"
