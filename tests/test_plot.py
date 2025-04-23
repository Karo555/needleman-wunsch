import os
from src.aligner.plot import plot_matrix


def test_plot_matrix(tmp_path):
    # Simple 2x2 matrix
    matrix = [[0, -1], [-1, 1]]
    out_file = tmp_path / "heatmap.png"
    fig = plot_matrix(matrix, str(out_file))
    assert fig is not None
    assert os.path.exists(str(out_file))
    # File should not be empty
    assert out_file.stat().st_size > 0
