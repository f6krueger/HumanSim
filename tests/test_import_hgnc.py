import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
import import_hgnc  # noqa: E402

TSV = ("hgnc_id\tsymbol\tname\tlocus_group\tlocus_type\tstatus\tlocation\n"
       "HGNC:1\tMSTN\tmyostatin\tprotein-coding gene\tgene with protein product\tApproved\t2q32.2\n"
       "HGNC:2\tOLDG\told gene\tprotein-coding gene\tgene with protein product\tEntry Withdrawn\t1p1\n"
       "HGNC:3\tMIR1\tmicroRNA 1\tnon-coding RNA\tRNA, micro\tApproved\t1q1\n")


def test_import_keeps_approved_genes_and_filters(tmp_path):
    src = tmp_path / "h.txt"
    src.write_text(TSV)
    assert [g["symbol"] for g in import_hgnc.read_hgnc(src)] == ["MIR1", "MSTN"]
    assert [g["symbol"] for g in import_hgnc.read_hgnc(src, protein_coding_only=True)] == ["MSTN"]
    out = tmp_path / "idx.yaml"
    assert import_hgnc.main([str(src), "--out", str(out)]) == 0
    assert {g["symbol"] for g in yaml.safe_load(out.read_text())["genes"]} == {"MIR1", "MSTN"}


def test_import_stops_on_missing_column(tmp_path):
    src = tmp_path / "bad.txt"
    src.write_text("symbol\tname\nA\tb\n")
    with pytest.raises(SystemExit):
        import_hgnc.read_hgnc(src)
