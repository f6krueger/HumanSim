import http.client
import json
import threading
from http.server import ThreadingHTTPServer

import pytest

from humansim.server import ApiError, App, make_handler


@pytest.fixture(scope="module")
def app():
    return App()


def test_meta_lists_everything_the_ui_needs(app):
    m = app.meta()
    assert {"parameters", "groups", "hormones", "drugs", "genes", "catalogue", "cycle_phases"} <= set(m)
    assert {d["id"] for d in m["drugs"]} >= {"acetaminophen", "ibuprofen", "midazolam"}
    json.dumps(m, allow_nan=False)


def test_resolve_edit_ripples_and_is_json_safe(app):
    r = app.resolve({"sex": "male", "solve": True, "edits": {"liver_mass": 2.5}})
    assert r["values"]["liver_mass"]["flags"] == ["edited"]
    assert r["values"]["resting_metabolic_rate"]["value"] > r["values"]["resting_metabolic_rate"]["baseline"]
    assert r["solver"]["converged"] and not r["hypothetical"]
    json.dumps(r, allow_nan=False)


def test_hypothetical_and_germline_are_reported(app):
    r = app.resolve({"solve": True, "gene_edits": [
        {"gene": "MSTN", "locus": "MSTN", "state": "loss_of_function", "mode": "germline"}]})
    assert r["hypothetical"] and r["germline"]
    r2 = app.resolve({"gene_edits": [{"gene": "X", "hypothesis": [{"param": "heart_mass", "op": "ratio", "value": 1.3}]}]})
    assert r2["hypothetical"] and "hypothetical" in r2["values"]["heart_mass"]["flags"]


def test_bad_input_gives_clear_errors(app):
    for bad in ({"sex": "other"}, {"age": "old"}, {"age": -1}, {"edits": {"height": "tall"}}, {"edits": {"nope": 1}},
                {"genotype": {"CYP2C9": "*1/*9"}}, {"hormones": {"body_mass": 80}},
                {"gene_edits": [{"gene": "X", "locus": "MSTN"}]},
                {"gene_edits": [{"gene": "X", "hypothesis": [{"param": "heart_mass", "op": "ratio", "value": "x"}]}]},
                {"sex": "male", "cycle_day": 5}, {"edits": {"body_mass": float("inf")}}):
        with pytest.raises(ApiError):
            app.resolve(bad)


def test_compare_pk_hpg_cycle(app):
    c = app.compare({"edits": {"liver_mass": 2.0}})
    assert c["male"]["values"]["height"]["value"] > c["female"]["values"]["height"]["value"]
    pk = app.pk({"drug": "ibuprofen", "dose_mg": 800, "infusion_h": 1, "hours": 24, "genotype": {"CYP2C9": "*1/*3"}})
    base = app.pk({"drug": "ibuprofen", "dose_mg": 800, "infusion_h": 1, "hours": 24})
    assert pk["half_life_h"] > base["half_life_h"] and "genotype_effect" in pk["flags"] and not base["hypothetical"]
    assert len(pk["t"]) <= 401 and pk["cmax_mg_l"] > 0
    for bad in ({"drug": "nope", "dose_mg": 1}, {"drug": "ibuprofen", "dose_mg": -5}, {"drug": "ibuprofen"},
                {"drug": "ibuprofen", "dose_mg": 10, "hours": 100000}, {"drug": "ibuprofen", "dose_mg": 10, "repeats": 500},
                {"drug": "ibuprofen", "dose_mg": 10, "genotype": {"CYP2C9": "*3/*3"}}):
        with pytest.raises(ApiError):
            app.pk(bad)
    h = app.hpg({"days": 100, "added_testosterone": 800, "dose_days": 50})
    assert min(h["LH"]) < 0.5 * h["reference"]["LH"] and h["assumed"]
    assert h["T"][-1] == pytest.approx(h["reference"]["T"], rel=0.02)
    with pytest.raises(ApiError):
        app.hpg({"days": -1})
    cy = app.cycle()
    assert len(cy["days"]) == 28 and max(cy["series"]["lh"]) == pytest.approx(cy["series"]["lh"][13])


def test_http_server_serves_page_api_and_refuses_foreign_hosts(app):
    srv = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(app))
    port = srv.server_address[1]
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        def call(method, path, body=None, host=None):
            c = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
            c.putrequest(method, path, skip_host=True)
            c.putheader("Host", host or f"127.0.0.1:{port}")
            data = json.dumps(body).encode() if body is not None else (b"" if method == "POST" else None)
            if data is not None:
                c.putheader("Content-Length", str(len(data)))
            c.endheaders(data)
            r = c.getresponse()
            return r.status, r.read()

        s, b = call("GET", "/")
        assert s == 200 and b"HumanSim" in b
        s, b = call("POST", "/api/resolve", {"edits": {"body_mass": 90}})
        assert s == 200 and json.loads(b)["values"]["body_mass"]["value"] == 90
        assert call("POST", "/api/resolve", {"sex": "x"})[0] == 400
        assert call("GET", "/nope")[0] == 404
        assert call("GET", "/", host="evil.example.com")[0] == 403
        c = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
        c.request("POST", "/api/resolve", body=b"{not json")
        assert c.getresponse().status == 400
    finally:
        srv.shutdown()
        srv.server_close()
