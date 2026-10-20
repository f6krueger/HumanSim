"""Local web server for the HumanSim UI (standard library only).

`python -m humansim serve` starts it on 127.0.0.1 and serves one page (web/index.html) plus a small
JSON API. It is a local tool, not a hosted service: it binds to the loopback address only, has no
accounts or storage, and every request is computed fresh from the library. The API functions below
are plain functions of a JSON payload, so they can be tested without a socket.
"""

from __future__ import annotations

import json
import math
import re
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from . import cycle as _cycle
from . import hpg as _hpg
from .context import Context, resolve
from .couplings import load_couplings
from .editing import Effect, GeneEdit, load_catalogue
from .genetics import load_genes, validate_genotype
from .modifiers import load_modifiers
from .pharmacology import Dose, load_drugs, simulate
from .registry import Registry
from .solver import solve

WEB_DIR = Path(__file__).parent / "web"
MAX_POINTS = 400


class ApiError(Exception):
    def __init__(self, message: str, status: int = 400):
        super().__init__(message)
        self.status = status


class App:
    """Loads the data once and answers API calls."""

    def __init__(self) -> None:
        self.reg = Registry.load()
        self.mods = load_modifiers()
        self.cps = load_couplings(self.reg)
        self.drugs = load_drugs()
        self.genes = load_genes()
        self.catalogue = load_catalogue()

    # ---- helpers -----------------------------------------------------------------------------
    @staticmethod
    def _num(x, what: str) -> float:
        try:
            v = float(x)
        except (TypeError, ValueError):
            raise ApiError(f"{what} must be a number") from None
        if not math.isfinite(v):
            raise ApiError(f"{what} must be finite")
        return v

    def context(self, p: dict) -> Context:
        sex = p.get("sex", "male")
        if sex not in ("male", "female"):
            raise ApiError("sex must be male or female")
        age: object = p.get("age", "adult")
        if age != "adult":
            age = self._num(age, "age")
            age = int(age) if float(age).is_integer() else age
            if age < 0:
                raise ApiError("age must be 0 or more")
        edits = {k: self._num(v, f"edit {k}") for k, v in (p.get("edits") or {}).items()}
        hormones = {k: self._num(v, f"hormone {k}") for k, v in (p.get("hormones") or {}).items()}
        genotype = dict(p.get("genotype") or {})
        try:
            validate_genotype(self.genes, genotype)
        except (KeyError, ValueError) as exc:
            raise ApiError(str(exc.args[0] if exc.args else exc)) from None
        gene_edits = []
        for i, g in enumerate(p.get("gene_edits") or []):
            try:
                hyp = tuple(Effect(param=h["param"], op=h["op"], value=float(h["value"]), tier="hypothetical",
                                   rationale="supplied by the user as a hypothesis; not data")
                            for h in g.get("hypothesis") or [])
                gene_edits.append(GeneEdit(gene=g["gene"], locus=g.get("locus") or None, state=g.get("state") or None,
                                           hypothesis=hyp, mode=g.get("mode", "somatic"),
                                           fraction=float(g.get("fraction", 1.0)), tissue=g.get("tissue") or None))
            except (KeyError, ValueError, TypeError) as exc:
                raise ApiError(f"gene edit {i + 1}: {exc.args[0] if exc.args else exc}") from None
        kwargs = dict(sex=sex, age=age, growth_z=self._num(p.get("z", 0.0), "z"), edits=edits,
                      genotype=genotype, gene_edits=tuple(gene_edits), hormone_state=hormones,
                      cycle_phase=p.get("cycle_phase") or None)
        day = p.get("cycle_day")
        if day not in (None, ""):
            if sex != "female":
                raise ApiError("cycle day applies to women")
            ctx = _cycle.cycle_context(self.reg, self._num(day, "cycle day"), age=age, edits=edits,
                                       genotype=genotype, gene_edits=tuple(gene_edits))
            ctx.hormone_state.update(hormones)
            return ctx
        return Context(**kwargs)

    def state(self, ctx: Context, do_solve: bool):
        try:
            return solve(self.reg, ctx, self.mods, self.cps) if do_solve else resolve(self.reg, ctx, self.mods)
        except (ValueError, KeyError) as exc:
            raise ApiError(str(exc.args[0] if exc.args else exc)) from None

    @staticmethod
    def _clean(v):
        return None if v is None or (isinstance(v, float) and not math.isfinite(v)) else v

    # ---- endpoints ---------------------------------------------------------------------------
    def meta(self) -> dict:
        params = [{"id": p.id, "name": p.name, "unit": p.unit, "group": p.group, "kind": p.kind}
                  for p in self.reg]
        return {
            "parameters": params, "groups": self.reg.groups(),
            "hormones": [p.id for p in self.reg.by_group("hormones")],
            "drugs": [{"id": d.id, "name": d.name, "genes": [e["gene"] for e in d.genotype_effects]}
                      for d in self.drugs.values()],
            "genes": [{"id": g.id, "name": g.name, "alleles": sorted(g.alleles),
                       "frequencies": {k: round(v, 4) for k, v in g.diplotype_frequencies("european").items()}}
                      for g in self.genes.values()],
            "catalogue": [{"id": loc.id, "gene": loc.gene, "name": loc.name,
                           "states": [{"state": s, "description": d,
                                       "effects": [{"param": e.param, "op": e.op, "value": e.value, "tier": e.tier,
                                                    "rationale": e.rationale} for e in effs]}
                                      for s, (d, effs) in loc.states.items()]}
                          for loc in self.catalogue.values()],
            "cycle_phases": ["follicular", "midcycle", "luteal", "postmenopausal"],
        }

    def resolve(self, p: dict) -> dict:
        ctx = self.context(p)
        st = self.state(ctx, bool(p.get("solve")))
        base = resolve(self.reg, Context(sex=ctx.sex, age=ctx.age, cycle_phase=ctx.cycle_phase,
                                         growth_z=ctx.growth_z), self.mods)
        rows = {}
        for prm in self.reg:
            v, b = st.values.get(prm.id), base.values.get(prm.id)
            if v is None and b is None:
                continue
            rows[prm.id] = {"value": self._clean(v), "baseline": self._clean(b),
                            "flags": sorted(st.flags.get(prm.id, ()))}
        flags = set().union(*st.flags.values()) if st.flags else set()
        return {
            "values": rows,
            "trace": [{"param": t.param, "stage": t.stage, "source": t.source, "before": self._clean(t.before),
                       "after": self._clean(t.after), "note": t.note} for t in st.trace],
            "solver": ({k: v for k, v in st.solver.items() if k != "skipped"} | {"skipped": len(st.solver["skipped"])}
                       if st.solver else None),
            "notes": [str(n) for n in st.notes], "phase": st.phase,
            "hypothetical": bool({"hypothetical", "downstream_of_hypothetical"} & flags),
            "germline": "germline" in flags,
        }

    def compare(self, p: dict) -> dict:
        out = {}
        for sex in ("male", "female"):
            out[sex] = self.resolve({**p, "sex": sex, "cycle_day": p.get("cycle_day") if sex == "female" else None})
        return out

    def pk(self, p: dict) -> dict:
        drug = self.drugs.get(p.get("drug"))
        if drug is None:
            raise ApiError(f"unknown drug {p.get('drug')!r}")
        ctx = self.context(p)
        st = self.state(ctx, bool(p.get("solve", True)))
        dose = self._num(p.get("dose_mg"), "dose")
        infusion = self._num(p.get("infusion_h", 0.0), "infusion time")
        hours = self._num(p.get("hours", 24.0), "hours")
        if dose <= 0 or infusion < 0 or not 0 < hours <= 500:
            raise ApiError("need dose > 0, infusion >= 0 and 0 < hours <= 500")
        repeats = int(self._num(p.get("repeats", 1), "repeats"))
        every = self._num(p.get("every_h", 8.0), "interval")
        if not 1 <= repeats <= 50 or (repeats > 1 and every <= 0):
            raise ApiError("repeats must be 1-50 with a positive interval")
        try:
            r = simulate(self.reg, st, drug, [Dose(dose, i * every, infusion) for i in range(repeats)], hours, self.mods)
        except ValueError as exc:
            raise ApiError(str(exc)) from None
        step = max(1, len(r.times_h) // MAX_POINTS)
        return {"t": r.times_h[::step], "c": r.conc_mg_l[::step], "clearance_l_h": r.clearance_l_h,
                "volume_l": r.volume_l, "half_life_h": r.half_life_h, "cmax_mg_l": r.cmax_mg_l,
                "auc": r.auc_mg_h_l(), "flags": sorted(r.flags), "notes": r.notes, "drug": drug.name,
                "hypothetical": "hypothetical_inputs" in r.flags}

    def hpg(self, p: dict) -> dict:
        days = self._num(p.get("days", 120), "days")
        raise_to = self._num(p.get("added_testosterone", 0.0), "added testosterone")
        dose_days = self._num(p.get("dose_days", 60), "dose days")
        if not 1 <= days <= 1000 or raise_to < 0 or dose_days < 0:
            raise ApiError("days 1-1000, added testosterone >= 0, dose days >= 0")
        tau = _hpg.ASSUMED["tau_t_days"]
        rate = raise_to / tau  # ng/dL per day that holds testosterone raised by about raise_to
        r = _hpg.simulate(self.reg, days, dose_rate=lambda t: rate if t < dose_days else 0.0, dt=0.02)
        step = max(1, len(r.trajectory.t) // MAX_POINTS)
        return {"t": r.trajectory.t[::step], "T": r.trajectory.series("T")[::step], "LH": r.trajectory.series("LH")[::step],
                "C": r.trajectory.series("C")[::step], "reference": r.reference, "assumed": True,
                "assumptions": _hpg.ASSUMED}

    def cycle(self, p: dict | None = None) -> dict:
        days = [float(d) for d in range(1, 29)]
        prof = _cycle.profile(self.reg, days)
        return {"days": days, "series": {h: [prof[d].get(h) for d in days] for h in _cycle.HORMONES},
                "units": {h: self.reg[h].unit for h in _cycle.HORMONES}, "interpolated": True}


ROUTES_GET = {"/api/meta": "meta", "/api/cycle": "cycle"}
ROUTES_POST = {"/api/resolve": "resolve", "/api/compare": "compare", "/api/pk": "pk", "/api/hpg": "hpg"}


def make_handler(app: App):
    index = (WEB_DIR / "index.html").read_bytes()

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *a):  # quiet
            pass

        def _send(self, status: int, body: bytes, ctype: str) -> None:
            self.send_response(status)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            self.wfile.write(body)

        def _json(self, status: int, obj) -> None:
            self._send(status, json.dumps(obj, allow_nan=False).encode(), "application/json")

        def _host_ok(self) -> bool:
            # refuse requests whose Host is not loopback (DNS-rebinding protection for a local tool)
            return re.fullmatch(r"(127\.0\.0\.1|localhost)(:\d+)?", self.headers.get("Host", "")) is not None

        def do_GET(self):
            if not self._host_ok():
                return self._json(403, {"error": "forbidden host"})
            path = self.path.split("?")[0]
            if path in ("/", "/index.html"):
                return self._send(200, index, "text/html; charset=utf-8")
            if path in ROUTES_GET:
                return self._json(200, getattr(app, ROUTES_GET[path])())
            self._json(404, {"error": "not found"})

        def do_POST(self):
            if not self._host_ok():
                return self._json(403, {"error": "forbidden host"})
            path = self.path.split("?")[0]
            if path not in ROUTES_POST:
                return self._json(404, {"error": "not found"})
            try:
                n = int(self.headers.get("Content-Length", "0"))
                if n > 1_000_000:
                    raise ApiError("request too large", 413)
                payload = json.loads(self.rfile.read(n) or b"{}")
                if not isinstance(payload, dict):
                    raise ApiError("body must be a JSON object")
                self._json(200, getattr(app, ROUTES_POST[path])(payload))
            except ApiError as exc:
                self._json(exc.status, {"error": str(exc)})
            except json.JSONDecodeError:
                self._json(400, {"error": "body is not valid JSON"})
            except Exception as exc:  # a bug, not the caller's fault: say so without leaking internals
                self._json(500, {"error": f"internal error: {type(exc).__name__}"})

    return Handler


def serve(port: int = 8765, open_browser: bool = False) -> None:
    app = App()
    srv = ThreadingHTTPServer(("127.0.0.1", port), make_handler(app))
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    print(f"HumanSim UI at {url}  (Ctrl+C to stop). Local only; nothing is stored.")
    if open_browser:
        import webbrowser
        threading.Timer(0.5, webbrowser.open, args=(url,)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()
