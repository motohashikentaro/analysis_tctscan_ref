"""Independent end-to-end checks of run-wide templates and held-out timing comparisons."""
import array
import csv
import json
import math
from pathlib import Path
import random
import statistics
import subprocess
import tempfile
import ROOT

PROJECT = Path(__file__).resolve().parent.parent
ARTIFACTS = Path(tempfile.mkdtemp(prefix="template_validation_", dir=PROJECT / "build"))


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def run(source, destination, flags=(), success=True):
    plot_flags = [] if "--plot-segments" in flags else ["--plot-segments", "1"]
    cmd = [str(PROJECT / "bin/peak_timing"), str(source), "--events", "99,7", "--segments", "0:79",
           "--baseline-window-ns", "48:51", "--signal-window-ns", "51.5:56.5", *plot_flags,
           "--output-dir", str(destination), *flags]
    p = subprocess.run(cmd, text=True, capture_output=True, cwd=PROJECT)
    (ARTIFACTS / (destination.name + ".log")).write_text(p.stdout + p.stderr)
    assert (p.returncode == 0) == success, p.stdout + p.stderr
    assert destination.exists() == success


def main():
    ROOT.gROOT.SetBatch(True)
    source = ARTIFACTS / "fixture.root"
    out = ROOT.TFile(str(source), "CREATE")
    info = ROOT.TTree("scan_info", "points")
    wave = ROOT.TTree("waveforms", "shots")
    event, segment = array.array("q", [0]), array.array("i", [0])
    count = array.array("i", [160])
    x, y, z = [array.array("d", [0]) for _ in range(3)]
    ns = array.array("i", [111] * 4)
    dt = array.array("d", [1e-10] * 4)
    t0 = array.array("d", [0] * 4)
    clipped = array.array("b", [0] * 4)
    volts = [ROOT.std.vector("float")() for _ in range(4)]
    for tree in (info, wave):
        tree.Branch("event_number", event, "event_number/L")
        for name, a in (("x_um", x), ("y_um", y), ("z_um", z)):
            tree.Branch(name, a, name + "/D")
    info.Branch("nsegments", count, "nsegments/I")
    wave.Branch("segment_number", segment, "segment_number/I")
    for name, a, typ in (("nsamples", ns, "I"), ("dt_s", dt, "D"), ("t0_s", t0, "D"), ("clipped", clipped, "O")):
        wave.Branch(name, a, f"{name}[4]/{typ}")
    for c in range(4):
        wave.Branch("voltage_c" + str(c + 1), volts[c])
    samples, truth = {}, {}
    rng = random.Random(719)
    for k, identifier in enumerate((42, 7, 99)):
        event[0], x[0] = identifier, 100 * k
        info.Fill()
        for n in range(count[0]):
            segment[0] = n
            for c in range(4):
                clipped[c] = int(n == 3)
                t0[c] = (47.8 + rng.uniform(-.05, .05)) * 1e-9
                amp = .26 + c * .01 if k < 2 else .06 + c * .005
                peak = 53 + .02 * k + .008 * c + rng.gauss(0, .01)
                volts[c].clear()
                for i in range(111):
                    t = (t0[c] + i * dt[c]) * 1e9
                    d = t - peak
                    shape = math.exp(-.5 * (d / (.2 if d < 0 else .65)) ** 2)
                    v = .004 - amp * shape + rng.gauss(0, .0015)
                    if n == 5:
                        d2 = t - peak - 1
                        v -= .10 * math.exp(-.5 * (d2 / (.2 if d2 < 0 else .65)) ** 2)
                    volts[c].push_back(v)
                samples[identifier, n, c + 1] = (t0[c], list(volts[c]))
                truth[identifier, n, c + 1] = peak - .2 * math.sqrt(2 * math.log(2))
            wave.Fill()
    info.Write(); wave.Write(); out.Close()
    train_flags = ["--template-min-count", "40", "--template-max-count", "100", "--template-per-position", "60"]
    first = ARTIFACTS / "first"
    run(source, first, train_flags)
    data = rows(first / "peak_timing.csv")
    assert len(data) == 2 * 40 * 4
    assert {int(r["segment_number"]) % 2 for r in data} == {1}
    uses = rows(first / "template_uses.csv")
    assert {int(r["segment_number"]) % 2 for r in uses} == {0}
    assert {int(r["event_number"]) for r in uses} == {42, 7}, "Training must include unrequested point 42 and exclude weak point 99"
    for c in range(1, 5):
        initial = [r for r in uses if r["channel"] == str(c) and r["stage"] == "0"]
        assert len(initial) == 100
        assert sum(r["event_number"] == "42" for r in initial) == 50
    errors = []
    for r in data:
        identity = tuple(int(r[n]) for n in ("event_number", "segment_number", "channel"))
        start, v = samples[identity]
        baseline = statistics.mean(v[i] for i in range(len(v)) if 48 <= (start + i * 1e-10) * 1e9 <= 51)
        peak = min((i for i in range(len(v)) if 51.5 <= (start + i * 1e-10) * 1e9 <= 56.5), key=v.__getitem__)
        assert math.isclose(float(r["sample_minimum_ns"]), (start + peak * 1e-10) * 1e9, abs_tol=1e-11)
        if r["linear_ok"] == "1":
            level = .5 * (baseline - v[peak])
            crossings = [i for i in range(peak) if (start + i * 1e-10) * 1e9 >= 51 and baseline - v[i] < level <= baseline - v[i + 1]]
            j = crossings[-1]
            crossing = (start + j * 1e-10) * 1e9 + (level - baseline + v[j]) / (v[j] - v[j + 1]) * .1
            assert math.isclose(float(r["linear_cfd50_ns"]), crossing, abs_tol=1e-10)
        if identity[1] == 3:
            assert r["fit_reason"] == "clipped"
        if identity[1] == 5:
            assert r["fit_reason"] == "high_chi2"
        if r["common"] == "1":
            errors.append(float(r["template_cfd50_ns"]) - truth[identity])
            assert math.isclose(float(r["peak_delta_ps"]), 1000 * (float(r["sample_minimum_ns"]) - float(r["template_peak_ns"])), abs_tol=1e-9)
    assert len(errors) > .9 * len(data)
    rmse = math.sqrt(statistics.mean(e * e for e in errors))
    assert rmse < .015, rmse
    for s in rows(first / "summary.csv"):
        common = [r for r in data if r["event_number"] == s["event_number"] and r["channel"] == s["channel"] and r["common"] == "1"]
        assert int(s["n_common"]) == len(common)
        key = {"linear_cfd50": "linear_cfd50_ns", "template_cfd50": "template_cfd50_ns", "sample_minimum": "sample_minimum_ns", "template_peak": "template_peak_ns"}.get(s["metric"])
        if key:
            times = [float(r[key]) for r in common]
            assert math.isclose(float(s["mean_ns"]), statistics.mean(times), abs_tol=1e-11)
            assert math.isclose(float(s["stddev_ns"]), statistics.stdev(times), abs_tol=1e-11)
    second = ARTIFACTS / "reloaded"
    run(source, second, ["--template-input", str(first / "templates")])
    assert (first / "peak_timing.csv").read_bytes() == (second / "peak_timing.csv").read_bytes()
    # Inadequate pool must fail rather than silently widening cuts or using odd shots.
    run(source, ARTIFACTS / "insufficient", ["--template-min-count", "1000", "--template-max-count", "1000"], success=False)
    run(source, ARTIFACTS / "even_plot", ["--plot-segments", "0"], success=False)
    run(source, ARTIFACTS / "override", ["--template-input", str(first / "templates"), "--template-min-count", "40"], success=False)
    before = (first / "peak_timing.csv").read_bytes()
    cmd = [str(PROJECT / "bin/peak_timing"), str(source), "--events", "7", "--output-dir", str(first)]
    assert subprocess.run(cmd, capture_output=True).returncode != 0
    assert (first / "peak_timing.csv").read_bytes() == before
    # A complete template with a different input identity must not be silently applied to another run.
    mismatch = ARTIFACTS / "mismatch_templates"
    mismatch.mkdir()
    for c in range(1, 5):
        txt = (first / f"templates/C{c}.template").read_text()
        (mismatch / f"C{c}.template").write_text(txt.replace("ROOT_MD5:", "OTHER_MD5:", 1))
    run(source, ARTIFACTS / "mismatched_run", ["--template-input", str(mismatch)], success=False)
    for result in (first, second):
        assert (result / "COMPLETE").exists()
        json.loads((result / "conditions.json").read_text())
        pngs = list(result.rglob("*.png"))
        assert len(pngs) == 24
        for p in pngs:
            assert p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    report = {"status": "passed", "artifacts": str(ARTIFACTS), "common": len(errors), "evaluated": len(data), "toa_rmse_ps": rmse * 1000,
              "checks": "run-wide balanced even training, independent linear CFD/summary, odd evaluation, time truth, clipped/double pulse rejection, persistence equality, insufficient pool, PNGs"}
    (PROJECT / "build/template_regression_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
