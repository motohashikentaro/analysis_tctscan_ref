"""End-to-end CFD map check on a non-square XY fixture, independent of C++ estimators."""
import array
import csv
import json
import math
from pathlib import Path
import statistics
import subprocess
import tempfile

import ROOT

PROJECT = Path(__file__).resolve().parent.parent
ARTIFACTS = Path(tempfile.mkdtemp(prefix="propagation_validation_", dir=PROJECT / "build"))


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    ROOT.gROOT.SetBatch(True)
    source = ARTIFACTS / "fixture.root"
    output = ARTIFACTS / "result"
    file = ROOT.TFile(str(source), "CREATE")
    info = ROOT.TTree("scan_info", "positions")
    wave = ROOT.TTree("waveforms", "shots")
    event, segment = array.array("q", [0]), array.array("i", [0])
    nsegments = array.array("i", [24])
    x, y, z = [array.array("d", [0]) for _ in range(3)]
    ns = array.array("i", [102] * 4)
    dt = array.array("d", [1e-10] * 4)
    t0 = array.array("d", [0] * 4)
    clipped = array.array("b", [0] * 4)
    voltage = [ROOT.std.vector("float")() for _ in range(4)]
    for tree in (info, wave):
        tree.Branch("event_number", event, "event_number/L")
        for name, a in (("x_um", x), ("y_um", y), ("z_um", z)):
            tree.Branch(name, a, name + "/D")
    info.Branch("nsegments", nsegments, "nsegments/I")
    wave.Branch("segment_number", segment, "segment_number/I")
    for name, a, fmt in (("nsamples", ns, "I"), ("dt_s", dt, "D"),
                         ("t0_s", t0, "D"), ("clipped", clipped, "O")):
        wave.Branch(name, a, name + "[4]/" + fmt)
    for c in range(4):
        wave.Branch("voltage_c" + str(c + 1), voltage[c])
    samples = {}
    for k, identifier in enumerate((42, 7, 99, 10, 3, 81)):
        event[0] = identifier
        x[0], y[0] = (-100, 50)[k % 2], (-200, 0, 300)[k // 2]
        info.Fill()
        leading = k % 4
        for n in range(24):
            segment[0] = n
            for c in range(4):
                t0[c] = (47.91 + n * .003 + c * .00007) * 1e-9
                clipped[c] = int(k == 0 and n == 0)
                amplitude = .030 if c == leading else .010
                if c == (leading + 1) % 4:
                    amplitude = .035 if n % 2 == 0 else .015
                if k == 5 or (k == 0 and n == 1):
                    amplitude = .001
                center = 53.4 + .02 * k + .04 * c + .005 * (n % 5 - 2)
                voltage[c].clear()
                times = [t0[c] * 1e9 + i * dt[c] * 1e9 for i in range(102)]
                for t in times:
                    voltage[c].push_back(-.04 * c - amplitude * max(0., 1 - abs(t - center) / .8))
                samples[identifier, n, c + 1] = (times, list(voltage[c]))
            wave.Fill()
    info.Write()
    wave.Write()
    ROOT.TNamed("source_trace_directory", "fixture_2dYX").Write()
    file.Close()

    command = [str(PROJECT / "bin/propagation_xy_maps"), str(source), "--output-dir", str(output)]
    result = subprocess.run(command, capture_output=True, text=True)
    (ARTIFACTS / "run.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stderr
    assert (output / "COMPLETE").exists()
    conditions = json.loads((output / "conditions.json").read_text())
    assert conditions["cfd_fraction"] == .5 and conditions["min_valid"] == 20
    rows = read_csv(output / "propagation_delay.csv")
    assert len(rows) == 6 * 24
    grouped = {}
    for row in rows:
        event_id, n, c = (int(row[key]) for key in ("event_number", "segment_number", "channel"))
        times, values = samples[event_id, n, c]
        if row["valid"] != "1":
            assert row["reason"] in ("clipped", "low_amplitude")
            continue
        baseline = statistics.mean(v for t, v in zip(times, values) if 48 <= t <= 51)
        peak = min((i for i, t in enumerate(times) if 51.5 <= t <= 56.5), key=values.__getitem__)
        amplitude = baseline - values[peak]
        level = .5 * amplitude
        brackets = [i for i in range(peak) if times[i] >= 51 and
                    baseline - values[i] < level <= baseline - values[i + 1]]
        j = brackets[-1]
        expected = times[j] + (level - (baseline - values[j])) / (values[j] - values[j + 1]) * (times[j + 1] - times[j])
        assert math.isclose(float(row["delay_ns"]), expected, abs_tol=1e-10)
        grouped.setdefault(event_id, []).append(expected)
    summaries = read_csv(output / "position_values.csv")
    root = ROOT.TFile(str(output / "maps.root"), "READ")
    assert root.Get("propagation_delay_xy").GetNbinsX() == 2
    assert root.Get("propagation_delay_xy").GetNbinsY() == 3
    identifiers = (42, 7, 99, 10, 3, 81)
    for row in summaries:
        identifier = int(row["event_number"])
        k = identifiers.index(identifier)
        c = int(row["leading_channel"])
        if k != 5:
            assert c == k % 4 + 1, "Fixed mean-amplitude leader must override per-shot leaders and raw minima"
        xbin, ybin = int(row["ix"]) + 1, int(row["iy"]) + 1
        selected = [r for r in rows if int(r["event_number"]) == identifier]
        assert {int(r["channel"]) for r in selected} == {c}, "Leader changed within position"
        values = grouped.get(identifier, [])
        assert int(row["n_valid"]) == len(values)
        assert math.isclose(float(row["valid_fraction"]), len(values) / 24)
        assert root.Get("leading_channel_xy").GetBinContent(xbin, ybin) == c
        valid = len(values) >= 20
        assert root.Get("valid_mask_xy").GetBinContent(xbin, ybin) == valid
        if values:
            for column, expected in (("mean_delay_ns", statistics.mean(values)),
                                      ("median_delay_ns", statistics.median(values)),
                                      ("stddev_ns", statistics.stdev(values)),
                                      ("stderr_ns", statistics.stdev(values) / math.sqrt(len(values)))):
                assert math.isclose(float(row[column]), expected, abs_tol=1e-10), column
            assert math.isclose(root.Get("propagation_delay_xy").GetBinContent(xbin, ybin), statistics.mean(values), abs_tol=1e-10)
        else:
            assert row["mean_delay_ns"] == "" and row["map_valid"] == "0"
    root.Close()
    assert sum(int(r["count"]) for r in read_csv(output / "failure_counts.csv")) == 26
    for png in output.glob("*.png"):
        assert png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert len(list(output.glob("*.png"))) == 10
    # Entirely masked maps must remain drawable and explicitly masked in ROOT.
    masked = ARTIFACTS / "masked"
    result = subprocess.run(command[:-1] + [str(masked), "--min-valid", "25"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert all(r["map_valid"] == "0" for r in read_csv(masked / "position_values.csv"))
    for flags in (("--cfd-fraction", "1"), ("--min-valid", "0"), ("--min-snr", "nan"),
                  ("--alignment", "sample-index")):
        rejected = ARTIFACTS / ("rejected_" + flags[0][2:])
        result = subprocess.run(command[:-1] + [str(rejected), *flags], capture_output=True)
        assert result.returncode != 0 and not rejected.exists()
    assert subprocess.run(command, capture_output=True).returncode != 0, "Existing output must not be overwritten"
    print("PASS propagation end-to-end: independent CFD/statistics, fixed amplitude leader, non-square XY, masks, cuts, metadata, PNG/ROOT, CLI and overwrite checks")
    print("Artifacts:", ARTIFACTS)


if __name__ == "__main__":
    main()
