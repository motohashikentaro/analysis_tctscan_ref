"""Compare the new applications with archived measurements and numerical results."""

import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parent.parent
OLD = ROOT.parent / "analysis_tctscan_old"
ARTIFACTS = Path(tempfile.mkdtemp(prefix="regression_", dir=ROOT / "build"))
REPORT = {"artifacts": str(ARTIFACTS), "checks": []}


def run(program, run_name, flags, destination, success=True):
    command = [str(ROOT / "bin" / program), str(ROOT / "data" / (run_name + ".root"))]
    command += flags + ["--output-dir", str(destination)]
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    log = ARTIFACTS / ("command_" + str(len(list(ARTIFACTS.glob("command_*.log")))) + ".log")
    log.write_text(" ".join(command) + "\n" + result.stdout + result.stderr)
    assert (result.returncode == 0) == success, str(log) + "\n" + result.stderr[-2000:]
    if success:
        assert (destination / "COMPLETE").is_file()
        metadata = json.loads((destination / "conditions.json").read_text())
        assert len(metadata["source_sha256"]) == 64
    else:
        assert not destination.exists(), "Failed analysis created a final output"


def rows(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def compare_csv(actual, expected, key, ignore=(), select=lambda row: True, fit=False):
    got = {tuple(row[column] for column in key): row for row in rows(actual)}
    wanted = {tuple(row[column] for column in key): row for row in rows(expected) if select(row)}
    assert got.keys() == wanted.keys(), f"Row identities differ: {actual}"
    maximum = {}
    for identity, row in got.items():
        for column, reference in wanted[identity].items():
            if column in ignore:
                continue
            assert column in row, f"Missing {column} in {actual}"
            value = row[column]
            if reference == value:
                continue
            try:
                number = float(value)
                expected_number = float(reference)
            except ValueError:
                raise AssertionError((actual.name, identity, column, value, reference))
            difference = abs(number - expected_number)
            maximum[column] = max(maximum.get(column, 0.0), difference)
            tolerance = 1e-10 if fit else 5e-13
            assert math.isclose(number, expected_number, rel_tol=tolerance, abs_tol=1e-18), (
                actual.name, identity, column, value, reference, difference
            )
    REPORT["checks"].append({"csv": str(actual), "reference": str(expected), "rows": len(got), "max_absolute_differences": maximum})


def check_pngs(directory):
    images = list(directory.rglob("*.png"))
    assert images
    for path in images:
        assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", str(path)
    REPORT["checks"].append({"png_directory": str(directory), "count": len(images)})


def check_conversion():
    import ROOT as root

    fixture = ARTIFACTS / "conversion"
    fixture.mkdir()
    base = fixture / "scan_fixture"
    traces = fixture / "scan_fixture_2dYX"
    traces.mkdir()
    original_run = "scan_20261021_154341"
    with (ROOT / "data" / (original_run + ".csv")).open(newline="") as source:
        table = list(csv.reader(source))
    header = next(row for row in table if row and row[0] == "index")
    wanted = [row for row in table if row and row[0] in {"25", "84"}]
    assert len(wanted) == 2
    with base.with_suffix(".csv").open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(header)
        writer.writerows(wanted)
    for event in (25, 84):
        for channel in range(1, 5):
            name = f"C{channel}--scan--{event:05d}.trc"
            (traces / name).symlink_to((ROOT / "data" / (original_run + "_2dYX") / name).resolve())
    converted = fixture / "converted.root"
    result = subprocess.run([str(ROOT / "bin/trc_to_root"), str(base), str(converted)], text=True, capture_output=True)
    (fixture / "conversion.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stderr
    original = root.TFile.Open(str(ROOT / "data" / (original_run + ".root")))
    rebuilt = root.TFile.Open(str(converted))
    original_wave = original.Get("waveforms")
    rebuilt_wave = rebuilt.Get("waveforms")
    original_info = original.Get("scan_info")
    rebuilt_info = rebuilt.Get("scan_info")
    assert rebuilt_info.GetEntries() == 2 and rebuilt_wave.GetEntries() == 4000
    raw_samples = 0
    for position, event in enumerate((25, 84)):
        original_info.GetEntry(event)
        rebuilt_info.GetEntry(position)
        for name in ("event_number", "x_um", "y_um", "z_um", "nsegments", "status", "stage_timestamp"):
            assert str(getattr(original_info, name)) == str(getattr(rebuilt_info, name)), name
        for name in ("nsamples", "gain_V_per_count", "vertical_offset_V", "dt_s", "comm_type", "nominal_bits", "scope_timestamp_bytes"):
            assert list(getattr(original_info, name)) == list(getattr(rebuilt_info, name)), name
        for segment in (0, 1, 999, 1000, 1999):
            original_wave.GetEntry(event * 2000 + segment)
            rebuilt_wave.GetEntry(position * 2000 + segment)
            for name in ("event_number", "segment_number", "x_um", "y_um", "z_um"):
                assert getattr(original_wave, name) == getattr(rebuilt_wave, name), name
            for name in ("dt_s", "t0_s", "trigger_elapsed_s", "clipped", "nsamples"):
                assert list(getattr(original_wave, name)) == list(getattr(rebuilt_wave, name)), name
            for channel in range(1, 5):
                for prefix in ("raw_c", "voltage_c"):
                    actual = list(getattr(rebuilt_wave, prefix + str(channel)))
                    assert actual == list(getattr(original_wave, prefix + str(channel))), prefix
                    raw_samples += len(actual)
    rebuilt.Close()
    original.Close()
    assert not converted.with_suffix(".root.partial").exists()
    REPORT["checks"].append({"conversion": "2 positions, 4000 shots, raw/voltage/timing/scaling metadata", "compared_samples": raw_samples})


def main():
    average_cases = {
        "scan_20261021_154341": [0, 84, 168],
        "scan_20261021_164141": [0, 80, 160],
    }
    for run_name, events in average_cases.items():
        destination = ARTIFACTS / ("average_" + run_name)
        flags = ["--events", ",".join(map(str, events)), "--segments", "0:1999"]
        run("average_waveforms", run_name, flags, destination)
        for event in events:
            reference = OLD / "results/average_waveforms" / run_name / "root_style" / f"event_{event:05d}_seg0000_n2000_trigger.csv"
            if not reference.exists():
                reference = reference.parent.parent / reference.name
            compare_csv(destination / "waveforms" / f"event_{event:05d}.csv", reference,
                        ("event_number", "channel", "sample_index"))
        check_pngs(destination)
        before = hashlib.sha256((destination / "averages.root").read_bytes()).hexdigest()
        result = subprocess.run([str(ROOT / "bin/average_waveforms"), str(ROOT / "data" / (run_name + ".root")),
                                 *flags, "--output-dir", str(destination)], capture_output=True)
        assert result.returncode != 0
        assert before == hashlib.sha256((destination / "averages.root").read_bytes()).hexdigest()
        print("PASS average", run_name, flush=True)

    windows = ["--baseline-window-ns", "48:51", "--signal-window-ns", "51.5:56.5"]
    for run_name in ("scan_20261021_164141", "scan_20261022_112522"):
        destination = ARTIFACTS / ("amplitude_" + run_name)
        run("scan_amplitude", run_name, windows, destination)
        compare_csv(destination / "amplitudes.csv", OLD / "results/scan_amplitude" / run_name / "amplitudes.csv",
                    ("event_number", "channel"))
        check_pngs(destination)
        print("PASS amplitude", run_name, flush=True)

    xy_run = "scan_20261021_154341"
    xy_output = ARTIFACTS / "xy"
    run("scan_xy_maps", xy_run, windows, xy_output)
    for filename, identity in (("channel_values.csv", ("event_number", "channel")), ("position_values.csv", ("event_number",))):
        compare_csv(xy_output / filename, OLD / "results/scan_xy_maps" / xy_run / filename, identity)
    check_pngs(xy_output)
    import ROOT as root
    sparse_output = ARTIFACTS / "sparse_sample_index"
    run("average_waveforms", xy_run, ["--events", "84", "--segments", "0,10,50", "--alignment", "sample-index"], sparse_output)
    input_file = root.TFile.Open(str(ROOT / "data" / (xy_run + ".root")))
    input_tree = input_file.Get("waveforms")
    samples = {}
    starts = {}
    for segment in (0, 10, 50):
        input_tree.GetEntry(84 * 2000 + segment)
        for channel in range(1, 5):
            values = list(getattr(input_tree, "voltage_c" + str(channel)))
            samples.setdefault(channel, []).append(values)
            starts.setdefault(channel, []).append(input_tree.t0_s[channel - 1])
    for row in rows(sparse_output / "waveforms/event_00084.csv"):
        channel = int(row["channel"])
        sample = int(row["sample_index"])
        expected_voltage = sum(waveform[sample] for waveform in samples[channel]) / 3
        expected_time = sum(starts[channel]) / 3 + sample * input_tree.dt_s[channel - 1]
        assert math.isclose(float(row["mean_voltage_V"]), expected_voltage, rel_tol=1e-12, abs_tol=1e-18)
        assert math.isclose(float(row["time_s"]), expected_time, rel_tol=1e-12, abs_tol=1e-22)
        assert row["n_used"] == "3"
    input_file.Close()
    check_pngs(sparse_output)
    REPORT["checks"].append({"sparse_sample_index": "independent raw-shot comparison passed"})
    reference_file = root.TFile.Open(str(OLD / "results/scan_xy_maps" / xy_run / "maps.root"))
    actual_file = root.TFile.Open(str(xy_output / "maps.root"))
    maps = 0
    for key in reference_file.GetListOfKeys():
        reference = reference_file.Get(key.GetName())
        if not reference.InheritsFrom("TH2"):
            continue
        actual = actual_file.Get(key.GetName())
        assert actual and reference.GetNcells() == actual.GetNcells()
        assert [reference.GetBinContent(i) for i in range(reference.GetNcells())] == [actual.GetBinContent(i) for i in range(actual.GetNcells())]
        for axis in ("GetXaxis", "GetYaxis"):
            expected_axis = getattr(reference, axis)()
            actual_axis = getattr(actual, axis)()
            assert [expected_axis.GetBinLowEdge(i) for i in range(1, expected_axis.GetNbins() + 2)] == [actual_axis.GetBinLowEdge(i) for i in range(1, actual_axis.GetNbins() + 2)]
        maps += 1
    assert maps == 11
    actual_file.Close()
    reference_file.Close()
    REPORT["checks"].append({"xy_root_maps": maps, "bins_and_edges": "exact match"})
    print("PASS XY CSV and 11 ROOT maps", flush=True)

    subprocess.run(["python3", "tests/regression_template.py"], cwd=ROOT, check=True)
    REPORT["checks"].append({"template_timing": "independent synthetic e2e; see build/template_regression_report.json"})

    negatives = [
        ("average_waveforms", xy_run, ["--events", "999999"]),
        ("average_waveforms", xy_run, ["--events", "0", "--segments", "0,0"]),
        ("average_waveforms", xy_run, ["--events", "0", "--segments", "2000"]),
        ("average_waveforms", xy_run, ["--events", "0", "--segments", "0", "--nsegments", "1"]),
        ("average_waveforms", xy_run, ["--events", "0", "--tag", "../bad"]),
        ("scan_amplitude", xy_run, windows),
        ("scan_xy_maps", "scan_20261021_164141", windows),
        ("peak_timing", xy_run, ["--events", "0", "--segments", "0", "--plot-segments", "1", *windows]),
    ]
    for index, (program, run_name, flags) in enumerate(negatives):
        run(program, run_name, flags, ARTIFACTS / ("rejected_" + str(index)), success=False)
    REPORT["checks"].append({"rejected_invalid_commands": len(negatives), "overwrite": "refused, saved ROOT unchanged"})
    check_conversion()
    REPORT["status"] = "passed"
    (ROOT / "build/regression_report.json").write_text(json.dumps(REPORT, indent=2) + "\n")
    print("PASS all regression checks; artifacts:", ARTIFACTS, flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        REPORT["status"] = "failed"
        (ROOT / "build/regression_report.json").write_text(json.dumps(REPORT, indent=2) + "\n")
        raise
