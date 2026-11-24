"""Check the diagnostic end to end using independent residual/chi2 accounting."""
import array
import bisect
import csv
import json
import math
from pathlib import Path
import random
import subprocess
import tempfile
import ROOT

PROJECT = Path(__file__).resolve().parents[2]
ARTIFACTS = Path(tempfile.mkdtemp(prefix='initial_fit_check_', dir=PROJECT / 'build'))


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    ROOT.gROOT.SetBatch(True)
    source = ARTIFACTS / 'fixture.root'
    output = ARTIFACTS / 'result'
    file = ROOT.TFile(str(source), 'CREATE')
    info, wave = ROOT.TTree('scan_info', ''), ROOT.TTree('waveforms', '')
    event, segment, count = array.array('q', [0]), array.array('i', [0]), array.array('i', [120])
    xyz = [array.array('d', [0]) for _ in range(3)]
    ns, dt, t0, clip = array.array('i', [111]*4), array.array('d', [1e-10]*4), array.array('d', [0]*4), array.array('b', [0]*4)
    voltage = [ROOT.std.vector('float')() for _ in range(4)]
    for tree in (info, wave):
        tree.Branch('event_number', event, 'event_number/L')
        for name, value in zip(('x_um', 'y_um', 'z_um'), xyz):
            tree.Branch(name, value, name+'/D')
    info.Branch('nsegments', count, 'nsegments/I')
    wave.Branch('segment_number', segment, 'segment_number/I')
    for name, value, typ in [('nsamples', ns, 'I'), ('dt_s', dt, 'D'), ('t0_s', t0, 'D'), ('clipped', clip, 'O')]:
        wave.Branch(name, value, f'{name}[4]/{typ}')
    for c in range(4):
        wave.Branch('voltage_c'+str(c+1), voltage[c])
    rng = random.Random(729)
    raw = {}
    for e in (7, 42):
        event[0] = e
        xyz[0][0] = 100 if e == 42 else 0
        info.Fill()
        for n in range(120):
            segment[0] = n
            start = 47.8+rng.uniform(-.05, .05)
            peak = 53+rng.gauss(0, .02)
            values = []
            for i in range(111):
                time = start+i*.1
                x = time-peak
                v = .004-.27*math.exp(-.5*(x/(.2 if x < 0 else .65))**2)+rng.gauss(0, .003)
                if n % 6 == 0:
                    x = time-peak-1
                    v -= .12*math.exp(-.5*(x/(.2 if x < 0 else .65))**2)
                values.append(v)
            for c in range(4):
                t0[c] = start*1e-9
                voltage[c].clear()
                for v in values:
                    voltage[c].push_back(v)
            raw[e, n] = (t0[0], list(voltage[0]))
            wave.Fill()
    info.Write(); wave.Write(); file.Close()
    command = [str(PROJECT/'build/check_initial_fit'), str(source), '--channels', '1', '--template-min-count', '20',
               '--output-dir', str(output)]
    result = subprocess.run(command, cwd=PROJECT, capture_output=True, text=True)
    (ARTIFACTS/'command.log').write_text(result.stdout+result.stderr)
    assert result.returncode == 0, result.stdout+result.stderr
    status = rows(output/'channel_status.csv')[0]
    fits = rows(output/'C1/first_fits.csv')
    assert int(status['n_selected']) == len(fits) == 120
    assert int(status['n_finite']) == 120
    assert int(status['n_pass_original_cut']) < 120, 'Rejected fits must survive in diagnostic output'
    assert any(r['reason'] == 'high_chi2' for r in fits)
    assert all(int(r['segment_number']) % 2 == 0 for r in fits)
    for r in fits:
        expected = float(r['sse_v2']) / (float(r['noise_rms_mv'])*1e-3)**2
        assert math.isclose(float(r['chi2']), expected, rel_tol=1e-12)
        assert math.isclose(float(r['chi2_ndf']), expected/int(r['ndf']), rel_tol=1e-12)
    bins = rows(output/'C1/chi2_bins.csv')
    assert sum(int(b['n']) for b in bins) == 120
    assert sum(int(b['n_pass_original_cut']) for b in bins) == int(status['n_pass_original_cut'])
    for b in bins:
        matches = [r for r in fits if float(b['lower_inclusive']) <= float(r['chi2_ndf']) < float(b['upper_exclusive'])]
        assert len(matches) == int(b['n'])
    examples = rows(output/'C1/examples.csv')
    ids = [(int(r['event_number']), int(r['segment_number'])) for r in examples]
    assert len(ids) == len(set(ids))
    for ex in examples:
        b = next(b for b in bins if b['bin'] == ex['bin'])
        matches = sorted((r for r in fits if float(b['lower_inclusive']) <= float(r['chi2_ndf']) < float(b['upper_exclusive'])),
                         key=lambda r: (float(r['chi2_ndf']), int(r['event_number']), int(r['segment_number'])))
        expected = matches[int(ex['rank_in_bin'])]
        assert (expected['event_number'], expected['segment_number']) == (ex['event_number'], ex['segment_number'])
        assert (output/ex['png']).is_file()
    template = rows(output/'C1/initial_template.csv')
    tx, ty = [float(r['time_ns']) for r in template], [float(r['normalized_signal']) for r in template]
    fit_by_id = {(int(r['event_number']), int(r['segment_number'])): r for r in fits}
    sums = {}
    for r in rows(output/'C1/example_samples.csv'):
        identity = int(r['event_number']), int(r['segment_number'])
        start, values = raw[identity]
        i = int(r['sample_index']);f = fit_by_id[identity]
        assert math.isclose(float(r['time_ns']), (start+i*1e-10)*1e9, abs_tol=1e-11)
        assert math.isclose(float(r['voltage_mv']), values[i]*1e3, abs_tol=1e-10)
        if r['in_fit_window'] == '1':
            x = float(r['time_ns'])-float(f['fit_cfd50_ns'])
            j = min(len(tx)-2, bisect.bisect_right(tx, x)-1)
            value = ty[j]+(x-tx[j])/(tx[j+1]-tx[j])*(ty[j+1]-ty[j])
            model = float(f['fit_baseline_mv'])-float(f['fit_amplitude_mv'])*value
            assert math.isclose(float(r['fit_mv']), model, abs_tol=1e-8)
            z = (values[i]*1e3-model)/float(f['noise_rms_mv'])
            assert math.isclose(float(r['chi2_contribution']), z*z, rel_tol=1e-7, abs_tol=1e-8)
            sums[identity] = sums.get(identity, 0)+z*z
    for identity, total in sums.items():
        assert math.isclose(total, float(fit_by_id[identity]['chi2']), rel_tol=1e-9)
    rf = ROOT.TFile.Open(str(output/'diagnostics.root'))
    assert rf.Get('C1/chi2_ndf_linear').GetEntries() == 120
    assert rf.Get('C1/chi2_log').GetEntries() == 120
    rf.Close()
    before = (output/'C1/first_fits.csv').read_bytes()
    assert subprocess.run(command, cwd=PROJECT, capture_output=True).returncode != 0
    assert before == (output/'C1/first_fits.csv').read_bytes()
    manifest = json.loads((output/'conditions.json').read_text())
    assert manifest['original_chi2_ndf_cut'] == 3
    assert (output/'COMPLETE').exists()
    for p in output.rglob('*.png'):
        assert p.read_bytes()[:8] == b'\x89PNG\r\n\x1a\n'
    report = {'status': 'passed', 'artifacts': str(ARTIFACTS), 'selected': 120, 'examples': len(examples),
              'checks': 'finite high-chi2 retained; independent bins, example ranks, sample/model/residual chi2, ROOT counts; overwrite protection'}
    (PROJECT/'build/check_initial_fit_test_report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
