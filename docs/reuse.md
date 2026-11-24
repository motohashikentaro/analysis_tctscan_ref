# 読み込み・平均波形・特徴量の再利用

共通処理は `tctscan` 名前空間の関数・クラスとして公開しています。新しい解析は `app/` に追加し、Makefileの `APPS` に実行ファイル名を追加すると、共通実装とリンクできます。

## 元波形を読む

```cpp
#include "rootfile_reader.h"
#include "analysis_support.h"
#include "waveform_analyzer.h"

tctscan::RootfileReader reader(input_path);

const auto& points = reader.GetScanPoints();
const auto& point = reader.GetScanPoint(84);
std::vector<int> segments = tctscan::SelectSegments(point.nsegments, {{}, 0, 500});

std::vector<tctscan::ShotData> shots = reader.ReadShots(84, segments);
tctscan::ShotData one_shot = reader.ReadShot(84, 10);
```

`event_number` は測定点番号、`segment_number` は0始まりのショット番号です。ROOTのentry番号との対応はReaderが管理します。複数点を使うときは、最初に `reader.IndexEvents({0, 84, 168});` を呼ぶと、一度の索引作成で準備できます。

`ReadShot` / `ReadShots` が返す波形は自分のデータを保持します。次の `GetEntry` で内容が入れ替わるROOTバッファへの参照ではありません。

`shot.waveforms[0]` がC1、`[3]` がC4です。`Waveform` は `time_s`、`voltage_v`、`dt_s`、`clipped` を持ちます。位置は `position.x_um` などで参照します。

## 平均して特徴量を求める

```cpp
tctscan::AverageCondition average_condition;
average_condition.alignment = tctscan::Alignment::kTrigger;

std::array<tctscan::AverageWaveform, tctscan::kNChannels> averages =
    tctscan::AverageWaveforms(shots, average_condition);

tctscan::FeatureCondition feature_condition;
feature_condition.baseline = {48.0, 51.0};
feature_condition.signal = {51.5, 56.5};

tctscan::WaveformFeature from_average =
    tctscan::ExtractWaveformFeature(averages[0].waveform, feature_condition);

tctscan::WaveformFeature from_shot =
    tctscan::ExtractWaveformFeature(one_shot.waveforms[0], feature_condition);
```

平均化の戻り値には、波形本体、使用数 `n_used`、clipped数 `n_clipped`、入力開始時刻の範囲が含まれます。`ShotData` の集合を渡す版は、測定点が混ざっていないことも検査します。1チャンネルだけの `std::vector<Waveform>` を受け取る版もあります。

`ExtractWaveformFeature` の入力は元波形・平均波形で共通です。戻り値には `amplitude_mv`、`baseline_mv`、`peak_voltage_mv`、`peak_time_ns` と品質フラグが含まれます。

## 平均波形を保存して再利用する

`average_waveforms` は `averages.root` を標準出力します。公開関数で直接保存・読み込みもできます。

```cpp
#include "rootfile_writer.h"

tctscan::AveragedEvent event;
event.point = point;
event.segments = segments;
event.waveforms = averages;
std::vector<tctscan::AveragedEvent> events{event};
tctscan::WriteAverageRootFile(events, source_directory, conditions_json, output_path);

std::vector<tctscan::AveragedEvent> restored_events =
    tctscan::ReadAverageRootFile(average_root_path);

auto feature = tctscan::ExtractWaveformFeature(restored_events[0].waveforms[0].waveform, feature_condition);
```

`AveragedEvent` は測定点、使用したショット番号、4chの平均波形をまとめた構造体です。ROOTでは `averages` TTreeに測定点×チャンネルで保存し、時間・電圧はdoubleで保持します。元TRCや全ショットのROOTを再読込せず、保存した平均波形から別の特徴量を求められます。

## TRCと測定点CSVを読む

```cpp
#include "trc_reader.h"
#include "scan_csv_reader.h"

tctscan::TraceData trace = tctscan::ReadTrc(trace_path, 1); // C1
std::vector<tctscan::ScanCsvRow> rows = tctscan::ReadScanCsv(csv_path);

int index = segment * trace.nsamples + sample;
double voltage_v = trace.raw[index] * trace.gain_v_per_count - trace.vertical_offset_v;
double time_s = trace.t0_s[segment] + sample * trace.dt_s;
```

TRCのバイナリ解釈とCSVの列解釈は別の関数です。変換全体は `ConvertScanToRoot` で実行できます。TRCの対応範囲・同期検査は旧実装を引き継いでいます。

## 描画・保存の分担

- `waveform_analyzer`：平均・ベースライン・波高。ROOTファイルを開かず、計算結果を返します。
- `peak_timing_analyzer`：1波形の最小時刻・フィットと品質判定。
- `scan_analyzer`：座標の整列、格子、チャンネル順位。
- `scan_analysis`：Readerと上の計算を組み合わせたスキャン解析。
- `plot_supporter`：色、余白、軸、文字の共通設定。
- `waveform_plotter` / `scan_plotter`：計算済みの結果をROOTで描画。
- `result_writer`：保存先、条件、CSVと正常終了処理。

計算条件は `struct` で受け取ります。計算関数はコマンドライン引数や保存ディレクトリを参照しません。

## テンプレートフィット（2026-11-24）

`template_timing.h` はROOTの読込・描画に依存しないAPIです。負極性の元サンプル波形を入力し、
`V(t)=b-A*T(t-tau)`、`T` のピーク=1、立ち上がりCFD50=0として扱います。

```cpp
#include "template_timing.h"

// 既存テンプレートの再利用。元の時間軸を渡し、位置ごとの遅延を引かない。
auto pulse_template = tctscan::LoadPulseTemplate("templates/C2.template");
tctscan::PropagationCondition quality; // amplitude >=8 mV, S/N >=5, fraction=0.5
const auto result = tctscan::FitTemplate(shot.waveforms[1], pulse_template, quality);
if (result.fit_ok) {
    double toa_ns = result.toa_ns;             // フィット関数のCFD50
    double peak_ns = result.peak_time_ns;       // toa + template.peak_time_ns
    double amplitude_v = result.amplitude_v;
    (void)toa_ns; (void)peak_ns; (void)amplitude_v;
}
```

`FitTemplate` はチャンネルのパルスが同じ形状を持つという仮定でのみ有効です。
`fit_ok=false` の場合も、最小サンプル、線形CFD、探索できたフィット、残差があれば保持します。
`high_chi2` の数値を成功値として集計しないでください。`linear.valid` は線形CFDの品質判定です。
ノイズがゼロの場合はchi2を定義できないため `undefined_noise` として失敗します。
誤差共分散や時刻の推定誤差は現時点で返しません。サンプル間のノイズ相関を補正していないchi2です。

新規作成は `MakeTemplateCandidate` → `SelectTemplateCandidates` → `BuildPulseTemplate` →
`SavePulseTemplate` です。`MakeTemplateCandidate` は偶数segmentか、振幅・S/N・飽和・窓・pol2を検査します。
`SelectTemplateCandidates` は測定点別の候補群をシャッフルし、ラウンドロビンで上限まで採用します。
メモリを抑えるため、呼出側は各位置の全候補を一様なreservoir samplingで最大100本に絞って渡します。
`peak_timing` がこのROOT入力からの一連の処理の実行例です。選択をラン全体で確定してから、フィット対象点を絞ります。

`PulseTemplate` は作成条件、入力識別子、チャンネル、各候補のショットID・更新段階・採否を保持します。
保存形式は `TCTSCAN_TEMPLATE_V1` のテキストです。読込は規格化・格子・ピーク・偶数ショット履歴を検査し、
外挿は許しません。アプリではさらにROOT内容のMD5・ch・特徴量窓の一致を検査します。
独自の呼出側でも入力識別子とchの照合は呼出側の責任です。

独立した分解能評価は作成に使わない奇数ショットで行います。API自体は任意のショットをフィットできるため、
作成ショットを評価へ混ぜない責任は呼出側にあります。偶数ショットに後からフィットした結果を独立評価とは扱いません。
ラン共通形状が有効かは位置・振幅別の残差・chi2・有効率で確認してください。

初回の残差を調べる場合は `BuildInitialPulseTemplate` を使えます。同じ候補・条件で通常ビルダーのstage=0だけを
実行し、更新前に返します。通常の `FitTemplate` はhigh_chi2でも求まった数値を保持するため、
初回分布は `fit_ok` で絞らず `std::isfinite(result.chi2_ndf)` で集め、数値なしの失敗を別集計してください。
診断用プログラム・出力・検証は `checks/template_fit/README.md` を参照してください。
