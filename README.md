# analysis_tctscan

TRCと測定点CSVをROOTへ変換し、平均波形、1D波高、XYマップ、個別ショットの時刻を解析します。

元のプロジェクトは `../analysis_tctscan_old/` に保存しています。新しい `data/` は実ディレクトリで、既存の各測定ファイル・トレースディレクトリをバックアップ側へのシンボリックリンクで参照します。新しい変換結果は新しい側へ保存できます。`review` は旧検証記録へのリンクです。参照元として `_old` も保持してください。

## ビルドと検証

```bash
conda activate root_env
make -j2
make test
make test-template
```

`root-config` が示すコンパイラ・C++規格・ROOTライブラリを使用します。ビルドはプロジェクトのルートで行います。

既存の実データ・数値結果との比較は次のコマンドです。

```bash
make regression
```

`../analysis_tctscan_old` と既存3Runを読み、新しい検証結果を `build/regression_*/` に保存します。要約は `build/regression_report.json`、各コマンドのログは検証ディレクトリに残ります。

## 使い方

平均波形：

```bash
./bin/average_waveforms data/scan_20261021_154341.root \
    --events 0,84,168 --segments 0:1999 --tag all_shots_trigger
```

1D波高：

```bash
./bin/scan_amplitude data/scan_20261021_164141.root \
    --baseline-window-ns 48:51 --signal-window-ns 51.5:56.5 \
    --tag standard
```

XYマップ：

```bash
./bin/scan_xy_maps data/scan_20261021_154341.root \
    --baseline-window-ns 48:51 --signal-window-ns 51.5:56.5 \
    --tag standard
```

leading chのpropagation delayマップ：

```bash
./bin/propagation_xy_maps data/scan_20261021_154341.root \
    --cfd-fraction 0.5 --min-amplitude-mv 8 --min-snr 5 --min-valid 20 \
    --tag cfd50
```

各位置で、トリガー基準で平均した4ch波形の `baseline - minimum` が最大のchを選びます。
そのchを位置ごとに固定し、個別ショットの立ち上がりCFD時刻を計算します。
トリガーからの到達遅延をpropagation delayとし、基準位置との差やch間の固定遅延補正は行いません。
既存のamplitudeマップの順位（最小電圧順）とは区別しています。同じamplitudeなら小さいch番号を優先します。

表示値は有効ショットの到達時刻の平均［ns］です。既定の品質条件は振幅8 mV以上・S/N 5以上、
clippedでない・ピークが信号窓端でない・CFD交点が存在することです。各閾値とCFD割合は上の引数で変更できます。
有効ショットが `--min-valid` 未満の位置は灰色で表示します。平均化はtrigger方式に固定します。
`--segments`、時間窓、`--tag`、`--output-dir` は既存解析と共通です。

出力は `results/<Run>/propagation_xy_maps/<tag>/` です。

- `propagation_delay.csv`：選んだchの各ショットの時刻、振幅、ノイズ、S/N、採否・理由。
- `position_values.csv`：leading/subleading ch、amplitude差、同値フラグ、件数、平均・中央値・標準偏差・標準誤差［ns］。
- `failure_counts.csv`：位置別の除外理由と件数。
- `maps.root` と10枚のPNG：平均・中央値の到達遅延、標準偏差・標準誤差［ps］、有効率・件数、leading ch・振幅・振幅差、有効マスク。
- `conditions.json`、`COMPLETE`：条件・定義・ソース識別情報と正常終了の記録。

ROOT内のマスクされた時間ビンは0を仮置きするため、必ず `valid_mask_xy` と併せて使用してください。
標準偏差・標準誤差には2件以上が必要です。`--min-valid 1` の場合も1件だけの位置ではこれらを描画しません。
CSVは件数が少なくても計算可能な統計値を保存し、`map_valid` で描画対象かを区別します。
leading chの切替による段差には読み出し系のch間遅延差も含まれます。

既知波形による検証は `make test`、非正方XY入力からCSV・ROOT・PNGまでの独立照合は `make test-propagation` で実行できます。

個別ショットの時刻解析（テンプレートCFD50と線形CFD50の比較）：

```bash
./bin/peak_timing data/scan_20261021_164141.root \
    --events 80,110,140 --segments 0:1999 --channels 2,3,4 \
    --plot-segments 1,11,51 --tag template_comparison
```

テンプレートは評価対象点によらずラン全体の偶数ショットから作り、比較は指定範囲の奇数ショットのみで行います。
`--plot-segments` にも奇数を指定してください。通常のサンプル・区分線形補間・テンプレート曲線、CFD50・ピーク位置、
残差を重ねて確認できます。ガウシアン解析は削除しました。既存の結果ファイルは変更しません。

作成条件はpol2振幅200–350 mV、S/N 20以上、非飽和などの品質条件です。
各位置・chから一様reservoir samplingで最大100本、位置間のラウンドロビンで最大2000本/chを選び、
初期作成および各更新・最終検査で最低100本必要です。足りなければ条件を自動緩和せず停止します。
`--template-amplitude-mv`、`--template-min-snr`、`--template-min-count`、`--template-max-count`、
`--template-per-position`、`--template-seed` で設定できます。詳細は `--help` を参照してください。

評価の品質条件は振幅8 mV・S/N 5以上、フィットchi2/ndf ≤3です。作成用の200 mV下限は評価には使いません。
線形CFD50とフィットCFD50は両方成功した同じショットで比較し、成功率・失敗理由も保存します。
到達時刻分布とピーク時刻分布は共通の原点だけを引いて表示し、手法間の平均差を保持します。
CSVの時刻は元のトリガー基準のns、時刻差は列名どおりpsまたはsamplesです。
分布の標準偏差は同期系の揺らぎを含み、センサー単体の時間分解能ではありません。

保存したテンプレートを同じランの別位置へ再利用する例：

```bash
./bin/peak_timing data/scan_20261021_164141.root \
    --events 60 --segments 0:1999 --channels 2,3,4 --plot-segments 1 \
    --template-input results/scan_20261021_164141/peak_timing/template_comparison/templates \
    --tag template_reuse
```

`--template-input` とテンプレート作成条件の同時指定は拒否します。ROOT内容・ch・時間窓が一致する必要があります。
初期条件で全chが作成成功するとは限りません。今回のランのC1は残差選別後の候補不足で停止することを確認しています。
関数の使い方は `docs/reuse.md`、数値定義は `docs/numerical_definitions.md` を参照してください。

時間窓の既定値はベースライン `48:51 ns`、信号 `51.5:56.5 ns` です。上の時間窓の引数は省略できます。別の値を指定すると、その窓だけを変更できます。

TRC変換の書式は `./bin/trc_to_root INPUT_BASE [OUTPUT.root]` です。`INPUT_BASE` にはCSVの `.csv` とトレースディレクトリの `_1dX` / `_2dYX` などを除いた共通名を指定します。既存のROOTファイルは上書きしません。既存3Runの変換済みROOTは `data/` から参照できます。

詳細な引数は各実行ファイルの `--help` で表示できます。従来の `--events`、`--segments`、`--segment-start`、`--nsegments`、`--alignment`、`--output-dir` も利用できます。`--segments` と個数による指定は併用できません。

## 結果

既定の保存先は `results/<Run名>/<解析名>/<条件ラベル>/` です。`--tag` が条件ラベルで、既定値は `standard` です。

```text
results/scan_20261021_154341/
    average_waveforms/all_shots_trigger/
        conditions.json
        summary.csv
        averages.root
        waveforms/event_00084.png
        waveforms/event_00084.csv
        COMPLETE
    scan_xy_maps/standard/
        conditions.json
        channel_values.csv
        position_values.csv
        maps.root
        amplitude_maps_xy.png
        ...
```

- 平均波形は測定点別のPNG・CSVと、全選択点をまとめた `averages.root` を出力します。
- 1D波高は `amplitudes.csv` と位置に対する4chのグラフを出力します。
- XY解析はチャンネル別・位置別CSV、11個のTH2Dを含むROOT、7枚のPNGを出力します。
- 時刻解析は `peak_timing.csv`、`summary.csv`、`failure_counts.csv` を出力します。
- `templates/` にテンプレートと格子CSV、`template_selection.csv` にラン全体の候補選別件数、`template_uses.csv` に選択候補の全更新履歴を保存します。
- `template_validation.csv` は評価点×ch×振幅帯の有効件数・残差・chi2です。`distributions/event_*/` はCFD50・ピークの分布比較と振幅／位相診断、`waveforms/event_*/` は元波形・フィットと残差図です。

`conditions.json` に入力、実際に使用した測定点・ショット番号、時間窓、平均方式、実行引数、ソースのSHA-256を保存します。正常終了時に結果ディレクトリを確定し、`COMPLETE` を置きます。同じ保存先が存在するときは停止するため、比較する条件には別のラベルを付けます。`--output-dir` では保存先全体を指定できます。

入力パスは実行位置を基準に解釈し、既定の結果保存先はビルドしたプロジェクトの `results/` を基準にします。

## コードを読む・再利用する

最初に [app/average_waveforms.cpp](app/average_waveforms.cpp)、[include/waveform_data.h](include/waveform_data.h)、[include/rootfile_reader.h](include/rootfile_reader.h)、[include/waveform_analyzer.h](include/waveform_analyzer.h) を読むと、読み込みから計算までを追えます。

| 場所 | 内容 |
|---|---|
| `include/` | データ構造、条件、公開する関数・クラス |
| `src/` | 読み込み・数値計算・描画・保存の実装 |
| `app/` | 各解析の条件指定と処理の組み合わせ |
| `bin/`, `build/` | 実行ファイル、ビルド・検証生成物 |
| `tools/` | バイナリ表示、ゲイン確認などの診断 |
| `tests/` | 既知波形、ROOT入出力、既存結果との比較 |

関数の使用例は [docs/reuse.md](docs/reuse.md)、計算の定義は [docs/numerical_definitions.md](docs/numerical_definitions.md) にまとめています。

今回の確認範囲と一致結果は [docs/validation.md](docs/validation.md) に記録しています。

診断ツールは入力ファイルを指定して使用できます。

```bash
python3 tools/vertical_gain_check.py data/scan_20261021_154341_2dYX/C1--scan--00000.trc
```

## 初回テンプレートフィットのチェック

条件の検討用に `checks/template_fit/` を用意しています。
`make check-template-fit` 後、`./build/check_initial_fit INPUT.root` で、更新前の初期テンプレートに対する
全選択候補のχ²分布と、χ²区間ごとの代表波形・フィット・残差を保存します。
χ²>3の波形も表示し、通常解析の閾値は変更しません。詳しくは `checks/template_fit/README.md` を参照してください。
