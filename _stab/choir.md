<style>
body {
  max-width: 704px;
  margin: auto;
  padding: 32px 8px;
}

canvas {
  /* image-rendering: pixelated; */
  display: inline-block;
  border-style: solid;
  border-width: 0.2px;
  border-color: #000000;
}

.numberInputLabel {
  display: inline-block;
  text-align: right;
  width: 160px;
  margin: 0 8px 4px 0;
  padding-right: 8px;
  border-right: solid 3px #606060;
}

.numberInputRange {
  display: inline-block;
  max-width: 220px;
  width: 50%;
  margin: 0 8px 0 0;
}

.numberInputNumber {
  display: inline-block;
  max-width: 100px;
  width: 20%;
  margin: 0;
}

code {
  overflow-x: scroll;
  overflow-y: hidden;
  white-space: pre;
}

.katex {
  font-size: 1.3em !important;
}
</style>

# コーラスの合成
## 倍音
|倍音|音量|
|-:|-:|
|1|0.5|
|2|1|
|3|0.5|
|4|1|
|8-9|0.1|
|20-24|0.05|

## PADsynth
1. 適当なウェーブテーブルを用意。
2. ウェーブテーブルの周波数を8~10倍上げる。
3. ウェーブテーブルの周波数成分を7~8左シフト。
4. BandWidthを40~200でレンダリング。
5. 3kHzでローパス。

- 基本周波数 220Hz
