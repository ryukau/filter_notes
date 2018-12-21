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

# ばね-ダンパ波動方程式のパラメータ設定
ばねの項 $k u$ とダンパの項 $a \dot{u}$ を加えて質量の係数 $m$ をつけた波動方程式です。

$$
m \ddot{u} + a \dot{u} + k u = c^2 \nabla^2 u
$$

実際の物理量を使ってパラメータ $m$ 、 $a$ 、 $k$ 、 $c$ を決めます。

## 単位
この文章では $\rm{[m/s^2]}$ のように単位を角カッコで囲みます。

||英語|単位|数式中|備考|
|-|-|-|-|-|
|長さ（メートル）|length (metre)|$\rm{m}$|$L$||
|時間（秒）|time (second)|$\rm{s}$|$t$||
|質量（キログラム）|mass (kilogram)|$\rm{kg}$|$m$||
|力（ニュートン）|force (Newton)|$\rm{N}$|$F$|$\rm{N=kg \cdot m/s^2}$|
|圧力（パスカル）|pressure (Pascal)|$\rm{Pa}$|$p$|$\rm{Pa=N/m^2}$|
|密度|density|$\rm{kg/m^3}$|$\rho$||
|粘度|viscosity|$\rm{Pa \cdot s}$|$\mu$|$\rm{Pa \cdot s=kg \cdot m^{-1} \cdot s^{-1}}$|
|体積弾性率|bulk modulus|$\rm{Pa}$|$B$||
|ヤング率|Young's modulus|$\rm{Pa}$|$E$||

## パラメータ設定
ばね-ダンパ波動方程式の左右の項は単位が力 $\rm{[N]}$ です。

3次元のとき u にかかる力の図。

前回までのシミュレーションでは空間を格子状に分割していました。質量 $m\,\rm{[kg]}$ は材料の密度 $\rho\,\rm{[kg/m^3]}$ と格子の各辺の長さ $\Delta_x\,\rm{[m]}$ から計算できます。

$$
m=\rho\Delta_x^3
$$

ダンパの項の係数 $a\,\rm{[kg/s]}$ は材料の粘度 $\mu\,\rm{[Pa \cdot s]}$ 、格子の表面積 $6 \Delta_x^2\,\rm{[m]}$ 、格子の辺の長さ $\Delta_x\,\rm{[m]}$ から計算できます。

$$
a = \frac{6 \Delta_x^2}{\Delta_x} \mu = 6 \Delta_x \mu
$$

ばね定数 $k\,\rm{[N/m]}$ はヤング率 $E\,\rm{[Pa]}$ 、格子の表面積 $6 \Delta_x^2\,\rm{[m]}$ 、格子の辺の長さ $\Delta_x\,\rm{[m]}$ から計算できます。

$$
k = \frac{6 \Delta_x^2}{\Delta_x} E = 6 \Delta_x E
$$

波の伝達速度 $c\,\rm{[m/s]}$ は[体積弾性率](http://hyperphysics.phy-astr.gsu.edu/hbase/permot3.html#c1) $B\,\rm{[Pa]}$ と密度 $\rho\,\rm{[kg/m^3]}$ から計算できます。

$$
c = \sqrt{\frac{B}{\rho}}
$$

Eの値については[Young's Modulus](http://homepages.uc.edu/~bortnelj/labs/Physics%201%20experiments/Youngs%20Modulus/Young's%20Modulus%20htm.htm)を参照。真鍮が90、銅が130、鉄が200程度。

- 弦の太さとシンバルの厚さ
- 長さ $L$ と断面積 $A$ のイメージ図

## 材料とパラメータ
$\Delta_x = 0.4 / 32$ 、 $A = \Delta_x^2$ 、 $L = 0.2 \cdot 10^{-3}$


|材料|$\rho$|$E$|$B$||$m$|$a$|$k$|$c$|
|-|-|-|-|-|-|-|-|-|
|鉄|7.8 g/cm^3|200 GPa||||||4500|
|真鍮||90|||||||
|水|||||||||
|空気|||||||||

## 実験
鉄の塊。

## 参考サイト
### 理論
- [Friction - Wikipedia](https://en.wikipedia.org/wiki/Friction)
- [Deformation (mechanics) - Wikipedia](https://en.wikipedia.org/wiki/Deformation_(mechanics))
- [Viscosity - Wikipedia](https://en.wikipedia.org/wiki/Viscosity)
- [Stiffness - Wikipedia](https://en.wikipedia.org/wiki/Stiffness)
- [Attenuation - Wikipedia](https://en.wikipedia.org/wiki/Attenuation)

### 値
- [Speed of Sound](http://hyperphysics.phy-astr.gsu.edu/hbase/Sound/souspe2.html)
- [Young's Modulus](http://homepages.uc.edu/~bortnelj/labs/Physics%201%20experiments/Youngs%20Modulus/Young's%20Modulus%20htm.htm)
- [Anatomy of a Cymbal | SABIAN Cymbals](https://www.sabian.com/en/pages/anatomy-of-a-cymbal)
