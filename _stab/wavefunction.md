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

input[type="button"] {
  display: inline-block;
  background-color: #ffffff;
  border: 2px solid #aaaaaa;
  font-size: 16px;
  height: 32px;
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

div[id=divWave] {
  margin-bottom: 32px;
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

# シュレディンガー方程式
$$
i \hbar \frac{\partial}{\partial t} |\Psi(\mathbf{r}, t)⟩
= \hat{H} |\Psi(\mathbf{r}, t)⟩
$$

$i$ は虚数単位。

$\hbar$ は [reduced Pranck constant](https://en.wikipedia.org/wiki/Planck_constant#Value) 。

$\Psi$ は波動関数。

$\mathbf{r}$ は位置ベクトル。

$t$ は時間。

$\hat{H}$ は Hamiltonian operator 。

## 波動関数

- [Wave Functions](http://quantum.phys.cmu.edu/CQT/chaps/cqt02.pdf)
- [quantum mechanics - What is a wave function in simple language? - Physics Stack Exchange](https://physics.stackexchange.com/questions/249239/what-is-a-wave-function-in-simple-language)

## bra-ket記法
braは行ベクトルを表す。

$$
⟨\phi| = [\phi_0, \phi_1, \phi_2, \dots]
$$

ketは列ベクトルを表す。

$$
|\psi⟩ = \begin{bmatrix}
\psi_0\\
\psi_1\\
\psi_2\\
\vdots
\end{bmatrix}
$$

ベクトル $\phi$ とベクトル $\psi$ の内積 $\langle\phi, \psi⟩$ はbra-ket記法で次のように書ける。

$$
⟨\phi, \psi⟩
= ⟨\phi|\psi⟩
$$

左右を入れ替えてket、braの順にすると直積になる。

$$
\phi\otimes\psi = |\psi⟩ ⟨\phi|
$$

- [Quantum Physics II, Lecture Notes 4 - MIT8_05F13_Chap_04.pdf](https://ocw.mit.edu/courses/physics/8-05-quantum-physics-ii-fall-2013/lecture-notes/MIT8_05F13_Chap_04.pdf)

## 複素共役と複素共役転置
複素共役 (complex conjugate) は複素数の虚部の正負を反転させた数。

$$
\begin{aligned}
z = a + ib\\
\overline{z} = a - ib
\end{aligned}
$$

複素共役転置 (conjugate transpose あるいは Hermitian transpose) 。

$$
(\mathbf{A}^{*})_{ij} = \overline{\mathbf{A}}_{ji}
$$
