# MultipleZetaValueHoffmanAlgebra.jl 技術解説・利用ガイド

本パッケージは、**多重ゼータ値 (MZV)** の代数構造、特に **Hoffman 代数** $\mathfrak{H} = \mathbb{Q}\langle x, y \rangle$ とそのインデックス表記 $\mathfrak{H}^1$ を計算機上で厳密に扱うための Julia ライブラリです。

多重ゼータ値の専門的な背景を持つユーザーに向けて、数学的定義と実装の対応を中心に解説します。

---

## 1. 代数構造の定義と表現

本パッケージでは、代数 $\mathfrak{H}, \mathfrak{H}^1, \mathfrak{H}^0$ の元を以下の型で表現します。

### 単項式 (Word)

| 数学的な対象 | Julia 型 | 内部表現 | 備考 |
| :--- | :--- | :--- | :--- |
| **Word** ($x, y$ の語) | `HoffmanWord` | `Tuple{Vararg{Int}}` | $x=0, y=1$ のビット列 |
| **Index** ($k_1, \dots, k_r$) | `IndexWord` | `Tuple{Vararg{Int}}` | 正整数の組 |

* **不変性**: これらの Word 型は **Immutable** です。一度生成された語の内容を変更することはできません。
* **生成例**: `HoffmanWord(1, 0)` は語 $yx$ を、`IndexWord(2, 1)` はインデックス $(2, 1)$ を表します。

### 線形結合 (Expression)

| 数学的な対象 | Julia 型 | 内部表現 |
| :--- | :--- | :--- |
| **有理係数線形結合** | `Hoffman`, `Index` | `Dict{Word, Rational{BigInt}}` |

* 内部的には、語をキー、有理数係数を値とする辞書で管理されています。
* 係数が 0 になった項は自動的に削除され、常に簡約された状態で保持されます。

---

## 2. 表記の相互変換と「向き」の定義

Hoffman 代数表記とインデックス表記の変換には、生成元 $z_k$ の定義（向き）が重要です。本パッケージではキーワード引数 `orientation` で制御します。

| 向き (`orientation`) | 定義 | $z_2$ の Hoffman 表現 | 収束条件（Admissibility） |
| :--- | :--- | :--- | :--- |
| **`:left` (デフォルト)** | $z_k = y x^{k-1}$ | `10` ($yx$) | 末尾が $x$ (0) / $k_r > 1$ |
| **`:right`** | $z_k = x^{k-1} y$ | `01` ($xy$) | 先頭が $x$ (0) / $k_1 > 1$ |

### 変換コード例

```julia
# インデックス (2, 1) を Hoffman 代数の元へ変換
h = Hoffman(IndexWord(2, 1), orientation=:left) # -> yx y (10 1)

# Hoffman 代数の元をインデックス表記へ変換
idx = Index(h, orientation=:left) # -> (2, 1)
```

---

## 3. 代数演算と積の構造

### 連結積 (Concatenation Product)

通常の乗算演算子 `*` は、非可換な語の連結として定義されています。

* `x * y` $\to$ $xy$
* `x^3` $\to$ $xxx$

### 各種の代数的な積

MZV の基本的な関係式を生む各種の積が実装されています。

1. **シャッフル積 (`shuffle_product`)**:
    * 記号: `ш`, `⨝`
    * 反復積分表示に基づく積。全 $\mathfrak{H}$ 上で定義されます。
2. **調和積/スタッフル積 (`stuffle_product`)**:
    * 記号: `∗`
    * 級数表示に基づく積。$\mathfrak{H}^1$ 上で定義されます。
3. **等号付きスタッフル積 (`star_stuffle_product`)**:
    * 記号: `⋆`
    * 等号付き多重ゼータ値 (MZSV) の積規則に従います。

---

## 4. 正規化 (Regularization)

収束しない多重級数や発散する反復積分から有限部分を取り出すための正規化が、代数的に定式化されています。

* **シャッフル正規化 (`reg_sh`)**: シャッフル積の意味での有限部分を返します。
* **調和正規化 (`reg_st`)**: スタッフル積の意味での有限部分を返します。
* **正規化多項式**: `shuffle_regularization_polynomial` 等を使用すると、$T = \zeta(1)$ と置いた $T$ の多項式として結果を得られます。

### $\rho$ 写像

スタッフル正規化の変数 $T$ とシャッフル正規化の変数 $T$ を結びつける $\mathbb{R}$-線形写像 $\rho: \mathbb{R}[T] \to \mathbb{R}[T]$ が実装されています。
定義: $\rho(e^{Tu}) = A(u)e^{Tu}$ ただし $A(u) = \exp\left(\sum_{n \ge 2} \frac{(-1)^n}{n} \zeta(n) u^n\right)$。

---

## 5. 導分関係式と双対性

### 導分 (`dell`, `∂`)

多重ゼータ値の線形関係式を生成する導分 $\partial_n$ が利用可能です。

* 定義例 (`:left`): $\partial_n(x) = y(x+y)^{n-1}x, \partial_n(y) = -y(x+y)^{n-1}x$。
* `dell(h, n)` により $\partial_n(h)$ を計算します。

### 双対性 (`dual`)

* **標準的な双対 (`dual`)**: 語を反転させ $x \leftrightarrow y$ を入れ替える反自己同型。
* **Hoffman 双対 (`Hoffman_dual`)**: 収束端の $y$ を固定し、他を入れ替える自己同型。
* **Landen 双対 (`Landen_dual`)**: ランデン結合に対応する変換。

---

## 6. 演算子システム (Operator System)

インデックス $\mathbf{k} = (k_1, \dots, k_r)$ に対して作用する演算子を組み合わせて、複雑な操作を簡潔に記述できます。

* **基本演算子**:
  * `up(n)` / `⬆️`: 先頭成分を $+n$ する。
  * `left(n)` / `⬅️`: 左側に `1` を $n$ 個追加する。
  * `τ`: 双対写像 $\tau$。
  * `∂(n)`: 導分 $\partial_n$。

* **合成と作用**:
    演算子同士を `*` で合成し、インデックスに `*` で作用させます。作用は **右から左** へ順次適用されます。

    ```julia
    op = left(1) * up(1)
    result = op * Index(2, 1)  # -> Index(1, 3, 1)
    ```

---

## 7. 便利な表示機能

* **`sortedprint(obj)`**: 線形結合の項を重さ（weight）や長さ、語の辞書順にソートして、整形表示します。
* **`upper_represent(true)`**: 演算子の表示において、累乗を上付き文字（例: $\partial_1^2$）で表示するように設定します。
* **`T`**: 多項式環の変数 $T$ は定数としてエクスポートされており、`3*T^2 + 2*T + 1` のように直接記述できます。

---
**参考文献**:

* 荒川恒男・金子昌信『多重ゼータ値入門』
* M. Hoffman, *The algebra of multiple harmonic series*, J. of Algebra
* K. Ihara, M. Kaneko, and D. Zagier, *Derivation and double shuffle relations for multiple zeta values*, Compos. Math.
