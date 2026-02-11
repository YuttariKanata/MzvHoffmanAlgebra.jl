# MultipleZetaValueHoffmanAlgebra.jl 開発者・高度専門家向けガイド

本ドキュメントは、**MultipleZetaValueHoffmanAlgebra.jl** のコードベースの拡張、最適化、あるいは理論的検証を行う開発者および数理専門家を対象としています。

---

## 1. 型システムと設計思想

本パッケージの核心は、多重ゼータ値の代数構造 $\mathfrak{H} = \mathbb{Q}\langle x, y \rangle$ および $\mathfrak{H}^1 = \mathbb{Q}\langle z_1, z_2, \dots \rangle$ を Julia の型システムへ写像することにあります。

### 階層構造

* **AbstractWord**: 単項式（語）の抽象型。内部表現はハッシュ効率とメモリ効率を考慮し `Tuple{Vararg{Int}}` を採用しています。
  * `HoffmanWord`: $x=0, y=1$ のビット列。
  * `IndexWord`: 正整数の列 $(k_1, \dots, k_r)$。
* **ZetaExpr**: 線形結合の抽象型。
  * `Hoffman`, `Index`: `Dict{Word, Rational{BigInt}}` をラップし、有理数係数による有限線形結合を表現します。
* **Poly{A}**: 多項式環 $R[T]$。係数 $A$ に `Hoffman` や `Index` を取ることで、正規化パラメーター $T$ を含む代数計算を可能にします。

### 不変性と破壊的操作

* **Immutability**: すべての `Word` 型は不変 (Immutable) です。
* **_add! ヘルパー**: 線形結合の構築時には、内部関数 `_add!(dict, word, coeff)` を使用してください。これは係数が 0 になった項を即座に削除し、辞書のスパース性を維持するための最適化です。

---

## 2. 低レイヤー演算（monomials.jl）

高レベルな積演算（`shuffle_product`, `stuffle_product`）の実装は、パフォーマンス上の理由から `monomials.jl` に分離された単項式レベルの演算に依存しています。

### シャッフル積の最適化

* **動的計画法 (DP)**: インデックス表記におけるシャッフル積 `monomial_sh_i` は、Hoffman 表記への変換を経由せず、直接ブロックベースの DP アルゴリズムで計算されます。これにより、重さや深さが大きい語の計算において、中間オブジェクトの生成を劇的に抑制しています。
* **Suffix-recursive**: 低レイヤーの実装では、Vector の `push!` や `append!` と相性の良い後方再帰的 (Suffix-recursive) な定義を採用しています。

---

## 3. 向き (Orientation) と正規化の厳密性

MZV の代数的定式化において、生成元 $z_k$ の定義（向き）の選択は計算結果に直結します。

### Orientation の統一

パッケージ全体で `:left` ($z_k = yx^{k-1}$) をデフォルトとしていますが、変換・正規化・導分の各関数は必ず `orientation` キーワード引数を受け取る設計になっています。

* **`:left`**: 発散成分は語の末尾（インデックスの最後）に現れます。
* **`:right`**: 発散成分は語の先頭（インデックスの最初）に現れます。

### 正規化カーネル (`_regularize`)

`regularization.jl` の `_regularize` 関数は、代数的な正規化（$y$ または $1$ を不定元 $T$ と見なす操作）を抽象化しています。

1. 語を「収束部分」と「発散端の $n$ 個の連続成分」に分解。
2. 代数的な積演算（shuffle または stuffle）の逆を適用し、再帰的に有限部分を抽出。

---

## 4. 演算子システム (AbstractOp)

インデックス代数 $\mathfrak{H}^1$ への作用素は、`AbstractOp` を継承する構造体群として実装されています。

* **Operator**: `Vector{AbstractOp}` のラッパーであり、積 `*` による合成はベクトルの連結に対応します。
* **clean 関数**: 演算子の合成時、隣接する同種の作用素の合算や、双対作用素 $\tau$ の対合性 ($\tau^2 = id$) による簡約を自動的に行います。
* **作用の優先順位**: `op * idx` の形式で作用させる際、演算子は **右から左** （数学的な写像の合成順序）へ順次適用されます。

---

## 5. 線形写像の実装パターン

新しい代数関係式や準同型写像を実装する場合、以下のテンプレートが推奨されます。

```julia
function my_linear_map(h::Hoffman)
    new_terms = Dict{HoffmanWord, Rational{BigInt}}()
    for (w, c) in h.terms
        # 単項式に対する作用を計算 (結果は Hoffman 型)
        image = apply_to_word(w) 
        # _add! を用いて結果をマージ
        for (iw, ic) in image.terms
            _add!(new_terms, iw, ic * c)
        end
    end
    return Hoffman(new_terms)
end
```

---

## 6. テストと品質管理

数学的ライブラリである性質上、微小なバグが巨大な計算結果の誤りを生むため、多層的なテストを要求します。

* **Property-based Testing**: `test_all.jl` に含まれる `rand_hoffman` 等を用いたランダムテストにより、加法の交換則や積の結合則、分配則を検証してください。
* **既知データ検証**: 歴史的に知られている $\rho$ 写像の値や導分関係式の結果を用いた回帰テスト（`test_datas.jl`）をパスする必要があります。
* **BigInt の強制**: 係数計算には必ず `Rational{BigInt}` を使用してください。MZV の関係式では、重さ 20 程度でも分子・分母が通常の 64bit 整数を容易に超過します。

---
**参照元:**

* DEVELOPER_GUIDE.md: 型定義
* monomials.jl: 単項式アルゴリズム
* regularization.jl: $\rho$ 写像の実装
* types.jl: 作用素階層構造
