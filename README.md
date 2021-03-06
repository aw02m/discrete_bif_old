# descrete_bif
離散力学系(差分方程式)の分岐解析ツールです．  
本プログラムは3部に分かれています:

1. pp (Phase portrait; 相平面描画ツール)
2. fix (Fixed point; 固定点計算及び特性定数の表示)
3. bif (Bifurcation; 分岐集合の計算)

## pp概要
ppは相平面をリアルタイムに描画します．(Original:Tetsushi Ueta)
固定点計算のための近似値取得に使用してください．  

### 動作環境
* numpy
* scipy
* matplotlib

### pp入力ファイル要素
* "xrange", "yrange" : x,y描画区間．
* "alpha" : 打点の透明度(alpha値)．0~1の範囲．
* "break" : 一回の描画に計算するステップ数．2000以上は重いかも．
* "x0" : プログラム実行時の初期値．
* "params" : 力学系のパラメタ．
* "dparams" : プログラム内でパラメタを変更する際の刻み幅．
* "explode" : 座標が発散したと判定するまでのノルム値．
* "func" : ここに差分方程式を書き込みます．pythonの文法に則って記述してください．ただしパラメタは`data.dict['params'][i]`で表記します．内部的にはこのフォイールドの文字列がpython内で`eval()`されます．

### 使用法
`$ python pp.py [json input file]`で実行．プログラム実行中は次のキーを押下することで各種パラメタ変更，ファイル出力等が行えます．
* `w` : `__ppout__.json`ファイルに現在の情報を保存します．ファイルがすでに存在している場合は上書きされます．
* `p` : 変化させるパラメタを指定します．現在の指定インデックスはターミナルに表示されます．
* `[矢印キー]` : 指定しているパラメタを"dparams"に従って変化させます．
* `[space]` : 描画をリセットします．
* その他の機能についてはProf. Uetaのオリジナルリポジトリを参照ください．

## fix概要
fixはppで取得した固定点情報をもとに，パラメタを変化させながらNewton法にて精度の良い固定点計算を行います．

### 動作環境
* Eigen3 : 線形代数ライブラリです．比較的新しい関数を使用しているため，gitリポジトリの最新バージョンを利用してください．Arch Linuxなひとは`# pamac build eigen-git`でインストールされます．
* nlohmann : jsonライブラリ．
* cmake : Makefileの自動生成に用います．面倒な人はMakefileを自力で書いてください．

### fix入力ファイル概要
* "fixed" : ppにてjsonファイルを出力した際に追加されます．この値が固定点計算の初期値として使用されます．
* "period" : 固定点の周期を指定します．ppにて周期を確認して正確な整数値を与えてください．
* "inc_param" : 変化させるパラメタを指定します．"params"のインデックスで指定してください．
* "delta_inc" : パラメタ変分量です．経験的に0.1~0.001の範囲で設定すると良いかと思います．
* "inc_iter" : 指定した回数固定点を計算します．(指定した計算回数に到達もしくは発散しないとプログラムは固定点を計算し続けます．)
* "max_iter" : Newton法の最大ステップを指定します．Newton法は通常数回の繰り返しで収束するため，10~32の整数値を指定します．
* "eps" : Newton法の収束判定を与えます．誤差が指定数値以下になった場合に計算が終了します．

### 使用法
`./main [input json file]`にて実行．コンパイルは`fix`ディレクトリに新たに`build`ディレクトリを作成し，その中で`cmake`するとMakefileが自動で作成されます．  
力学系の写像及びその微分は`ds_func.cpp`に記述してください．
計算に成功すると固定点座標，パラメタ値，特性定数，特性定数のノルム・偏角が出力されます．特性定数のノルムが1に近いもの(分岐点)をピックアップしてbifプログラムに渡してください．

## bif概要
bifはfixで求めた分岐点の近似値をもとに分岐集合を計算します．NS,G,PDにそれぞれプログラムが別れていますが，中身はほとんど一緒です．なお，NSに関してはアルゴリズムの特性上すべての局所分岐を計算できます．複素数を扱うため遅くなりますが，Eigen側でGPU演算を有効にしていると(標準で有効)気になりません．

### bif入力ファイル概要
* "var_param" : 分岐計算にて変数として扱うパラメタを指定します．
* "dif_strip" : 数値微分を用いる際の微小量を指定します．本バージョンでは2階微分に一部数値微分が使用されています．

### 使用法
fixと同じです．新たにbuildディレクトリを作成し，その中でcmakeしてMakefileを生成してください．
力学系の写像及びその微分は`ds_func.cpp`に記述してください．
指定回数の計算が終了すると`out`ファイルに計算結果が保存されます．gnuplot(`plot.gp`として内包)やpythonのmatplotlibを使うなどしてグラフにしてください．

## TO DO
* bifプログラムの統合
* 2階微分
* オブジェクト指向化（C++は関数よりメソッドのほうが早いらしい）
