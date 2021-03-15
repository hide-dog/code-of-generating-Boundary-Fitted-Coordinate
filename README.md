# code-of-generating-Boundary-Fitted-Coordinate  
agg.jl  
格子数や読み込むデータセットを変更し、実行することで、境界適合格子を生成できます  
格子数、出力ファイル名、読み込みファイル名、O型格子の半径、O型格子の中心x座標、データタイプを選択し実行します。  
読み込みに対応しているデータセットは2つ用意しており、NACA0012タイプとNACA64A010タイプがあります。前者は上面と下面が別々に記載されたデータセット。後者は後端から下面を通り、一周するように記載されたデータセットです。  
詳しくは中身を見てみてください。

agg_v2.jl  
物体近傍の格子を細かくするため、格子を中央に寄せることができます。
壁面の第一格子幅を決めることで、拡大率を計算し生成します。  
  
pre.jl  
ホームページで紹介していたコードを使用する場合は，こちらのpre.jlを使用ください．

plot.py  
格子確認用  

reference  
データの出典元をこちらに記載しています  

This software is released under the MIT License, see LICENSE.
