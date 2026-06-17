# a

在N维空间中 有正方形网格，

每个网格包含大量粒子M{质量mass,速度\vec {vel}}，服从各向同性的正态分布，粒子在每个网格中均匀方便。

已知大量粒子的统计学量:

* 总质量 $Mass = \sum _m^M m.mass$

* 平均速度 $\vec {Vel} = \frac{1}{Mass} \sum _m^M m.mass*m.\vec {vel}$

* 平均速度方差 $VelVarSq = \frac{1}{Mass} \sum _m^M m.mass*(m.\vec {vel} - \vec {Vel})^2$

* 总动量 $Momentum = \sum _m^M m.mass*m.\vec {vel} = Mass * \vec {Vel}$

* 总能量 $Energy = 0.5 \sum _m^M m.mass*m.\vec {vel}^2$

网格大小为 $Volume$，

对于一个边界$edge_dir_vec$，其方向表示边界的方向，长度表示边界的大小（面积）

经过$dt$时间，一些粒子穿过了边界，对于每个粒子m，其穿过边界的概率为$max(vel \cdot edge_dir_vec,0) \cdot dt \cdot edge_area / Volume$

求所有穿过的粒子的 总质量，总动量，总能量