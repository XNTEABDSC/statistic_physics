对于随机变量X~N(\mu,\sigma^2)，请给出E(x*I(x>0)), E(x^2*I(x>0)), E(x^3*I(x>0))


在N维正方形网格中，每个格子的体积为 $\text{Vol} = e^N$，边的（N-1维）面积为 $\text{Edge} = e^{N-1}$，法方向单位向量为 $\mathbf{n} = \text{DirVec}$。设粒子速度分布为各向同性正态分布，格子总质量为 $\text{Mass}$，平均速度向量为 $\mathbf{v}_{\text{mean}}$，速度方差为 $v_{\text{sq}} = \frac{1}{\text{Mass}} \sum dm \cdot (\mathbf{v} - \mathbf{v}_{\text{mean}})^2$。

定义：
$$
\mu = \mathbf{v}_{\text{mean}} \cdot \mathbf{n}, \quad \sigma = \sqrt{\frac{v_{\text{sq}}}{N}}, \quad a = \frac{\mu}{\sigma}.
$$
令 $\phi(a) = \frac{1}{\sqrt{2\pi}} e^{-a^2/2}$ 为标准正态概率密度函数，$\Phi(a) = \int_{-\infty}^a \phi(x) dx$ 为标准正态累积分布函数。

则在时间 $dt$ 内通过该边的粒子总质量、总动量和总能量分别为：

$$
\boxed{
\begin{aligned}
\Delta M &= \text{Mass} \cdot \frac{dt \cdot \text{Edge}}{\text{Vol}} \cdot \left( \mu \Phi(a) + \sigma \phi(a) \right), \\
\Delta \mathbf{P} &= \text{Mass} \cdot \frac{dt \cdot \text{Edge}}{\text{Vol}} \cdot \left[ \left( \mu \Phi(a) + \sigma \phi(a) \right) \mathbf{v}_{\text{mean}} + \sigma^2 \Phi(a) \mathbf{n} \right], \\
\Delta E &= \text{Mass} \cdot \frac{dt \cdot \text{Edge}}{\text{Vol}} \cdot \frac{1}{2} \left[ \mu \Phi(a) \left( |\mathbf{v}_{\text{mean}}|^2 + (N+2)\sigma^2 \right) + \sigma \phi(a) \left( |\mathbf{v}_{\text{mean}}|^2 + (N+1)\sigma^2 \right) \right].
\end{aligned}}
$$

其中，$\Delta M$、$\Delta \mathbf{P}$ 和 $\Delta E$ 分别表示通过边的质量、动量和能量。

考虑一个 $N$ 维空间中的正方形网格单元，边长为 $e$，体积为 $\mathrm{Vol} = e^N$。单元内粒子总质量为 $\mathrm{Mass}$，平均速度为 $\boldsymbol{v}_{\text{mean}}$，速度方差为 $v_{\text{sq}} = \frac{1}{\mathrm{Mass}} \sum \mathrm{d} mass \cdot (\boldsymbol{v} - \boldsymbol{v}_{\text{mean}})^2$。速度分布为各向同性的正态分布，即 $\boldsymbol{v} \sim \mathcal{N}(\boldsymbol{v}_{\text{mean}}, \sigma^2 \mathbf{I})$，其中 $\sigma^2 = v_{\text{sq}} / N$。

设单元的一个面为 $N-1$ 维超平面，面积为 $\mathrm{Edge} = e^{N-1}$，单位法向量为 $\boldsymbol{n} = \mathrm{DirVec}$。在时间 $\mathrm{d}t$ 内，粒子通过该面的概率为 $\max(\boldsymbol{v} \cdot \boldsymbol{n}, 0) \cdot \mathrm{d}t \cdot \mathrm{Edge} / \mathrm{Vol}$。通过积分速度分布，得到通过该面的粒子总质量、总动量和总能量如下：

定义：
- $u_n = \boldsymbol{v}_{\text{mean}} \cdot \boldsymbol{n}$（平均速度在法方向的分量），
- $\sigma = \sqrt{v_{\text{sq}} / N}$，
- $\xi = u_n / \sigma$，
- $\phi(x) = \frac{1}{\sqrt{2\pi}} e^{-x^2/2}$（标准正态密度函数），
- $\Phi(x) = \int_{-\infty}^x \phi(t) \, \mathrm{d}t$（标准正态累积分布函数），
- $G_1 = u_n \Phi(\xi) + \sigma \phi(\xi)$。

则：
- 通过的质量：
  $$
  \Delta M = \mathrm{d}t \cdot \frac{\mathrm{Edge}}{\mathrm{Vol}} \cdot \mathrm{Mass} \cdot G_1.
  $$
- 通过的动量：
  $$
  \Delta \boldsymbol{P} = \mathrm{d}t \cdot \frac{\mathrm{Edge}}{\mathrm{Vol}} \cdot \mathrm{Mass} \cdot \left( G_1 \boldsymbol{v}_{\text{mean}} + \sigma^2 \Phi(\xi) \boldsymbol{n} \right).
  $$
- 通过的能量：
  $$
  \Delta E = \frac{1}{2} \mathrm{d}t \cdot \frac{\mathrm{Edge}}{\mathrm{Vol}} \cdot \mathrm{Mass} \cdot \left( \sigma^2 \Phi(\xi) + G_1 \left( (N-1)\sigma^2 + |\boldsymbol{v}_{\text{mean}}|^2 \right) \right).
  $$

其中，$\mathrm{Edge}/\mathrm{Vol} = 1/e$。这些表达式基于速度分布为各向同性正态分布的假设，并包含了所有阶矩的贡献。