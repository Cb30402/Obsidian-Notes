---
tags:
  - 量子力学
---
WKB 近似考虑的问题情景是: 能量为 $E$ 的粒子穿过势能 $V(x)$ 的区域, 这里的 $V(x)$ 是随 $x$ 变化而"缓慢"变化的.
在定态薛定谔方程的章节, 我们已经有了一些结论. 当 $E>V$ 时, 波函数形式为 $$\psi(x)=Ae^{\pm ikx},\quad k=\sqrt{2m(E-V)}/\hbar$$ 而当 $E<V$ 时, 波函数形式为 $$\psi(x)=Ae^{\pm \kappa x},\quad \kappa=\sqrt{2m(V-E)}/\hbar$$ 但在经典拐点 $(E=V)$ 附近, 问题变得棘手.
我们先讨论拐点以外的"经典区域".

#### **9.1 "经典"区域**
薛定谔方程为 $$- \frac{\hbar^{2}}{2m} \frac{\mathrm{d}^{2}\psi}{\mathrm{d}x^{2}}+V(x)\psi=E \psi$$ 可以改写为 $$\frac{\mathrm{d}^{2}\psi}{\mathrm{d}x^{2}}=- \frac{p^{2}}{\hbar^{2}} \psi$$ 其中 $$p(x)=\sqrt{2m[E-V(x)]}$$ 可以说是经典意义上的动量函数.
不妨考虑 $E>V(x)$ 的情况, 此时 $p(x)$ 为实数, 我们考虑对波函数进行拆分: 拆分为振幅 $A(x)$ 和相位 $\phi(x)$ 两部分, 即 $$\psi(x)=A(x)e^{i \phi(x)}$$ 这一形式的波函数代入改写后的薛定谔方程中得到 $$A''+2iA'\phi'+iA \phi''-A(\phi')^{2}=- \frac{p^{2}}{\hbar^{2}} A$$ 这个方程可以根据实部虚部分离为两个方程, 即为 $$A''-A(\phi')^{2}=- \frac{p^{2}}{\hbar^{2}}A\quad 与\quad  2A'\phi'+A \phi''=0$$ 后一个式子是虚部对应的, 我们迅速注意到它等同于 $$(A^{2}\phi')'=0$$ 这说明 $A^{2}\phi'$ 某种意义上是一个守恒量, 我们可以得到一个关系为 $$A=\frac{C}{\sqrt{|\phi'|}}$$ 将这一关系代回实部对应的方程之前, 我们做**近似**(这里是 WKB 的近似) 为: $A''/A\to0$, 也就是振幅 $A$ 的变化非常缓慢, 这样我们代入整理得到实部的方程变为 $$(\phi')^{2}=\frac{p^{2}}{\hbar^{2}}$$ 这可以解得 $$\phi(x)=\pm \frac{1}{\hbar} \int p(x)\mathrm{d}x$$ 我们处理一下系数, 将 $\hbar$ 一并纳入 $C$ 之中, 得到最终波函数形式为 $$\psi(x)\approx 
\frac{C}{\sqrt{p(x)}} e^{\pm \frac{i}{\hbar}\int p(x)\mathrm{d}x}$$ 如果 $E<V$, 只是 $p(x)$ 变为虚数, 得到的结果为 $$\psi(x)\approx \frac{C}{\sqrt{ |p(x)|}} e^{\pm \frac{1}{\hbar}\int |p(x)|\mathrm{d}x}$$
#### **9.2 隧道效应**

考虑粒子穿过顶部不平的势垒 $V(x)$, 势垒左端波函数形式为 $$\psi(x)=Ae^{ikx}+Be^{-ikx}$$
右端波函数形式为 $$\psi(x)=Fe^{ikx}$$
根据散射理论, 透射率为 $$T=\frac{|F|^{2}}{|A|^2}$$
而在势垒区域, 由 WKB 近似给出波函数形式为 $$\psi(x)\approx \frac{C}{\sqrt{|p(x)| }} e^{\frac{1}{\hbar} \int_{0}^{x}|p(x')|\mathrm{d}x'}+\frac{D}{\sqrt{|p(x)| }} e^{-\frac{1}{\hbar} \int_{0}^{x}|p(x')|\mathrm{d}x'}$$ 势垒又高又宽时, 指数增长项的系数 $C$ 应当很小, 这时入射透射的相对振幅取决于指数衰减项, 也就是 $$\frac{|F|}{|A|}\sim e^{-\frac{1}{\hbar} \int_{0}^{x}|p(x')|\mathrm{d}x'}$$ 那么透射率自然为 $$T\approx e^{-2\gamma},\quad \gamma=-\frac{1}{\hbar}\int_{0}^{a}|p(x)|\mathrm{d}x$$
#### **9.3 连接公式**

这一节我们讨论拐点处的 WKB 近似的形式. 拐点处两侧波函数的形式在 WKB 近似下为 $$\psi(x)\approx\left\{\begin{matrix}
 \frac{1}{\sqrt{p(x)}}  \Big[ Be^{\frac{i}{\hbar} \int_{x}^{0}}p(x')\mathrm{d}x'+Ce^{-\frac{i}{\hbar} \int_{x}^{0}}p(x')\mathrm{d}x' \Big],\quad x< 0\\
\frac{1}{\sqrt{|p(x)|}} De^{\frac{1}{\hbar} \int_{0}^{x}|p(x')|\mathrm{d}x'},\quad x> 0
\end{matrix}\right.$$ 我们的任务是找到系数 $B,C,D$ 之间的关系, 从而将两端的公式连接, 以实现经典拐点处的"平滑过渡".
在拐点(我们事先将拐点移动到了 $x=0$)处, 考虑对势能做**近似**为 $$V(x)\approx E+V'(0)x$$ 将这个线性化的势能代入到薛定谔方程中, 化简得到 $$\frac{\mathrm{d}^{2}\psi_{p}}{\mathrm{d}x^{2}}=\alpha^{3}x \psi_{p},\quad \alpha=\Big[ \frac{2m}{\hbar^{2}}V'(0) \Big]^{1/3}$$
做一步变量代换 $z=\alpha x$, 于是得到很简洁的形式为 $$\frac{\mathrm{d}^{2}\psi_{p}}{\mathrm{d}z^{2}}=z \psi_{p}$$ 这个微分方程是 Airy 方程, 它的解为 Airy 函数, 通解形式为 $$\psi_{p}(x)=a\mathrm{Ai}(\alpha x)+b\mathrm{Bi}(\alpha x)$$
在前面变量代换的基础上, 现在我们的动量函数为 $$p(x)\approx \hbar \alpha^{3/2}\sqrt{-x}$$
在右端区域 $z>0$, 积分项整理出 $$\int_{0}^{x}|p(x')|\mathrm{d}x'= \frac{2}{3} \hbar (\alpha x)^{3/2}$$ 从而 WKB 近似波函数写为 $$\psi(x)\approx \frac{D}{\sqrt{\hbar}\alpha^{3/4}x^{1/4}}e^{- \frac{2}{3}(\alpha x)^{3/2}}$$ 但同时我们可以写出 $z$ 很大条件下的 Airy 函数的近似为 $$\psi_{p}(x)\approx \frac{a}{2\sqrt{\pi}(\alpha x)^{1/4}}e^{- \frac{2}{3}(\alpha x)^{3/2}}+\frac{b}{\sqrt{\pi}(\alpha x)^{1/4}}e^{\frac{2}{3}(\alpha x)^{3/2}}$$ 比较 $\psi$ 和 $\psi_{p}$, 我们迅速得到 $$a=\sqrt{\frac{4\pi}{\alpha\hbar}}D,\quad b=0$$
类似地, 我们在左端区域能整理出积分项为 $$\int_{x}^{0}p(x')\mathrm{d}x'\approx \frac{2}{3}\hbar(-\alpha x)^{3/2}$$ WKB 近似给出 $$\psi(x)\approx \frac{1}{\sqrt{\hbar}\alpha^{3/4}(-x)^{1/4}}\Big[ Be^{i \frac{2}{3}(-\alpha)^{3/2}}+Ce^{-i \frac{2}{3}(-\alpha x)^{3/2}} \Big]$$ 同时我们写出 $-z$ 很大时的 Airy 函数的近似为 $$\begin{aligned}
\psi_{p}(x) & \approx \frac{a}{\sqrt{\pi}(-\alpha x)^{1 / 4}} \sin \left[\frac{2}{3}(-\alpha x)^{3 / 2}+\frac{\pi}{4}\right] \\
& =\frac{a}{\sqrt{\pi}(-\alpha x)^{1 / 4}} \frac{1}{2 i}\left[e^{i \pi / 4} e^{i \frac{2}{3}(-\alpha x)^{3 / 2}}-e^{-i \pi / 4} e^{-i \frac{2}{3}(-\alpha x)^{3 / 2}}\right]
\end{aligned}$$ 依然比较 $\psi$ 与 $\psi_{p}$, 系数上的关系为 $$\frac{a}{2i \sqrt{\pi }}e^{i \pi/4}=\frac{B}{\sqrt{\hbar \alpha}}\quad,\frac{-a}{2i \sqrt{\pi }}e^{-i \pi/4}=\frac{C}{\sqrt{\hbar \alpha}}$$
与前面的系数关系联系一下, 我们可以得到 $$B=-ie^{i \pi}D,\quad C=ie^{-i\pi/4
}D$$