# HomoclinicContinuation_COCO
Three seperate continuation schemes for the continuation and bifurcation analysis of homoclinic orbits in ordinary differential equations built on the COCO Continuation Toolbox in MATLAB. More precisely, we consider homoclinic connections of equations in the form

$$\frac{dx}{dt} = f(x,\mu),$$

with state variables $x \in \mathbb{R}^n$ and parameters $\mu \in \mathbb{R}^2$. We largly follow the approx of (Kuznetsov, Champneys, 1994), and consider a boudnary value problem (BVP) in the following form

<div style="background-color: white; padding: 15px; border: 1px solid #ddd; border-radius: 5px;">
  ![Equation](https://latex.codecogs.com/svg.image?\bg{white}\begin{aligned}&f(\mathbf{u}_0,\mu)=0,\\&\frac{d\mathbf{u}(t)}{dt}-T&space;f(\mathbf{u}(t),\lambda)=0,\\&\int^1_0\frac{d\tilde{\mathbf{u}}(t)}{dt}\mathbf{u}(t)\,dt=0,\\&J(\mathbf{w}_u,\lambda)\mathbf{v}_{u,i}=\lambda_{u,i}\mathbf{v}_{u,i},\quad&space;i=1,\dots,n_s,\\&J(\mathbf{w}_s,\lambda)\mathbf{v}_{s,i}=\lambda_{s,i}\mathbf{v}_{s,i},\quad&space;i=1,\dots,n_u,\\&\mathbf{v}^*_{ui}\mathbf{v}_{ui}=1,\quad&space;i=1,\dots,n_u,\\&\mathbf{v}^*_{si}\mathbf{v}_{si}=1,\quad&space;i=1,\dots,n_s,\\&L_s(\mathbf{w}_s-\mathbf{u}_0)=0,\\&L_u(\mathbf{w}_u-\mathbf{u}_0)=0.\end{aligned})
</div>

