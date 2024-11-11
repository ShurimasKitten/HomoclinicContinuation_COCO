# HomoclinicContinuation_COCO
Homoclinic connections are trajectories $\mathbf u(t)$ in the phase space of a dynamical system that leave a saddle equilibrium $\mathbf u_0$ and eventually return to the same equilibrium as time tends to positive and negative infinity. Such connections are important for understanding the dynamics of a system, as they often organize nearby bifurcations, making their computation a starting point in unraveling the bifurcation diagram. 

<div style="border: 1px solid black; padding: 10px; display: inline-block;">
  <img src="Images/NormalHomo.png" alt="Description of the image">
</div>

This repository presents a continuation scheme for the continuation and bifurcation analysis of homoclinic connections in ordinary differential equations, utilizing the COCO Continuation Toolbox in MATLAB, largely following (Kuznetsov, Champneys, 1994). We do not aim to provide a detailed description of the resulting dynamics but instead refer the reader to "Elements of Applied Bifurcation Theory" by Yuri Kuznetsov as an entry point to the literature. 

We begin by formulating the boundary value problem (BVP) for path-following a homoclinic trajectory with respect to the parameters $\mu$. We then describe the detection of codimension-two bifurcations, followed by working examples of the continuation of homoclinic connections in a four dimension climate model. In these examples, we demonstrate the detection of a codimension-two Belyakov point and a codimension-two resonance point with real eignvalues. In the appendices, we give the description of this climate model, along with the technical details regarding the detection of codimension-two bifurcation

## The boundary value problem
We consider differential equations equation of the form

$$\frac{d\mathbf u}{dt} = f(\mathbf u(t),\mu),$$

with state variables $\mathbf u \in \mathbb{R}^n$ and bifurcation parameters $\mu \in \mathbb{R}^2$. A homoclinic orbit with $\mathbf u(t)$ for $t \xrightarrow{} \pm \infty$ may be continued with respect to the parameters $\mu$, and will generically trace out a codimension-one curve in the parameter space. To this end, we implement a boundary value problem (BVP), which employs projection boundary conditions and an integral phase condition, while simultaneously tracking the equilibrium point $\mathbf u_0$. The BVP is given by 

<div style="display: flex; justify-content: center; align-items: center; background-color: white; padding: 15px; border: 1px solid #ddd; border-radius: 5px;">
  
$$f(\mathbf u_0,\mu)=0$$  

$$J(\mathbf u_0, \mu) \mathbf v_{u,j} = \lambda_{u,j}\mathbf v_{u,j}, \quad i = 1,\dots,n_u$$  

$$J(\mathbf u_0, \mu) \mathbf v_{s,i} = \lambda_{s,i}\mathbf v_{s,i}, \quad i = 1, \dots, n_s$$  

$$\mathbf v^\dagger_{u,i}\mathbf v_{u,j} = 1, \quad j = 1, \dots, n_u$$  

$$\mathbf v^\dagger_{s,i}\mathbf v_{s,i} = 1, \quad i = 1, \dots, n_s$$  

$$\int^1_0 \frac{d\tilde{\mathbf u}(t)}{dt} \mathbf u(t) dt = 0$$  

$$L_s(\mu) \cdot (\mathbf w_s - \mathbf u_0) = 0$$  

$$L_u(\mu) \cdot (\mathbf w_u - \mathbf u_0) = 0$$

</div>

where $J$ is the jacobian of $f$, $\mathbf u(t)$ represents the homoclinic solution at the current continuation step, and $\tilde{\mathbf u}(t)$ represents it at the previous step. The endpoints of the homoclinic connection are given by $\mathbf w_{s}$ and $\mathbf w_{u}$, which lie in the stable and unstable linear eignspaces of $\mathbf u_0$, respectively. These eigenspaces are spanned by the eignvectors $\mathbf v_{s,i}$ nad $\mathbf v_{u,j}$. 

The projection operators, $L_s$ and $L_u$, are reconstructed at each continuation step to ensure that they vary continuously with the parameters $\mu\in\mathbb R^2,$ following the approach of (Kuznetsov, Champneys, 1994). More precisely, we solve the linear system

$$U_s(\mu)\left(V(\mu)^T V(\tilde{\mu})\right) = V(\tilde{\mu})V(\tilde{\mu})^T,$$

for $U_s$, where $V(\mu)$ is a matrix whose columns span the orthogonal complement of the unstable (linear) eignspace $E_u^\perp(\mu)$ at the current parameter value $\mu$, and similarly columns of $V(\tilde\mu)$ span the eignspace $E_u^\perp(\tilde\mu)$ at the previous parameter value $\tilde\mu$. We also solve the linear system

$$U_u(\mu)\left(W(\mu)^T W(\tilde{\mu})\right) = W(\tilde{\mu})W(\tilde{\mu})^T,$$

where the columns of $W(\mu)$ and $W(\tilde\mu)$ belong to the span of the orthogonal complement of the stable eigenspaces $E_s^\perp(\mu)$ and $E_s^\perp(\tilde\mu)$, respectively. The desired projection operators are then given by 

$$L_s(\mu) = U_s(\mu)V(\alpha)^T,$$

$$L_u(\mu) = U_u(\mu)W(\alpha)^T,$$

which both vary smoothly with the parameters $\mu$.

## Codimension-two homoclinic bifurcations
A homoclinic connection can become degenerate at isolated codimension-two points along a homoclinic bifurcation curve. These degeneracies act as organisation centers for the surrounding parameter space. We monitor the type of equilibrium involved in the homoclinic connection and detect several codimension-two homoclinic bifurcations. The classification of equilibria is given in the following table.

<div align="center">
  
| EQtype ($\pm{}$) | Name                                      |
|:-------:|:-----------------------------------------:|
| 1     | Saddle with real leading eigenvalues                            |
| 2     | Saddle-focus with 2D stable manifold                      |
| 3     | Saddle-focus with 2D unstable manifold    |
| 4     | Bi-focus |
  
</div>

Similarly, the detected codimension-two homoclinic bifurcations are classified as follows.

<div align="center">
  
| Symbol  | Name                                      |
|:-------:|:-----------------------------------------:|
| NSS     | Neutral saddle                            |
| NSF     | Neutral saddle-focus                      |
| DRS     | Double real leading stable eigenvalues    |
| DRU     | Double real leading unstable eigenvalues  |
| TLS     | Three leading stable eigenvalues          |
| TLR     | Three leading unstable eigenvalues        |
| NDS     | Neutrally Divergent Saddle-Focus (Stable) |
| NDU     | Neutrally Divergent Saddle-Focus (Unstable)|
| OFS     | Orbit flip w.r.t stable manifold          |
| OFU     | Orbit flip w.r.t unstable manifold        |
| H       | Shilnikov-Hopf bifurcation                |
| S       | Non-central homoclinic saddle-node bifurcation |
| RES     | Zero of the saddle value                  |
| EqType     | Change in the type of equilibrium                  |

</div>

In addition, when an orbit flip (OFS or OFU) is detected, its type - A, B, or C - is printed for further classification. 

## Working example
The boundary value problem (BVP) is demonstrated using a four-dimensional climate model as a representative example. We begin by performing one-parameter continuation and branching a periodic solution from a Hopf bifurcation point  
```markdown
# Load settings
# IMPORTANT: Inside 'loadDefaultSettings()', you must identify yout COCO compatible function handle and a secondary ODE handle. 
[probSettings, thmEq, thmPO, thmHB, thmSN, thmHom, thmSNPst, thmSNPun, thmPDst] = loadDefaultSettings();

# Update settings and run one-parameter continuation
probSettings.contSettings.h0 = 1e-2;
probSettings.contSettings.PtMX = [1000 1000];
probSettings.contSettings.h_max = 2e-2;
run1Dcont(probSettings, 'EQ_run1', [-2.39e-3, -2.6, 0.015], 'mu', [-6e-3 0.0]);

# Collect the Hopf points
HB_labs = coco_bd_labs('EQ_run1', 'HB');

# Update settings and branch off the second Hopf point. We turn bifurcation detection 'off'.
probSettings.corrSettings.TOL = 1e-4;
probSettings.collSettings.NTST = 150;
probSettings.contSettings.PtMX = [0 1000];
probSettings.contSettings.h0 = 1e-2;
probSettings.contSettings.h_max = 2e2;
PO_hb2po(probSettings, 'EQ_run1', HB_labs(2), 'PO_run1', 'off') 

# Plot time-series near homoclinic at LAB=80
figure(1)
hold on
[s, d] = po_read_solution('PO_run1', 80);
plot(s.tbp(:,1), s.xbp(:,1))

# Initilise homoclinic continuation. NOTE: We need to set mesh adaoption off (NAdapt = 0) due to a bug. 
probSettings.corrSettings.TOL = 1e-4;
probSettings.collSettings.NTST = 150;
probSettings.contSettings.PtMX = [1000 1000];
probSettings.contSettings.h0 = 1e-2;
probSettings.contSettings.h_max = 2e-2;
probSettings.collSettings.NAdapt = 0;
prob = proj_isol2hom(fnPOi, 90, homSet);
coco(prob, 'Hom_run1', [], 1, {'mu', 'eta', 'RES', 'isSF'})
```


# Known bugs and TODO

- Known bugs that I pinky promise to fix:
  - COCOs mesh-adaption has to be turned off (i.e. NAdapt = 0). There seems to be a conflict with updating the phase condition before each continuation step. 

- TODO:
  - Computation of the adjoint problem.
  - Following from the adjoint problem, the test functions for inclination flips.
  - Orientability index.
  - More efficient continuation of eigenspaces, such as employing the Riccati equation approach as implemented in MATCONT.

# Appendix

## Appendix A: Detection of codimension-two homoclinic bifurcations 
```
      ▄▄▄▄▄    ▄▄▄▄▄
    ▄█▀▀▀▀▀█▄ ▄█▀▀▀▀█▄
    █░░░░░░░█ █░░░░░░█
    █░░░░░░░█ █░░░░░░█
    █░░░░░░░█ █░░░░░░█
    ▀█▄▄▄▄▄█▀ ▀█▄▄▄▄█▀
  ▄▄█▄▄▄▄▄█▄▄▄█▄▄▄▄▄█▄▄▄▄
▄█░░█░░░░█░░░█░░░░█░░░█░█▄
█░░█░░░░░█░░█░░░░░█░░░█░░█
█░█░░░░░░░█░█░░░░░█░░░█░░█
▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀
```

## Appendix B: A four-dimensional climate model 
Examples of homoclinic continuation is given using a four-dimensional differential equation motivated by climate science. The model takes the form

$$\frac{dx_N}{dt} = -\delta_N(x_N-1) + \frac{\Psi}{2}(x_E-x_B) + \frac{|\Psi|}{2}(x_E+x_B-2x_N) + W(x_E-x_N) - \mathbf{K}_N(x_N-x_B),$$
 
$$\frac{dy_N}{dt} = \mu_N + \frac{\Psi}{2}(y_E-y_B) + \frac{|\Psi|}{2}(y_E+y_B-2y_N) + W(y_E-y_N) - \mathbf{K}_N(y_N-y_B),$$

$$\frac{dy_E}{dt} = \frac{1}{\widetilde V_E}\left(\mu_E + \frac{\Psi}{2}(y_B-y_N) + \frac{|\Psi|}{2}(y_B+y_N-2y_E) - W(y_E-y_N) - \mathbf{K}_E(y_E-y_B)\right),$$

$$\frac{dx_B}{dt} = \frac{1}{\widetilde V_B}\left(\frac{\Psi}{2}(x_N-x_E) + \frac{|\Psi|}{2}(x_N+x_E-2x_B) + \mathbf{K}_N(x_N-x_B) + \mathbf{K}_E(x_E-x_B)\right),$$

where the so-called convective exchange functions are given by

$$\mathbf{K}_N = \kappa_d + \frac{1}{2}(\kappa_c^N - \kappa_d)\left(1 + \tanh\left(\frac{(y_N - y_B) - (x_N - x_B) - \eta}{\varepsilon}\right)\right), $$

$$\mathbf{K}_E = \kappa_d +  \frac{1}{2}(\kappa_c^E - \kappa_d)\left(1 + \tanh\left(\frac{(y_E - y_B) - (x_E - x_B) - \eta}{\varepsilon}\right)\right),$$

and the advective strength by

$$\Psi = \alpha_T\left(T^a_N-T_0)(y_N - x_N - (y_E-x_E)\right).$$

The remainder of the variables appear in the MATLAB code and a comprehensive description of its dynamics will appear in future publication. 




