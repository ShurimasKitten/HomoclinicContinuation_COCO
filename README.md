# HomoclinicContinuation_COCO
This repository presents a continuation scheme for the continuation and bifurcation analysis of homoclinic connections in ordinary differential equations, utilizing the COCO Continuation Toolbox in MATLAB and largely following (Kuznetsov, Champneys, 1994).

Homoclinic connections are trajectories in the phase space of a dynamical system that leave a saddle equilibrium point and then return to the same equilibrium as time tends to positive and negative infinity. Such connections play a critical role in understanding a system's dynamics, including chaotic behavior, since they often organize nearby bifurcations. Understanding the structure of homoclinic connections in a parameter space is a starting point for unraveling the bifurcation diagram.

We consider equations of the form

$$\frac{dx}{dt} = f(x,\mu),$$

with state variables $x \in \mathbb{R}^n$ and parameters $\mu \in \mathbb{R}^2$. 

## Boundary value problem

The first continuation scheme is the most stable and uses projection boundary conditions on the equilibrium eigenspaces and an integral phase condition to close the problem. The full boundary value problem becomes:

<div style="background-color: white; padding: 15px; border: 1px solid #ddd; border-radius: 5px;">
  
$$f(\mathbf u_0,\mu)=0,$$  

$$\int^1_0 \frac{d\tilde{\mathbf u}(t)}{dt} \mathbf u(t)\,dt = 0,$$  

$$J(\mathbf w_u,\lambda)\mathbf v_{u,i} = \lambda_{u,i}\mathbf v_{u,i}, \quad i = 1,\dots,n_u,$$  

$$J(\mathbf w_s, \lambda)\mathbf v_{s,i} = \lambda_{s,i}\mathbf v_{s,i}, \quad i = 1, \dots, n_s,$$  

$$\mathbf v^\dagger_{u,i}\mathbf v_{u,i} = 1, \quad i = 1, \dots, n_u,$$  

$$\mathbf v^\dagger_{s,i}\mathbf v_{s,i} = 1, \quad i = 1, \dots, n_s,$$  

$$L_s(\mathbf w_s - \mathbf u_0) = 0,$$  

$$L_u(\mathbf w_u - \mathbf u_0) = 0.$$

</div>

The projection operators $L_s$ and $L_u$ are reconstructed at each continuation step and ensured to vary continiously with the parameters $\mu\in\mathbb R^2,$ following the approach of (Kuznetsov, Champneys, 1994). 


## A four-dimensional climate model
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


## MATLAB useadge
The boundary value problem (BVP) is demonstrated using a four-dimensional climate model as a representative example. We begin by performing one-parameter continuation and branching a periodic solution from a Hopf bifurcation point  
```markdown
# Load settings
[probSettings, thmEq, thmPO, thmHB, thmSN, thmHom, thmSNPst, thmSNPun, thmPDst] = loadDefaultSettings();

# Update settings and run one-parameter continuation
probSettings.contSettings.h0 = 1e-2;
probSettings.contSettings.PtMX = [1000 1000];
probSettings.contSettings.h_max = 2e-2;
run1Dcont(temp1D, 'EQ_run1', [-2.39e-3, -2.6, 0.015], 'mu', [-6e-3 0.0]);

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
prob = proj_isol2hom(fnPOi, 90, homSet);
coco(prob, 'Hom_run1', [], 1, {'mu', 'eta', 'RES', 'isSF'})





