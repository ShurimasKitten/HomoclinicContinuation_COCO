function f = AutoFuncODE(x,p,t)
   %%% State values
    y_N = x(1);
    y_E = x(2);
    x_B = x(3);
    x_N = x(4);

    m = p(1);
    eta = p(2);
    k2 = p(3);
    ep = p(4);
 
     %%% Pars
    SV=10.^6;              
    sig=2.1.*10.^4.*SV; 
    alphaT=1.7  .*10.^(-4);                                                                    
    alphaS=0.8  .*10.^(-3);                                                                 
    TaN=7;                                                                                 
    TaL=25              ;                                                                
    T0=2.65                ;                                                             
    S0 =35;                                                                                   
    SB =34.538;                                                                             
    y_B = alphaS ./alphaT   .* (SB - S0) ./(TaN - T0);   
    VN=7.2106  .*10.^15     ;                                                                              
    VL=6.3515   .* 10.^16       ;                                                                         
    VB=18.*VN     ;
    g=3.17  .*10^-8;      
    phi=alphaT  .*(TaN-T0);                                                                   
    x_E=(TaL-T0) ./(TaN-T0);    
    % x_N = (TaN - T0) ./(TaN - T0);
    TildeVE=VL ./VN;                                                                         
    TildeVB=VB ./VN;                                                                         
    W = 5.456.*SV ./sig;    
    delta_N = g.*VN./sig; 

    %%% Advective functions
    ps=phi.*((y_N-x_N)-(y_E-x_E));   
    Psi = ps.*tanh(ps./ep);

    %%% Convective functions
    k2n = k2;
    k2l = k2./2.0;
    k1  = 10.^-2.*SV./sig;
    kappaN = k1+0.5.*(k2n-k1)   .*(1+tanh(((y_N-y_B) - (x_N-x_B) - eta) ./ep));                                        
    kappaE = k1+0.5.*(k2l-k1)   .*(1+tanh(((y_E-y_B) - (x_E-x_B) - eta) ./ep));       

    %%% System equations
    % y_N
    ynp = m + ps./2.*(y_E-y_B) + Psi./2.*(y_E+y_B-2.*y_N) + W.*(y_E-y_N) - kappaN.*(y_N-y_B);

    % Y_E
    ylp = 1./TildeVE.*(-m + ps./2.*(y_B-y_N) + Psi./2.*(y_B+y_N-2.*y_E) - W.*(y_E-y_N) - kappaE.*(y_E-y_B));
    
    % x_B
    xbp = 1./TildeVB.*(ps./2.*(x_N-x_E) + Psi./2.*(x_N+x_E-2.*x_B) + kappaN.*(x_N-x_B) + kappaE.*(x_E-x_B));
  
    % x_N
    xnp = -delta_N.*(x_N - 1) + ps./2.*(x_E-x_B) + Psi./2.*(x_E + x_B-2.*x_N) + W.*(x_E-x_N) - kappaN.*(x_N-x_B);


    f = [ynp;ylp;xbp; xnp];
end



