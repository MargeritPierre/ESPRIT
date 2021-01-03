# ESPRIT: a generalized algorithm for MATLAB

*A function implementing a generalized signal pole estimation algorithm called ESPRIT.* 

The orginal ESPRIT method can be found in `ROY, Richard et KAILATH, Thomas. ESPRIT-estimation of signal parameters via rotational invariance techniques. IEEE Transactions on acoustics, speech, and signal processing, 1989, vol. 37, no 7, p. 984-995`. 

The theory of the present implementation can be found (in French) in the chapter IV of `MARGERIT, Pierre. Caractérisation large bande du comportement dynamique linéaire des structures hétérogènes viscoélastiques anisotropes: application à la table d'harmonie du piano. 2018. Thèse de doctorat. Paris Est.`.

### This ESPRIT algorithm implements:
- Multidimensionnal signals (RD-ESPRIT)
- Standard `exp` and Paired `cos` components
- Least-Squares or Total Least-Squares regression (LS/TLS-ESPRIT)
- Multiple Invariance or Multiple Resolution (MI-ESPRIT)
- Signal order criterion: MDL, ESTER, SAMOS
- Stabilization Diagrams
- Signal subspace tracking
- Signal decimation
