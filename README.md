# Lattice DFT
An python script implementation of the Density Functional Theory (cDFT) for Lattice fluids in 3D geometries.

The cDFT is the extension of the equation of state to treat inhomogeneous fluids. For a fluid with temperature T, total volume V, and chemical potential $\mu$ specified, the grand potential, $\Omega$, is written as

$$\Omega[\rho(\boldsymbol{r})] = F[\rho (\boldsymbol{r})] +  \int_{V} [ V^{(\text{ext})}(\boldsymbol{r}) - \mu ]\rho(\boldsymbol{r}) d\boldsymbol{r}$$

[^1] 

# Examples


# Pre-requisites

- Numpy, Scipy, Matplotlib, 
- SciencePlot
```bash
pip install SciencePlots
```

----
# References

[^1]: [Monson, P. A. 2012. “Understanding Adsorption/Desorption Hysteresis for Fluids in Mesoporous Materials Using Simple Molecular Models and Classical Density Functional Theory.” Microporous and Mesoporous Materials 160: 47–66.](https://doi.org/10.1016/j.micromeso.2012.04.043.)
