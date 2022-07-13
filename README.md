# Diffuse-field-reciprocity-for-hydrodynamics
An efficient reciprocal approach for hydrodynamic response estimation in random spreading seas

## Abstract 

Estimating the hydrodynamic response of an offshore structure in a random spreading sea can lead to large computational costs. In this paper the actual spreading sea is replaced by an idealized diffuse wave field and the diffuse field reciprocity (DFR) relationship is derived analytically and verified against diffraction analysis for offshore application. The DFR approach is very efficient because only the added damping matrix is required to compute the cross spectrum of the wave loading and if normalised to the peak amplitude of a spreading sea, an upper bound response can be obtained using the reciprocal approach. And this is demonstrated using a spar type floating wind turbine. Given that the hydrodynamic coefficients are routine outputs for offshore structural design, engineers would obtain the upper bound response without additional computational cost using this new approach. 

## Method - Diffuse field reciprocity (DFR)

The “diffuse field reciprocity (DFR)” technique was first derived from the vibro-acoustics field, and it can be used to express the cross spectrum of the wave loading in terms of the resistive impedance to radiation in a reciprocal manner:

$\huge{\mathop{\mathbb{E}}[\mathbf{f}_b\mathbf{f}_b^\text{T}]=\frac{4 E(\omega)d\omega}{\pi\omega n(\omega)}\text{Im}(\mathbf{D}\_{\text{dir}}(\omega))}$

In the offshore application, the DFR relationship can be found as:

$\huge{\mathop{\mathbb{E}}[\mathbf{f}_b\mathbf{f}_b^\text{T}] \propto S\_{\eta\eta}\mathbf{C}\_p}$ 

where $S\_{\eta\eta}$ is the wave power spectrum and $\mathbf{C}\_p$ is the potential damping matrix. The potential damping is the result of waves travelling away due to radiation in a potential flow.

## Application to a spar type floating wind turbine 

<img src="/files/FroceComp.png" height="250" width="600">

<img src="/files/ResComp.png" height="250" width="600">



