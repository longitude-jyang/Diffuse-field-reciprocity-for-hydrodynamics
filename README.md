# Diffuse-field-reciprocity-for-hydrodynamics
An efficient reciprocal approach for hydrodynamic response estimation in random spreading seas

## Abstract 

Estimating the hydrodynamic response of an offshore structure in a random spreading sea can lead to large computational costs. In this paper the actual spreading sea is replaced by an idealized diffuse wave field and the diffuse field reciprocity (DFR) relationship is derived analytically and verified against diffraction analysis for offshore application. The DFR approach is very efficient because only the added damping matrix is required to compute the cross spectrum of the wave loading and if normalised to the peak amplitude of a spreading sea, an upper bound response can be obtained using the reciprocal approach. And this is demonstrated using a spar type floating wind turbine. Given that the hydrodynamic coefficients are routine outputs for offshore structural design, engineers would obtain the upper bound response without additional computational cost using this new approach. 

## Motivation - An efficient method for response estimation in a spreading sea 

Realistic sea states contain waves that are dependent on frequency as well as direction. The response of an offshore structure can be sensitive to both the frequency and the heading of an incident ocean wave, and the most realistic way of ensuring that this effect is captured in a numerical simulation is to model the wave environment as a random spreading sea. For a linear system subject to a random wave excitation, the covariance matrix of the structure response can be found as:

$$\huge{\mathbf{\sigma}^2 =\int\mathbf{H}\[\int\mathbf{f}\_b\mathbf{f}\_b^\text{T} D(\theta)d\theta\]\mathbf{H}^{*T}(\omega)S\_{\eta\eta}(\omega)d\omega}$$

The key assumption here is that the spreading function D(θ) can be a viewed as a measure for the fraction of sea states, out of a large number of ensemble, that are coming from direction θ. With this view, the expression inside the bracket is the expected cross spectrum of the blocked force by ensemble average and the covariance matrix of the hydrodynamic response can be more compactly expressed as: 

$$\huge{\mathbf{\sigma}^2 =\int\mathbf{H} \mathop{\mathbb{E}}[\mathbf{f}\_b\mathbf{f}\_b^\text{T}]\mathbf{H}^{*T}(\omega)S\_{\eta\eta}(\omega)d\omega}$$

At first sight, it might appear that we haven’t gained any computational saving by converting the expression. However, if the spreading sea is replaced by an idealised diffuse wave field, the DFR technique can be used to turn the expectation into a simple algebra expression. 

## Method - Diffuse field reciprocity (DFR)

The “diffuse field reciprocity (DFR)” technique was first derived from the vibro-acoustics field, and it can be used to express the cross spectrum of the wave loading in terms of the resistive impedance to radiation in a reciprocal manner:

$$\huge{\mathop{\mathbb{E}}[\mathbf{f}_b\mathbf{f}_b^\text{T}]=\frac{4 E(\omega)d\omega}{\pi\omega n(\omega)}\text{Im}(\mathbf{D}\_{\text{dir}}(\omega))}$$

In the offshore application, the DFR relationship can be found as:

$$\huge{\mathop{\mathbb{E}}[\mathbf{f}_b\mathbf{f}_b^\text{T}] \propto S\_{\eta\eta}\mathbf{C}\_p}$$

where $S\_{\eta\eta}$ is the wave power spectrum and $\mathbf{C}\_p$ is the potential damping matrix. The potential damping is the result of waves travelling away due to radiation in a potential flow.

## Application to a spar type floating wind turbine 

<img src="/files/FroceComp.png" height="250" width="600">

<img src="/files/ResComp.png" height="250" width="600">



