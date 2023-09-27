# PhotomaskGen_HKUST
Software developed in 2021 for non-graphical generation of greyscale photomasks in .gdsii format.  Developed by Dr. Christopher Salmean at the Hong Kong University of Science and Techology.

This code was used in the generation of photomasks for my [published article](https://www.sciencedirect.com/science/article/pii/S0017931023006129). A more detailed explanation of the theoretical basis can be found in the linked article.  

To summarise, gradient microstructures may be generated using a single photolithographic step, by locally altering the penetration depth of the incident light. This is done by breaking the shape down into constituent pixels, each smaller than the maximum resolution of the mask aligner, and controlling the light throughput of each pixel by filling the pixel to a specified extent. This necessitates the individual placement of around 250 million elements, so a Python code was written for the procedural, hierarchical generation of the desired photmasks.  

![1-s2 0-S0017931023006129-gr2_lrg](https://github.com/csalmean/PhotomaskGen_HKUST/assets/133036780/ceea24b4-c833-4d8c-a672-fea34f6b87a3)

In this example, substrates with ramp-, thorn-, NACA0050-, dimple-, and torus-shaped gradient microstructures were fabricated by this means:  

![1-s2 0-S0017931023006129-gr5](https://github.com/csalmean/PhotomaskGen_HKUST/assets/133036780/11d8d17c-0eaf-4d0c-bd55-0dd113fa029f)
