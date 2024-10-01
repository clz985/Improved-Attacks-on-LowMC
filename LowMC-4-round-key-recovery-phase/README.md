The code is used to verify the correctness of attack on LowMC-4-round-key-recovery-phase

We generate valid affine layers and save them to file Affine_layers.txt, so you can load it directly from this file.

You need to update the path of Affine_layers.txt in your code so you can run  code LowMC-4round.sage correctly.

You can also generate affine layers by yourself with function create_affine_layers in LowMC-4round.sage and you need to check if affine layers are valid.

If not, you need to call this function until it produces valid affine layers.

In LowMC-4round.sage, we set secret key as a vector with only first 3 bits are 0 and the rest are 1.

So you can check whether variables x0 to x128 are the same as secret key.