# Tropical cryptography

This Git repository contains our implementation of some protocols of public key cryptography as part of our master's degree dissertation *State of the art of tropical cryptography, 2021*.

This Git repository contains our implementation of the following protocols of tropical cryptography:
- *Tropical cryptography I, 2014* by D. Grigoriev and V. Shpilrain: both protocols.
  - `tropical_stickel.sage` for the first protocol.
  - `tropical_public_key_encryption` for the second protocol *(not added yet)*.
- *Tropical cryptography II, 2018* by D. Grigoriev and V. Shpilrain: both protocols.
  - `tropical_semidirect_product.sage` for both *(not added yet)*.
- *Modifying the tropical version of Stickel's key exchange protocol, 2019* by A. Muanalifah and S. Sergeev: both protocols.
  - `tropical_linde_delapuente.sage` for the protocol based on Linde-De La Puente matrices *(not added yet)*.
  - `tropical_jones.sage` for the protocol based on Jones matrices *(not added yet)*.
  - `model.lp` and `tropical_jones.py` for the ASP implementation used in `tropical_jones.sage`.

It contains the implementation of two standard protocols, on which tropical schemes are based:
- The Diffie-Hellman public key exchange protocol over <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{Z}/n\mathbb{Z}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{Z}/n\mathbb{Z}" title="\mathbb{Z}/n\mathbb{Z}" /></a> in `diffie_hellman.sage` *(not added yet)*.
- *A new method for exchanging secret keys, 2005* by E. Stickel: general protocol and protocol based on a semiring of matrices as suggested by the author in `stickel.sage` *(not added yet)*

Finally, the implementation of attacks on the above tropical protocols:
- *Analysis of a key exchange protocol based on tropical matrix algebra, 2018* by M. Kotov and A. Ushakov: both attacks.
  - `tropical_stickel.sage` for the heuristic attack *(not added yet)*.
  - `tropical_simplex.sage` for the full attack based on the simplex algortihm. This implementation belongs to M. Kotov and A. Ushakov. Their original implementation in Gap is available at https://github.com/mkotov/tropical *(not added yet)*.
- *Remarks on a tropical key exchange system, 2020* by D. Rudy and C. Monico: `tropical_semidirect_product_attack.sage` *(not added yet)*.

A last file `tropical_algebra.sage` contains the implementation of general functions essential for most of our implementations.


A. Herledan Le Merdy and C. Foucault
