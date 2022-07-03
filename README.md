# Tropical cryptography

This Git repository contains our implementation of some protocols of public key cryptography as part of our master's degree report *State of the art of tropical cryptography, 2021*, `State of the art of tropical cryptography.pdf`.

It contains our implementation of the following protocols of tropical cryptography:
- *Tropical cryptography I, 2014* by D. Grigoriev and V. Shpilrain: both protocols.
  - `tropical_stickel.sage` for the first protocol.
  - `tropical_public_key_encryption.sage` for the second protocol.
- *Tropical cryptography II, 2018* by D. Grigoriev and V. Shpilrain: first protocol.
  - `tropical_semidirect_product.sage` for the first protocol, since the second protocol is not considered in our report.
- *Modifying the tropical version of Stickel's key exchange protocol, 2019* by A. Muanalifah and S. Sergeev: both protocols.
  - `tropical_linde.sage` for the protocol based on Linde-De La Puente matrices.
  - `tropical_jones.sage` for the protocol based on Jones matrices.
  - `model.lp` and `tropical_jones.py` for the ASP implementation used in `tropical_jones.sage`. The environment Clingo (https://github.com/potassco/clingo) is required to run these files.

It contains the implementation of two standard protocols, on which tropical schemes are based:
- *A new method for exchanging secret keys, 2005* by E. Stickel:
  - the general theoretical protocol in `stickel_theoretical_scheme.sage`.
  - the protocol based on a semiring of matrices in `stickel_matrix_scheme.sage`.

Finally, the implementation of attacks on the above tropical protocols:
- *Analysis of a key exchange protocol based on tropical matrix algebra, 2018* by M. Kotov and A. Ushakov: both attacks.
  - `tropical_stickel.sage` for the heuristic attack.
  - `tropical_simplex.sage` for the full attack based on the simplex algortihm. This implementation belongs to M. Kotov and A. Ushakov. Their original implementation in Gap is available at https://github.com/mkotov/tropical.
- *Cryptanalysis of Stickel's Key Exchange Scheme 2008* by V. Shpilrain: `stickel_matrix_scheme.sage` *(not added yet)*.
- *Modifying the tropical version of Stickel's key exchange protocol, 2019* by A. Muanalifah and S. Sergeev:
    - `tropical_linde.sage` for the heuristic attacks on the Linde-De La Puente protocol.
- *Remarks on a tropical key exchange system, 2020* by D. Rudy and C. Monico: `tropical_semidirect_product_attack.sage`.
- *A closer look at the tropical cryptography, 2020* by S. Isaac and D. Kahrobaei: their implementation is available at https://github.com/steveisaac/TropicalCryptography.

A last file `tropical_class.sage` contains the implementation of general functions required for most of our programs.


A. Herledan Le Merdy and C. Foucault
