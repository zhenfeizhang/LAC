# LAC
Lattice based pke and ke scheme

## LAC-v2
This directory contains implementations using OpenSSL.

## LAC-x86
This directory contains implementations without third-party libraries.

## LAC-v3a and LAC-v3b
This director contains implementations with the new version of LAC with improved parameters and efficiency. Main modifications:
1:The decryption error rate is decreased below 2^{-l}, where l=128,192,256 is the security level. 
2:Light version, LAC-light, without BCH code for embedded processors such as Cortex M4.
3:Extended version LAC-v3b with modulus q=256, which is more efficient than LAC-v3a with q=251
4:Constant-time BCH decode, secret/error sample, accumulation based polynomial multiplicaton.
