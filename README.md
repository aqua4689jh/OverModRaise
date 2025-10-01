# OverModRaise: Reducing Modulus Consumption of CKKS Bootstrapping

## Build and Run Instructions

```bash
# Compile
g++ -g proposal_main.cpp -o filename -I.

# Launch under gdb
gdb ./filename
run

# If you encounter a stack overflow or segmentation fault due to limited stack size
ulimit -s unlimited
```

## Reproducing the throughput results in Table 4 of the paper

```bash
# execute the following eight functions in proposal_main.cpp one at a time
# logN = 16
evalroundplus_test_S2C_first_adaptive_S2C<16, 37, 42, 4, 5, 58>();
erpluspar_12_S2C_first_adaptive_S2C<16, 38, 42, 4, 5, 16, 26, 30, 58>();
erpluspar_23_S2C_first_adaptive_S2C<16, 38, 42, 4, 5, 16, 26, 30, 58>();
erpluspar_all_together_bootstrap_test_S2C_first_adaptive_S2C<16, 38, 42, 4, 5, 16, 26, 30, 58>();

# logN = 15
evalroundplus_test_S2C_first_adaptive_S2C<15, 27, 30, 7, 7, 49>();
erpluspar_12_S2C_first_adaptive_S2C<15, 28, 30, 7, 7, 10, 20, 23, 49>();
erpluspar_23_S2C_first_adaptive_S2C<15, 28, 30, 7, 7, 10, 20, 23, 49>();
erpluspar_all_together_bootstrap_test_S2C_first_adaptive_S2C<15, 28, 30, 7, 7, 10, 20, 23, 49>();
```

## Reproducing the throughput results in Table 7 of the paper

```bash
# execute the following eight functions in proposal_main.cpp one at a time
# logN = 16
plain_test_S2C_first_adaptive_S2C<16, 58, 42, 4, 5, 58>();
plain_12_S2C_first_adaptive_S2C<16, 59, 42, 4, 5, 16, 26, 8, 58>();
plain_23_S2C_first_adaptive_S2C<16, 59, 42, 4, 5, 16, 26, 8, 58>();
plain_all_together_bootstrap_test_S2C_first_adaptive_S2C<16, 59, 42, 4, 5, 16, 26, 8, 58>();

# logN = 15
plain_test_S2C_first_adaptive_S2C<15, 49, 30, 7, 7, 49>();
plain_12_S2C_first_adaptive_S2C<15, 50, 30, 7, 7, 10, 20, 8, 49>();
plain_23_S2C_first_adaptive_S2C<15, 50, 30, 7, 7, 10, 20, 8, 49>();	
plain_all_together_bootstrap_test_S2C_first_adaptive_S2C<15, 50, 30, 7, 7, 10, 20, 8, 49>();
```
