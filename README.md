# Post-quantum Oblivious Transfers - Kyber Branch
This is the emp-ot implementation of [Masny-Rindal](https://eprint.iacr.org/2019/706) 1-out-of-2 OT protocol based on CRYSTALs-Kyber. We adapt this implementation from the former [libOTe implementation](https://github.com/osu-crypto/libOTe) of the same protocol.

The other Saber-based implementations can be found in the `master` branch.

## Dependency
Our implementation uses [emp-tool](https://github.com/emp-toolkit/emp-tool). Here is how to install it, more details can be found in their [README](https://github.com/emp-toolkit/emp-tool/blob/master/README.md).

1. `wget https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/scripts/install.py`
2. `python install.py --deps --tool `
    1. You can use `--ot=[release]` to install a particular branch or release
    2. By default it will build for Release. `-DCMAKE_BUILD_TYPE=[Release|Debug]` option is also available.
    3. No sudo? Change [`CMAKE_INSTALL_PREFIX`](https://cmake.org/cmake/help/v2.8.8/cmake.html#variable%3aCMAKE_INSTALL_PREFIX).

## Build
After installing the dependencies, you can build the project with the following commands:

```shell
cmake .
make -j
./run ./bin/test_mr_kyber
```

where `test_saber` is our test program for the SABER-based OT, running our Saber-based OTs (including NPOTs, Simplest OTs, and MR OTs) and some CDH-based classical OTs. The output should be like:

```
connected
connected
128 COOTs (average over 1000 runs):	12.921 ms
128 COOTs (average over 1000 runs):	12.928 ms
128 NPOTs (average over 1000 runs):	14.866 ms
128 NPOTs (average over 1000 runs):	14.898 ms
128 MR Kyber OTs (average over 1000 runs):	17.684 ms
128 MR Kyber OTs (average over 1000 runs):	17.681 ms
```

Credited to [emp-ot](https://github.com/emp-toolkit/emp-ot) for the OT implementation,
 one can test these OTs on two different machines in a real network setting, more details to use their features can be found in their [README](https://github.com/emp-toolkit/emp-ot/blob/master/README.md).