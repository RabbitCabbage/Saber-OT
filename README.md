# Post-quantum Oblivious Transfers
Implementations of post-quantum oblivious transfers under [emp-ot](https://github.com/emp-toolkit/emp-ot) framework, including:
- Our new Naor-Pinkas-like 1-out-of-2 OT protocol based on [SABER](https://github.com/KULeuven-COSIC/SABER);
- A Simplest-like 1-out-of-2 OT protocol based on SABER;
- [Masny-Rindal](https://eprint.iacr.org/2019/706) 1-out-of-2 OT protocol based on SABER;
- Masny-Rindal 1-out-of-2 OT protocol based on [CRYSTALs-Kyber](https://pq-crystals.org/kyber/index.shtml).

The first two protocols are our new protocols inspired by [Naor-Pinkas OT](https://dl.acm.org/doi/pdf/10.5555/365411.365502) and [Simplest OT](https://eprint.iacr.org/2015/267.pdf). In our paper, we mainly focus on the first Naor-Pinkas-like construction. The other two are following the [Masny-Rindal](https://eprint.iacr.org/2019/706) construction from uniform key agreement, as baselines. The last Kyber-based implementation is in the branch `kyber`.

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
./run ./bin/test_saber
```

where `test_saber` is our test program for the SABER-based OT, running our Saber-based OTs (including NPOTs, Simplest OTs, and MR OTs) and some CDH-based classical OTs. The output should be like:

```
connected
connected
128 COOTs (average over 1000 runs):	13.921 ms
128 COOTs (average over 1000 runs):	13.91 ms
128 NPOTs (average over 1000 runs):	14.667 ms
128 NPOTs (average over 1000 runs):	14.678 ms
128 NP Saber OTs (average over 1000 runs):	8.221 ms
128 NP Saber OTs (average over 1000 runs):	8.235 ms
128 Simplest Saber OTs (average over 1000 runs):	6.349 ms
128 Simplest Saber OTs (average over 1000 runs):	6.36 ms
128 MRSaber OTs (average over 1000 runs):	18.568 ms
128 MRSaber OTs (average over 1000 runs):	18.587 ms
```

Credited to [emp-ot](https://github.com/emp-toolkit/emp-ot) for the OT implementation,
 one can test these OTs on two different machines in a real network setting, more details to use their features can be found in their [README](https://github.com/emp-toolkit/emp-ot/blob/master/README.md).