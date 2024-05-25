# Saber OT
Implementation adapted from [emp-ot](https://github.com/emp-toolkit/emp-ot) and [SABER](https://github.com/KULeuven-COSIC/SABER).

Run the test with following shell commands:
```shell
cmake .
make -j
./run ./bin/test_saber
```
Sample outputs:
```
connected
connected
128 COOTs:      Tests passed.   8 ms
128 COOTs:      Tests passed.   9 ms
128 NPOTs:      Tests passed.   9 ms
128 NPOTs:      Tests passed.   9 ms
128 NP Saber OTs(1 secret):     Tests passed.   5 ms
128 NP Saber OTs(1 secret):     Tests passed.   6 ms
128 Simplest Saber OTs: Tests passed.   5 ms
128 Simplest Saber OTs: Tests passed.   5 ms
```