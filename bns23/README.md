# 2d Poisson equation

## Single locale version

Code taken from https://chapel-lang.org/blog/posts/bns2

```shell
chpl nsPoisson.chpl --fast
./nsPoisson --sourceMag=500 --makePlots=true
open initial.png solution.png
```

## Distributed version

See the `nsPoissonStencil.chpl` file.

See https://chapel-lang.org/blog/posts/bns3
